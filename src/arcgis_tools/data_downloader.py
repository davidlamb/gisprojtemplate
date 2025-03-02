
from pathlib import Path
from enum import Enum
import requests
import arcpy
import zipfile
from . static_tools import helper_functions

class apiname(Enum):
    """
    API name
    """
    RIDB_API = "ridb_api"
    ALTFUEL_API = "ev_station"
    FRA_API = "fra"
    TNM_API = "tnm_api"
    WARNING = 5
    ERROR = 6


class downloader(object):

    """
    Description:
        Base Class for downloaders
        All url variables should be URL encoded (e.g. space = %20)
        all classes should have a download data overwrite that obtains the data. They may include special functions to post-process the data.
    Args:
        baseFCPath(Path): path to boundary feature. Entire feature class is used, not a single feature
        outDir (Path): files will be downloaded to this director. If the file exists, it will not be downloaded.
        timeout (int): How long before the server responds is the request cancelled.
        checkFileExists(bool): may not be implemented.
    """

    def __init__(self,baseFCPath:Path,outDir:Path,timeout:int=60,checkFileExists:bool=True) -> None:
        self.baseFCPath = baseFCPath
        self.outDir = outDir
        self.timeout = timeout
        self.checkFileExists = checkFileExists
        self.rootURL = ""
        self.files = []

    def extent_as_string_latlong_url(self)->str:
        wgs84 = arcpy.SpatialReference(4326)
        desc = arcpy.Describe(str(self.baseFCPath))
        ext = desc.extent
        ll = arcpy.PointGeometry(arcpy.Point(ext.XMin,ext.YMin),desc.spatialReference)
        ur = arcpy.PointGeometry(arcpy.Point(ext.XMax,ext.YMax),desc.spatialReference)
        llprj = ll.projectAs(wgs84)
        urprj = ur.projectAs(wgs84)
        return "{},{},{},{}".format(llprj.firstPoint.X,llprj.firstPoint.Y,urprj.firstPoint.X,urprj.firstPoint.Y)

    def extent_as_string(fc:Path)->str:
        desc = arcpy.Describe(str(fc))
        rect = "{} {} {} {}".format(desc.extent.XMin,desc.extent.YMin,desc.extent.XMax,desc.extent.YMax)
        return rect
    def buildURL(self):
        pass

    def downloadData(self):
        pass

class tnm_elevation_downloader(downloader):

    def __init__(self, baseFCPath: Path, outDir: Path, timeout: int = 60, checkFileExists: bool = True) -> None:
        """
        Change URLVars["datasets"] to a different DEM dataset. The default is 1/3 arc-second.
        1/9 arc second: National%20Elevation%20Dataset%20%28NED%29%201%2F9%20arc-second
        1 arc second: National%20Elevation%20Dataset%20%28NED%29%201%20arc-second
        
        """
        super().__init__(baseFCPath, outDir, timeout, checkFileExists)
        self.rootURL = "https://tnmaccess.nationalmap.gov/api/v1/products?"
        self.URLVars = {"datasets":"National%20Elevation%20Dataset%20(NED)%201/3%20arc-second"}
        self.URLVars["bbox"] = self.extent_as_string_latlong_url()

    def buildURL(self):
        turl = self.rootURL
        vars = ["{}={}".format(k,v) for k,v in self.URLVars.items()]
        turl += "&".join(vars)
        return turl
    
    def downloadData(self):
        filelist = []
        api_url = self.buildURL()
        response = requests.get(api_url,timeout=self.timeout)
        result = response.json()

        if 'message' in result:
            if result['message'] == 'Endpoint request timed out':
                print("timed out, trying again...")
                response = requests.get(api_url,timeout=self.timeout*2)
                result = response.json()
        #raise Exception("Error getting the DEM from the api.")
        if 'total' not in result:
            raise Exception("Error getting the DEM from the api.")

        print("Total DEM Files that may need to be downloaded: {}".format(result["total"]))

        for item in result["items"]:
            downloadurl = item["downloadURL"]
            name = downloadurl.split("/")[-1]
            if ".tif" in name:
                outfile = self.outDir / name
                if outfile.exists() == False:
                    imgresponse = requests.get(downloadurl)
                    with open(outfile,"wb") as file:
                        file.write(imgresponse.content)
                filelist.append(str(outfile))
            else:
                print("Results are not tif files")
                return filelist
        self.files = filelist

    def mosaic_and_append(self,outws:Path,mosaic_name:str):
        """
        mosaics raster dataset, clips to extent

        Args:
            rasterFiles (list): list of raster paths
            outws (Path): output workspace

        """
        wstype = arcpy.Describe(str(outws)).workspaceType

        desc = arcpy.Describe(str(self.baseFCPath))

        clip_name = mosaic_name + "clp"

        rect = "{} {} {} {}".format(desc.extent.XMin,desc.extent.YMin,desc.extent.XMax,desc.extent.YMax)

        if wstype == "FileSystem":
            clip_name+=".tif"
            mosaic_name+=".tif"

        clipout = outws / clip_name

        if arcpy.Exists(str(clipout)):
            arcpy.Delete_management(str(clipout))

        with arcpy.EnvManager(scratchWorkspace=str(outws), workspace=str(outws)):
            print("Mosaicing Rasters")
            res = arcpy.management.MosaicToNewRaster(self.files, str(outws), mosaic_name , desc.spatialReference, "32_BIT_FLOAT", None, 1, "LAST", "FIRST").getOutput(0)
            print("Clipping to study area Rasters")
            cres = arcpy.management.Clip(res, rect, out_raster=str(clipout), in_template_dataset=str(self.baseFCPath), clipping_geometry= "NONE", maintain_clipping_extent="NO_MAINTAIN_EXTENT").getOutput(0)
            

class tnm_hydrography_huc4(downloader):
    def __init__(self, baseFCPath: Path, outDir: Path, timeout: int = 60, checkFileExists: bool = True) -> None:
        super().__init__(baseFCPath, outDir, timeout, checkFileExists)
        self.rootURL = "https://tnmaccess.nationalmap.gov/api/v1/products?"
        self.URLVars = {"prodExtents":"HU-4%20Subregion","prodFormats":"FileGDB"}
        self.URLVars["bbox"] = self.extent_as_string_latlong_url()
        #https://tnmaccess.nationalmap.gov/api/v1/products?bbox=-113.05825604854886,34.128395659780594,-111.74383373528363,35.15719297261837&prodExtents=HU-4%20Subregion&prodFormats=FileGDB


    def buildURL(self):
        turl = self.rootURL
        vars = ["{}={}".format(k,v) for k,v in self.URLVars.items()]
        turl += "&".join(vars)
        return turl
    
    def downloadData(self):
        filelist = []
        api_url = self.buildURL()
        response = requests.get(api_url,timeout=self.timeout)
        result = response.json()
        if 'message' in result:
            if result['message'] == 'Endpoint request timed out':
                print("timed out, trying again...")
                response = requests.get(api_url,timeout=self.timeout*2)
                result = response.json()
        #raise Exception("Error getting the DEM from the api.")
        if 'total' not in result:
            raise Exception("Error getting the HUC-4 file geodatabase from the api.")

        print("Total Files that may need to be downloaded: {}".format(result["total"]))

        for item in result["items"]:
            if item["extent"] == "HU-4 Subregion" and item["format"] == "FileGDB":
                downloadurl = item["downloadURL"]
                name = downloadurl.split("/")[-1]
                if ".zip" in name:
                    outfile = self.outDir / name
                    if outfile.exists() == False:
                        fresponse = requests.get(downloadurl)
                        with open(outfile,"wb") as file:
                            file.write(fresponse.content)
                    filelist.append(str(outfile))

        self.files = filelist

    def extract_and_merge_flow_lines(self,outws:Path,mergeName:str,waterBody:bool=True,whereClause:str=None,deleteZipFile:bool=True):
        """
        Description:
            Unzips Hydrography geodatabase. Merges the Hydrography/NHDFlowline feature class to a new one. Use the where clause to limit which lines are merged
        Args:
            outws(Path): File geodatabase where the merged data will be saved
            mergeName(str): name of the output feature class
            waterBody(bool): merge waterbodies too
            whereClause(str): filter by this for the NHDFlowline
            deleteZipFile(bool): delete zip file after extracting.
        """
        if len(self.files)==0:
            print("No files...")
            return None
        extractedFiles = []
        print("Extracting zip files.")
        for x in self.files:
            x = Path(x)
            f = x.name.replace('.zip','.gdb')
            extractedFiles.append(self.outDir/f)
            archive = zipfile.ZipFile(x)
            for file in archive.namelist():
                if file.startswith('{}/'.format(f)):
                    archive.extract(file, self.outDir)
                    
            del archive
            if deleteZipFile == True:
                x.unlink()
        waterbodies = []
        files = [str(x/"Hydrography/NHDFlowline") for x in extractedFiles]

        if waterBody == True:
            waterbodies = [str(x/"Hydrography/NHDWaterbody") for x in extractedFiles]
        
        sr = arcpy.Describe(files[0]).spatialReference

        outline = helper_functions.drop_add_featureclass_with_template(outws,mergeName,files[0],"POLYLINE",sr)
        print("Merging flow lines...")
        if whereClause is not None:
            print("Using where clause")
            arcpy.management.Append(files, str(outline), "TEST", None, '', whereClause)
        else:
            print("no where clause")
            arcpy.management.Append(files, str(outline), "TEST", None, '')

        if len(waterbodies)>0:
            print("Merging water bodies...")
            sr = arcpy.Describe(waterbodies[0]).spatialReference
            outpoly= helper_functions.drop_add_featureclass_with_template(outws,f"{mergeName}_wb",waterbodies[0],"POLYGON",sr)
            arcpy.management.Append(waterbodies, str(outpoly), "TEST", None, '')
