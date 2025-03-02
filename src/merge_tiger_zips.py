from pathlib import Path
import arcpy
from arcgis_tools.static_tools import helper_functions
from zipfile import ZipFile



ROOT_DIR = Path(__file__).parents[0]
BASE_DIR = Path(__file__).parents[1]
WORKING_DIR = Path(r"C:\working")
INPUTS = WORKING_DIR / "streetmap_na/tiger"
OUTPUTS = WORKING_DIR / "streetmap_na/tiger"

FIELDS = ["SHAPE@","LINEARID","FULLNAME","RTTYP","MTFCC"]
CREATE_FIELDS = {"LINEARID": ['TEXT', 22], "FULLNAME": ['TEXT', 100], "RTTYP": ['TEXT', 1], "MTFCC": ['TEXT', 5]}

def main(tigerFolder):
    for f in tigerFolder.glob("*.zip"):
        print(f.name)
        with ZipFile(f, 'r') as zObject:
            tempOutFolder = OUTPUTS / f.name.replace(".zip","")
            zObject.extractall(path=tempOutFolder)

        break


if __name__ == "__main__":

    fgdb = helper_functions.drop_add_fgdb(OUTPUTS,"tiger_routes.gdb")
    sr = helper_functions.get_wgs84_sr()
    arcpy.management.CreateFeatureDataset(str(fgdb), "routes", sr)
    fgdb = fgdb / "routes"
    helper_functions.drop_add_featureclass(fgdb,"all_routes","POLYLINE",sr)
    main(INPUTS)

    print('\ndone')
