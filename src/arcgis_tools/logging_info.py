import logging
import arcpy
import datetime
from pathlib import Path
from enum import Enum


class logTypes(Enum):
    """
    Logging output types
    """
    PYLOGGING = 1
    ARCPYMESSAGE = 2
    ALLLOGGING = 3
    INFORMATION = 4
    WARNING = 5
    ERROR = 6

class arcpyMessenger(object):
    def __init__(self):
        self.sendMessage = True

    def info(message):
        dt = datetime.datetime.now()
        arcpy.AddMessage("@ {0} INFO  Message: {1}".format(dt,message))

    def warning(message):
        dt = datetime.datetime.now()
        arcpy.AddWarning("@ {0} WARNING Message: {1}".format(dt,message))

    def error(message):
        dt = datetime.datetime.now()
        arcpy.AddError("@ {0} ERROR Message: {1}".format(dt,message))

class custLogger(object):

    def __init__(self,loggingType:logTypes,logPath:Path=None,logFileName:str=None):
        self.all_logs = []
        self.LOGGINGTYPE = loggingType
        
        if self.LOGGINGTYPE == loggingType.PYLOGGING or self.LOGGINGTYPE == loggingType.ALLLOGGING:
            logger = logging.getLogger('Custom Logger')
            format = logging.Formatter("@ %(asctime)s %(levelname)s Message: %(message)s")
            logger.setLevel(logging.INFO)
            shandler = logging.StreamHandler()
            shandler.setFormatter(format)
            shandler.setLevel(logging.INFO)
            logger.addHandler(shandler)
            if logPath:
                file_output = logPath / logFileName
                handler = logging.FileHandler(str(file_output))
                handler.setFormatter(format)
                handler.setLevel(logging.INFO)
                logger.addHandler(handler)
            
            self.all_logs.append(logger)


        elif self.LOGGINGTYPE == loggingType.ARCPYMESSAGE or self.LOGGINGTYPE == loggingType.ALLLOGGING:
            logger = arcpyMessenger()
            self.all_logs.append(logger)      


    def send_message(self,message:str,message_type:logTypes):
        for logger in self.all_logs:
            if message_type == logTypes.INFORMATION:
                logger.info(message)
            if message_type == logTypes.WARNING:
                logger.warning(message)
            if message_type == logTypes.ERROR:
                logger.error(message)
