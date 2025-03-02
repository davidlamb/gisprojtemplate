import logging
import datetime
message="test"
logging.basicConfig(format="%(levelname)s @ %(asctime)s;  Message: %(message)s",level="INFO")
logging.info("Test INFO")
logging.warning("OMG")
logging.error("error")
dt = datetime.datetime.now()
print("INFO @ {0}; Message: {1}".format(dt,message))