import sys
from ._version import __version__

if sys.version_info[0] < 3 and sys.version_info[1] < 7:
    raise Exception('Python 3.4 or greater is required to use this package.')