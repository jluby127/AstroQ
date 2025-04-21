import os
from kpfcc import driver
print("Importing all of KPF-CC modules")

__version__='2.1.0'
_ROOT = os.path.abspath(os.path.dirname(__file__))
DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)),'data')
