## Python script to test running eco functions ##
## Run python scripts using Python 3.5 or newer ##
from ctypes import cdll
lib = cdll.LoadLibrary('../../libeco.so')
print(lib)
success = lib.blah()
