from numpy cimport ndarray
from libcpp cimport bool
import numpy as np
import os

DEF S_STR_LEN = 20;
DEF ERR_LEN = 1024;

include "WaterAdiabatClimate.pyx"

# utils
cdef pystring2cstring(str pystring):
  # add a null c char, and convert to byes
  cdef bytes cstring = (pystring+'\0').encode('utf-8')
  return cstring

class ClimaException(Exception):
    pass
