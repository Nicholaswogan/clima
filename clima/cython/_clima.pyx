from numpy cimport ndarray
from libcpp cimport bool
import numpy as np
import os

DEF S_STR_LEN = 20;
DEF ERR_LEN = 1024;

include "futils.pyx"
include "OpticalProperties.pyx"
include "ClimaRadtranWrk.pyx"
include "Radtran.pyx"
include "WaterAdiabatClimate.pyx"

# utils
cdef pystring2cstring(str pystring):
  # add a null c char, and convert to byes
  cdef bytes cstring = (pystring+'\0').encode('utf-8')
  return cstring

cdef c2stringarr(ndarray c_str_arr, int str_len, int arr_len):  
  bs = c_str_arr[:-1].tobytes()
  return [bs[i:i+str_len].decode().strip() for i in range(0, str_len*arr_len, str_len)]

class ClimaException(Exception):
    pass
