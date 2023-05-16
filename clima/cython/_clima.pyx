from numpy cimport ndarray
from libc.stdint cimport uintptr_t
from libcpp cimport bool
import numpy as np
import ctypes as ct
import os

DEF S_STR_LEN = 20;
DEF ERR_LEN = 1024;

include "futils.pyx"
include "OpticalProperties.pyx"
include "ClimaRadtranWrk.pyx"
include "Radtran.pyx"
include "AdiabatClimate.pyx"

# version
cdef extern void clima_version_get(char *version_c)

def _clima_version():
  cdef char version_c[100+1]
  clima_version_get(version_c)
  return version_c.decode("utf-8").strip()
  
__version__ = _clima_version()

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
