import os
os.environ['OMP_NUM_THREADS'] = '1'
 
from ._clima import WaterAdiabatClimate, ClimaException
from ._clima import rebin # rebin routine from futils