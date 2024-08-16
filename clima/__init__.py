import os
os.environ['OMP_NUM_THREADS'] = '1'
 
from ._clima import AdiabatClimate, ClimaException
from ._clima import rebin, rebin_with_errors, grid_at_resolution, make_bins # rebin routine from futils
from ._clima import __version__