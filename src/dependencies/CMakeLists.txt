
# fortran-yaml-c
CPMAddPackage(
  NAME fortran-yaml-c
  VERSION 0.2.1
  GITHUB_REPOSITORY "Nicholaswogan/fortran-yaml-c"
  GIT_TAG "v0.2.1"
  EXCLUDE_FROM_ALL ON
)

# h5fortran
CPMAddPackage(
  NAME h5fortran
  VERSION 4.5.0
  OPTIONS
    "CMAKE_POSITION_INDEPENDENT_CODE ON"
  GITHUB_REPOSITORY "geospace-code/h5fortran"
  GIT_TAG "v4.5.0"
  EXCLUDE_FROM_ALL ON
)

# futils
CPMAddPackage(
  NAME futils
  VERSION 0.1.7
  GITHUB_REPOSITORY "Nicholaswogan/futils"
  GIT_TAG "v0.1.7"
  EXCLUDE_FROM_ALL ON
)

# finterp
CPMAddPackage(
  NAME finterp
  VERSION 1.3.0
  GITHUB_REPOSITORY "jacobwilliams/finterp"
  GIT_TAG "1.3.0"
  DOWNLOAD_ONLY ON
)

add_library(finterp ${finterp_SOURCE_DIR}/src/linear_interpolation_module.F90 )

CPMAddPackage(
  NAME dop853
  VERSION 1.2.0
  GITHUB_REPOSITORY "jacobwilliams/dop853"
  GIT_TAG "1.2.0"
  DOWNLOAD_ONLY ON
)

add_library(dop853 
  ${dop853_SOURCE_DIR}/src/dop853_constants.f90
  ${dop853_SOURCE_DIR}/src/dop853_module.f90
)

CPMAddPackage(
  NAME minpack
  VERSION 2.0.0
  GITHUB_REPOSITORY "fortran-lang/minpack"
  GIT_TAG "v2.0.0-rc.1"
  DOWNLOAD_ONLY ON
)

add_library(minpack
  ${minpack_SOURCE_DIR}/src/minpack.f90
)

# mrgrnk is from ORDERPACK 2.0: http://www.fortran-2000.com/rank/#1.1
add_library(mrgrnk mrgrnk_mod.f90) 
