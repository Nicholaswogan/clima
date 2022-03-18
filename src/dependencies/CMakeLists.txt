
# fortran-yaml-c
CPMAddPackage(
  NAME fortran-yaml-c
  VERSION 0.1.2
  GITHUB_REPOSITORY "Nicholaswogan/fortran-yaml-c"
  GIT_TAG "v0.1.2"
  EXCLUDE_FROM_ALL ON
)

# h5fortran
CPMAddPackage(
  NAME h5fortran
  VERSION 4.4.4
  GITHUB_REPOSITORY "geospace-code/h5fortran"
  GIT_TAG "v4.4.4"
  EXCLUDE_FROM_ALL ON
)

# futils
CPMAddPackage(
  NAME futils
  VERSION 0.1.0
  URL "https://github.com/Nicholaswogan/futils/archive/7fc7a57eaee65ad6a825e7fed584eb8a0269d020.zip"
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

add_library(finterp ${finterp_SOURCE_DIR}/src/linear_interpolation_module.F90)
