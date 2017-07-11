# --------------------------------------------------------------------------- #
# Author:       Denis Gagnon               <denis.gagnon@emt.inrs.ca>         #
# Date:         2017-06-07                                                    #
# Description:  CMake module to find strattocalculator.                       #
# --------------------------------------------------------------------------- #

# -- LibFindMacros for convenience
# https://github.com/Tronic/cmake-modules/blob/master/LibFindMacros.cmake
include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(strattocalculator_PKGCONF strattocalculator)

# Include dir
find_path(strattocalculator_INCLUDE_DIR
  NAMES strattocalculator
  PATHS ${strattocalculator_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(strattocalculator_LIBRARY
    NAMES libStrattoCalculator StrattoCalculator libStrattoCalculator.so libStrattoCalculator.dylib
    PATHS ${strattocalculator_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
set(strattocalculator_PROCESS_INCLUDE strattocalculator_INCLUDE_DIR)
set(strattocalculator_PROCESS_LIB strattocalculator_LIBRARY)
libfind_process(strattocalculator)