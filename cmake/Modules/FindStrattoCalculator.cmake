# --------------------------------------------------------------------------- #
# Author:       Denis Gagnon               <denis.gagnon@emt.inrs.ca>         #
# Date:         2017-06-07                                                    #
# Description:  CMake module to find StrattoCalculator.                       #
# ----

# -- LibFindMacros for convenience
# https://github.com/Tronic/cmake-modules/blob/master/LibFindMacros.cmake
include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(StrattoCalculator_PKGCONF StrattoCalculator)

# Include dir
find_path(StrattoCalculator_INCLUDE_DIR
  NAMES StrattoCalculator
  PATHS ${StrattoCalculator_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(StrattoCalculator_LIBRARY
    NAMES libStrattoCalculator StrattoCalculator libStrattoCalculator.so libStrattoCalculator.dylib
    PATHS ${StrattoCalculator_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
set(StrattoCalculator_PROCESS_INCLUDE StrattoCalculator_INCLUDE_DIR)
set(StrattoCalculator_PROCESS_LIB StrattoCalculator_LIBRARY)
libfind_process(StrattoCalculator)