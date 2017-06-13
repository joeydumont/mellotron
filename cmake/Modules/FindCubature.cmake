# --------------------------------------------------------------------------- #
# Author:       Joey Dumont                   <joey.dumont@gmail.com>         #
# Date:         2016-10-26                                                    #
# Description:  Attempts to find the Cubature library.                        #
# License:      CC0                                                           #
#               <https://creativecommons.org/publicdomain/zero/1.0/>          #
#                                                                             #
# This CMake module attempts to find the Cubature library. Once it is done,   #
# the following will be defined:                                              #
#     - CUBATURE_FOUND: if system has Cubature and everything is parsed;      #
#     - CUBATURE_INCLUDE_DIR: folder where headers are.                       #
#     - CUBATURE_LIBRARY: library against which to link.                      #
# ----------------------------------------------------------------------------#

# -- LibFindMacros for convenience
# https://github.com/Tronic/cmake-modules/blob/master/LibFindMacros.cmake
include(LibFindMacros)

# -- We use pkg-config to give information about the
libfind_pkg_check_modules(CUBATURE_PKGCONF Cubature)

find_path(CUBATURE_INCLUDE_DIR
  NAMES Cubature
  PATHS ${CUBATURE_PKGCONF_INCLUDE_DIRS}
  HINT external/Cubature/include
)

# Finally the library itself
find_library(CUBATURE_LIBRARY 
    NAMES Cubature libCubature
    PATHS ${CUBATURE_PKGCONF_LIBRARY_DIRS}
    HINT external/Cubature/
)
  
# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(CUBATURE_PROCESS_INCLUDE CUBATURE_INCLUDE_DIR)
set(CUBATURE_PROCESS_LIB CUBATURE_LIBRARY)
libfind_process(CUBATURE)
