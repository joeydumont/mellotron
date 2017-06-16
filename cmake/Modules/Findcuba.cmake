# --------------------------------------------------------------------------- #
# Author:       Denis Gagnon                  <denis.gagnon@emt.inrs.ca>      #
# Date created: 2017-06-13                                                    #
# Description:  Attempts to find the cuba library.                            #
# ----------------------------------------------------------------------------#

# -- LibFindMacros for convenience
# https://github.com/Tronic/cmake-modules/blob/master/LibFindMacros.cmake
include(LibFindMacros)

# -- We use pkg-config to give information about the
libfind_pkg_check_modules(cuba_PKGCONF Cubature)

find_path(cuba_INCLUDE_DIR
  NAMES cuba.h
  PATHS ${cuba_PKGCONF_INCLUDE_DIRS}
  HINT external/Cuba
)

# Finally the library itself
find_library(cuba_LIBRARY 
    NAMES cuba libcuba
    PATHS ${CUBATURE_PKGCONF_LIBRARY_DIRS}
    HINT external/Cuba
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(cuba_PROCESS_INCLUDE cuba_INCLUDE_DIR)
set(cuba_PROCESS_LIB cuba_LIBRARY)
libfind_process(cuba)
