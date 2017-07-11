# --------------------------------------------------------------------------- #
# Author:       Denis Gagnon               <denis.gagnon@emt.inrs.ca>         #
# Date:         2017-06-07                                                    #
# Description:  CMake module to find meshpi.                                  #
# ----

# -- LibFindMacros for convenience
# https://github.com/Tronic/cmake-modules/blob/master/LibFindMacros.cmake
include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(meshpi_PKGCONF meshpi)

# Include dir
find_path(meshpi_INCLUDE_DIR
  NAMES meshpi
  PATHS ${meshpi_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(meshpi_LIBRARY
    NAMES libMeshPI MeshPI libMeshPI.so libMeshPI.dylib
    PATHS ${meshpi_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
set(meshpi_PROCESS_INCLUDE meshpi_INCLUDE_DIR)
set(meshpi_PROCESS_LIB meshpi_LIBRARY)
libfind_process(meshpi)