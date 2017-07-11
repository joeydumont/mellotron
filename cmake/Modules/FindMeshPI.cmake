# --------------------------------------------------------------------------- #
# Author:       Denis Gagnon               <denis.gagnon@emt.inrs.ca>         #
# Date:         2017-06-07                                                    #
# Description:  CMake module to find MeshPI.                                  #
# ----

# -- LibFindMacros for convenience
# https://github.com/Tronic/cmake-modules/blob/master/LibFindMacros.cmake
include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(MeshPI_PKGCONF MeshPI)

# Include dir
find_path(MeshPI_INCLUDE_DIR
  NAMES MeshPI
  PATHS ${MeshPI_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(MeshPI_LIBRARY
    NAMES libMeshPI MeshPI libMeshPI.so libMeshPI.dylib
    PATHS ${MeshPI_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
set(MeshPI_PROCESS_INCLUDE MeshPI_INCLUDE_DIR)
set(MeshPI_PROCESS_LIB MeshPI_LIBRARY)
libfind_process(MeshPI)