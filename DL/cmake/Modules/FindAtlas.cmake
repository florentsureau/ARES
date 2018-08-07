# Find the Atlas (and Lapack) libraries
#
# The following variables are optionally searched for defaults
#  Atlas_ROOT_DIR:            Base directory where all Atlas components are found
#
# The following are set after configuration is done:
#  Atlas_FOUND
#  Atlas_INCLUDE_DIRS
#  Atlas_LIBRARIES
#  Atlas_LIBRARYRARY_DIRS

set(Atlas_INCLUDE_SEARCH_PATHS
  $ENV{Atlas_ROOT_DIR}
  $ENV{Atlas_ROOT_DIR}/include
  /usr/include/atlas
  /usr/include/atlas-base
)

set(Atlas_LIB_SEARCH_PATHS
  $ENV{Atlas_ROOT_DIR}
  $ENV{Atlas_ROOT_DIR}/lib
  /usr/lib/atlas
  /usr/lib/atlas-base
)

find_path(Atlas_CBLAS_INCLUDE_DIR   NAMES cblas.h   PATHS ${Atlas_INCLUDE_SEARCH_PATHS} NO_DEFAULT_PATH)
find_path(Atlas_CLAPACK_INCLUDE_DIR NAMES clapack.h PATHS ${Atlas_INCLUDE_SEARCH_PATHS} NO_DEFAULT_PATH)
find_path(Atlas_LAPACKE_INCLUDE_DIR NAMES lapacke.h PATHS ${Atlas_INCLUDE_SEARCH_PATHS} NO_DEFAULT_PATH)

find_library(Atlas_CBLAS_LIBRARY NAMES  ptcblas_r ptcblas cblas_r cblas PATHS ${Atlas_LIB_SEARCH_PATHS} NO_DEFAULT_PATH)
find_library(Atlas_BLAS_LIBRARY NAMES   atlas_r   atlas                 PATHS ${Atlas_LIB_SEARCH_PATHS} NO_DEFAULT_PATH)
find_library(Atlas_FBLAS_LIBRARY NAMES   f77blas  af77blas f77blas_atlas       PATHS ${Atlas_LIB_SEARCH_PATHS} NO_DEFAULT_PATH)
find_library(Atlas_LAPACK_LIBRARY NAMES alapack_r alapack lapack_atlas lapack  PATHS ${Atlas_LIB_SEARCH_PATHS} NO_DEFAULT_PATH)
find_library(Atlas_LAPACKE_LIBRARY NAMES alapacke_r alapacke lapacke_atlas lapacke  PATHS ${Atlas_LIB_SEARCH_PATHS} NO_DEFAULT_PATH)

set(LOOKED_FOR
  Atlas_CBLAS_INCLUDE_DIR
  Atlas_CLAPACK_INCLUDE_DIR
  Atlas_LAPACKE_INCLUDE_DIR

  Atlas_CBLAS_LIBRARY
  Atlas_FBLAS_LIBRARY
  Atlas_BLAS_LIBRARY
  Atlas_LAPACK_LIBRARY
  Atlas_LAPACKE_LIBRARY
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Atlas DEFAULT_MSG ${LOOKED_FOR})

if(ATLAS_FOUND)
  set(Atlas_INCLUDE_DIR ${Atlas_CBLAS_INCLUDE_DIR} ${Atlas_CLAPACK_INCLUDE_DIR} ${Atlas_LAPACKE_INCLUDE_DIR})
  set(Atlas_LIBRARIES ${Atlas_LAPACKE_LIBRARY} ${Atlas_LAPACK_LIBRARY} ${Atlas_CBLAS_LIBRARY} ${Atlas_FBLAS_LIBRARY} ${Atlas_BLAS_LIBRARY})
  mark_as_advanced(${LOOKED_FOR})
  message(STATUS "ATLAS_INCLUDE=${Atlas_INCLUDE_DIR}")
  message(STATUS "ATLAS_LIBRARIES=${Atlas_LIBRARIES}")
  message(STATUS "Found Atlas (include: ${Atlas_CBLAS_INCLUDE_DIR}, library: ${Atlas_BLAS_LIBRARY})")
endif(ATLAS_FOUND)
