cmake_minimum_required(VERSION 2.6)
project(ISAP_DL CXX Fortran)

## SYSTEM TYPE
if(WIN32)
    message(STATUS "WIN32 platform")
endif(WIN32)
if(APPLE)
    message(STATUS "APPLE platform")
endif(APPLE)
if(NOT WIN32)
    include(FindPkgConfig)
endif(NOT WIN32)

## CMAKE OPTIONS
option(WITH_ENABLE_GPROF                                "Enable profiling output"   OFF)
option(WITH_USING_ARMADILLO                             "Enable Armadillo if available" ON)
option(WITH_TIMEINFOMANAGER                             "Enable Timing" OFF)
option(WITH_ENABLE_VERBOSE_DEBUG_IN_SparseCoding        "Enable debug verbosity in SparseCoding"   OFF)
option(WITH_ENABLE_VERBOSE_DEBUG_IN_DictionaryLearning  "Enable debug verbosity in DictionaryLearning"   OFF)
option(WITH_ENABLE_VERBOSE_INFO_IN_MatrixDecomposition  "Enable Ingo verbosity in MatrixDecomposition"   OFF)

## APPLY CMAKE OPTIONS
if(WITH_ENABLE_GPROF)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -pg")
endif()
if(WITH_TIMEINFOMANAGER)
    add_definitions(-D ENABLE_TIMEINFOMANAGER)
endif()
if(WITH_ENABLE_VERBOSE_DEBUG_IN_SparseCoding)
    add_definitions(-D ENABLE_VERBOSE_DEBUG_IN_SparseCoding)
endif()
if(WITH_ENABLE_VERBOSE_DEBUG_IN_DictionaryLearning)
    add_definitions(-D ENABLE_VERBOSE_DEBUG_IN_DictionaryLearning)
endif()
if(WITH_ENABLE_VERBOSE_INFO_IN_MatrixDecomposition)
    add_definitions(-D ENABLE_VERBOSE_INFO_IN_MatrixDecomposition)
endif()

# List source files
file(GLOB_RECURSE LOCAL_SOURCES src/*.c* src/libdl/[a-Z]*.c* src/libtools/*.c* src/libtime/[a-Z]*.c* )
# List header files (if exist)
file(GLOB_RECURSE LOCAL_HEADERS src/*.h* src/libdl/*.h* src/libtools/*.h* src/libtime/[a-Z]*.h* )

# to run in debug mode, call cmake with command:
#       cmake -DCMAKE_BUILD_TYPE=Debug

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if(DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
else(DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
endif(DEFINED CMAKE_BUILD_TYPE)


if(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message("-- GNU Compiler Found")
    set(FORTRAN_FOUND 1)
else(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message("==================================================")
    message("NO GNU FORTRAN COMPILER. Few binaries will not be compiled.")
    message("==================================================")
endif(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER_ID MATCHES "GNU")


#####################################
#FIND EXTERNAL PACKAGES
# Find libcfitsio
# Locate cfitsio using pkg-config
if(WIN32)
    find_package(CFITSIO REQUIRED cfitsio)
else(WIN32)
    pkg_check_modules(CFITSIO REQUIRED cfitsio)
endif(WIN32)
include_directories(${CFITSIO_INCLUDE_DIRS})
link_directories(${CFITSIO_LIBRARY_DIRS})

# Locate GSL using pkg-config
if(WIN32)
    find_package(GSL REQUIRED gsl)
else(WIN32)
    pkg_check_modules(GSL REQUIRED gsl)
endif(WIN32)
if (GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    link_directories(${GSL_LIBRARY_DIRS})
endif(GSL_FOUND)

# Locate the Armadillo library # (RECOMMANDED)
if(WITH_USING_ARMADILLO)
    find_package(Armadillo REQUIRED CONFIG HINTS $ENV{ARMAPATH})
    message(STATUS "Found Armadillo ? ${Armadillo_DIR} ${Armadillo_CONFIG}")
    set(ARMADILLO_FOUND $ARMADILLO_FOUND)
    #if (Armadillo_FOUND)
    #    message(STATUS "Found Armadillo : ${ARMADILLO_LIBRARIES}")
    #    include_directories(${ARMADILLO_INCLUDE_DIRS})
    #endif (Armadillo_FOUND)
else(WITH_USING_ARMADILLO)
    message(STATUS "Disabled Armadillo detection using WITH_USING_ARMADILLO")
    set(ARMADILLO_FOUND 0)
endif(WITH_USING_ARMADILLO)


find_package(Atlas) #(RECOMMANDED IF NO ARMADILLO, SHOULD HAVE LAPACKE COMPILED)
if (ATLAS_FOUND)
    include_directories(${Atlas_INCLUDE_DIR})
    message(STATUS "ATLAS INCLUDE:${Atlas_INCLUDE_DIR}")
endif (ATLAS_FOUND)

# tclap
include_directories($ENV{TCLAP_INCLUDE_DIRS})

#####################################
# Compilation flags
if(APPLE)
   set(CMAKE_CXX_FLAGS "-DMACOS -Wno-write-strings -fPIC -fopenmp -std=c++11")
   set(CMAKE_CXX_FLAGS_RELEASE "-g0 -O2  -fomit-frame-pointer -DNDEBUG")
   set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb3 -Wall -DENABLE_VERBOSE_DEBUG")
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -ggdb3 -Wall")
   set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-fopenmp -fno-inline -femit-class-debug-always -fno-default-inline  -fno-builtin -ansi -fkeep-inline-functions -fno-function-cse  -fno-peephole ")
   set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "-fopenmp")
   if (FORTRAN_FOUND)
        set(CMAKE_Fortran_FLAGS "-DGFORTRAN -fno-second-underscore -fopenmp -fPIC -frecursive")
        set(CMAKE_Fortran_FLAGS_DEBUG "-fbounds-check -g")
        set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp")
        set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -fbounds-check -ggdb3 -fopenmp")
   	set(CMAKE_EXE_LINKER_FLAGS "-lgfortran")
   endif(FORTRAN_FOUND)
else(APPLE)
    set(CMAKE_CXX_FLAGS "-Wno-write-strings -Wl,--allow-multiple-definition -fopenmp -std=c++11 -fpic -fopenmp")
    set(CMAKE_CXX_FLAGS_RELEASE "-g0 -ffast-math -O2 -fomit-frame-pointer -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-ggdb3 -Wall -g -DENABLE_VERBOSE_DEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -ggdb3 -Wall")
    set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "-fopenmp")
    if (FORTRAN_FOUND)
        set(CMAKE_Fortran_FLAGS "-DGFORTRAN -fno-second-underscore -fPIC -frecursive -fopenmp")
        set(CMAKE_Fortran_FLAGS_DEBUG "-fbounds-check")
        set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
        set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -fbounds-check")
   	set(CMAKE_EXE_LINKER_FLAGS "-lgfortran")
    endif(FORTRAN_FOUND)
endif(APPLE)

#ADD FLAGS IF RECOMMANDED LIBRARIES FOUND
if(ARMADILLO_FOUND)
    add_definitions(-D_ARMA_DEC)
    message(STATUS "Add-def _ARMA_DEC to ${PROJECT_NAME}")
endif(ARMADILLO_FOUND)
if(ATLAS_FOUND)
    add_definitions(-D_ATLAS_DEC)
    message(STATUS "Add-def _ATLAS_DEC to ${PROJECT_NAME}")
endif(ATLAS_FOUND)


#
# Build tools library
#
FILE(GLOB src_tools "${PROJECT_SOURCE_DIR}/src/libtools/*.cc" "${PROJECT_SOURCE_DIR}/src/libtools/*.c*" "${PROJECT_SOURCE_DIR}/src/libtools/*.h*")
include_directories("${PROJECT_SOURCE_DIR}/src/libtools")
add_library(tools STATIC ${src_tools})
target_link_libraries(tools ${CFITSIO_LIBRARIES})

#
# Build libtime library
#
FILE(GLOB src_time "${PROJECT_SOURCE_DIR}/src/libtime/*.c*" "${PROJECT_SOURCE_DIR}/src/libtime/*.h*")
include_directories("${PROJECT_SOURCE_DIR}/src/libtime")
add_library(libtime STATIC ${src_time})

#
# Build dl library
#
set(DL_EXT_LINK "")
if(WIN32)
    if(ARMADILLO_FOUND)
        message(STATUS "adding link ${ARMADILLO_LIBRARIES}")
        set(DL_EXT_LINK )
    endif(ARMADILLO_FOUND)
else(WIN32)
    if(ARMADILLO_FOUND)
        set(DL_EXT_LINK ${ARMADILLO_LIBRARIES})
    endif(ARMADILLO_FOUND)
    if(ATLAS_FOUND)
        set(DL_EXT_LINK ${DL_EXT_LINK} ${Atlas_LIBRARIES})
        if (FORTRAN_FOUND)
            set(DL_EXT_LINK ${DL_EXT_LINK} gfortran)
        endif(FORTRAN_FOUND)
    endif(ATLAS_FOUND)
endif(WIN32)
message(STATUS "DL_EXT=${DL_EXT_LINK}")
FILE(GLOB src_dl  "${PROJECT_SOURCE_DIR}/src/libdl/*.cpp")
include_directories("${PROJECT_SOURCE_DIR}/src/libdl")
add_library(dl STATIC ${src_dl})
target_link_libraries(dl armadillo tools libtime ${GSL_LIBRARIES} ${DL_EXT_LINK})

#
# Compile and link executables
#
add_executable(DL_test_matdec  ${PROJECT_SOURCE_DIR}/src/DL_test_matrix_decomposition.cpp )
add_executable(DL_learning  ${PROJECT_SOURCE_DIR}/src/DL_learning.cpp )
add_executable(sparse_decomp  ${PROJECT_SOURCE_DIR}/src/Sparse_decomposition.cpp)

target_link_libraries(DL_test_matdec dl)
target_link_libraries(DL_learning dl)
target_link_libraries(sparse_decomp dl)

set(BINDL    DL_test_matdec DL_learning sparse_decomp)

#
# Install (by default in the project directory)
#

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

#install headers
FILE(GLOB inc_lib "${PROJECT_SOURCE_DIR}/src/libtools/*.h*" "${PROJECT_SOURCE_DIR}/src/libdl/*.h*" "${PROJECT_SOURCE_DIR}/src/TimeInfoManager/*.h*")
INSTALL(FILES ${inc_lib} DESTINATION include)

# Install libraries
INSTALL(TARGETS dl tools DESTINATION lib)

# install executables
INSTALL(TARGETS ${BINDL} DESTINATION bin)
