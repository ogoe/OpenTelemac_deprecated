cmake_minimum_required(VERSION 2.6)
#use Fortran
enable_language (Fortran)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME STREQUAL "gfortran")
# gfortran
#  set (CMAKE_Fortran_FLAGS_RELEASE "-fconvert=big-endian -frecord-marker=4 -funroll-all-loops -fno-f2c -O3")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-fconvert=big-endian -frecord-marker=4 -fno-f2c -O0 -g")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m64")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -m64 -O0 -g")
elseif (Fortran_COMPILER_NAME STREQUAL "g95")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
else (Fortran_COMPILER_NAME STREQUAL "gfortran")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME STREQUAL "gfortran")

#Project-data
project(telemac Fortran)
SET(TELEMAC_VERSION_MAJOR 6)
SET(TELEMAC_VERSION_MINOR 0)
SET(TELEMAC_VERSION v${TELEMAC_VERSION_MAJOR}p${TELEMAC_VERSION_MINOR})

set(CMAKE_INSTALL_PREFIX $ENV{PROJECT})

add_subdirectory(bief/bief_${TELEMAC_VERSION}/sources)
add_subdirectory(damocles/damo_${TELEMAC_VERSION}/sources)
add_subdirectory(paravoid/paravoid_${TELEMAC_VERSION}/sources)
add_subdirectory(mumpsvoid/mumpsvoid_${TELEMAC_VERSION}/sources)
add_subdirectory(special/special_${TELEMAC_VERSION}/sources)
add_subdirectory(artemis/arte_${TELEMAC_VERSION}/sources)
add_subdirectory(sisyphe/sisyphe_${TELEMAC_VERSION}/sources)
add_subdirectory(tomawac/toma_${TELEMAC_VERSION}/sources)
add_subdirectory(telemac2d/tel2d_${TELEMAC_VERSION}/sources)
add_subdirectory(bin/sel2vtk_dir)

#Dependencies
add_dependencies(special${TELEMAC_VERSION} bief${TELEMAC_VERSION} )

add_dependencies(sisyphe${TELEMAC_VERSION}
                 bief${TELEMAC_VERSION} damo${TELEMAC_VERSION} paravoid${TELEMAC_VERSION}
                 mumpsvoid${TELEMAC_VERSION} special${TELEMAC_VERSION} )

add_dependencies(sisyphe${TELEMAC_VERSION}.exe sisyphe${TELEMAC_VERSION} )

add_dependencies(telemac2d${TELEMAC_VERSION} 
                 bief${TELEMAC_VERSION} damo${TELEMAC_VERSION} paravoid${TELEMAC_VERSION}
                 mumpsvoid${TELEMAC_VERSION} special${TELEMAC_VERSION} sisyphe${TELEMAC_VERSION})

add_dependencies(telemac2d${TELEMAC_VERSION}.exe telemac2d${TELEMAC_VERSION} )

add_dependencies(tomawac${TELEMAC_VERSION}
                 bief${TELEMAC_VERSION} damo${TELEMAC_VERSION} paravoid${TELEMAC_VERSION}
                 mumpsvoid${TELEMAC_VERSION} special${TELEMAC_VERSION} )

add_dependencies(tomawac${TELEMAC_VERSION}.exe tomawac${TELEMAC_VERSION} )

add_dependencies(artemis${TELEMAC_VERSION}
                 bief${TELEMAC_VERSION} damo${TELEMAC_VERSION} paravoid${TELEMAC_VERSION}
                 mumpsvoid${TELEMAC_VERSION} special${TELEMAC_VERSION} )

add_dependencies(artemis${TELEMAC_VERSION}.exe artemis${TELEMAC_VERSION} )

