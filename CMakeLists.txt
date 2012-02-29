#CMake main-buildfile for opentelemac
#Author: Oliver Goethel (Email: coding at ogoethel.de)
#
cmake_minimum_required(VERSION 2.8)
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

find_package(MPI)
# if (MPI_Fortran_FOUND)
# Message ("MPI found")
# endif (MPI_Fortran_FOUND)
#
#Project-data
project(telemac Fortran)
SET(TELEMAC_VERSION_MAJOR 6)
SET(TELEMAC_VERSION_MINOR 1)
SET(TELEMAC_VERSION v${TELEMAC_VERSION_MAJOR}p${TELEMAC_VERSION_MINOR})

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
#
set(CMAKE_INCLUDE_CURRENT_DIR On)

add_subdirectory(bief/sources)
add_subdirectory(damocles/sources)
add_subdirectory(paravoid/sources)
add_subdirectory(mumpsvoid/sources)
add_subdirectory(special/sources)

if (MPI_Fortran_FOUND)
add_subdirectory(parallel/sources)
endif (MPI_Fortran_FOUND)

add_subdirectory(artemis/sources)
add_subdirectory(postel3d/sources)
add_subdirectory(stbtel/sources)
add_subdirectory(spartacus2d/sources)
add_subdirectory(sisyphe/sources)
add_subdirectory(tomawac/sources)
add_subdirectory(telemac2d/sources)
add_subdirectory(telemac3d/sources)
add_subdirectory(bin/sel2vtk_dir)
add_subdirectory(bin/conviee)


#Dependencies
add_dependencies(special${TELEMAC_VERSION} bief${TELEMAC_VERSION} )

add_dependencies(sisyphe${TELEMAC_VERSION}
                 bief${TELEMAC_VERSION} damo${TELEMAC_VERSION} paravoid${TELEMAC_VERSION}
                 mumpsvoid${TELEMAC_VERSION} special${TELEMAC_VERSION} )

add_dependencies(sisyphe${TELEMAC_VERSION}.exe sisyphe${TELEMAC_VERSION} )

add_dependencies(telemac2d${TELEMAC_VERSION} 
                 bief${TELEMAC_VERSION} damo${TELEMAC_VERSION} paravoid${TELEMAC_VERSION}
                 mumpsvoid${TELEMAC_VERSION} special${TELEMAC_VERSION} sisyphe${TELEMAC_VERSION}
                 tomawac${TELEMAC_VERSION} )

add_dependencies(telemac2d${TELEMAC_VERSION}.exe telemac2d${TELEMAC_VERSION} )

add_dependencies(telemac3d${TELEMAC_VERSION} 
                 bief${TELEMAC_VERSION} damo${TELEMAC_VERSION} paravoid${TELEMAC_VERSION}
                 mumpsvoid${TELEMAC_VERSION} special${TELEMAC_VERSION} sisyphe${TELEMAC_VERSION}
                 telemac2d${TELEMAC_VERSION} tomawac${TELEMAC_VERSION} )

add_dependencies(telemac3d${TELEMAC_VERSION}.exe telemac3d${TELEMAC_VERSION} )

add_dependencies(tomawac${TELEMAC_VERSION}
                 bief${TELEMAC_VERSION} damo${TELEMAC_VERSION} paravoid${TELEMAC_VERSION}
                 mumpsvoid${TELEMAC_VERSION} special${TELEMAC_VERSION} )

add_dependencies(tomawac${TELEMAC_VERSION}.exe tomawac${TELEMAC_VERSION} )

add_dependencies(artemis${TELEMAC_VERSION}
                 bief${TELEMAC_VERSION} damo${TELEMAC_VERSION} paravoid${TELEMAC_VERSION}
                 mumpsvoid${TELEMAC_VERSION} special${TELEMAC_VERSION} )

add_dependencies(artemis${TELEMAC_VERSION}.exe artemis${TELEMAC_VERSION} )

add_dependencies(spartacus2d${TELEMAC_VERSION}
                 bief${TELEMAC_VERSION} damo${TELEMAC_VERSION} paravoid${TELEMAC_VERSION}
                 mumpsvoid${TELEMAC_VERSION} special${TELEMAC_VERSION} )

add_dependencies(spartacus2d${TELEMAC_VERSION}.exe spartacus2d${TELEMAC_VERSION} )

add_dependencies(postel3d${TELEMAC_VERSION}
                 bief${TELEMAC_VERSION} damo${TELEMAC_VERSION} paravoid${TELEMAC_VERSION}
                 special${TELEMAC_VERSION} )

add_dependencies(postel3d${TELEMAC_VERSION}.exe postel3d${TELEMAC_VERSION} )

add_dependencies(stbtel${TELEMAC_VERSION}
                 bief${TELEMAC_VERSION} damo${TELEMAC_VERSION} paravoid${TELEMAC_VERSION}
                 special${TELEMAC_VERSION} )

add_dependencies(stbtel${TELEMAC_VERSION}.exe stbtel${TELEMAC_VERSION} )



#Parallel dependencies
if (MPI_Fortran_FOUND)
  add_dependencies(artemis${TELEMAC_VERSION}_MP.exe artemis${TELEMAC_VERSION} parallel${TELEMAC_VERSION})
  add_dependencies(sisyphe${TELEMAC_VERSION}_MP.exe sisyphe${TELEMAC_VERSION} parallel${TELEMAC_VERSION})
  add_dependencies(telemac2d${TELEMAC_VERSION}_MP.exe telemac2d${TELEMAC_VERSION} parallel${TELEMAC_VERSION})
  add_dependencies(telemac3d${TELEMAC_VERSION}_MP.exe telemac3d${TELEMAC_VERSION} parallel${TELEMAC_VERSION})
  add_dependencies(tomawac${TELEMAC_VERSION}_MP.exe tomawac${TELEMAC_VERSION} parallel${TELEMAC_VERSION})
  add_dependencies(spartacus2d${TELEMAC_VERSION}_MP.exe spartacus2d${TELEMAC_VERSION} parallel${TELEMAC_VERSION})
endif (MPI_Fortran_FOUND)
