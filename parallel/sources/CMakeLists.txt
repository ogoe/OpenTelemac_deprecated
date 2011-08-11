#CMake buildfile for telemac/parallel
#Author: Oliver Goethel (Email: coding at ogoethel.de)
#
#
FILE(GLOB CTR_SOURCES parallel/sources *.f)
#
file(GLOB to_remove partel.f gretel_autop.f gredelhyd_autop.f gredelmet_autop.f 
                    gredelpts_autop.f gredelseg_autop.f diffsel.f)
#
list(REMOVE_ITEM CTR_SOURCES ${to_remove}) 
#define the library
add_library(parallel${TELEMAC_VERSION} STATIC ${CTR_SOURCES})
#
set_target_properties(parallel${TELEMAC_VERSION} PROPERTIES PREFIX "")
#
#target_link_libraries(parallel${TELEMAC_VERSION} ${MPI_Fortran_LIBRARIES})
#
include_directories(${MPI_Fortran_INCLUDE_PATH})
#
if(MPI_COMPILE_FLAGS) 
  set_target_properties(parallel${TELEMAC_VERSION} PROPERTIES 
    COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}") 
endif() 

if(MPI_LINK_FLAGS) 
  set_target_properties(parallel${TELEMAC_VERSION} PROPERTIES 
    LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}") 
endif() 
#
install(TARGETS parallel${TELEMAC_VERSION}
            RUNTIME DESTINATION parallel$/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION parallel$/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION parallel$/${Fortran_COMPILER_NAME})
#
#Parallel executables
#partel is only buiöd if libmetis is found in one of the given diretories
#
add_executable(gretel.exe gretel_autop.f)
#
install(TARGETS gretel.exe
            RUNTIME DESTINATION parallel$/${Fortran_COMPILER_NAME})
#
FIND_LIBRARY(METIS_LIBRARY metis
    /usr/local/lib
    /usr/lib
    ${CMAKE_SOURCE_DIR}/parallel
  )

IF(METIS_LIBRARY)
#
  add_executable(partel.exe partel.f)
#
  target_link_libraries(partel.exe ${METIS_LIBRARY})
#
  install(TARGETS partel.exe
            RUNTIME DESTINATION parallel$/${Fortran_COMPILER_NAME})
#
ENDIF(METIS_LIBRARY)