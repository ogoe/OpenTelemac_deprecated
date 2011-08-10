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


