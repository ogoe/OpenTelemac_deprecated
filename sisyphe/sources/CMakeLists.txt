#CMake buildfile for telemac/sisyphe
#Author: Oliver Goethel (Email: coding at ogoethel.de)
#
include_directories("${PROJECT_BINARY_DIR}/bief/sources")
#
FILE(GLOB CTR_SOURCES sisyphe/sources *.f)
#
#
# file(GLOB to_remove bedload_interact_CV.f)
# #
# list(REMOVE_ITEM CTR_SOURCES ${to_remove}) 
#
#define the library
add_library(sisyphe${TELEMAC_VERSION} STATIC ${CTR_SOURCES})


set_target_properties(sisyphe${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS sisyphe${TELEMAC_VERSION}
            RUNTIME DESTINATION sisyphe/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION sisyphe/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION sisyphe/${Fortran_COMPILER_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/declarations_sisyphe.mod ${CMAKE_CURRENT_BINARY_DIR}/interface_sisyphe_bedload.mod
              ${CMAKE_CURRENT_BINARY_DIR}/interface_sisyphe.mod ${CMAKE_CURRENT_BINARY_DIR}/interface_sisyphe_suspension.mod
        DESTINATION sisyphe/${Fortran_COMPILER_NAME})


link_directories("${PROJECT_BINARY_DIR}/special/sources"
                 "${PROJECT_BINARY_DIR}/paravoid/sources"
                 "${PROJECT_BINARY_DIR}/mumpsvoid/sources"
                 "${PROJECT_BINARY_DIR}/bief/sources" 
                 "${PROJECT_BINARY_DIR}/damocles/sources"
                 "${PROJECT_BINARY_DIR}/sisyphe/sources")

add_executable(sisyphe${TELEMAC_VERSION}.exe homere_sisyphe.f)
target_link_libraries(sisyphe${TELEMAC_VERSION}.exe sisyphe${TELEMAC_VERSION}
                                                    bief${TELEMAC_VERSION}
                                                    damo${TELEMAC_VERSION}
                                                    mumpsvoid${TELEMAC_VERSION}
                                                    paravoid${TELEMAC_VERSION}
                                                    special${TELEMAC_VERSION} )

install(TARGETS sisyphe${TELEMAC_VERSION}.exe
            RUNTIME DESTINATION sisyphe/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION sisyphe/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION sisyphe/${Fortran_COMPILER_NAME})
#
#Parallel:
if (MPI_Fortran_FOUND)

add_executable(sisyphe${TELEMAC_VERSION}_MP.exe homere_sisyphe.f)
target_link_libraries(sisyphe${TELEMAC_VERSION}_MP.exe sisyphe${TELEMAC_VERSION}
                                                    bief${TELEMAC_VERSION}
                                                    damo${TELEMAC_VERSION}
                                                    mumpsvoid${TELEMAC_VERSION}
                                                    parallel${TELEMAC_VERSION}
                                                    special${TELEMAC_VERSION}
                                                    ${MPI_Fortran_LIBRARIES}    )

install(TARGETS sisyphe${TELEMAC_VERSION}_MP.exe
            RUNTIME DESTINATION sisyphe/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION sisyphe/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION sisyphe/${Fortran_COMPILER_NAME})

endif (MPI_Fortran_FOUND)
