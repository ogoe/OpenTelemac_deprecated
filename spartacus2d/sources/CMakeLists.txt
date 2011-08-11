#CMake buildfile for telemac/spartacus2d
#Author: Oliver Goethel (Email: coding at ogoethel.de)
#
include_directories("${PROJECT_BINARY_DIR}/bief/sources")
#
FILE(GLOB CTR_SOURCES spartacus2d/sources *.f)
#
#
#
#define the library
add_library(spartacus2d${TELEMAC_VERSION} STATIC ${CTR_SOURCES})


set_target_properties(spartacus2d${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS spartacus2d${TELEMAC_VERSION}
            RUNTIME DESTINATION spartacus2d/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION spartacus2d/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION spartacus2d/${Fortran_COMPILER_NAME})


link_directories("${PROJECT_BINARY_DIR}/special/sources"
                 "${PROJECT_BINARY_DIR}/paravoid/sources"
                 "${PROJECT_BINARY_DIR}/bief/sources" 
                 "${PROJECT_BINARY_DIR}/damocles/sources"
                 "${PROJECT_BINARY_DIR}/spartacus2d/sources")

add_executable(spartacus2d${TELEMAC_VERSION}.exe spartacus2d.f)
target_link_libraries(spartacus2d${TELEMAC_VERSION}.exe spartacus2d${TELEMAC_VERSION}
                                                    bief${TELEMAC_VERSION}
                                                    damo${TELEMAC_VERSION}
                                                    mumpsvoid${TELEMAC_VERSION}
                                                    paravoid${TELEMAC_VERSION}
                                                    special${TELEMAC_VERSION} )

install(TARGETS spartacus2d${TELEMAC_VERSION}.exe
            RUNTIME DESTINATION spartacus2d/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION spartacus2d/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION spartacus2d/${Fortran_COMPILER_NAME})
#
#Parallel:
if (MPI_Fortran_FOUND)

add_executable(spartacus2d${TELEMAC_VERSION}_MP.exe spartacus2d.f)
target_link_libraries(spartacus2d${TELEMAC_VERSION}_MP.exe spartacus2d${TELEMAC_VERSION}
                                                    bief${TELEMAC_VERSION}
                                                    damo${TELEMAC_VERSION}
                                                    mumpsvoid${TELEMAC_VERSION}
                                                    parallel${TELEMAC_VERSION}
                                                    special${TELEMAC_VERSION}
                                                    ${MPI_Fortran_LIBRARIES}    )

install(TARGETS spartacus2d${TELEMAC_VERSION}_MP.exe
            RUNTIME DESTINATION spartacus2d/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION spartacus2d/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION spartacus2d/${Fortran_COMPILER_NAME})

endif (MPI_Fortran_FOUND)
