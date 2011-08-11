#CMake buildfile for telemac/telemac3d
#Author: Oliver Goethel (Email: coding at ogoethel.de)
#
include_directories("${PROJECT_BINARY_DIR}/bief/sources"
                    "${PROJECT_BINARY_DIR}/special/sources"
                    "${PROJECT_BINARY_DIR}/sisyphe/sources"
                    "${PROJECT_BINARY_DIR}/tomawac/sources"
                    "${PROJECT_BINARY_DIR}/telemac2d/sources" )
#
FILE(GLOB CTR_SOURCES telemac3d/sources *.f)
#define the library
add_library(telemac3d${TELEMAC_VERSION} STATIC ${CTR_SOURCES})


set_target_properties(telemac3d${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS telemac3d${TELEMAC_VERSION}
            RUNTIME DESTINATION telemac3d/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION telemac3d/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION telemac3d/${Fortran_COMPILER_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/declarations_telemac3d.mod ${CMAKE_CURRENT_BINARY_DIR}/interface_telemac3d.mod
        DESTINATION telemac3d/${Fortran_COMPILER_NAME})


link_directories("${PROJECT_BINARY_DIR}/special/sources"
                 "${PROJECT_BINARY_DIR}/paravoid/sources"
                 "${PROJECT_BINARY_DIR}/mumpsvoid/sources"
                 "${PROJECT_BINARY_DIR}/bief/sources" 
                 "${PROJECT_BINARY_DIR}/damocles/sources"
                 "${PROJECT_BINARY_DIR}/sisyphe/sources"
                 "${PROJECT_BINARY_DIR}/tomawac/sources"
                 "${PROJECT_BINARY_DIR}/telemac2d/sources"
                 "${PROJECT_BINARY_DIR}/telemac3d/sources")

add_executable(telemac3d${TELEMAC_VERSION}.exe homere_telemac3d.f)
target_link_libraries(telemac3d${TELEMAC_VERSION}.exe telemac3d${TELEMAC_VERSION}
                                                    sisyphe${TELEMAC_VERSION}
                                                    telemac2d${TELEMAC_VERSION}
                                                    tomawac${TELEMAC_VERSION}
                                                    bief${TELEMAC_VERSION}
                                                    damo${TELEMAC_VERSION}
                                                    mumpsvoid${TELEMAC_VERSION}
                                                    paravoid${TELEMAC_VERSION}
                                                    special${TELEMAC_VERSION} )

install(TARGETS telemac3d${TELEMAC_VERSION}.exe
            RUNTIME DESTINATION telemac3d/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION telemac3d/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION telemac3d/${Fortran_COMPILER_NAME})
#
#Parallel:
if (MPI_Fortran_FOUND)

add_executable(telemac3d${TELEMAC_VERSION}_MP.exe homere_telemac3d.f)
target_link_libraries(telemac3d${TELEMAC_VERSION}_MP.exe telemac3d${TELEMAC_VERSION}
                                                    sisyphe${TELEMAC_VERSION}
                                                    telemac2d${TELEMAC_VERSION}
                                                    tomawac${TELEMAC_VERSION}
                                                    bief${TELEMAC_VERSION}
                                                    damo${TELEMAC_VERSION}
                                                    mumpsvoid${TELEMAC_VERSION}
                                                    parallel${TELEMAC_VERSION}
                                                    special${TELEMAC_VERSION}
                                                    ${MPI_Fortran_LIBRARIES}    )

install(TARGETS telemac3d${TELEMAC_VERSION}_MP.exe
            RUNTIME DESTINATION telemac3d/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION telemac3d/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION telemac3d/${Fortran_COMPILER_NAME})

endif (MPI_Fortran_FOUND)
