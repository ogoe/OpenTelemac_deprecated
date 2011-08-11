#CMake buildfile for telemac/postel3d
#Author: Oliver Goethel (Email: coding at ogoethel.de)
#
include_directories("${PROJECT_BINARY_DIR}/bief/sources")
#
FILE(GLOB CTR_SOURCES postel3d/sources *.f)
#
#
#define the library
add_library(postel3d${TELEMAC_VERSION} STATIC ${CTR_SOURCES})


set_target_properties(postel3d${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS postel3d${TELEMAC_VERSION}
            RUNTIME DESTINATION postel3d/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION postel3d/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION postel3d/${Fortran_COMPILER_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/declarations_postel3d.mod
        DESTINATION postel3d/${Fortran_COMPILER_NAME})


link_directories("${PROJECT_BINARY_DIR}/special/sources"
                 "${PROJECT_BINARY_DIR}/paravoid/sources"
                 "${PROJECT_BINARY_DIR}/bief/sources" 
                 "${PROJECT_BINARY_DIR}/damocles/sources"
                 "${PROJECT_BINARY_DIR}/postel3d/sources")

add_executable(postel3d${TELEMAC_VERSION}.exe homere_postel3d.f)
target_link_libraries(postel3d${TELEMAC_VERSION}.exe postel3d${TELEMAC_VERSION}
                                                    bief${TELEMAC_VERSION}
                                                    damo${TELEMAC_VERSION}
                                                    paravoid${TELEMAC_VERSION}
                                                    special${TELEMAC_VERSION} )

install(TARGETS postel3d${TELEMAC_VERSION}.exe
            RUNTIME DESTINATION postel3d/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION postel3d/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION postel3d/${Fortran_COMPILER_NAME})
