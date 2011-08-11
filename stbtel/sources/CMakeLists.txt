#CMake buildfile for telemac/stbtel
#Author: Oliver Goethel (Email: coding at ogoethel.de)
#
include_directories("${PROJECT_BINARY_DIR}/bief/sources")
#
FILE(GLOB CTR_SOURCES stbtel/sources *.f)
#
#
#define the library
add_library(stbtel${TELEMAC_VERSION} STATIC ${CTR_SOURCES})


set_target_properties(stbtel${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS stbtel${TELEMAC_VERSION}
            RUNTIME DESTINATION stbtel/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION stbtel/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION stbtel/${Fortran_COMPILER_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/declarations_stbtel.mod
        DESTINATION stbtel/${Fortran_COMPILER_NAME})


link_directories("${PROJECT_BINARY_DIR}/special/sources"
                 "${PROJECT_BINARY_DIR}/paravoid/sources"
                 "${PROJECT_BINARY_DIR}/bief/sources" 
                 "${PROJECT_BINARY_DIR}/damocles/sources"
                 "${PROJECT_BINARY_DIR}/stbtel/sources")

add_executable(stbtel${TELEMAC_VERSION}.exe homere_stbtel.f)
target_link_libraries(stbtel${TELEMAC_VERSION}.exe stbtel${TELEMAC_VERSION}
                                                    bief${TELEMAC_VERSION}
                                                    damo${TELEMAC_VERSION}
                                                    paravoid${TELEMAC_VERSION}
                                                    special${TELEMAC_VERSION} )

install(TARGETS stbtel${TELEMAC_VERSION}.exe
            RUNTIME DESTINATION stbtel/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION stbtel/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION stbtel/${Fortran_COMPILER_NAME})

