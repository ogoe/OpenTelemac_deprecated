#CMake buildfile for telemac/bief
#Author: Oliver Goethel (Email: coding at ogoethel.de)
#
#define the library
FILE(GLOB CTR_SOURCES bief/sources *.f)

add_library(bief${TELEMAC_VERSION} STATIC ${CTR_SOURCES})

set_target_properties(bief${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS bief${TELEMAC_VERSION}
            RUNTIME DESTINATION bief/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION bief/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION bief/${Fortran_COMPILER_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/bief.mod ${CMAKE_CURRENT_BINARY_DIR}/bief_def.mod 
                 ${CMAKE_CURRENT_BINARY_DIR}/declarations_telemac.mod ${CMAKE_CURRENT_BINARY_DIR}/m_med.mod ${CMAKE_CURRENT_BINARY_DIR}/streamline.mod
            DESTINATION bief/${Fortran_COMPILER_NAME})



#include_directories("${CMAKE_SOURCE_DIR}")
#SET(LIBRARY_OUTPUT_NAME archiv)
