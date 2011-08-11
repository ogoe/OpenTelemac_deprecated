#CMake buildfile for telemac/damocles
#Author: Oliver Goethel (Email: coding at ogoethel.de)
#
FILE(GLOB CTR_SOURCES damocles/sources *.f)
#define the library
add_library(damo${TELEMAC_VERSION} STATIC ${CTR_SOURCES})

set_target_properties(damo${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS damo${TELEMAC_VERSION}
            RUNTIME DESTINATION damocles/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION damocles/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION damocles/${Fortran_COMPILER_NAME})


