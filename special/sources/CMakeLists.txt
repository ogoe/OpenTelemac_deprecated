#
include_directories("${PROJECT_BINARY_DIR}/bief/sources")
#
#FILE->GLOB won't work with special, as there are architecture dependend versions
#FILE(GLOB CTR_SOURCES special/sources *.f)
#define the library
set(CTR_SOURCES plante.f gracestop.f)
add_library(special${TELEMAC_VERSION} STATIC ${CTR_SOURCES})

set_target_properties(special${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS special${TELEMAC_VERSION}
            RUNTIME DESTINATION special/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION special/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION special/${Fortran_COMPILER_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gracestop.mod 
              DESTINATION special/${Fortran_COMPILER_NAME})



#include_directories("${CMAKE_SOURCE_DIR}")
#SET(LIBRARY_OUTPUT_NAME archiv)
