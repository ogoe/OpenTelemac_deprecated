set(CMAKE_INCLUDE_CURRENT_DIR On)
include_directories("${PROJECT_BINARY_DIR}/bief/bief_${TELEMAC_VERSION}/sources")
#define the library
add_library(special${TELEMAC_VERSION} STATIC gracestop.f     plante.f )

set_target_properties(special${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS special${TELEMAC_VERSION}
            RUNTIME DESTINATION special/special_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION special/special_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION special/special_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gracestop.mod 
              DESTINATION special/special_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME})



#include_directories("${CMAKE_SOURCE_DIR}")
#SET(LIBRARY_OUTPUT_NAME archiv)
