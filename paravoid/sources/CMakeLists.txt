#
FILE(GLOB CTR_SOURCES paravoid/sources *.f)
#define the library
add_library(paravoid${TELEMAC_VERSION} STATIC ${CTR_SOURCES})

set_target_properties(paravoid${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS paravoid${TELEMAC_VERSION}
            RUNTIME DESTINATION paravoid/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION paravoid/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION paravoid/${Fortran_COMPILER_NAME})


