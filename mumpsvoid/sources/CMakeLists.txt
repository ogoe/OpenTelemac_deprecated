#
FILE(GLOB CTR_SOURCES mumpsvoid/sources *.f)
#define the library
add_library(mumpsvoid${TELEMAC_VERSION} STATIC ${CTR_SOURCES})

set_target_properties(mumpsvoid${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS mumpsvoid${TELEMAC_VERSION}
            RUNTIME DESTINATION mumpsvoid/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION mumpsvoid/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION mumpsvoid/${Fortran_COMPILER_NAME})


