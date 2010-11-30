#define the library
add_library(mumpsvoid${TELEMAC_VERSION} STATIC solve_mumps.f pre4_mumps.f)

set_target_properties(mumpsvoid${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS mumpsvoid${TELEMAC_VERSION}
            RUNTIME DESTINATION mumpsvoid/mumpsvoid_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION mumpsvoid/mumpsvoid_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION mumpsvoid/mumpsvoid_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME})


