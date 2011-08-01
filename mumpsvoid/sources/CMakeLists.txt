#define the library
add_library(mumpsvoid${TELEMAC_VERSION} STATIC solve_mumps.f pre4_mumps.f)

set_target_properties(mumpsvoid${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS mumpsvoid${TELEMAC_VERSION}
            RUNTIME DESTINATION mumpsvoid/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION mumpsvoid/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION mumpsvoid/${Fortran_COMPILER_NAME})


