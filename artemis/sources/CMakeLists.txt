#
include_directories("${PROJECT_BINARY_DIR}/bief/sources" 
                    "${PROJECT_BINARY_DIR}/special/sources")
#
FILE(GLOB CTR_SOURCES artemis/sources *.f)
## parallel file with mpi.h
file(GLOB to_remove build_global_front.f)
#
list(REMOVE_ITEM CTR_SOURCES ${to_remove}) 
#
#define the library
add_library(artemis${TELEMAC_VERSION} STATIC ${CTR_SOURCES})


set_target_properties(artemis${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS artemis${TELEMAC_VERSION}
            RUNTIME DESTINATION artemis/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION artemis/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION artemis/${Fortran_COMPILER_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/declarations_artemis.mod ${CMAKE_CURRENT_BINARY_DIR}/interface_artemis.mod
        DESTINATION artemis/${Fortran_COMPILER_NAME})


link_directories("${PROJECT_BINARY_DIR}/special/sources"
                 "${PROJECT_BINARY_DIR}/paravoid/sources"
                 "${PROJECT_BINARY_DIR}/mumpsvoid/sources"
                 "${PROJECT_BINARY_DIR}/bief/sources" 
                 "${PROJECT_BINARY_DIR}/damocles/sources"
                 "${PROJECT_BINARY_DIR}/artemis/sources")

add_executable(artemis${TELEMAC_VERSION}.exe homere_artemis.f)
target_link_libraries(artemis${TELEMAC_VERSION}.exe artemis${TELEMAC_VERSION}
                                                    bief${TELEMAC_VERSION}
                                                    damo${TELEMAC_VERSION}
                                                    mumpsvoid${TELEMAC_VERSION}
                                                    paravoid${TELEMAC_VERSION}
                                                    special${TELEMAC_VERSION} )

install(TARGETS artemis${TELEMAC_VERSION}.exe
            RUNTIME DESTINATION artemis/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION artemis/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION artemis/${Fortran_COMPILER_NAME})

