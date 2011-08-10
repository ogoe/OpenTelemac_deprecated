#
include_directories("${PROJECT_BINARY_DIR}/bief/sources")
#
FILE(GLOB CTR_SOURCES tomawac/sources *.f)
#define the library
add_library(tomawac${TELEMAC_VERSION} STATIC ${CTR_SOURCES})


set_target_properties(tomawac${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS tomawac${TELEMAC_VERSION}
            RUNTIME DESTINATION tomawac/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION tomawac/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION tomawac/${Fortran_COMPILER_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/declarations_tomawac.mod
              ${CMAKE_CURRENT_BINARY_DIR}/tomawac_mpi.mod ${CMAKE_CURRENT_BINARY_DIR}/tomawac_mpi_tools.mod
        DESTINATION tomawac/${Fortran_COMPILER_NAME})


link_directories("${PROJECT_BINARY_DIR}/special/sources"
                 "${PROJECT_BINARY_DIR}/paravoid/sources"
                 "${PROJECT_BINARY_DIR}/mumpsvoid/sources"
                 "${PROJECT_BINARY_DIR}/bief/sources" 
                 "${PROJECT_BINARY_DIR}/damocles/sources"
                 "${PROJECT_BINARY_DIR}/tomawac/sources")

add_executable(tomawac${TELEMAC_VERSION}.exe homere_tomawac.f)
target_link_libraries(tomawac${TELEMAC_VERSION}.exe tomawac${TELEMAC_VERSION}
                                                    bief${TELEMAC_VERSION}
                                                    damo${TELEMAC_VERSION}
                                                    mumpsvoid${TELEMAC_VERSION}
                                                    paravoid${TELEMAC_VERSION}
                                                    special${TELEMAC_VERSION} )

install(TARGETS tomawac${TELEMAC_VERSION}.exe
            RUNTIME DESTINATION tomawac/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION tomawac/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION tomawac/${Fortran_COMPILER_NAME})

#
#Parallel:
if (MPI_Fortran_FOUND)

add_executable(tomawac${TELEMAC_VERSION}_MP.exe homere_tomawac.f)
target_link_libraries(tomawac${TELEMAC_VERSION}_MP.exe tomawac${TELEMAC_VERSION}
                                                    bief${TELEMAC_VERSION}
                                                    damo${TELEMAC_VERSION}
                                                    mumpsvoid${TELEMAC_VERSION}
                                                    parallel${TELEMAC_VERSION}
                                                    special${TELEMAC_VERSION}
                                                    ${MPI_Fortran_LIBRARIES}         )

install(TARGETS tomawac${TELEMAC_VERSION}_MP.exe
            RUNTIME DESTINATION tomawac/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION tomawac/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION tomawac/${Fortran_COMPILER_NAME})

endif (MPI_Fortran_FOUND)