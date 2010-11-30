#
set(CMAKE_INCLUDE_CURRENT_DIR On)
include_directories("${PROJECT_BINARY_DIR}/bief/bief_${TELEMAC_VERSION}/sources" 
                    "${PROJECT_BINARY_DIR}/special/special_${TELEMAC_VERSION}/sources")
#
#define the library
add_library(artemis${TELEMAC_VERSION} STATIC declarations_artemis.f interface_artemis.f 
artemis.f        berkho.f         borh.f           calcfw.f         
calcmn.f         calcqb.f         calctm.f         calcue.f         
calre2.f         calres.f         cntpre.f         condih.f         
corfon.f         dirale.f         dismoy.f         entart.f         
fwspec.f         homere_artemis.f lecdon_artemis.f leclim_artemis.f 
masque_artemis.f nomvar_artemis.f perale.f  phbor.f   point_artemis.f  
radia1.f         radia2.f         spd.f            spe.f            
utimp.f)


set_target_properties(artemis${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS artemis${TELEMAC_VERSION}
            RUNTIME DESTINATION artemis/arte_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION artemis/arte_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION artemis/arte_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/declarations_artemis.mod ${CMAKE_CURRENT_BINARY_DIR}/interface_artemis.mod
        DESTINATION artemis/arte_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME})


link_directories("${PROJECT_BINARY_DIR}/special/special_${TELEMAC_VERSION}/sources"
                 "${PROJECT_BINARY_DIR}/paravoid/paravoid_${TELEMAC_VERSION}/sources"
                 "${PROJECT_BINARY_DIR}/mumpsvoid/mumpsvoid_${TELEMAC_VERSION}/sources"
                 "${PROJECT_BINARY_DIR}/bief/bief_${TELEMAC_VERSION}/sources" 
                 "${PROJECT_BINARY_DIR}/damocles/damo_${TELEMAC_VERSION}/sources"
                 "${PROJECT_BINARY_DIR}/artemis/arte_${TELEMAC_VERSION}/sources")

add_executable(artemis${TELEMAC_VERSION}.exe homere_artemis.f)
target_link_libraries(artemis${TELEMAC_VERSION}.exe artemis${TELEMAC_VERSION}
                                                    bief${TELEMAC_VERSION}
                                                    damo${TELEMAC_VERSION}
                                                    mumpsvoid${TELEMAC_VERSION}
                                                    paravoid${TELEMAC_VERSION}
                                                    special${TELEMAC_VERSION} )

install(TARGETS artemis${TELEMAC_VERSION}.exe
            RUNTIME DESTINATION artemis/arte_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION artemis/arte_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION artemis/arte_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME})

