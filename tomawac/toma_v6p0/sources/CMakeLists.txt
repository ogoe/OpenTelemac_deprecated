#
set(CMAKE_INCLUDE_CURRENT_DIR On)
include_directories("${PROJECT_BINARY_DIR}/bief/bief_${TELEMAC_VERSION}/sources")
#
#define the library
add_library(tomawac${TELEMAC_VERSION} STATIC declarations_tomawac.f tomawac_mpi.f tomawac_mpi_tools.f
anacos.f anamar.f anaven.f angles.f bornes.f 
condiw.f conw4d.f conwac.f corfon.f 
cormar.f couuti.f cqueue.f 
delfra.f dump2d.f 
ecrete.f ecrspe.f 
foncro.f fpread.f frem01.f frem02.f fremoy.f frepic.f fsprd1.f 
fsprd2.f fsprd3.f 
gammln.f 
homere_tomawac.f 
impr.f inip4d.f iniphy.f inipie.f  initab.f 
iniven.f intang.f inter4d.f interp_tomawac.f 
kerbou.f kmoyen.f 
lecdoi.f lecdon.f lecdon_tomawac.f lecham.f leclim_tomawac.f lecsui.f 
lecuti.f limwac.f  preleo_mpi.f 
maruti.f nomvar_tomawac.f noudon.f noumar.f 
p_read2.f p_writ2.f pied4d.f pieds.f point_tomawac.f preleo.f prenl1.f 
prepro.f preqt2.f propa.f 
qbbj78.f qbrek1.f qbrek2.f qbrek3.f qbrek4.f qfrot1.f qgauss.f qmout1.f qnlin1.f 
qtria1.f qtria2.f qwind1.f qwind2.f 
radiat.f 
semimp.f sor3d.f speini.f spejon.f spetma.f spread.f stress.f 
tautot.f temp.f tetmoy.f totnrj.f transf.f 
ustar1.f ustar2.f 
vartel.f venuti.f vitfon.f 
wac.f wnscou.f wpower.f )


set_target_properties(tomawac${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS tomawac${TELEMAC_VERSION}
            RUNTIME DESTINATION tomawac/toma_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION tomawac/toma_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION tomawac/toma_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/declarations_tomawac.mod
              ${CMAKE_CURRENT_BINARY_DIR}/tomawac_mpi.mod ${CMAKE_CURRENT_BINARY_DIR}/tomawac_mpi_tools.mod
        DESTINATION tomawac/toma_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME})


link_directories("${PROJECT_BINARY_DIR}/special/special_${TELEMAC_VERSION}/sources"
                 "${PROJECT_BINARY_DIR}/paravoid/paravoid_${TELEMAC_VERSION}/sources"
                 "${PROJECT_BINARY_DIR}/mumpsvoid/mumpsvoid_${TELEMAC_VERSION}/sources"
                 "${PROJECT_BINARY_DIR}/bief/bief_${TELEMAC_VERSION}/sources" 
                 "${PROJECT_BINARY_DIR}/damocles/damo_${TELEMAC_VERSION}/sources"
                 "${PROJECT_BINARY_DIR}/tomawac/toma_${TELEMAC_VERSION}/sources")

add_executable(tomawac${TELEMAC_VERSION}.exe homere_tomawac.f)
target_link_libraries(tomawac${TELEMAC_VERSION}.exe tomawac${TELEMAC_VERSION}
                                                    bief${TELEMAC_VERSION}
                                                    damo${TELEMAC_VERSION}
                                                    mumpsvoid${TELEMAC_VERSION}
                                                    paravoid${TELEMAC_VERSION}
                                                    special${TELEMAC_VERSION} )

install(TARGETS tomawac${TELEMAC_VERSION}.exe
            RUNTIME DESTINATION tomawac/toma_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION tomawac/toma_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION tomawac/toma_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME})

