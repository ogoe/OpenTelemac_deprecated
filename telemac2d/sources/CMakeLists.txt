#
set(CMAKE_INCLUDE_CURRENT_DIR On)
include_directories("${PROJECT_BINARY_DIR}/bief/sources"
                    "${PROJECT_BINARY_DIR}/special/sources"
                    "${PROJECT_BINARY_DIR}/sisyphe/sources" )
#
#define the library
add_library(telemac2d${TELEMAC_VERSION} STATIC friction_def.f
declarations_telemac2d.f  interface_telemac2d.f  
m_coupling_estel3d.f 
homere_telemac2d.f                     
lecdon_telemac2d.f                     
akepin.f assignstr.f astro.f  bilan.f  
bilant1.f bord.f    caldt.f  
carafr.f cdl.f      cdlproj.f celerite.f chpcon.f clhuvt.f  
clsing.f cltrac.f   coefmat.f coefro.f complim.f condin.f  condin_adj.f  
corfon.f corpor.f   cornor.f  correction_depth_2d.f  corstr.f  
corvis.f cosake.f   cost_function.f cubeequation.f 
debimp.f debsce.f   def_zones.f  
difsou.f disper.f   dragcoeff.f dragfo.f entete.f exlim.f fdnrst.f filter_h.f flot.f 
flucin.f flucint.f fluhyd.f 
fluroe.f flusec_telemac2d.f flusew.f   flusrc.f flutrac.f fluxe.f   
fluxpr_telemac2d.f fricti.f friction.f friction_bubble.f 
friction_lindner.f  friction_calc.f friction_quad.f 
friction_unif.f friction_zones.f friction_user.f 
friction_read.f friction_scan.f friction_init.f friction_choice.f 
gestio.f  gradnod.f gradnodt.f gradz.f 
homere_adj_t2d.f hpropa.f interpol.f  
href.f   incide.f   initstr.f intemp.f isitok.f  
kepscl.f kepsil.f   kepsin.f lagran.f  
lecsip.f lecsng.f   loiden.f loinoy.f  maj.f majtrac.f 
marast.f masbas2d.f maskob.f   matbou.f mesures.f meteo.f  metgra.f 
newstr.f nomvar_telemac2d.f output_telemac2d.f point_telemac2d.f point_adj_t2d.f  
oil_spill.f poros.f poro11.f 
prebor.f preres_telemac2d.f  propag.f propag_adj.f propin_telemac2d.f  
prosou.f q.f  read_fic_curves.f read_fic_sources.f read_sections_telemac2d.f 
reinit.f rescue.f resolu.f  read_fic_frliq.f 
rotne0.f siphon.f sl.f smagor.f smoothing_flux.f smtrac.f 
sorflo.f source_telemac2d.f  spectre.f sta_dis_cur.f steady.f tel4del.f telemac2d.f  testeur.f 
thomps.f tr.f trsce.f  utimp_telemac2d.f 
valro.f  vfcfl.f vistur.f vit.f     
volfin.f vusce.f vvsce.f wrihyd.f zeropsi.f zerophi.f )


set_target_properties(telemac2d${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS telemac2d${TELEMAC_VERSION}
            RUNTIME DESTINATION telemac2d/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION telemac2d/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION telemac2d/${Fortran_COMPILER_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/declarations_telemac2d.mod ${CMAKE_CURRENT_BINARY_DIR}/interface_telemac2d.mod
              ${CMAKE_CURRENT_BINARY_DIR}/friction_def.mod ${CMAKE_CURRENT_BINARY_DIR}/m_coupling_estel3d.mod
        DESTINATION telemac2d/${Fortran_COMPILER_NAME})


link_directories("${PROJECT_BINARY_DIR}/special/sources"
                 "${PROJECT_BINARY_DIR}/paravoid/sources"
                 "${PROJECT_BINARY_DIR}/mumpsvoid/sources"
                 "${PROJECT_BINARY_DIR}/bief/sources" 
                 "${PROJECT_BINARY_DIR}/damocles/sources"
                 "${PROJECT_BINARY_DIR}/sisyphe/sources"
                 "${PROJECT_BINARY_DIR}/telemac2d/sources")

add_executable(telemac2d${TELEMAC_VERSION}.exe homere_telemac2d.f)
target_link_libraries(telemac2d${TELEMAC_VERSION}.exe telemac2d${TELEMAC_VERSION}
                                                    sisyphe${TELEMAC_VERSION}
                                                    bief${TELEMAC_VERSION}
                                                    damo${TELEMAC_VERSION}
                                                    mumpsvoid${TELEMAC_VERSION}
                                                    paravoid${TELEMAC_VERSION}
                                                    special${TELEMAC_VERSION} )

install(TARGETS telemac2d${TELEMAC_VERSION}.exe
            RUNTIME DESTINATION telemac2d/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION telemac2d/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION telemac2d/${Fortran_COMPILER_NAME})

