#
set(CMAKE_INCLUDE_CURRENT_DIR On)
include_directories("${PROJECT_BINARY_DIR}/bief/bief_${TELEMAC_VERSION}/sources")
#
#define the library
add_library(sisyphe${TELEMAC_VERSION} STATIC declarations_sisyphe.f         
interface_sisyphe_bedload.f    
interface_sisyphe_suspension.f 
interface_sisyphe.f            
homere_sisyphe.f               
bedload_main.f             qsform.f
sis_arret.f                bedload_meyer.f            rescue_sisyphe.f
bedload_bailard.f          diricl.f     ride.f
bedload_bijker.f           entete_sisyphe.f
bedload_nerbed_vf.f        gf_user.f                  suspension_bilan.f
bedload_calcdw.f           suspension_computation.f  
bedload_dibwat.f           init_avai.f                suspension_dispersion.f 
bedload_diffin.f           bedload_solidischarge.f    init_compo.f
bedload_cheng.f            layer.f                    suspension_erosion.f   
suspension_erosion_coh.f   suspension_flux_mixte.f 
bedload_direction.f        bedload_solvs_fe.f         init_constant.f
suspension_depot.f
bedload_solvs_vf.f         integ.f                    suspension_conv.f
bedload_effpnt.f           bedload_soulsby.f          suspension_bijker.f
suspension_fredsoe.f
bedload_einst.f            flusec_sisyphe.f           fluxpr_sisyphe.f
bedload_engel.f            bedload_vanrijn.f          suspension_listing.f
bedload_engel_old.f        bedload_seccurrent.f       bilan_sisyphe.f
suspension_main.f
bedload_evol.f             calcuw.f                   lecdon_sisyphe.f   
bedload_formula.f          leclis.f                   suspension_rouse.f
tobw_sisyphe.f
coefro_sisyphe.f           maskab.f                         tob_sisyphe.f
bedload_hiding_factor.f    condim_sisyphe.f           mean_grain_size.f
ride_VR.f
bedload_hunz_meyer.f       condim_susp.f              noerod.f
suspension_evol.f
condis_sisyphe.f           nomvar_sisyphe.f           vitchu_sisyphe.f
dredgesim_interface.f
bedload_interact.f         conlit.f                   point_sisyphe.f  
corstr_sisyphe.f           predes.f                         maxslope.f 
init_transport.f           init_mixte.f               init_zero.f
init_sediment.f
neighbourpoints.f	point_height.f		point_volume.f	
calc_neighb_grad.f	calc_critical_bottom_angle.f	neighbourelements.f
sandslide.f
sisyphe.f                 read_sections_sisyphe.f rescue_sisyphe_notperma.f   tassement.f)


set_target_properties(sisyphe${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS sisyphe${TELEMAC_VERSION}
            RUNTIME DESTINATION sisyphe/sisyphe_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION sisyphe/sisyphe_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION sisyphe/sisyphe_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/declarations_sisyphe.mod ${CMAKE_CURRENT_BINARY_DIR}/interface_sisyphe_bedload.mod
              ${CMAKE_CURRENT_BINARY_DIR}/interface_sisyphe.mod ${CMAKE_CURRENT_BINARY_DIR}/interface_sisyphe_suspension.mod
        DESTINATION sisyphe/sisyphe_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME})


link_directories("${PROJECT_BINARY_DIR}/special/special_${TELEMAC_VERSION}/sources"
                 "${PROJECT_BINARY_DIR}/paravoid/paravoid_${TELEMAC_VERSION}/sources"
                 "${PROJECT_BINARY_DIR}/mumpsvoid/mumpsvoid_${TELEMAC_VERSION}/sources"
                 "${PROJECT_BINARY_DIR}/bief/bief_${TELEMAC_VERSION}/sources" 
                 "${PROJECT_BINARY_DIR}/damocles/damo_${TELEMAC_VERSION}/sources"
                 "${PROJECT_BINARY_DIR}/sisyphe/sispyhe_${TELEMAC_VERSION}/sources")

add_executable(sisyphe${TELEMAC_VERSION}.exe homere_sisyphe.f)
target_link_libraries(sisyphe${TELEMAC_VERSION}.exe sisyphe${TELEMAC_VERSION}
                                                    bief${TELEMAC_VERSION}
                                                    damo${TELEMAC_VERSION}
                                                    mumpsvoid${TELEMAC_VERSION}
                                                    paravoid${TELEMAC_VERSION}
                                                    special${TELEMAC_VERSION} )

install(TARGETS sisyphe${TELEMAC_VERSION}.exe
            RUNTIME DESTINATION sisyphe/sisyphe_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION sisyphe/sisyphe_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION sisyphe/sisyphe_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME})

