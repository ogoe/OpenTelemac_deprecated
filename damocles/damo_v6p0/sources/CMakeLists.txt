#define the library
add_library(damo${TELEMAC_VERSION} STATIC aidelu.f carlu.f classe.f cmd.f damoc.f damocle.f damocles.f
dico.f    influ.f intlu.f  loglu.f longlu.f majus.f
minus.f   next.f  precar.f prev.f  preval.f realu.f )

set_target_properties(damo${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS damo${TELEMAC_VERSION}
            RUNTIME DESTINATION damocles/damo_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION damocles/damo_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION damocles/damo_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME})


