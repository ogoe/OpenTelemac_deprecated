#define the library
add_library(paravoid${TELEMAC_VERSION} STATIC errpvm.f get_mpi_parameters.f  p_dmax.f p_dmin.f  p_dsum.f p_exit.f 
p_imax.f p_imin.f p_init.f p_iread.f p_iread_c.f p_isum.f p_iwrit.f p_iwrit_c.f p_lsum.f p_mail.f 
p_mpi_address.f p_mpi_address2.f p_mpi_address3.f 
p_mpi_alltoall.f p_mpi_alltoallv.f p_mpi_type_commit.f p_mpi_type_extent.f p_mpi_type_free.f 
p_mpi_type_lb.f p_mpi_type_struct.f p_mpi_type_ub.f 
p_read.f p_sync.f p_wait_paraco.f p_writ.f build_global_bound.f 
p_mpi_alltoallv_toma1.f p_mpi_alltoallv_toma2.f p_mpi_alltoallv_toma3.f p_mpi_alltoallv_toma4.f )

set_target_properties(paravoid${TELEMAC_VERSION} PROPERTIES PREFIX "")

install(TARGETS paravoid${TELEMAC_VERSION}
            RUNTIME DESTINATION paravoid/paravoid_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            LIBRARY DESTINATION paravoid/paravoid_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME}
            ARCHIVE DESTINATION paravoid/paravoid_${TELEMAC_VERSION}/${Fortran_COMPILER_NAME})


