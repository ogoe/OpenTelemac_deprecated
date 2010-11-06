C                       ******************************
                        SUBROUTINE DREDGESIM_INTERFACE
C                       ******************************
C
     *(OPTION)
C
C***********************************************************************
C SISYPHE VERSION 6.0      02/02/2009 J-M HERVOUET (LNHE) 01 30 87 80 18
C                                                                                               
C COPYRIGHT EDF      
C***********************************************************************
C
C  FUNCTION: THIS IS THE INTERFACE TO DREDGESIM, CONTAINING ALL 
C            DEPENDENCIES TO DREDGESIM LIBRARIES 
C
C  FOR REAL INTERFACING WITH DREDGESIM, COMMENTS "!" MUST BE REMOVED
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   OPTION       | -->| 1 : INITIALISATION (CALLED IN SISYPHE)
C |                |    | 2 : CALLED EVERY TIME STEP (FROM 
C |                |    |     BEDLOAD_POSTTREATMENT)
C |                |    | 3 : END  (CALLED IN SISYPHE)
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C PROGRAMME APPELANT : SISYPHE, BEDLOAD_POSTTREATMENT
C PROGRAMMES APPELES : 
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC, ONLY : EXTENS
      USE DECLARATIONS_SISYPHE, ONLY : DREDGESIM,DT,NPOIN,NSICLA,DZF_GF,
     *                                 ZFCL_C,AVAIL,MESH,SIS_FILES
!     USE p_sisyphe_ui, ONLY : init_and_setup_ds,clear_sisydredge
!     USE p_dredgesim_ui, ONLY : stop_dredgesim,clear_dredgesim
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: OPTION
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER(LEN=250) :: DREDGEINP,SEDGEO
      DOUBLE PRECISION,ALLOCATABLE :: AVAI_GF(:,:)
      INTEGER I,J
C
C-----------------------------------------------------------------------
C
      IF(OPTION.EQ.1) THEN
C
C     INITIALISATION
C
        DREDGEINP = ''
        sedgeo = ''
        IF(NCSIZE.GT.1) then
C         INPUT FILE FOR DREDGESIM
          DREDGEINP = trim('SISMAF'//extens(ncsize-1,ipid))
          sedgeo = trim('SISGEO'//extens(ncsize-1,ipid))
        ELSE
          dredgeinp = 'SISMAF'
          sedgeo = 'SISGEO'
        ENDIF
!       CALL init_and_setup_ds(SIS_FILES(SISMAF)%LU,dredgeinp,
!                              SIS_FILES(SISGEO)%LU,
!                              sedgeo, 
!    &                         ncsize,ipid,mesh%ilmax,
!    &                         mesh%ikp%i,mesh%ikm%i,mesh%nachb%i,
!    &                         mesh%indpu%i,mesh%nhp%i,mesh%nhm%i)
C
      ELSEIF(OPTION.EQ.2) THEN
C
C     CALL FROM WITHIN BEDLOAD_POSTTREATMENT
C
C       allocating avai_gf
        ALLOCATE(AVAI_GF(NPOIN,NSICLA))
C       INITIALIZATION OF THE ELEVATION TO ADD
        CALL OS('X=0     ',X=DZF_GF)
!       CALL run_dredgesim(DT)
!       avai_gf = get_sm_node_sediment_fraction()
        DO J = 1, NPOIN
          IF(DZF_GF%R(J).GT.0.D0) THEN
            DO I = 1, NSICLA
              ZFCL_C%ADR(I)%P%R(J) = ZFCL_C%ADR(I)%P%R(J) +
     &                               DZF_GF%R(J)*AVAI_GF(J,I)
            ENDDO
          ELSE
            DO I = 1, NSICLA
              ZFCL_C%ADR(I)%P%R(J) = ZFCL_C%ADR(I)%P%R(J) +
     &                               DZF_GF%R(J)*AVAIL(J,1,I)   
            ENDDO
          ENDIF
        ENDDO
        DEALLOCATE(avai_gf)
C
      ELSEIF(OPTION.EQ.3) THEN
C
C     CLOSING
C
!       CALL stop_dredgesim()
!       CALL clear_dredgesim()
!       CALL clear_sisydredge()
C
      ELSE
C
C     ERROR
C
        IF(LNG.EQ.1) WRITE(LU,*) 'MAUVAISE OPTION POUR DREDGESIM'
        IF(LNG.EQ.2) WRITE(LU,*) 'BAD OPTION FOR DREDGESIM'
        CALL PLANTE(1)
        STOP        
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
