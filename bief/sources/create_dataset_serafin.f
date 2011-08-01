C                       *********************************
                        SUBROUTINE CREATE_DATASET_SERAFIN
C                       *********************************
C
     *(NFIC,TITRE,NVAR,NOMVAR,OUTVAR)
C
C***********************************************************************
C BIEF VERSION 6.0                                        05/09/2008  RN
C                            
C***********************************************************************
C
C FUNCTION :
C CREATE A DATA SET FOR A GIVEN FILE FORMAT IN THE FILE WITH THE LOGICAL
C UNIT NFILE. THE TITLE OF THE DATASET IS GIVEN AS A 72 CHARACTER
C STRING.
C THE TABLE NOMVAR CONTAINS ALL POSSIBLE VARIABLES TO OUTPUT (IE THE
C NAME OF ALL VARIABLES IN THE OUTPUT BLOCK). THE LOGICAL OUTVAR
C INDICATES FOR EACH VARIABLES WHETHER IT WILL BE WRITTEN OR NOT TO THE
C DATA FILE.
C
C***********************************************************************
C
      IMPLICIT NONE
!
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER                          , INTENT(IN) :: NFIC
      CHARACTER(LEN=72)                , INTENT(IN) :: TITRE
      INTEGER                          , INTENT(IN) :: NVAR
      CHARACTER(LEN=32),DIMENSION(NVAR), INTENT(IN) :: NOMVAR
      LOGICAL          ,DIMENSION(NVAR), INTENT(IN) :: OUTVAR
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER           :: NSOR ! NUMBER OF VARIABLES TO WRITE
      CHARACTER(LEN=80) :: TITSEL
      DOUBLE PRECISION XBID(2)
      INTEGER IB(10),ISTAT,I,IBID(1)
      CHARACTER*2 CBID
!
!***********************************************************************
!     IF(DEBUG) CALL PROC_BEGIN('CREATE_DATASET_SERAFIN')
!***********************************************************************
!
      REWIND NFIC
!
!   LEC/ECR 1   : NOM DU FICHIER GEOMETRIQUE.
!
      TITSEL = TITRE // 'SERAFIN '
      CALL ECRI2(XBID,IBID,TITSEL,80,'CH',NFIC,'STD',ISTAT)
!
!   LEC/ECR 2   : NOMBRE DE FONCTIONS DE DISCRETISATION 1 ET 2
! NOTA : ON N'UTILISE PAS CETTE FONCTIONNALITE DE SERAFIN. TOUTES LES
! VARIABLES SONT DE LA MEME DISCRETISATION (NODALE). 
!
      IB(1)=0
      IB(2)=0
      DO I=1,NVAR
        IF(OUTVAR(I)) IB(1) = IB(1) + 1
      ENDDO
      CALL ECRI2(XBID,IB,CBID,2,'I ',NFIC,'STD',ISTAT)
      NSOR =  IB(1)  +  IB(2)
!
!   LEC/ECR 3 : NOMS ET UNITES DES VARIABLES
!
      IF(NVAR.GE.1) THEN
       DO I=1,NVAR
         IF(OUTVAR(I)) THEN
          CALL ECRI2(XBID,IBID,NOMVAR(I)(1:32),32,'CH',NFIC,'STD',ISTAT)
         ENDIF
       ENDDO
      ENDIF
!
!***********************************************************************
!     IF(DEBUG) CALL PROC_END('CREATE_DATASET_SERAFIN')
!***********************************************************************
!
      RETURN
      END
