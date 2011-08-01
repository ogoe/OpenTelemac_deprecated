C                       *****************************
                        SUBROUTINE WRITE_MESH_SERAFIN
C                       *****************************
C
     *(NFIC,MESH,NPLAN,DATE,TIME,I_ORIG,J_ORIG,FFORMAT)
C
C***********************************************************************
C BIEF VERSION 6.0           01/04/2009          R NEBAUER (LNHE) 
C***********************************************************************
C
C FONCTION : WRITES HEADER DATA IN A SERAFIN FORMAT FILE
C            THEN ONLY THE RECORDS OF RESULTS WILL BE ADDED
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   NFIC         | -->| LOGICAL UNIT OF FILE
C |   MESH         | -->| MESH STRUCTURE
C |   NPLAN        | -->| NUMBER OF PLANES IN 3D, 1 IN 2D
C |   DATE, TIME   | -->| DATE AND TIME : YEAR MONTH DAY HOUR MIN SEC
C |   I_ORIG,J_ORIG| -->| COORDINATES OF ORIGIN OF MESH (IN METRES) 
C |   FFORMAT      | -->| FILE FORMAT
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER        ,       INTENT(IN) :: NFIC   ! LU FILE
      TYPE(BIEF_MESH),       INTENT(IN) :: MESH   ! MESH STRUCTURE
      INTEGER        ,       INTENT(IN) :: NPLAN  ! NUMBER OF PRISM LAYERS
      INTEGER, DIMENSION(3), INTENT(IN) :: DATE   ! DATE INFO
      INTEGER, DIMENSION(3), INTENT(IN) :: TIME   ! TIME INFO
      INTEGER        ,       INTENT(IN) :: I_ORIG ! X ORIGIN OF MESH
      INTEGER        ,       INTENT(IN) :: J_ORIG ! Y ORIGIN OF MESH
      CHARACTER(LEN=8) ,     INTENT(IN) :: FFORMAT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER              :: NDP      ! NODES PER ELEMENT
      INTEGER              :: IELEM    ! ELEMENT TYPE
      DOUBLE PRECISION     :: XBID(2)
      INTEGER              :: IBID(1)
      CHARACTER*2          :: CBID
      INTEGER              :: IB(10)
      INTEGER              :: ISTAT    ! RETURN VALUE
      INTEGER              :: I        ! LOOP COUNTER
      INTEGER              :: ERR      ! ERROR CODE
      INTEGER, ALLOCATABLE :: IPOBO(:) ! BORDER NODE ID'S
      INTEGER, ALLOCATABLE :: IKLES(:) ! GLOBAL CONNECTIVITY
      LOGICAL              :: YA_IPOBO
      LOGICAL              :: YA_IKLES
      CHARACTER*2 RF
      INTEGER, DIMENSION(:), POINTER :: IKLE
      INTEGER :: NELEM
!
!***********************************************************************
!     IF(DEBUG) CALL PROC_BEGIN('WRITE_MESH_SERAFIN')
!***********************************************************************
!
      IF(FFORMAT.EQ.'SERAFIND') THEN
        RF = 'R8'
      ELSE
        RF = 'R4'
      END IF
!
      SELECT CASE (MESH%TYPELM )
        CASE (10) ! TRIANGLES
           NDP = 3
        CASE (20) ! QUADRANGLES
           NDP = 4
        CASE (40) ! PRISMES
           NDP = 6
        CASE (50) ! PRISMES DECOUPES EN TETRAS
           NDP = 4
        CASE DEFAULT
          IF(LNG.EQ.1) THEN
            WRITE(LU,*)'TYPE ELEMENT ERRONE POUR SERAFIN :',MESH%TYPELM
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*)'BAD ELEMENT TYPE FOR SERAFIN :',MESH%TYPELM
          ENDIF
          CALL PLANTE(1)
          STOP
      END SELECT
!
      YA_IPOBO = .FALSE.
      YA_IKLES = .FALSE.
!
!   LEC/ECR 4   : LISTE DE 10 PARAMETRES ENTIERS
!
        IB(1) = 1
        DO I = 2,10
          IB(I) = 0
        ENDDO  
!
!       ORIGIN COORDINATES IN METRES
!
        IB(3)=I_ORIG
        IB(4)=J_ORIG
C
C       NUMBER OF PLANES IN 3D
C
C       FOLLOWING LINE WILL GIVE NPLAN=0 IN 2D
C       THIS IS WHAT IS EXPECTED IN GRETEL
C       TO REMOVE IF(NPLAN.GT.1) CHANGE GRETEL ACCORDINGLY.
        IF(NPLAN.GT.1) IB(7)=NPLAN
!
!PARA   MARQUAGE POUR ANNONCER LA LECTURE DE KNOLG
        IF(NCSIZE.GT.1) THEN
          IB(8)=MESH%NPTFR
          IB(9)=NPTIR
        ENDIF
!PARAFIN
        NELEM = MESH%NELEM
!
!   Y-A-T-IL PASSAGE DE LA DATE ?
        IF(DATE(1)+DATE(2)+DATE(3)+TIME(1)+TIME(2)+TIME(3).NE.0) THEN
         IB(10) = 1
        ENDIF
!   ECRITURE DU TABLEAU DE 10 PARAMETRES
        CALL ECRI2(XBID,IB,CBID,10,'I ',NFIC,'STD',ISTAT)
!   PASSAGE DE LA DATE
        IF(IB(10).EQ.1) THEN
          IB(1)=DATE(1)
          IB(2)=DATE(2)
          IB(3)=DATE(3)
          IB(4)=TIME(1)
          IB(5)=TIME(2)
          IB(6)=TIME(3)
          CALL ECRI2(XBID,IB,CBID,6,'I ',NFIC,'STD',ISTAT)
        ENDIF
!
!   LEC/ECR 5 : 4 ENTIERS
!
      IF(MESH%TYPELM.NE.50) THEN
        IB(1) = NELEM
      ELSE
!       TETRAHEDRONS REGROUPED INTO PRISMS
        IB(1)=NELEM/3
      ENDIF
      IB(2) = MESH%NPOIN
!
      IF(MESH%TYPELM.NE.50) THEN
        IB(3) = NDP
      ELSE
!       TETRAHEDRONS REGROUPED INTO PRISMS
        IB(3) = 6
      ENDIF
      IB(4) = 1
      CALL ECRI2(XBID,IB,CBID,4,'I ',NFIC,'STD',ISTAT)
!
!   LEC/ECR 6 : IKLE
!
      IF(MESH%TYPELM.NE.50) THEN
        ALLOCATE(IKLES(NELEM*NDP),STAT=ERR)
      ELSE
!       TETRAHEDRONS REGROUPED INTO PRISMS
        ALLOCATE(IKLES(NELEM*2)  ,STAT=ERR)
      ENDIF
      IF(ERR.NE.0) STOP 'ECRGEO : ALLOCATION DE IKLES'
      YA_IKLES = .TRUE.
!     INVERSION DE IKLE EN IKLES POUR SELAFIN
      IKLE => MESH%IKLE%I
      IF(MESH%TYPELM.NE.50) THEN
        DO I = 1,NDP
          DO IELEM = 1,NELEM
            IKLES((IELEM-1)*NDP+I) = IKLE((I-1)*NELEM+IELEM)
          ENDDO
        ENDDO
      ELSE
!     TETRAHEDRONS REGROUPED INTO PRISMS
        DO IELEM  = 1,NELEM/3
          IKLES((IELEM-1)*6+1) = IKLE(      IELEM)
          IKLES((IELEM-1)*6+2) = IKLE(NELEM+IELEM)
          IKLES((IELEM-1)*6+3) = IKLE(NELEM+IELEM)
          IKLES((IELEM-1)*6+4) = IKLE(      IELEM)+MESH%NPOIN/NPLAN
          IKLES((IELEM-1)*6+5) = IKLE(NELEM+IELEM)+MESH%NPOIN/NPLAN
          IKLES((IELEM-1)*6+6) = IKLE(NELEM+IELEM)+MESH%NPOIN/NPLAN
        ENDDO
      ENDIF
!
      IF(MESH%TYPELM.NE.50) THEN
      CALL ECRI2(XBID,IKLES,CBID,NELEM*NDP,'I ',NFIC,'STD',ISTAT)
      ELSE
!     TETRAHEDRONS REGROUPED INTO PRISMS
      CALL ECRI2(XBID,IKLES,CBID,NELEM*2,'I ',NFIC,'STD',ISTAT)
      ENDIF
!
!   LEC/ECR 7 : IPOBO (CAS DES FICHIERS SANS PARALLELISME)
!
      IF(IB(8).EQ.0.AND.IB(9).EQ.0) THEN
!
        ALLOCATE(IPOBO(MESH%NPOIN),STAT=ERR)
        IF(ERR.NE.0) STOP 'ECRGEO : ALLOCATION DE IPOBO'
        YA_IPOBO = .TRUE.
        DO I=1,MESH%NPOIN
         IPOBO(I) = 0
        END DO  
!       ONLY LATERAL BOUNDARY POINTS WITH PRISMS
        DO I =1,MESH%NPTFR
         IPOBO(MESH%NBOR%I(I)) = I
        END DO
        CALL ECRI2(XBID,IPOBO,CBID,MESH%NPOIN,'I ',NFIC,'STD',ISTAT)
!
      ENDIF
!
      IF(IB(8).NE.0.OR.IB(9).NE.0) THEN
!
!   LEC/ECR  7.1 KNOLG (SEULEMENT EN CAS DE PARALLELISME)
!
      CALL ECRI2(XBID,MESH%KNOLG%I,CBID,MESH%NPOIN,'I ',NFIC,'STD',
     *           ISTAT)
!
      ENDIF
!
!   LEC/ECR  8 ET 9 : X ET Y  COORDONNEES DES POINTS DU MAILLAGE

      CALL ECRI2(MESH%X%R ,IBID,CBID,MESH%NPOIN,RF,NFIC,'STD',ISTAT)
      CALL ECRI2(MESH%Y%R ,IBID,CBID,MESH%NPOIN,RF,NFIC,'STD',ISTAT)
!
      IF(YA_IPOBO) DEALLOCATE(IPOBO)
      IF(YA_IKLES) DEALLOCATE(IKLES)

!***********************************************************************
!     IF(DEBUG) CALL PROC_END('WRITE_MESH_SERAFIN')
!***********************************************************************

      RETURN
      END
