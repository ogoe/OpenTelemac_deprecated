C                       *************************
                        SUBROUTINE FLUSEC_SISYPHE
C                       *************************
C
     *(U,V,H,QSXC,QSYC,CHARR,QSXS,QSYS,SUSP,
     * IKLE,NELMAX,NELEM,X,Y,DT,NCP,CTRLSC,INFO,TPS,KNOGL)
C
C***********************************************************************
C  SISYPHE VERSION 5.7 27/12/06       J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C  FONCTION : CALCUL DES FLUX A TRAVERS DES SECTIONS DE CONTROLE
C             ET CUMUL DE CES FLUX POUR OBTENIR LES VOLUMES OSCILLANTS.
C
C             MAILLAGES DE DIMENSION 2 ET HAUTEUR D'EAU CONSIDEREE
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   U,V          | -->| CHAMP DE VITESSE
C |   H            | -->| HAUTEUR D'EAU
C |   CS           | -->| BLOC DE TRACEURS EN SUSPENSION
C |   SUSP         | -->| LOGIQUE INDIQUANT DE PRENDRE EN COMPTE
C |                |    | LE BLOC DE TRACEURS EN SUSPENSION
C |   IKLE         | -->| TABLEAUX DE CONNECTIVITE LOCAL-GLOBAL
C |   XEL,YEL      | -->| COORDONNEES DES POINTS PAR ELEMENT
C |   NELMAX       | -->| NOMBRE MAXIMUM D'ELEMENTS.
C |   NELEM        | -->| NOMBRE D'ELEMENTS.
C |   X,Y          | -->| COORDONNEES DES POINTS DU MAILLAGE
C |   DT           | -->| PAS DE TEMPS.
C |   NCP          | -->| TWO TIMES THE NUMBER OF CONTROL SECTIONS
C |   CTRLSC       | -->| DONNEES SUR LES SECTIONS DE CONTROLE.
C |   INFO         | -->| SI OUI : IMPRESSIONS.
C |   TPS          | -->| TEMPS
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR :
C
C SOUS-PROGRAMMES APPELES :
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_SISYPHE, ONLY: CHAIN
      USE INTERFACE_SISYPHE, EX_FLUSEC_SISYPHE => FLUSEC_SISYPHE
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)          :: NELMAX,NELEM,NCP
      INTEGER, INTENT(IN)          :: IKLE(NELMAX,*) 
      INTEGER, INTENT(IN)          :: CTRLSC(NCP),KNOGL(*) 
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*),TPS,DT
       LOGICAL, INTENT(IN)          :: INFO,SUSP,CHARR
      TYPE(BIEF_OBJ), INTENT(IN)   :: U,V,H,QSXC,QSYC,QSXS,QSYS
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NSEMAX,ERR
      PARAMETER(NSEMAX=50)
C
      INTEGER IELEM,I1,I2,I3,ELBEST,IGBEST,ILBEST
      INTEGER ILPREC,ISEG,ISEC,NSEC,PT,DEP,ARR
C
      DOUBLE PRECISION DIST,DIST1,DIST2,DIST3
      DOUBLE PRECISION H1,H2,X1,Y1,X2,Y2,UN1,UN2,NX,NY,SUR6
C
      DOUBLE PRECISION, ALLOCATABLE :: FLX(:),VOLNEG(:),VOLPOS(:)
      DOUBLE PRECISION, ALLOCATABLE :: FLXS(:),FLXC(:)
      DOUBLE PRECISION, ALLOCATABLE :: VOLNEGS(:),VOLPOSS(:)
      DOUBLE PRECISION, ALLOCATABLE :: VOLNEGC(:),VOLPOSC(:)
      INTEGER, ALLOCATABLE :: NSEG(:),LISTE(:,:,:)
C
      LOGICAL DEJA
      DATA DEJA/.FALSE./
      LOGICAL :: CLASSIC=.FALSE. 
C
      SAVE LISTE,DEJA,NSEG,VOLNEG,VOLPOS,FLX,FLXS,VOLNEGS,VOLPOSS
      SAVE FLXC,VOLNEGC,VOLPOSC,CLASSIC
C
C-----------------------------------------------------------------------
C
      WRITE(lu,*) '-> entering FLUSEC_SISYPHE'
      WRITE(lu,*) 'NCP: ',NCP
      WRITE(lu,*) 'CTRLSC: ',CTRLSC(:)

      SUR6 = 1.D0/6.D0
      NSEC = NCP/2
C
C  RECHERCHE DES CHEMINS QUI JOIGNENT LES COUPLES DE POINTS :
C
      IF(.NOT.DEJA) THEN
C
C     ALLOCATION DYNAMIQUE DE FLX, VOLNEG, VOLPOS, ETC.
C
      ALLOCATE(FLX(NSEC)           ,STAT=ERR)
      ALLOCATE(VOLNEG(NSEC)        ,STAT=ERR)
      ALLOCATE(VOLPOS(NSEC)        ,STAT=ERR)
      ALLOCATE(NSEG(NCP)           ,STAT=ERR)
      ALLOCATE(LISTE(NCP,NSEMAX,2) ,STAT=ERR)
C     S FOR SUSPENSION, C FOR BEDLOAD (CHARRIAGE...)
      ALLOCATE(FLXS(NSEC)          ,STAT=ERR)
      ALLOCATE(VOLNEGS(NSEC)       ,STAT=ERR)
      ALLOCATE(VOLPOSS(NSEC)       ,STAT=ERR)
      ALLOCATE(FLXC(NSEC)          ,STAT=ERR)
      ALLOCATE(VOLNEGC(NSEC)       ,STAT=ERR)
      ALLOCATE(VOLPOSC(NSEC)       ,STAT=ERR)
C
      IF(ERR.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,100) ERR
        IF(LNG.EQ.2) WRITE(LU,200) ERR
100     FORMAT(1X,'FLUSEC : ERREUR A L''ALLOCATION DE MEMOIRE : ',/,1X,
     *            'CODE D''ERREUR : ',1I6)
200     FORMAT(1X,'FLUSEC: ERROR DURING ALLOCATION OF MEMORY: ',/,1X,
     *            'ERROR CODE: ',1I6)
      ENDIF
C
      IF (.NOT.ALLOCATED(CHAIN)) CLASSIC=.TRUE.
C
      DO ISEC=1,NSEC
        FLX(ISEC)=0.0D0
        VOLNEG(ISEC) =0.D0
        VOLPOS(ISEC) =0.D0
        VOLNEGS(ISEC)=0.D0
        VOLPOSS(ISEC)=0.D0
        VOLNEGC(ISEC)=0.D0
        VOLPOSC(ISEC)=0.D0
      ENDDO
C
      DO 60 ISEC =1,NSEC
C
C!jaj #### in the serial case, or "classical" in parallel, 
C     follow the algorithm of finding segment chains
C 
C NOTE: if you change the algorithm, change it in PARTEL as well
C
        IF (NCSIZE.LE.1 .OR. CLASSIC) THEN
C
        DEP = CTRLSC(1+2*(ISEC-1))
        ARR = CTRLSC(2+2*(ISEC-1))
        IF(NCSIZE.GT.1) THEN
          DEP=KNOGL(DEP)
          ARR=KNOGL(ARR)
          IF(DEP.EQ.0.AND.ARR.EQ.0) THEN
            NSEG(ISEC)=0
            GO TO 60
          ENDIF
          IF((DEP.EQ.0.AND.ARR.NE.0).OR.(DEP.NE.0.AND.ARR.EQ.0)) THEN
            NSEG(ISEC)=-1
            GO TO 60
          ENDIF
        ENDIF
        PT = DEP
        ISEG = 0
        DIST=(X(DEP)-X(ARR))**2+(Y(DEP)-Y(ARR))**2
10      CONTINUE
C
        DO 20 IELEM =1,NELEM
C
          I1 = IKLE(IELEM,1)
          I2 = IKLE(IELEM,2)
          I3 = IKLE(IELEM,3)
C         SI L'ELEMENT CONTIENT LE POINT COURANT :
          IF(PT.EQ.I1.OR.PT.EQ.I2.OR.PT.EQ.I3) THEN
            DIST1 = (X(I1)-X(ARR))**2 + (Y(I1)-Y(ARR))**2
            DIST2 = (X(I2)-X(ARR))**2 + (Y(I2)-Y(ARR))**2
            DIST3 = (X(I3)-X(ARR))**2 + (Y(I3)-Y(ARR))**2
            IF(DIST1.LT.DIST) THEN
              DIST = DIST1
              ELBEST = IELEM
              IGBEST = I1
              ILBEST = 1
              IF(I1.EQ.PT) ILPREC = 1
              IF(I2.EQ.PT) ILPREC = 2
              IF(I3.EQ.PT) ILPREC = 3
            ENDIF
            IF(DIST2.LT.DIST) THEN
              DIST = DIST2
              ELBEST = IELEM
              IGBEST = I2
              ILBEST = 2
              IF(I1.EQ.PT) ILPREC = 1
              IF(I2.EQ.PT) ILPREC = 2
              IF(I3.EQ.PT) ILPREC = 3
            ENDIF
            IF(DIST3.LT.DIST) THEN
              DIST = DIST3
              ELBEST = IELEM
              IGBEST = I3
              ILBEST = 3
              IF(I1.EQ.PT) ILPREC = 1
              IF(I2.EQ.PT) ILPREC = 2
              IF(I3.EQ.PT) ILPREC = 3
            ENDIF
          ENDIF
C
20      CONTINUE
C
        IF(IGBEST.EQ.PT) THEN
          IF(LNG.EQ.1) WRITE(LU,32)
          IF(LNG.EQ.2) WRITE(LU,33)
32        FORMAT(1X,'FLUSEC : BLOCAGE DE L''ALGORITHME')
33        FORMAT(1X,'FLUSEC : ALGORITHM FAILED')
          CALL PLANTE(1)
          STOP
        ELSE
          PT = IGBEST
          ISEG = ISEG + 1
          IF(ISEG.GT.NSEMAX) THEN
            IF(LNG.EQ.1) THEN
              WRITE(LU,*) 'FLUSEC : TROP DE SEGMENTS DANS UNE'
              WRITE(LU,*) '         SECTION. AUGMENTER NSEMAX'
            ENDIF
            IF(LNG.EQ.2) THEN
              WRITE(LU,*) 'FLUSEC: TOO MANY SEGMENTS IN A   '
              WRITE(LU,*) '        SECTION. INCREASE  NSEMAX'
            ENDIF
            STOP
          ENDIF
          LISTE(ISEC,ISEG,1) = IKLE(ELBEST,ILPREC)
          LISTE(ISEC,ISEG,2) = IKLE(ELBEST,ILBEST)
          IF(IGBEST.NE.ARR) GO TO 10
        ENDIF
C
        NSEG(ISEC) = ISEG
!
!jaj #### this part to be done in the parallel case; fill LISTE 
! with ready segment chains provided by PARTEL: see read_sections
! note: future optimisation - use CHAIN structure in the whole routine 
!
        ELSE
C
          NSEG(ISEC) = CHAIN(ISEC)%NSEG
          LISTE(ISEC,:,:)=0

          DO ISEG=1,NSEG(ISEC)
            LISTE(ISEC,ISEG,1) = CHAIN(ISEC)%LISTE(ISEG,1)
            LISTE(ISEC,ISEG,2) = CHAIN(ISEC)%LISTE(ISEG,2)
          END DO 

          WRITE(lu,*) 'chain@sisyphe -> liste@sisyphe:'
          WRITE(lu,*) 'isec,nseg(isec): ',isec,nseg(isec)
          DO ISEG=1,NSEG(ISEC)
            WRITE(lu,*) LISTE(ISEC,ISEG,:) 
          END DO 

        ENDIF 
C
60    CONTINUE
C
C     IF(.NOT.DEJA) THEN
      ENDIF
C
C-----------------------------------------------------------------------
C
      DEJA = .TRUE.
C
C-----------------------------------------------------------------------
C
      DO ISEC = 1 , NSEC
C
      FLX(ISEC)  = 0.D0
      FLXS(ISEC) = 0.D0
      FLXC(ISEC) = 0.D0
C
      IF(NSEG(ISEC).GE.1) THEN
C
C     COMPUTING THE FLUX DIRECTLY, REGARDLESS OF THE WEAK FORM
C     OF THE IMPERMEABILITY CONDITION
C
      DO ISEG = 1 , NSEG(ISEC)          
        I1 = LISTE(ISEC,ISEG,1)
        I2 = LISTE(ISEC,ISEG,2)
        X1 = X(I1)
        X2 = X(I2)
        Y1 = Y(I1)
        Y2 = Y(I2)
        H1 = H%R(I1)
        H2 = H%R(I2)
        NX = Y1-Y2
        NY = X2-X1
        UN1= U%R(I1)*NX + V%R(I1)*NY
        UN2= U%R(I2)*NX + V%R(I2)*NY
        FLX(ISEC) = FLX(ISEC) + ((H1+H2)*(UN1+UN2)+H2*UN2+H1*UN1)*SUR6
        IF(SUSP) THEN
          UN1= QSXS%R(I1)*NX + QSYS%R(I1)*NY
          UN2= QSXS%R(I2)*NX + QSYS%R(I2)*NY
          FLXS(ISEC) = FLXS(ISEC) + 0.5D0*(UN1+UN2)
        ENDIF
        IF(CHARR) THEN
          UN1= QSXC%R(I1)*NX + QSYC%R(I1)*NY
          UN2= QSXC%R(I2)*NX + QSYC%R(I2)*NY
          FLXC(ISEC) = FLXC(ISEC) + 0.5D0*(UN1+UN2)
        ENDIF
      ENDDO
C
      IF(FLX(ISEC).GT.0.D0) THEN
        VOLPOS(ISEC) = VOLPOS(ISEC) + FLX(ISEC)*DT
      ELSE
        VOLNEG(ISEC) = VOLNEG(ISEC) + FLX(ISEC)*DT
      ENDIF
C
      IF(SUSP) THEN
        IF(FLXS(ISEC).GT.0.D0) THEN
          VOLPOSS(ISEC) = VOLPOSS(ISEC) + FLXS(ISEC)*DT
        ELSE
          VOLNEGS(ISEC) = VOLNEGS(ISEC) + FLXS(ISEC)*DT
        ENDIF
      ENDIF
C
      IF(CHARR) THEN
        IF(FLXC(ISEC).GT.0.D0) THEN
          VOLPOSC(ISEC) = VOLPOSC(ISEC) + FLXC(ISEC)*DT
        ELSE
          VOLNEGC(ISEC) = VOLNEGC(ISEC) + FLXC(ISEC)*DT
        ENDIF
      ENDIF
C
C     IF(NSEG(ISEC).GT.1)...
      ENDIF
C
      ENDDO
C
C-----------------------------------------------------------------------
C
C     PRINTING THE RESULTS / !jaj here allreduces for values
C
      CALL FLUXPR_SISYPHE(NSEC,CTRLSC,FLX,VOLNEG,VOLPOS,INFO,TPS,
     *                    NSEG,NCSIZE,
     *                    FLXS,VOLNEGS,VOLPOSS,SUSP,
     *                    FLXC,VOLNEGC,VOLPOSC,CHARR)
C
C-----------------------------------------------------------------------
C
      WRITE(lu,*) '-> leaving FLUSEC_SISYPHE'
      RETURN
      END
