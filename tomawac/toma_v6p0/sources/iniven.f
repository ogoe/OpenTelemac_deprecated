C                       *****************
                        SUBROUTINE INIVEN
C                       *****************
C
     *(UV,VV,X,Y,NPOIN,NVEN, BINVEN,NBOR,NPTFR,AT,DDC,TV1,TV2,
     * NP,XRELV,YRELV,U1,V1,U2,V2,INDIC,NPMAX,ITR01)
C
C***********************************************************************
C  TOMAWAC VERSION 1.0    01/02/95        F.MARCOS     (LNH) 30 87 72 66
C***********************************************************************
C
C   FONCTION : CE SOUS-PROGRAMME PROJETE LA VALEUR DES VENTS
C              SUR LE MAILLAGE DE CALCUL ET INTERPOLE
C              A L'INSTANT INITIAL
C        (INSPIRE ENTRE AUTRES DE LA ROUTINE FOND DE TELEMAC 2D)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    UV,VV       !<-- !  VENT AUX NOEUDS DU MAILLAGE                 !
C !    X,Y         ! -->!  COORDONNEES DU MAILLAGE                     !
C !    NPOIN       ! -->!  NOMBRE DE POINTS DU MAILLAGE                !
C !    NVEN-F      ! -->!  NUMERO D'UNITE LOGIQUE DU FICHIER DES VENTS !
C !    BINVEN      ! -->!  BINAIRE DU FICHIER DES VENTS (SI INDIC>=2)  !
C !    NBOR        ! -->!  NUMEROTATION DES POINTS FRONTIERE           !
C !    NPTFR       ! -->!  NOMBRE DE  POINTS FRONTIERE                 !
C !    AT          ! -->!  TEMPS                                       !
C !    DDC         ! -->!  DATE DU DEBUT DU CALCUL                     !
C !    TV1         !<-->!  TEMPS DU CHAMPS DE VENT 1                   !
C !    TV2         !<-->!  TEMPS DU CHAMPS DE VENT 2                   !
C !    NP          !<-->!  NOMBRE DE POINTS RELEVES                    !
C !    XRELV       !<-->!  TABLEAU DES ABSCISSES DES POINTS RELEVES    !
C !    YRELV       !<-->!  TABLEAU DES ORDONNEES DES POINTS RELEVES    !
C !    UR,VR       !<-->!  TABLEAU DES COURANTS RELEVES                !
C !    U1,V1,U2,V2 !<-->!  VENT AUX NOEUDS DU MAILLAGE DU VENT         !
C !    INDIC       ! -->!  TYPE DE FORMAT DE LECTURE                   !
C !    NPMAX       ! -->!  NOMBRE DE POINTS RELEVES MAXIMUM            !
C !    ITR01       ! -->!  TABLEAU DE TRAVAIL ENTIER                   !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PREWAC
C
C SOUS-PROGRAMME APPELE : FASP
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TOMAWAC ,ONLY : MESH
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
      INTEGER NP,NVEN,NPOIN,NPTFR,INDIC,NCOL,NLIG,BID,I,J
      INTEGER NVAR,NI,ISTAT,IB(10),ITR01(*)
C
      INTEGER NPMAX,NBOR(NPTFR,2)
C
      DOUBLE PRECISION X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION UV(NPOIN),VV(NPOIN)
      DOUBLE PRECISION XRELV(NPMAX),YRELV(NPMAX)
      DOUBLE PRECISION U1(NPMAX),V1(NPMAX),U2(NPMAX),V2(NPMAX)
      DOUBLE PRECISION XMAX,XMIN,YMAX,YMIN,DX,DY,AT,TV1,TV2
      DOUBLE PRECISION DDC,DAT1,DAT2,COEF,Z(1),ATT, ATB(1)
C
      CHARACTER*3 BINVEN,C
      CHARACTER*72 TITCAS
      CHARACTER*32 TEXTE(10)
C
      DOUBLE PRECISION, ALLOCATABLE :: UR(:),VR(:)
      REAL, ALLOCATABLE :: W(:)
      ALLOCATE(W(MAX(NPMAX,72)))
C
C-----------------------------------------------------------------------
C        LECTURE DES POINTS RELEVES SUR UNITE LOGIQUE NVEN
C-----------------------------------------------------------------------
C
      IF(INDIC.EQ.1) THEN
C
      REWIND NVEN
C
C     ------------------------------------------------------------------
C     FORMAT WAM DIFFERENCES FINIES + INTERPOLATION AUX POINTS
C                DU MAILLAGE
C     ------------------------------------------------------------------
C
       READ(NVEN,10,END=100,ERR=100)
     * NCOL,NLIG,YMIN,YMAX,XMIN,XMAX,BID,BID
       DX=(XMAX-XMIN)/REAL(NCOL-1)
       DY=(YMAX-YMIN)/REAL(NLIG-1)
       NP=NCOL*NLIG
       ALLOCATE(UR(1:NPMAX),VR(1:NPMAX))
       WRITE(LU,*) '---------------------------------------------------'
       WRITE(LU,*) 'INIVEN : LECTURE DU FICHIER DE VENT'
       WRITE(LU,*) '         NOMBRE DE LIGNES   : ',NLIG
       WRITE(LU,*) '         NOMBRE DE COLONNES :',NCOL
       WRITE(LU,*) '         ABSCISSE OU LONGITUDE MINIMALE : ',XMIN
       WRITE(LU,*) '         ABSCISSE OU LONGITUDE MAXIMALE : ',XMAX
       WRITE(LU,*) '         ORDONNEE OU LATITUDE MINIMALE  : ',YMIN
       WRITE(LU,*) '         ORDONNEE OU LATITUDE MAXIMALE  : ',YMAX
       IF (NP.GT.NPMAX) THEN
        WRITE(LU,*) '**************************************************'
        WRITE(LU,*) ' LA DIMENSION PREVUE PAR DEFAUT POUR LE TABLEAU   '
        WRITE(LU,*) ' DE VENT :',NPMAX,' EST TROP FAIBLE POUR CONTENIR'
        WRITE(LU,*) ' CONTENIR LA TOTALITE DES DONNEES :',NCOL*NLIG
        WRITE(LU,*) '**************************************************'
        CALL PLANTE(0)
       ENDIF
       READ(NVEN,*) DAT1
       CALL TEMP(TV1,DAT1,DDC)
       IF (TV1.GT.AT) THEN
        WRITE(LU,*) '******************************************'
        WRITE(LU,*) ' LE PREMIER ENREGISTREMENT DU FICHIER DES'
        WRITE(LU,*) ' VENTS : ',DAT1,' EST POSTERIEUR AU TEMPS'
        WRITE(LU,*) ' DU DEBUT DU CALCUL',DDC
        WRITE(LU,*) '******************************************'
        CALL PLANTE(0)
       ENDIF
C
       DO 50 I=1,NCOL
          DO 40 J=1,NLIG
                XRELV((I-1)*NLIG+J)=XMIN+DX*(I-1)
                YRELV((I-1)*NLIG+J)=YMIN+DY*(J-1)
40        CONTINUE
50     CONTINUE
C
90    CONTINUE
      READ(NVEN,*,END=100,ERR=100)
      READ(NVEN,20,END=100,ERR=100)
     *       (UR(I),I=1,NP)
          READ(NVEN,*)
         READ(NVEN,20,END=100,ERR=100)
     *       (VR(I),I=1,NP)
          CALL OV( 'X=C     ' , U1 , Y , Z , 0.D0 , NPMAX)
          CALL OV( 'X=C     ' , V1 , Y , Z , 0.D0 , NPMAX)
          CALL FASP(X,Y,U1,NPOIN,XRELV,YRELV,UR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,0.D0)
          CALL FASP(X,Y,V1,NPOIN,XRELV,YRELV,VR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,0.D0)
C
           READ(NVEN,*) DAT2
           CALL TEMP(TV2,DAT2,DDC)
           IF (TV2.LT.AT) THEN
              TV1=TV2
              GOTO 90
           ENDIF
C
          READ(NVEN,*,END=100,ERR=100)
          READ(NVEN,20,END=100,ERR=100)
     *      (UR(I),I=1,NP)
        READ(NVEN,*,END=100,ERR=100)
        READ(NVEN,20,END=100,ERR=100)
     *      (VR(I),I=1,NP)
          CALL FASP(X,Y,U2,NPOIN,XRELV,YRELV,UR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,0.D0)
          CALL FASP(X,Y,V2,NPOIN,XRELV,YRELV,VR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,0.D0)
C
C
      ELSEIF (INDIC.EQ.2) THEN
C
       REWIND NVEN
C
C     ------------------------------------------------------------------
C     FORMAT TELEMAC, DISCRETISATION SUR LE MEME MAILLAGE
C         LES VARIABLES 1 ET 2 SONT LES COMPOSANTES X ET Y DU VENT
C     ------------------------------------------------------------------
C
C      LECTURE DU TITRE
C
       CALL LIT(X,W,IB,TITCAS,72,'CH',NVEN,BINVEN,ISTAT)
C
C      LECTURE DU NOMBRE DE VARIABLES ET DE LEURS NOMS
C
       CALL LIT(X,W,IB,C,2,'I ',NVEN,BINVEN,ISTAT)
       NVAR=IB(1)
       DO 80 I=1,NVAR
          CALL LIT(X,W,IB,TEXTE(I),32,'CH',NVEN,BINVEN,ISTAT)
80     CONTINUE
C
C      VARIABLES FORMAT ET GEOMETRIE
C
       CALL LIT(X,W,IB,C,10,'I ',NVEN,BINVEN,ISTAT)
       CALL LIT(X,W,IB,C, 4,'I ',NVEN,BINVEN,ISTAT)
       NP=IB(2)
       NI=IB(1)*IB(3)
       WRITE(LU,*) '--------------------------------------------'
       WRITE(LU,*) 'INIVEN : LECTURE DU FICHIER TELEMAC'
       WRITE(LU,*) '         TITRE DU CAS LU : ',TITCAS
       WRITE(LU,*) '         NOMBRE DE POINTS   : ',NP
       WRITE(LU,*) '--------------------------------------------'
       IF (NP.NE.NPOIN) THEN
       WRITE(LU,*) '***************************************************'
       WRITE(LU,*)'VOUS UTILISEZ UN FORMAT DE LECTURE TELEMAC SUPPOSANT'
       WRITE(LU,*)' L''EGALITE DES MAILLAGES LUS ET UTILISES. CECI EST '
       WRITE(LU,*)' IMPOSSIBLE CAR LE NOMBRE DE POINTS LUS :',NP,' EST'
       WRITE(LU,*)' DIFFERENT DU NOMBRE DE POINTS DU MAILLAGE :',NPOIN
       WRITE(LU,*) '***************************************************'
       CALL PLANTE(1)
       ENDIF
       ALLOCATE(UR(1:NPMAX),VR(1:NPMAX))
       CALL LIT(X,W,ITR01,C,NI,'I ',NVEN,BINVEN,ISTAT)
       CALL LIT(X,W,ITR01,C,NP,'I ',NVEN,BINVEN,ISTAT)
C
C      X ET Y
C
       CALL LIT(XRELV,W,IB,C,NP,'R4',NVEN,BINVEN,ISTAT)
       CALL LIT(YRELV,W,IB,C,NP,'R4',NVEN,BINVEN,ISTAT)
C
C      PAS DE TEMPS ET VARIABLES
C
       CALL LIT(ATB,W,IB,C,1,'R4',NVEN,BINVEN,ISTAT)
       ATT=ATB(1)*1.D2
       CALL TEMP(TV1,ATT,DDC)
       IF (TV1.GT.AT) THEN
          WRITE(LU,*) 'ERREUR DEMARAGE LECTURE',TV1,AT
        CALL PLANTE(0)
          ENDIF
110       CONTINUE
        CALL LIT(U1,W,IB,C,NP,'R4',NVEN,BINVEN,ISTAT)
        CALL LIT(V1,W,IB,C,NP,'R4',NVEN,BINVEN,ISTAT)
C
        CALL LIT(ATB,W,IB,C,1,'R4',NVEN,BINVEN,ISTAT)
        ATT=ATB(1)*1.D2
        CALL TEMP(TV2,ATT,DDC)
        IF (TV2.LT.AT) THEN
           TV1=TV2
           GOTO 110
          ENDIF
          CALL LIT(U2,W,IB,C,NP,'R4',NVEN,BINVEN,ISTAT)
          CALL LIT(V2,W,IB,C,NP,'R4',NVEN,BINVEN,ISTAT)
        WRITE(LU,*) 'TVENT1:',TV1
        WRITE(LU,*) 'TVENT2:',TV2
C
C
      ELSEIF (INDIC.EQ.3) THEN
C
        ALLOCATE(UR(1:NPMAX),VR(1:NPMAX))
        CALL VENUTI
     * (X,Y,NPOIN,NVEN,BINVEN,NBOR,NPTFR,AT,DDC,TV1,TV2,
     *  NP,XRELV,YRELV,UR,VR,U1,V1,U2,V2,NPMAX)
C
      ELSE
C
        WRITE(LU,*) '***********************************************'
        WRITE(LU,*) 'INIVEN : INDICATEUR DE FORMAT INCONNU   '
        WRITE(LU,*) '         POUR LE FICHIER DES VENTS :',INDIC
        WRITE(LU,*) '***********************************************'
        CALL PLANTE(1)
      ENDIF
C
C-----------------------------------------------------------------------
C   INTERPOLATION
C-----------------------------------------------------------------------
C
        COEF=(AT-TV1)/(TV2-TV1)
        DO 120 I=1,NPOIN
           UV(I)=(U2(I)-U1(I))*COEF+U1(I)
           VV(I)=(V2(I)-V1(I))*COEF+V1(I)
120     CONTINUE
C
C-----------------------------------------------------------------------
C
      DEALLOCATE(UR,VR)
C
C     FORMATS
C
10    FORMAT (2I4,4F9.3,2I2)
20    FORMAT (10F6.2)
C
      RETURN
C
C     EN CAS DE PROBLEME DE LECTURE ...
C
100   CONTINUE
      WRITE(LU,*) '*********************************************'
      WRITE(LU,*) '  ERREUR A LA LECTURE DU FICHIER DE VENT     '
      WRITE(LU,*) '      OU FIN DE FICHIER PREMATUREE           '
      WRITE(LU,*) '*********************************************'
      CALL PLANTE(1)
C
C-----------------------------------------------------------------------
C
      DEALLOCATE(W)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
