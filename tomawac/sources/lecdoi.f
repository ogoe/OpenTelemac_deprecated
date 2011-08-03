C                       *****************
                        SUBROUTINE LECDOI
C                       *****************
C
     *( UD , VD  , X  , Y  , NPOIN2, NDON , BINDON, NBOR , NPTFR, 
     *  AT , DDC , TV1, TV2, NP   , XRELV, YRELV , UR   , VR   ,
     *  TRA, U1  , V1 , U2 , V2   , INDIC, NPMAX , CHDON, NVAR )
C
C***********************************************************************
C  TOMAWAC VERSION 5.0
C***********************************************************************
C
C   FONCTION : CE SOUS-PROGRAMME PROJETE LA VALEUR DES DONNEES
C              DE COURANT OU DE VENT
C              SUR LE MAILLAGE DE CALCUL ET INTERPOLE
C              A L'INSTANT INITIAL
C        (INSPIRE ENTRE AUTRES DE LA ROUTINE FOND DE TELEMAC 2D)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    UD,VD       !<-- !  DONNEE AUX NOEUDS DU MAILLAGE               !
C !    X,Y         ! -->!  COORDONNEES DU MAILLAGE                     !
C !    NPOIN2      ! -->!  NOMBRE DE POINTS DU MAILLAGE                !
C !    NDON        ! -->!  NUMERO D'UNITE LOGIQUE DU FICHIER DE DONNEES!
C !    BINDON      ! -->!  BINAIRE DU FICHIER DES DONNEES              !
C !    NBOR        ! -->!  NUMEROTATION DES POINTS FRONTIERE           !
C !    NPTFR       ! -->!  NOMBRE DE  POINTS FRONTIERE                 !
C !    AT          ! -->!  TEMPS                                       !
C !    DDC         ! -->!  DATE DU DEBUT DU CALCUL                     !
C !    TV1         !<-->!  TEMPS DU CHAMPS DE DONNEES 1                !
C !    TV2         !<-->!  TEMPS DU CHAMPS DE DONNEES 2                !
C !    NP          !<-->!  NOMBRE DE POINTS RELEVES                    !
C !    XRELV       !<-->!  TABLEAU DES ABSCISSES DES POINTS RELEVES    !
C !    YRELV       !<-->!  TABLEAU DES ORDONNEES DES POINTS RELEVES    !
C !    UR,VR       !<-->!  TABLEAU DES DONNEES RELEVEES                !
C !    U1,V1,U2,V2 !<-->!  DONNEES AUX NOEUDS DU MAILLAGE A TV1 ET TV2 !
C !    INDIC       ! -->!  TYPE DE FORMAT DE LECTURE                   !
C !    NPMAX       ! -->!  NOMBRE DE POINTS RELEVES MAXIMUM            !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : CONDIW
C
C SOUS-PROGRAMME APPELE : COUUTI, VENUTI, FASP
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TOMAWAC ,ONLY : MESH
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER NP,NDON,NPOIN2,NPTFR,INDIC,NCOL,NLIG,BID,I,J
      INTEGER NVAR,ISTAT,IB(10),ID(2)
C
      INTEGER NPMAX,NBOR(NPTFR,2)
C
      DOUBLE PRECISION X(NPOIN2)    , Y(NPOIN2)
      DOUBLE PRECISION UD(NPOIN2)   , VD(NPOIN2)
      DOUBLE PRECISION XRELV(NPMAX) , YRELV(NPMAX)
      DOUBLE PRECISION UR(NPMAX)    , VR(NPMAX), TRA(NPMAX)
      DOUBLE PRECISION U1(NPOIN2)   , V1(NPOIN2)
      DOUBLE PRECISION U2(NPOIN2)   , V2(NPOIN2)
      DOUBLE PRECISION XMAX,XMIN,YMAX,YMIN,DX,DY,AT,TV1,TV2
      DOUBLE PRECISION DDC,DAT1,DAT2,COEF,Z(1),ATT, ATB(1)
C
      CHARACTER*3  BINDON,C
      CHARACTER*7  CHDON
      CHARACTER*72 TITCAS
      CHARACTER*32 TEXTE(10)
C
      REAL, ALLOCATABLE :: W(:)
      ALLOCATE(W(NPMAX))
C
C
C-----------------------------------------------------------------------
C        LECTURE DES POINTS RELEVES SUR UNITE LOGIQUE NDON
C-----------------------------------------------------------------------
C
      IF (INDIC.EQ.1) THEN
C
C      -----------------------------------------------------------------
C      FORMAT WAM DIFFERENCES FINIES + INTERPOLATION AUX POINTS
C                 DU MAILLAGE
C      -----------------------------------------------------------------
C
       REWIND NDON
C
       READ(NDON,10,END=100,ERR=100)
     *      NCOL,NLIG,YMIN,YMAX,XMIN,XMAX,BID,BID
       DX=(XMAX-XMIN)/REAL(NCOL-1)
       DY=(YMAX-YMIN)/REAL(NLIG-1)
       NP=NCOL*NLIG
       IF(LNG.EQ.1) THEN
        WRITE(LU,*) '--------------------------------------------------'
        WRITE(LU,*) 'LECDOI : LECTURE DU FICHIER DE ',CHDON
        WRITE(LU,*) '         NOMBRE DE LIGNES   : ',NLIG
        WRITE(LU,*) '         NOMBRE DE COLONNES :',NCOL
        WRITE(LU,*) '         ABSCISSE OU LONGITUDE MINIMALE : ',XMIN
        WRITE(LU,*) '         ABSCISSE OU LONGITUDE MAXIMALE : ',XMAX
        WRITE(LU,*) '         ORDONNEE OU LATITUDE MINIMALE  : ',YMIN
        WRITE(LU,*) '         ORDONNEE OU LATITUDE MAXIMALE  : ',YMAX
        IF (NP.GT.NPMAX) THEN
         WRITE(LU,*) '*************************************************'
         WRITE(LU,*) ' LA DIMENSION PREVUE PAR DEFAUT POUR LE TABLEAU '
         WRITE(LU,*) ' DE DONNEES :',NPMAX,' EST TROP FAIBLE POUR '
         WRITE(LU,*) ' CONTENIR LA TOTALITE DES DONNEES :',NP
         WRITE(LU,*) '*************************************************'
         CALL PLANTE(0)
        ENDIF
       ELSE
        WRITE(LU,*) '--------------------------------------------------'
        WRITE(LU,*)'LECDOI : READING OF THE ',CHDON,' DATA FILE '
        WRITE(LU,*)'         NUMBER OF LINES   : ',NLIG
        WRITE(LU,*)'         NUMBER OF COLUMNS : ',NCOL
        WRITE(LU,*)'         MINIMAL ABSCISSAE : ',XMIN
        WRITE(LU,*)'         MAXIMAL ABSCISSAE : ',XMAX
        WRITE(LU,*)'         MINIMAL ORDINATES : ',YMIN
        WRITE(LU,*)'         MAXIMAL ORDINATES : ',YMAX
        IF (NP.GT.NPMAX) THEN
         WRITE(LU,*) '*************************************************'
         WRITE(LU,*) ' THE DEFAULT DIMENSION ALLOWED FOR THE ARRAY OF '
         WRITE(LU,*) ' DATA :',NPMAX,' IS TOO LOW TO HOLD'
         WRITE(LU,*) ' ALL THE DATA :',NP
         WRITE(LU,*) '*************************************************'
         CALL PLANTE(0)
        ENDIF
       ENDIF
C      LECTURE DE LA DATE DU PREMIER ENREGISTREMENT DE DONNEES
       READ(NDON,*) DAT1
       CALL TEMP(TV1,DAT1,DDC)
       IF (TV1.GT.AT) THEN
        WRITE(LU,*) '*************************************************'
        IF(LNG.EQ.1) THEN
         WRITE(LU,*) ' LE PREMIER ENREGISTREMENT DU FICHIER DE ',CHDON
         WRITE(LU,*) '   ',DAT1,' EST POSTERIEUR AU TEMPS '
         WRITE(LU,*) '   DU DEBUT DU CALCUL',DDC
        ELSE
         WRITE(LU,*) ' THE FIRST RECORDING OF THE ',CHDON,' FILE '
         WRITE(LU,*) '   ',DAT1,' IS OLDER THAN THE DEGINNING '
         WRITE(LU,*) '   OF THE COMPUTATION',DDC
        ENDIF
        WRITE(LU,*) '*************************************************'
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
90     CONTINUE
       READ(NDON,*,END=100,ERR=100)
       READ(NDON,20,END=100,ERR=100)
     *       (UR(I),I=1,NP)
       READ(NDON,*)
       READ(NDON,20,END=100,ERR=100)
     *       (VR(I),I=1,NP)
       CALL OV( 'X=C     ' , U1 , Y , Z , 0.D0 , NPOIN2)
       CALL OV( 'X=C     ' , V1 , Y , Z , 0.D0 , NPOIN2)       
C
       READ(NDON,*) DAT2
       CALL TEMP(TV2,DAT2,DDC)
       IF (TV2.LT.AT) THEN
         TV1=TV2
         GOTO 90
       ENDIF
       CALL FASP(X,Y,U1,NPOIN2,XRELV,YRELV,UR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,0.D0)
       CALL FASP(X,Y,V1,NPOIN2,XRELV,YRELV,VR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,0.D0)
C
       READ(NDON,*,END=100,ERR=100)
       READ(NDON,20,END=100,ERR=100) (UR(I),I=1,NP)
       READ(NDON,*,END=100,ERR=100)
       READ(NDON,20,END=100,ERR=100)
     *      (VR(I),I=1,NP)
       CALL FASP(X,Y,U2,NPOIN2,XRELV,YRELV,UR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,0.D0)
       CALL FASP(X,Y,V2,NPOIN2,XRELV,YRELV,VR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,0.D0)
C
C
      ELSEIF(INDIC.EQ.3) THEN
C
C      -----------------------------------------------------------------
C      FORMAT TELEMAC,
C      LES VARIABLES 1 ET 2 SONT LES COMPOSANTES X ET Y DU VENT
C      -----------------------------------------------------------------
C
       REWIND NDON
       ID(1)=1
       ID(2)=2
C
C      LECTURE DU TITRE
C
       CALL LIT(X,W,IB,TITCAS,72,'CH',NDON,BINDON,ISTAT)
C
C      LECTURE DU NOMBRE DE VARIABLES ET DE LEURS NOMS
C
       CALL LIT(X,W,IB,C,2,'I ',NDON,BINDON,ISTAT)
       NVAR=IB(1)
       DO 80 I=1,NVAR
          CALL LIT(X,W,IB,TEXTE(I),32,'CH',NDON,BINDON,ISTAT)
80     CONTINUE
C
C      VARIABLES FORMAT ET GEOMETRIE
C
       CALL LIT(X,W,IB,C,10,'I ',NDON,BINDON,ISTAT)
       IF(IB(10).EQ.1) THEN
         CALL LIT(X,W,IB,C, 6,'I ',NDON,BINDON,ISTAT)
       ENDIF
       CALL LIT(X,W,IB,C, 4,'I ',NDON,BINDON,ISTAT)
       NP=IB(2)
       WRITE(LU,*) '--------------------------------------------'
       IF (LNG.EQ.1) THEN
        WRITE(LU,*)'LECDOI : LECTURE DU FICHIER TELEMAC'
        WRITE(LU,*) '         TITRE DU CAS LU : ',TITCAS
        WRITE(LU,*)'         NOMBRE DE POINTS   : ',NP
       ELSE
        WRITE(LU,*)'LECDOI : READING OF TELEMAC DATA FILE '
        WRITE(LU,*) '         FILE TITLE : ',TITCAS
        WRITE(LU,*)'         NUMBER OF POINTS   : ',NP
       ENDIF
       WRITE(LU,*) '--------------------------------------------'
       IF (NP.GT.NPMAX) THEN
        WRITE(LU,*) '**************************************************'
        IF(LNG.EQ.1) THEN
         WRITE(LU,*)
     *             ' LA DIMENSION PREVUE PAR DEFAUT POUR LE TABLEAU DE'
         WRITE(LU,*)
     *             ' DONNEES :',NPMAX,' EST TROP FAIBLE POUR CONTENIR'
         WRITE(LU,*) ' LA TOTALITE DES DONNEES :',NCOL*NLIG
        ELSE
         WRITE(LU,*) ' THE DEFAULT DIMENSION ALLOWED FOR THE ARRAY OF '
         WRITE(LU,*) ' DATA :',NPMAX,' IS TOO LOW TO HOLD'
         WRITE(LU,*) ' ALL THE DATA :',NP
        ENDIF
        WRITE(LU,*) '**************************************************'
        CALL PLANTE(0)
       ENDIF
C      TABLEAU ENTIER IKLE
       READ(NDON)
C      TABLEAU ENTIER IPOBO
       READ(NDON)
C
C      X ET Y
C
       CALL LIT(XRELV,W,IB,C,NP,'R4',NDON,BINDON,ISTAT)
       CALL LIT(YRELV,W,IB,C,NP,'R4',NDON,BINDON,ISTAT)
C
C      PAS DE TEMPS ET VARIABLES
C
       CALL LIT(ATB,W,IB,C,1,'R4',NDON,BINDON,ISTAT)
       IF(CHDON(1:1).EQ.'C') THEN
        TV1=ATB(1)
       ELSE
        ATT=ATB(1)*1.D2
        CALL TEMP(TV1,ATT,DDC)
       ENDIF
       IF (TV1.GT.AT) THEN
        WRITE(LU,*) '*************************************************'
        IF(LNG.EQ.1) THEN
         WRITE(LU,*) ' LE PREMIER ENREGISTREMENT DU FICHIER DE ',CHDON
         WRITE(LU,*) '   ',ATT,' EST POSTERIEUR AU TEMPS '
         WRITE(LU,*) '   DU DEBUT DU CALCUL',DDC
        ELSE
         WRITE(LU,*) ' THE FIRST RECORDING OF THE ',CHDON,' FILE '
         WRITE(LU,*) '   ',ATT,' IS OLDER THAN THE BEGINNING '
         WRITE(LU,*) '   OF THE COMPUTATION',DDC
        ENDIF
        WRITE(LU,*) '*************************************************'
        CALL PLANTE(0)
       ENDIF
C
110    CONTINUE
       DO I =1,NVAR
        IF(I.EQ.ID(1)) THEN
         CALL LIT(UR,W,IB,C,NP,'R4',NDON,BINDON,ISTAT)
        ELSEIF(I.EQ.ID(2)) THEN
         CALL LIT(VR,W,IB,C,NP,'R4',NDON,BINDON,ISTAT)
        ELSE
         READ(NDON)
        ENDIF
       ENDDO
C
       CALL LIT(ATB,W,IB,C,1,'R4',NDON,BINDON,ISTAT)
       IF(CHDON(1:1).EQ.'C') THEN
        TV2=ATB(1)
       ELSE
        ATT=ATB(1)*1.D2
        CALL TEMP(TV2,ATT,DDC)
       ENDIF
       IF (TV2.LT.AT) THEN
        TV1=TV2
        GOTO 110
       ENDIF
C
C      ITERPOLATION SPATIALE DES DONNEES AU TEMPS TV1
       CALL FASP(X,Y,U1,NPOIN2,XRELV,YRELV,UR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,1.D-6)
       CALL FASP(X,Y,V1,NPOIN2,XRELV,YRELV,VR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,1.D-6)
C
       DO I =1,NVAR
        IF(I.EQ.ID(1)) THEN
         CALL LIT(UR,W,IB,C,NP,'R4',NDON,BINDON,ISTAT)
        ELSEIF(I.EQ.ID(2)) THEN
         CALL LIT(VR,W,IB,C,NP,'R4',NDON,BINDON,ISTAT)
        ELSE
         READ(NDON)
        ENDIF
       ENDDO
       WRITE(LU,*) 'T',CHDON,'1:',TV1
       WRITE(LU,*) 'T',CHDON,'2:',TV2
C
C      ITERPOLATION SPATIALE DES DONNEES AU TEMPS TV2
       CALL FASP(X,Y,U2,NPOIN2,XRELV,YRELV,UR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,1.D-6)
       CALL FASP(X,Y,V2,NPOIN2,XRELV,YRELV,VR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,1.D-6)
C
C
      ELSEIF (INDIC.EQ.4) THEN
C       LECTURE D'UN FORMAT DEFINI PAR L'UTILISATEUR
        IF(CHDON(1:1).EQ.'C') THEN
C         LECTURE D'UN CHAMP DE COURANT
              CALL COUUTI
     *    (X,Y,NPOIN2,NDON,BINDON,NBOR,NPTFR,AT,DDC,TV1,TV2,
     *     NP,XRELV,YRELV,UR,VR,U1,V1,U2,V2,NPMAX)
        ELSEIF(CHDON(1:1).EQ.'V' .OR. CHDON(1:1).EQ.'W') THEN
C         LECTURE D'UN CHAMP DE VENT
          CALL VENUTI
     *    (X,Y,NPOIN2,NDON,BINDON,NBOR,NPTFR,AT,DDC,TV1,TV2,
     *     NP,XRELV,YRELV,UR,VR,U1,V1,U2,V2,NPMAX)
        ELSE
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'LE TYPE DE DONNEES A LIRE EST INCONNU'
          ELSE
            WRITE(LU,*) 'UNKNOWN DATA'
          ENDIF
            CALL PLANTE(0)
        ENDIF
C
      ELSE
        WRITE(LU,*) '************************************************'
        IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'LECDOI : INDICATEUR DE FORMAT INCONNU : ',INDIC
        WRITE(LU,*) '         POUR LE FICHIER DE ',CHDON
        ELSE
          WRITE(LU,*)'LECDOI : UNKNOWN INDICATOR OF FORMAT : ',INDIC
          WRITE(LU,*)'         FOR THE ',CHDON,' DATA FILE '
        ENDIF
        WRITE(LU,*) '************************************************'
        CALL PLANTE(0)
      ENDIF
C
C-----------------------------------------------------------------------
C   INTERPOLATION TEMPORELLE DES DONNEES
C-----------------------------------------------------------------------
C
      COEF=(AT-TV1)/(TV2-TV1)
      DO 120 I=1,NPOIN2
         UD(I)=(U2(I)-U1(I))*COEF+U1(I)
         VD(I)=(V2(I)-V1(I))*COEF+V1(I)
120   CONTINUE
C
C-----------------------------------------------------------------------
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
      WRITE(LU,*)'*********************************************'
      IF (LNG.EQ.1) THEN
         WRITE(LU,*)'  ERREUR A LA LECTURE DU FICHIER DE DONNEES  '
         WRITE(LU,*)'      OU FIN DE FICHIER PREMATUREE           '
      ELSE
         WRITE(LU,*)'  ERROR WHILE READING DATA FILE '
         WRITE(LU,*)'    OR UNEXPECTED END OF FILE           '
      ENDIF
      WRITE(LU,*)'*********************************************'
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
