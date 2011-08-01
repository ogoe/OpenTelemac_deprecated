C                       *****************
                        SUBROUTINE LECHAM
C                       *****************
C
     *( ZM , DZHDT, X    , Y     , NPOIN2, NDON , BINDON, NBOR  , NPTFR, 
     *  AT , DDC  , TM1  , TM2   , NP   , XRELV , YRELV , ZR   ,
     *  Z1 , Z2   , INDIM, NPMAX , IDHMA, NVAR  )
C
C***********************************************************************
C  TOMAWAC VERSION 5.0
C***********************************************************************
C
C   FONCTION : CE SOUS-PROGRAMME PROJETE LA VALEUR DES DONNEES
C              DE MAREE
C              SUR LE MAILLAGE DE CALCUL ET INTERPOLE
C              A L'INSTANT INITIAL
C        (INSPIRE ENTRE AUTRES DE LA ROUTINE FOND DE TELEMAC 2D)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    ZM          !<-- !  DONNEE AUX NOEUDS DU MAILLAGE               !
C !    X,Y         ! -->!  COORDONNEES DU MAILLAGE                     !
C !    NPOIN2      ! -->!  NOMBRE DE POINTS DU MAILLAGE                !
C !    NDON        ! -->!  NUMERO D'UNITE LOGIQUE DU FICHIER DE DONNEES!
C !    BINDON      ! -->!  BINAIRE DU FICHIER DES DONNEES              !
C !    NBOR        ! -->!  NUMEROTATION DES POINTS FRONTIERE           !
C !    NPTFR       ! -->!  NOMBRE DE  POINTS FRONTIERE                 !
C !    AT          ! -->!  TEMPS                                       !
C !    DDC         ! -->!  DATE DU DEBUT DU CALCUL                     !
C !    TM1         !<-->!  TEMPS DU CHAMPS DE DONNEES 1                !
C !    TM2         !<-->!  TEMPS DU CHAMPS DE DONNEES 2                !
C !    NP          !<-->!  NOMBRE DE POINTS RELEVES                    !
C !    XRELV       !<-->!  TABLEAU DES ABSCISSES DES POINTS RELEVES    !
C !    YRELV       !<-->!  TABLEAU DES ORDONNEES DES POINTS RELEVES    !
C !    ZR          !<-->!  TABLEAU DES DONNEES RELEVEES                !
C !    Z1,Z2       !<-->!  DONNEES AUX NOEUDS DU MAILLAGE A TM1 ET TM2 !
C !    INDIM       ! -->!  TYPE DE FORMAT DE LECTURE                   !
C !    NPMAX       ! -->!  NOMBRE DE POINTS RELEVES MAXIMUM            !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : CONDIW
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
      INTEGER NP,NDON,NPOIN2,NPTFR,INDIM,NCOL,NLIG,BID,I,J
      INTEGER NVAR,NI,ISTAT,IB(10),IDHMA
C
      INTEGER NPMAX,NBOR(NPTFR,2)
C
      DOUBLE PRECISION X(NPOIN2)    , Y(NPOIN2)
      DOUBLE PRECISION ZM(NPOIN2)   , DZHDT(NPOIN2)
      DOUBLE PRECISION XRELV(NPMAX), YRELV(NPMAX)
      DOUBLE PRECISION ZR(NPMAX)   , Z1(NPMAX)   , Z2(NPMAX)
      DOUBLE PRECISION XMAX,XMIN,YMAX,YMIN,DX,DY,AT,TM1,TM2
      DOUBLE PRECISION DDC,DAT1,DAT2,Z(1),ATT, ATB(1)
      DOUBLE PRECISION COE1, COE2
      REAL, ALLOCATABLE :: W(:)
C
      CHARACTER*3  BINDON,C
      CHARACTER*72 TITCAS
      CHARACTER*32 TEXTE(10)
C
      ALLOCATE(W(MAX(NPMAX,72)))
C
C-----------------------------------------------------------------------
C        LECTURE DES POINTS RELEVES SUR UNITE LOGIQUE NDON
C-----------------------------------------------------------------------
C
      IF (INDIM.EQ.1) THEN
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
        WRITE(LU,*) 'LECHAM : LECTURE DU FICHIER DE LA MAREE '
        WRITE(LU,*) '         NOMBRE DE LIGNES   : ',NLIG
        WRITE(LU,*) '         NOMBRE DE COLONNES : ',NCOL
        WRITE(LU,*) '         ABSCISSE OU LONGITUDE MINIMALE : ',XMIN
        WRITE(LU,*) '         ABSCISSE OU LONGITUDE MAXIMALE : ',XMAX
        WRITE(LU,*) '         ORDONNEE OU LATITUDE MINIMALE  : ',YMIN
        WRITE(LU,*) '         ORDONNEE OU LATITUDE MAXIMALE  : ',YMAX
        IF(NP.GT.NPMAX) THEN
         WRITE(LU,*) '*************************************************'
         WRITE(LU,*) ' LA DIMENSION PREVUE PAR DEFAUT POUR LE TABLEAU '
         WRITE(LU,*) ' DE DONNEES :',NPMAX,' EST TROP FAIBLE POUR '
         WRITE(LU,*) ' CONTENIR LA TOTALITE DES DONNEES :',NP
         WRITE(LU,*) '*************************************************'
         CALL PLANTE(1)
        ENDIF
       ELSE
        WRITE(LU,*) '--------------------------------------------------'
        WRITE(LU,*)'LECHAM : READING OF THE TIDE DATA FILE '
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
       CALL TEMP(TM1,DAT1,DDC)
       IF (TM1.GT.AT) THEN
        WRITE(LU,*) '*************************************************'
        IF(LNG.EQ.1) THEN
         WRITE(LU,*) 'LE PREMIER ENREGISTREMENT DU FICHIER DE LA MAREE'
         WRITE(LU,*) '   ',DAT1,' EST POSTERIEUR AU TEMPS '
         WRITE(LU,*) '   DU DEBUT DU CALCUL',DDC
        ELSE
         WRITE(LU,*) ' THE FIRST RECORDING OF THE TIDE DATA FILE '
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
       READ(NDON,20,END=100,ERR=100) (ZR(I),I=1,NP)
       CALL OV( 'X=C     ' , Z1 , Y , Z , 0.D0 , NPOIN2)      
C
       READ(NDON,*) DAT2
       CALL TEMP(TM2,DAT2,DDC)
       IF (TM2.LE.AT) THEN
         TM1=TM2
         GOTO 90
       ENDIF
       CALL FASP(X,Y,Z1,NPOIN2,XRELV,YRELV,ZR,NP,NBOR,MESH%KP1BOR%I,
     *                                                   NPTFR,0.D0)
C
       READ(NDON,*,END=100,ERR=100)
       READ(NDON,20,END=100,ERR=100) (ZR(I),I=1,NP)
       CALL FASP(X,Y,Z2,NPOIN2,XRELV,YRELV,ZR,NP,NBOR,MESH%KP1BOR%I,
     *                                                   NPTFR,0.D0)
C
C
      ELSEIF (INDIM.EQ.2) THEN
C
C      -----------------------------------------------------------------
C      FORMAT TELEMAC
C      -----------------------------------------------------------------
C
       REWIND NDON
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
       CALL LIT(X,W,IB,C, 4,'I ',NDON,BINDON,ISTAT)
       NP=IB(2)
       NI=IB(1)*IB(3)
       WRITE(LU,*) '--------------------------------------------'
       IF (LNG.EQ.1) THEN
        WRITE(LU,*)'LECHAM : LECTURE DU FICHIER TELEMAC'
        WRITE(LU,*) '         TITRE DU CAS LU : ',TITCAS
        WRITE(LU,*)'         NOMBRE DE POINTS   : ',NP
       ELSE
        WRITE(LU,*)'LECHAM : READING OF TELEMAC DATA FILE '
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
       TM1=ATB(1)
CFBG       ATT=ATB(1)*1.D2
CFBG       CALL TEMP(TM1,ATT,DDC)
       IF (TM1.GT.AT) THEN
        WRITE(LU,*) '*************************************************'
        IF(LNG.EQ.1) THEN
         WRITE(LU,*) ' LE PREMIER ENREGISTREMENT DU FICHIER '
         WRITE(LU,*) '   CONTENANT LE NIVEAU DE LA MAREE    '
!         WRITE(LU,*) '   ',ATT,' EST POSTERIEUR AU TEMPS '
!         WRITE(LU,*) '   DU DEBUT DU CALCUL',DDC
         WRITE(LU,*) '   ',AT,' EST POSTERIEUR AU TEMPS '
         WRITE(LU,*) '   DU DEBUT DU CALCUL',TM1
        ELSE
         WRITE(LU,*) ' THE FIRST RECORDING OF THE TIDE LEVEL FILE '
         WRITE(LU,*) '   ',ATT,' DATES LATER THAN THE DEGINNING '
         WRITE(LU,*) '   OF THE COMPUTATION',DDC
        ENDIF
        WRITE(LU,*) '*************************************************'
        CALL PLANTE(0)
       ENDIF
C
110    CONTINUE
       DO I =1,NVAR
        IF(I.EQ.IDHMA) THEN
         CALL LIT(Z1,W,IB,C,NP,'R4',NDON,BINDON,ISTAT)
        ELSE
         READ(NDON)
        ENDIF
       ENDDO
C

       CALL LIT(ATB,W,IB,C,1,'R4',NDON,BINDON,ISTAT)
       TM2=ATB(1)
CFBG       ATT=ATB(1)*1.D2
CFBG       CALL TEMP(TM2,ATT,DDC)
       IF (TM2.LE.AT) THEN
        TM1=TM2
        GOTO 110
       ENDIF
        CALL FASP(X,Y,ZR,NPOIN2,XRELV,YRELV,Z1,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,1.D-6)
!       write(*,*)'Z1 RELEVE PT 7176',Z1(7176)

        CALL OV( 'X=Y     ' , Z1 , ZR , Z , 0.D0 , NPOIN2)
       DO I =1,NVAR
        IF(I.EQ.IDHMA) THEN
          CALL LIT(Z2,W,IB,C,NP,'R4',NDON,BINDON,ISTAT)
        ELSE
         READ(NDON)
        ENDIF
       ENDDO

       WRITE(LU,*) 'TMAREE1:',TM1
       WRITE(LU,*) 'TMAREE2:',TM2

C
C       ITERPOLATION SPATIALE DES DONNEES
!        CALL FASP(X,Y,ZR,NPOIN2,XRELV,YRELV,Z1,NP,NBOR,MESH%KP1BOR%I,
!     *                                                    NPTFR,0.D0)
!       IF (X(2080)==577499.125d0) write(*,*)'Z1 RELEVE PT 2080',ZR(2080)
!       IF (X(1307)==577499.125d0) write(*,*)'Z1 RELEVE PT 1307',ZR(1307)
!        CALL OV( 'X=Y     ' , Z1 , ZR , Z , 0.D0 , NPOIN2)
        CALL FASP(X,Y,ZR,NPOIN2,XRELV,YRELV,Z2,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,1.D-6)
!        write(*,*)'Z2 RELEVE PT 7176',Z2(7176)
!        IF (X(2080)==577499.125d0) write(*,*)'Z2 RELEVE PT 2080',ZR(2080)
!        IF (X(1307)==577499.125d0) write(*,*)'Z2 RELEVE PT 1307',ZR(1307)
        CALL OV( 'X=Y     ' , Z2 , ZR , Z , 0.D0 , NPOIN2)
C
C
      ELSEIF (INDIM.EQ.3) THEN
C       LECTURE D'UN FORMAT DEFINI PAR L'UTILISATEUR
              CALL MARUTI
     *    (X,Y,NPOIN2,NDON,BINDON,NBOR,NPTFR,AT,DDC,TM1,TM2,
     *     NP,XRELV,YRELV,ZR,Z1,Z2,NPMAX)
C
      ELSE
        WRITE(LU,*) '************************************************'
        IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'LECHAM : INDICATEUR DE FORMAT INCONNU : ',INDIM
        WRITE(LU,*) '         POUR LE FICHIER DU NIVEAU DE LA MAREE '
        ELSE
          WRITE(LU,*)'LECHAM : UNKNOWN INDICATOR OF FORMAT : ',INDIM
          WRITE(LU,*)'         FOR THE TIDE LEVEL DATA FILE '
        ENDIF
        WRITE(LU,*) '************************************************'
        CALL PLANTE(0)
      ENDIF
C
C-----------------------------------------------------------------------
C   INTERPOLATION TEMPORELLE DES DONNEES
C   ET GRADIENT TEMPOREL DE LA MAREE
C-----------------------------------------------------------------------
C
      COE1=(TM2-TM1)
      IF (COE1.LT.1.D-4) THEN
         WRITE(LU,*) '****************************************'
         IF(LNG.EQ.1) THEN
           WRITE(LU,*) ' DEUX TEMPS IDENTIQUES                '
           WRITE(LU,*) ' DANS LE FICHIER DES HAUTEURS DE MAREE'
         ELSE
           WRITE(LU,*) ' TWO IDENTICAL TIMES IN THE TIDAL FILE'
         ENDIF
         WRITE(LU,*) '****************************************'
         CALL PLANTE(0)
      ENDIF
      COE2=(AT-TM1)/COE1
      DO 120 I=1,NPOIN2
         ATT     = (Z2(I)-Z1(I))
         ZM(I)   = ATT*COE2+Z1(I)
         DZHDT(I)= ATT/COE1
120   CONTINUE
C
C
C-----------------------------------------------------------------------
C
      DEALLOCATE(W)
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
      WRITE(LU,*)'**********************************************'
      IF (LNG.EQ.1) THEN
         WRITE(LU,*)'  ERREUR A LA LECTURE DU FICHIER DE LA MAREE '
         WRITE(LU,*)'      OU FIN DE FICHIER PREMATUREE           '
      ELSE
         WRITE(LU,*)'  ERROR WHILE READING TIDAL DATA '
         WRITE(LU,*)'    OR UNEXPECTED END OF FILE           '
      ENDIF
      WRITE(LU,*)'**********************************************'
      CALL PLANTE(0)
C
      RETURN
      END
