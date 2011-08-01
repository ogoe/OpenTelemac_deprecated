C                       *****************
                        SUBROUTINE NOUMAR
C                       *****************
C
     *(ZM , DZHDT, X  , Y  , NPOIN2, NDON , BINDON, NBOR , NPTFR,
     * AT , DDC  , TM1, TM2, NP    , XRELV, YRELV , ZR   , TRA  ,
     * Z1 , Z2   , INDIM, IDHMA , NVHMA )
C
C***********************************************************************
C  TOMAWAC VERSION 5.0
C***********************************************************************
C
C   FONCTION : CE SOUS-PROGRAMME CALCULE LA VALEUR DE LA MAREE
C              A L'INSTANT COURANT
C              SUR LE MAILLAGE DE CALCUL
C             (INSPIRE DE LA ROUTINE FOND DE TELEMAC 2D)
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
C !    BINDON      ! -->!  BINAIRE DU FICHIER DE DONNEES               !
C !    NBOR        ! -->!  NUMEROTATION DES POINTS FRONTIERE           !
C !    NPTFR       ! -->!  NOMBRE DE  POINTS FRONTIERE                 !
C !    AT          ! -->!  TEMPS                                       !
C !    DDC         ! -->!  DATE DU DEBUT DU CALCUL                     !
C !    TM1         !<-->!  TEMPS DU CHAMPS DE DONNEES 1                !
C !    TM2         !<-->!  TEMPS DU CHAMPS DE DONNEES 2                !
C !    NP          !<-->!  NOMBRE DE POINTS DU MAILLAGE DES DONNEES    !
C !    XRELV       !<-- !  TABLEAU DES ABSCISSES DES POINTS RELEVES    !
C !    YRELV       !<-- !  TABLEAU DES ORDONNEES DES POINTS RELEVES    !
C !    ZR          !<-->!  TABLEAU DES COURANTS RELEVES                !
C !    Z1,Z2       !<-->!  DONNEES AUX NOEUDS DU MAILLAGE              !
C !    INDIM       ! -->!  TYPE DE FORMAT DE LECTURE                   !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : SEMIMP
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
      INTEGER NP,NDON,NPOIN2,NPTFR,INDIM,I,ISTAT,IW(1)
C
      INTEGER NBOR(NPTFR,2), IDHMA,NVHMA
C
      DOUBLE PRECISION X(NPOIN2) ,Y(NPOIN2)
      DOUBLE PRECISION ZM(NPOIN2),ZR(NP), DZHDT(NPOIN2)
      DOUBLE PRECISION Z1(NPOIN2),Z2(NPOIN2)
      DOUBLE PRECISION XRELV(NP),YRELV(NP),TRA(NP)
      DOUBLE PRECISION AT,TM1,TM2
      DOUBLE PRECISION DDC,DAT2,DAT2B(1),Z(1),C,COE1,COE2,ATT
C
      CHARACTER*3 BINDON, C1
C
C-----------------------------------------------------------------------
C
      REAL, ALLOCATABLE :: W(:)
      ALLOCATE(W(NP))
C
C-----------------------------------------------------------------------
C
      IF (AT.GE.TM2) THEN
C
C       ----------------------------------------------------------------
C        ON CHANGE D'ENREGISTREMENT : 2->1 ET ON LIT UN NOUVEAU 2
C       ----------------------------------------------------------------
        TM1=TM2
        CALL OV('X=Y     ', Z1 , Z2 , Z , C , NPOIN2)
C
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) '   NOUMAR : LECTURE D''UN NOUVEL ENREGISTREMENT'
          WRITE(LU,*) '            DE LA HAUTEUR DE LA MAREE          '
        ELSE
          WRITE(LU,*) '   NOUMAR : READING A NEW RECORDING '
          WRITE(LU,*) '            OF THE TIDE LEVEL       '
        ENDIF
C
        IF (INDIM.EQ.1) THEN
C
C     ------------------------------------------------------------------
C          FICHIER DIFFERENCES FINIES FORMATTE DU TYPE WAM CYCLE 4
C     ------------------------------------------------------------------
 90        CONTINUE
C          LECTURE : DATE DE L'ENREGISTREMENT
           READ(NDON,*,END=100,ERR=100) DAT2
           CALL TEMP(TM2,DAT2,DDC)
C          LECTURE : DONNEES
           READ(NDON,*,END=100,ERR=100)
           READ(NDON,20,END=100,ERR=100) (ZR(I),I=1,NP)
C
           IF (TM2.LE.AT) THEN
             IF(LNG.EQ.1) THEN
               WRITE(LU,*) ' NOUMAR : ON SAUTE 1 ENREGISTREMENT ..'
             ELSE
               WRITE(LU,*) ' NOUMAR : JUMP OF 1 RECORDED DATA SERIES'
             ENDIF
             TM1=TM2
             CALL FASP(X,Y,Z1,NPOIN2,XRELV,YRELV,ZR,NP,NBOR,
     *                                         MESH%KP1BOR%I,NPTFR,0.D0)
             GOTO 90
           ENDIF
C
           CALL FASP(X,Y,Z2,NPOIN2,XRELV,YRELV,ZR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,0.D0)
C
        ELSEIF (INDIM.EQ.2) THEN
C
C     ------------------------------------------------------------------
C       FICHIER SELAFIN DU TYPE TELEMAC
C     ------------------------------------------------------------------
C
 95     CONTINUE
C       LECTURE : DATE DE L'ENREGISTREMENT
        CALL LIT(DAT2B,W,IW,C1,1,'R4',NDON,BINDON,ISTAT)
        TM2=DAT2B(1)
C       LECTURE : DONNEES
        DO I =1,NVHMA
          IF(I.EQ.IDHMA) THEN
            CALL LIT(ZR,W,IW,C1,NP,'R4',NDON,BINDON,ISTAT)
          ELSE
            READ(NDON)
          ENDIF
        ENDDO
C
        IF (TM2.LE.AT) THEN
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) ' NOUMAR : ON SAUTE 1 ENREGISTREMENT ..'
          ELSE
            WRITE(LU,*) ' NOUMAR : JUMP OF 1 RECORDED DATA SERIES'
          ENDIF
          TM1=TM2
C         INTERPOLATION SPATIALE DES DONNEES AU TEMPS 1
          CALL FASP(X,Y,Z1,NPOIN2,XRELV,YRELV,ZR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,1.D-6)
!     *                                                    NPTFR,0.D0)
          GOTO 95
        ENDIF
C
        WRITE(LU,*) 'TMENT1=',TM1
        WRITE(LU,*) 'TMENT2=',TM2
C       INTERPOLATION SPATIALE DES DONNEES AU TEMPS 2
        CALL FASP(X,Y,Z2,NPOIN2,XRELV,YRELV,ZR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,1.D-6)
!     *                                                    NPTFR,0.D0)
C
        ELSEIF (INDIM.EQ.3) THEN
C
C     ------------------------------------------------------------------
C        LECTURE D'UN FORMAT DEFINI PAR L'UTILISATEUR
C     ------------------------------------------------------------------
C
          CALL MARUTI
     *    (X,Y,NPOIN2,NDON,BINDON,NBOR,NPTFR,AT,DDC,TM1,TM2,
     *     NP,XRELV,YRELV,ZR,Z1,Z2,NP)
C
C
        ELSE
C
        WRITE(LU,*) '************************************************'
        IF(LNG.EQ.1) THEN
         WRITE(LU,*) 'NOUMAR : INDICATEUR DE FORMAT INCONNU : ',INDIM
        ELSE
          WRITE(LU,*)'NOUMAR : UNKNOWN INDICATOR OF FORMAT : ',INDIM
        ENDIF
        WRITE(LU,*) '************************************************'
        CALL PLANTE(1)
        ENDIF
C
      ENDIF
C
C     -------------------------------------------------
C       INTERPOLATION TEMPORELLE DES DONNEES
C        ET GRADIENT TEMPOREL DE LA MAREE
C     -------------------------------------------------
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
      DO I=1,NPOIN2
         ATT     = (Z2(I)-Z1(I))
         ZM(I)   = ATT*COE2+Z1(I)
         DZHDT(I)= ATT/COE1
      ENDDO
C
C-----------------------------------------------------------------------
C
C     FORMATS
C
20    FORMAT (10F6.2)
C
      DEALLOCATE(W)
      RETURN
C
C     EN CAS DE PROBLEME DE LECTURE ...
C
100   CONTINUE
      WRITE(LU,*)'*********************************************'
      IF (LNG.EQ.1) THEN
         WRITE(LU,*)'  ERREUR A LA LECTURE DU FICHIER DE MAREE '
         WRITE(LU,*)'      OU FIN DE FICHIER PREMATUREE        '
      ELSE
         WRITE(LU,*)'  ERROR WHILE READING DATA FILE '
         WRITE(LU,*)'    OR UNEXPECTED END OF FILE           '
      ENDIF
      WRITE(LU,*)'*********************************************'
      CALL PLANTE(0)
C
      RETURN
      END
