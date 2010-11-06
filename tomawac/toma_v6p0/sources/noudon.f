C                       *****************
                        SUBROUTINE NOUDON
C                       *****************
C
     *(UV , VV , X  , Y  , NPOIN, NDON , BINDON, NBOR, NPTFR,
     * AT , DDC, TV1, TV2, NP   , XRELV, YRELV , UR  , VR   ,
     * TRA, U1 , V1 , U2 , V2   , INDIC, CHDON , NVAR)
C
C***********************************************************************
C  TOMAWAC VERSION 5.0
C***********************************************************************
C
C   FONCTION : CE SOUS-PROGRAMME CALCULE LA VALEUR DU VENT OU
C              DU COURANT A L'INSTANT COURANT
C              SUR LE MAILLAGE DE CALCUL
C             (INSPIRE DE LA ROUTINE FOND DE TELEMAC 2D)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    UV,VV       !<-- !  DONNEE AUX NOEUDS DU MAILLAGE               !
C !    X,Y         ! -->!  COORDONNEES DU MAILLAGE                     !
C !    NPOIN       ! -->!  NOMBRE DE POINTS DU MAILLAGE                !
C !    NDON        ! -->!  NUMERO D'UNITE LOGIQUE DU FICHIER DE DONNEES!
C !    BINDON      ! -->!  BINAIRE DU FICHIER DE DONNEES               !
C !    NBOR        ! -->!  NUMEROTATION DES POINTS FRONTIERE           !
C !    NPTFR       ! -->!  NOMBRE DE  POINTS FRONTIERE                 !
C !    AT          ! -->!  TEMPS                                       !
C !    DDC         ! -->!  DATE DU DEBUT DU CALCUL                     !
C !    TV1         !<-->!  TEMPS DU CHAMPS DE DONNEES 1                !
C !    TV2         !<-->!  TEMPS DU CHAMPS DE DONNEES 2                !
C !    NP          !<-->!  NOMBRE DE POINTS DU MAILLAGE DES DONNEES    !
C !    XRELV       !<-- !  TABLEAU DES ABSCISSES DES POINTS RELEVES    !
C !    YRELV       !<-- !  TABLEAU DES ORDONNEES DES POINTS RELEVES    !
C !    UR,VR       !<-->!  TABLEAU DES COURANTS RELEVES                !
C !    U1,V1,U2,V2 !<-->!  DONNEES AUX NOEUDS DU MAILLAGE              !
C !    INDIC       ! -->!  TYPE DE FORMAT DE LECTURE                   !
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
      INTEGER NP,NDON,NPOIN,NPTFR,INDIC,I,ISTAT,NVAR,IW(1)
C
      INTEGER NBOR(NPTFR,2),ID(2)
C
      DOUBLE PRECISION X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION UV(NPOIN),VV(NPOIN),UR(NP),VR(NP)
      DOUBLE PRECISION U1(NPOIN),V1(NPOIN),U2(NPOIN),V2(NPOIN)
      DOUBLE PRECISION XRELV(NP),YRELV(NP),TRA(NP)
      DOUBLE PRECISION AT,TV1,TV2
      DOUBLE PRECISION DDC,DAT2,DAT2B(1),Z(1),C,COEF
C
      CHARACTER*3 BINDON, C1
      CHARACTER*7 CHDON
C
C-----------------------------------------------------------------------
C
      REAL, ALLOCATABLE :: W(:)
      ALLOCATE(W(NP))
C
C-----------------------------------------------------------------------
C
      IF (AT.GT.TV2) THEN
C
C       ----------------------------------------------------------------
C        ON CHANGE D'ENREGISTREMENT : 2->1 ET ON LIT UN NOUVEAU 2
C       ----------------------------------------------------------------
        TV1=TV2
        CALL OV('X=Y     ', U1 , U2 , Z , C , NPOIN)
        CALL OV('X=Y     ', V1 , V2 , Z , C , NPOIN)
C
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) '   NOUDON : LECTURE D''UN NOUVEL ENREGISTREMENT'
        ELSE
          WRITE(LU,*) '   NOUDON : READING A NEW RECORDING'
        ENDIF
C
        IF (INDIC.EQ.1) THEN
C
C     ------------------------------------------------------------------
C          FICHIER DIFFERENCES FINIES FORMATTE DU TYPE WAM CYCLE 4
C     ------------------------------------------------------------------
 90        CONTINUE
C          LECTURE : DATE DE L'ENREGISTREMENT
           READ(NDON,*,END=100,ERR=100) DAT2
           CALL TEMP(TV2,DAT2,DDC)
C          LECTURE : DONNEES
           READ(NDON,*,END=100,ERR=100)
           READ(NDON,20,END=100,ERR=100) (UR(I),I=1,NP)
           READ(NDON,*,END=100,ERR=100)
           READ(NDON,20,END=100,ERR=100) (VR(I),I=1,NP)
C
           IF (TV2.LT.AT) THEN
             IF(LNG.EQ.1) THEN
               WRITE(LU,*) ' NOUDON : ON SAUTE 1 ENREGISTREMENT ..'
             ELSE
               WRITE(LU,*) ' NOUDON : JUMP OF 1 RECORDED DATA SERIES'
             ENDIF
             TV1=TV2
             CALL FASP(X,Y,U1,NPOIN,XRELV,YRELV,UR,NP,NBOR,
     *                                       MESH%KP1BOR%I,NPTFR,0.D0)
             CALL FASP(X,Y,V1,NPOIN,XRELV,YRELV,VR,NP,NBOR,
     *                                       MESH%KP1BOR%I,NPTFR,0.D0)
             GOTO 90
           ENDIF
C
           CALL FASP(X,Y,U2,NPOIN,XRELV,YRELV,UR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,0.D0)
           CALL FASP(X,Y,V2,NPOIN,XRELV,YRELV,VR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,0.D0)
C
        ELSEIF (INDIC.EQ.3) THEN
C
C     ------------------------------------------------------------------
C       FICHIER SELAFIN DU TYPE TELEMAC
C     ------------------------------------------------------------------
C
        ID(1)=1
        ID(2)=2
 95     CONTINUE
C       LECTURE : DATE DE L'ENREGISTREMENT
        CALL LIT(DAT2B,W,IW,C1,1,'R4',NDON,BINDON,ISTAT)
        IF(CHDON(1:1).EQ.'C') THEN
         TV2=DAT2B(1)
        ELSE
         DAT2=DAT2B(1)*1.D2
         CALL TEMP(TV2,DAT2,DDC)
        ENDIF
C       LECTURE : DONNEES
       DO I =1,NVAR
        IF(I.EQ.ID(1)) THEN
         CALL LIT(UR,W,IW,C1,NP,'R4',NDON,BINDON,ISTAT)
        ELSEIF(I.EQ.ID(2)) THEN
         CALL LIT(VR,W,IW,C1,NP,'R4',NDON,BINDON,ISTAT)
        ELSE
         READ(NDON)
        ENDIF
       ENDDO
C
        IF (TV2.LT.AT) THEN
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) ' NOUDON : ON SAUTE 1 ENREGISTREMENT ..'
          ELSE
            WRITE(LU,*) ' NOUDON : JUMP OF 1 RECORDED DATA SERIES'
          ENDIF
          TV1=TV2
C         INTERPOLATION SPATIALE DES DONNEES
          CALL FASP(X,Y,U1,NPOIN,XRELV,YRELV,UR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,0.D0)
          CALL FASP(X,Y,V1,NPOIN,XRELV,YRELV,VR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,0.D0)
          GOTO 95
        ENDIF
C
        WRITE(LU,*) 'T',CHDON,'1:',TV1
        WRITE(LU,*) 'T',CHDON,'2:',TV2
C
C       INTERPOLATION SPATIALE DES DONNEES
        CALL FASP(X,Y,U2,NPOIN,XRELV,YRELV,UR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,0.D0)
        CALL FASP(X,Y,V2,NPOIN,XRELV,YRELV,VR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,0.D0)
C
        ELSEIF (INDIC.EQ.4) THEN
C
C     ------------------------------------------------------------------
C        LECTURE D'UN FORMAT DEFINI PAR L'UTILISATEUR
C     ------------------------------------------------------------------
C
          IF(CHDON(1:1).EQ.'C') THEN
          CALL COUUTI
     *    (X,Y,NPOIN,NDON,BINDON,NBOR,NPTFR,AT,DDC,TV1,TV2,
     *     NP,XRELV,YRELV,UR,VR,U1,V1,U2,V2,NP)
          ELSE
          CALL VENUTI
     *    (X,Y,NPOIN,NDON,BINDON,NBOR,NPTFR,AT,DDC,TV1,TV2,
     *     NP,XRELV,YRELV,UR,VR,U1,V1,U2,V2,NP)
          ENDIF
C
C
        ELSE
C
        WRITE(LU,*) '************************************************'
        IF(LNG.EQ.1) THEN
         WRITE(LU,*) 'NOUDON : INDICATEUR DE FORMAT INCONNU : ',INDIC
        ELSE
          WRITE(LU,*)'NOUDON : UNKNOWN INDICATOR OF FORMAT : ',INDIC
        ENDIF
        WRITE(LU,*) '************************************************'
        CALL PLANTE(1)
        ENDIF
C
      ENDIF
C
C       --------------------------------------------------------------
C          INTERPOLATION
C       --------------------------------------------------------------
C
      COEF=(AT-TV1)/(TV2-TV1)
      DO 60 I=1,NPOIN
         UV(I)=(U2(I)-U1(I))*COEF+U1(I)
         VV(I)=(V2(I)-V1(I))*COEF+V1(I)
60    CONTINUE
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
         WRITE(LU,*)'  ERREUR A LA LECTURE DU FICHIER DE DONNEES  '
         WRITE(LU,*)'      OU FIN DE FICHIER PREMATUREE           '
      ELSE
         WRITE(LU,*)'  ERROR WHILE READING DATA FILE '
         WRITE(LU,*)'    OR UNEXPECTED END OF FILE           '
      ENDIF
      WRITE(LU,*)'*********************************************'
      DEALLOCATE(W)
      CALL PLANTE(1)
C
      RETURN
      END
