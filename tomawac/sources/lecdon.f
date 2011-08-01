C                       *****************
                        SUBROUTINE LECDON
C                       *****************
C
     *( U , V , X, Y, NPOIN2, NDON, BINDON, NBOR, NPTFR, XRELV, YRELV,
     *  UR, VR, TRA01,TRA02,TRA03,IDTEL,NPTT,DONTEL,COURAN,INDIC,NPMAX,
     *  CHDON)
C
C***********************************************************************
C  COWADIS VERSION 1.0    01/02/95        F.MARCOS     (LNH) 30 87 72 66
C***********************************************************************
C
C   FONCTION : CE SOUS-PROGRAMME PROJETE LA VALEUR DES COURANTS OU VENTS
C              SUR LE MAILLAGE DE CALCUL
C        (INSPIRE DE LA ROUTINE FOND DE TELEMAC 2D)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    U,V         !<-- !  COURANT OU VENT AUX NOEUDS DU MAILLAGE      !
C !    X,Y         ! -->!  COORDONNEES DU MAILLAGE                     !
C !    NPOIN2      ! -->!  NOMBRE DE POINTS DU MAILLAGE                !
C !    NDON        ! -->!  NUMERO D'UNITE LOGIQUE DU FICHIER           !
C !    BINDON      ! -->!  BINAIRE DU FICHIER DES DONNEES  (INDIC>2)   !
C !    NBOR        ! -->!  NUMEROTATION DES POINTS FRONTIERE           !
C !    NPTFR       ! -->!  NOMBRE DE  POINTS FRONTIERE                 !
C !    XRELV       !<-->!  TABLEAU DES ABSCISSES DES POINTS RELEVES    !
C !    YRELV       !<-->!  TABLEAU DES ORDONNEES DES POINTS RELEVES    !
C !    UR,VR       !<-->!  TABLEAU DES COURANTS RELEVES                !
C !    TRA01       !<-->!  TABLEAU DE TRAVAIL                          !
C !    TRA02       !<-->!  TABLEAU DE TRAVAIL                          !
C !    TRA03       !<-->!  TABLEAU DE TRAVAIL                          !
C !    IDTEL       ! -->!  RANG DE LA VARIABLE TELEMAC A RECUPERER     !
C !    DONTEL      ! -->!  LOGIQUE INDIQUANT SI ON RECUPERE            !
C !                !    !  UNE VARIABLE TELEMAC                        !
C !    COURAN      ! -->!  LOGIQUE INDIQUANT LA PRESENCE DE DONNEES    !
C !    INDIC       ! -->!  TYPE DE FORMAT DE LECTURE                   !
C !    NPMAX       ! -->!  NOMBRE DE POINTS RELEVES MAXIMUM            !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : WAC
C
C SOUS-PROGRAMME APPELE : COUUTI,FASP
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
      INTEGER NP,NDON,NPOIN2,NPTFR,INDIC,NCOL,NLIG,BID,I,J
C
      INTEGER NPMAX,NBOR(NPTFR,2)
C
      DOUBLE PRECISION X(NPOIN2),Y(NPOIN2),U(NPOIN2),V(NPOIN2)
      DOUBLE PRECISION TRA01(NPMAX),TRA02(NPMAX),TRA03(NPOIN2)
      DOUBLE PRECISION XRELV(NPMAX),YRELV(NPMAX),UR(NPMAX),VR(NPMAX)
      DOUBLE PRECISION XMAX,XMIN,YMAX,YMIN,DX,DY,ATT,BDX(2)
C
      INTEGER NVAR,IB(10),ISTAT,NPTT,ID(3),IDTEL
C
      CHARACTER*3  BINDON,C
      CHARACTER*7  CHDON
      CHARACTER*32 TEXTE(20)
      CHARACTER*72 TITCAS
C
      LOGICAL DONTEL,COURAN
C
C-----------------------------------------------------------------------      
C
      REAL, ALLOCATABLE :: W(:)
      ALLOCATE(W(NPMAX))
C
C-----------------------------------------------------------------------
C        LECTURE DES POINTS RELEVES SUR UNITE LOGIQUE NDON
C-----------------------------------------------------------------------
C
      IF (INDIC.EQ.1) THEN
C
C     ------------------------------------------------------------------
C     FORMAT DE TYPE WAM - DIFFERENCES FINIES
C     ------------------------------------------------------------------
C
         READ(NDON,10,END=100,ERR=100)
     *      NCOL,NLIG,YMIN,YMAX,XMIN,XMAX,BID,BID
         DX=(XMAX-XMIN)/REAL(NCOL-1)
         DY=(YMAX-YMIN)/REAL(NLIG-1)
         NP=NCOL*NLIG
         WRITE(LU,*)
     *    '-----------------------------------------------------'
         IF (LNG.EQ.1) THEN
            WRITE(LU,*)'LECDON : LECTURE DU FICHIER DE DONNEES'
            WRITE(LU,*)'         NOMBRE DE LIGNES   : ',NLIG
            WRITE(LU,*)'         NOMBRE DE COLONNES : ',NCOL
            WRITE(LU,*)'         ABSCISSE MINIMALE : ',XMIN
            WRITE(LU,*)'         ABSCISSE MAXIMALE : ',XMAX
            WRITE(LU,*)'         ORDONNEE MINIMALE : ',YMIN
            WRITE(LU,*)'         ORDONNEE MAXIMALE : ',YMAX
         ELSE
            WRITE(LU,*)'LECDON : READING OF THE DATA FILE '
            WRITE(LU,*)'         NUMBER OF LINES   : ',NLIG
            WRITE(LU,*)'         NUMBER OF COLUMNS : ',NCOL
            WRITE(LU,*)'         MINIMAL ABSCISSAE : ',XMIN
            WRITE(LU,*)'         MAXIMAL ABSCISSAE : ',XMAX
            WRITE(LU,*)'         MINIMAL ORDINATES : ',YMIN
            WRITE(LU,*)'         MAXIMAL ORDINATES : ',YMAX
         ENDIF 
         WRITE(LU,*)
     *    '-----------------------------------------------------'
         IF (NP.GT.NPMAX) THEN
          WRITE(LU,*)
     *     '*****************************************************'
          IF (LNG.EQ.1) THEN
             WRITE(LU,*)
     *        ' LA DIMENSION PREVUE PAR DEFAUT POUR LE TABLEAU '
             WRITE(LU,*)
     *        ' DE DONNEES :',NPMAX,' EST TROP FAIBLE POUR '
             WRITE(LU,*)
     *        ' CONTENIR LA TOTALITE DES DONNEES :',NP
          ELSE
             WRITE(LU,*)
     *        ' THE DEFAULT DIMENSION ALLOWED FOR THE ARRAY OF '
             WRITE(LU,*)
     *        ' DATA :',NPMAX,' IS TOO LOW TO HOLD'
             WRITE(LU,*)
     *        ' ALL THE DATA :',NP
          ENDIF
          WRITE(LU,*)
     *        '*****************************************************'
          CALL PLANTE(0)
         ENDIF
         READ(NDON,*)
         READ(NDON,20,END=100,ERR=100)
     *      (UR(I),I=1,NCOL*NLIG)
         READ(NDON,*)
         READ(NDON,20,END=100,ERR=100)
     *      (VR(I),I=1,NCOL*NLIG)
         DO 30 I=1,NCOL
             DO 40 J=1,NLIG
                XRELV((I-1)*NLIG+J)=XMIN+DX*(I-1)
                YRELV((I-1)*NLIG+J)=YMIN+DY*(J-1)
40       CONTINUE
30       CONTINUE
C
C
      ELSEIF (INDIC.EQ.2) THEN
C 
C     ------------------------------------------------------------------
C     FORMAT DE TYPE SINUSX - SEMIS DE POINTS
C     ------------------------------------------------------------------
C
          DO 50 I=1,100000
            READ(NDON,*,END=60,ERR=100)
            NP=I
50        CONTINUE
60        CONTINUE
          WRITE(LU,*)
     *     '-----------------------------------------------------'
          IF (LNG.EQ.1) THEN
             WRITE(LU,*)'LECDON : LECTURE DU FICHIER DE DONNEES'
             WRITE(LU,*)'         NOMBRE DE POINTS   : ',NP
          ELSE
             WRITE(LU,*)'LECDON : READING OF THE DATA FILE '
             WRITE(LU,*)'         NUMBER OF POINTS   : ',NP
          ENDIF   
          WRITE(LU,*)
     *     '-----------------------------------------------------'
          IF (NP.GT.NPMAX) THEN
           WRITE(LU,*)
     *     '*****************************************************'
           IF (LNG.EQ.1) THEN
             WRITE(LU,*)
     *        ' LA DIMENSION PREVUE PAR DEFAUT POUR LE TABLEAU '
             WRITE(LU,*)
     *        ' DE DONNEES :',NPMAX,' EST TROP FAIBLE POUR '
             WRITE(LU,*)
     *        ' CONTENIR LA TOTALITE DES DONNEES :',NP
           ELSE
             WRITE(LU,*)
     *        ' THE DEFAULT DIMENSION ALLOWED FOR THE ARRAY OF '
             WRITE(LU,*)
     *        ' DATA :',NPMAX,' IS TOO LOW TO HOLD'
             WRITE(LU,*)
     *        ' ALL THE DATA :',NP
           ENDIF
           WRITE(LU,*)
     *        '*****************************************************'
           CALL PLANTE(0)
          ENDIF
          REWIND NDON
          DO 70 I=1,NP
             READ(NDON,*,ERR=100) XRELV(I),YRELV(I),UR(I),VR(I)
70        CONTINUE
C
C
      ELSEIF (INDIC.EQ.3) THEN
C 
C     ------------------------------------------------------------------
C     FORMAT DE TYPE TELEMAC - MAILLAGE EVENTUELLEMENT DIFFERENT 
C          (BINAIRE)                 DU MAILLAGE COWADIS
C     ------------------------------------------------------------------
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
80        CONTINUE
C
C      VARIABLES FORMAT ET GEOMETRIE
C
          READ(NDON)
          CALL LIT(X,W,IB,C,4,'I ',NDON,BINDON,ISTAT)
          NP=IB(2)
          WRITE(LU,*)
     *        '-----------------------------------------------------'
          IF (LNG.EQ.1) THEN
             WRITE(LU,*)'LECDON : LECTURE DU FICHIER TELEMAC'
             WRITE(LU,*)'         NOMBRE DE POINTS   : ',NP
          ELSE
             WRITE(LU,*)'LECDON : READING OF TELEMAC DATA FILE '
             WRITE(LU,*)'         NUMBER OF POINTS   : ',NP
          ENDIF   
          IF (NP.GT.NPMAX) THEN
           WRITE(LU,*)
     *        '*****************************************************'
           IF (LNG.EQ.1) THEN
             WRITE(LU,*)
     *        ' LA DIMENSION PREVUE PAR DEFAUT POUR LE TABLEAU '
             WRITE(LU,*)
     *        ' DE DONNEES :',NPMAX,' EST TROP FAIBLE POUR '
             WRITE(LU,*)
     *        ' CONTENIR LA TOTALITE DES DONNEES :',NP
           ELSE
             WRITE(LU,*)
     *        ' THE DEFAULT DIMENSION ALLOWED FOR THE ARRAY OF '
             WRITE(LU,*)
     *        ' DATA :',NPMAX,' IS TOO LOW TO HOLD'
             WRITE(LU,*)
     *        ' ALL THE DATA :',NP
           ENDIF
           WRITE(LU,*)
     *        '*****************************************************'
           CALL PLANTE(0)
           ENDIF
          READ(NDON)
          READ(NDON)
C
C      X ET Y
C
          CALL LIT(XRELV,W,IB,C,NP,'R4',NDON,BINDON,ISTAT)
          CALL LIT(YRELV,W,IB,C,NP,'R4',NDON,BINDON,ISTAT)
C
C      PAS DE TEMPS ET VARIABLES
C
          DO 110 J=1,(NPTT-1)*(NVAR+1)
             READ(NDON)
110       CONTINUE
C
          IF (DONTEL) ID(3)=IDTEL
          IF (COURAN) THEN
             ID(1)=1
             ID(2)=2
          ENDIF
C
          CALL LIT(BDX(1),W,IB,C,1,'R4',NDON,BINDON,ISTAT)
          ATT=BDX(1)
          DO 90 I=1,NVAR
            IF (I.EQ.ID(1)) THEN
               CALL LIT(UR,W,IB,C,NP,'R4',NDON,BINDON,ISTAT)
            ELSEIF (I.EQ.ID(2)) THEN
               CALL LIT(VR,W,IB,C,NP,'R4',NDON,BINDON,ISTAT)
            ELSEIF ((I.EQ.ID(3)).AND.(DONTEL)) THEN
               CALL LIT(TRA01,W,IB,C,NP,'R4',NDON,BINDON,ISTAT)
            ELSE
               READ(NDON)
            ENDIF
90        CONTINUE
C
C      IMPRESSIONS SUR LE LISTING
C
          IF (LNG.EQ.1) THEN
            WRITE(LU,*)'         TITRE DU CAS TELEMAC : '
            WRITE(LU,*)'           ',TITCAS
            WRITE(LU,*)'         TEMPS DE TELEMAC : ',ATT
            WRITE(LU,*)'         VARIABLES DE TELEMAC RETENUE(S) : '
          ELSE
            WRITE(LU,*)'         TITLE OF TELEMAC CASE : '
            WRITE(LU,*)'           ',TITCAS
            WRITE(LU,*)'         TIME OF TELEMAC RECORD : ',ATT
            WRITE(LU,*)'         VARIABLE(S) OF TELEMAC READ : '
          ENDIF  
          IF (COURAN) THEN
              WRITE(LU,*)'           ',TEXTE(ID(1))
              WRITE(LU,*)'           ',TEXTE(ID(2))
          ENDIF
          IF (DONTEL)
     *            WRITE(LU,*)'           ',TEXTE(ID(3))
          WRITE(LU,*)
     *        '-----------------------------------------------------'
C
C
      ELSEIF (INDIC.EQ.4) THEN
C
C     ------------------------------------------------------------------
C        LECTURE D'UN FORMAT DEFINI PAR L'UTILISATEUR
C     ------------------------------------------------------------------
        IF(CHDON(1:1).EQ.'C') THEN
C         LECTURE D'UN CHAMP DE COURANT
              CALL COUUTI
     *    (X,Y,NPOIN2,NDON,BINDON,NBOR,NPTFR,0.,0.,0.,0.,
     *     NP,XRELV,YRELV,UR,VR,TRA03,TRA03,TRA03,TRA03,NPMAX)
        ELSEIF(CHDON(1:1).EQ.'V' .OR. CHDON(1:1).EQ.'W') THEN
C         LECTURE D'UN CHAMP DE VENT
          CALL VENUTI
     *    (X,Y,NPOIN2,NDON,BINDON,NBOR,NPTFR,0.,0.,0.,0.,
     *     NP,XRELV,YRELV,UR,VR,TRA03,TRA03,TRA03,TRA03,NPMAX)
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
        WRITE(LU,*)'***********************************************'
        IF (LNG.EQ.1) THEN
          WRITE(LU,*)'LECDON : INDICATEUR DE FORMAT INCONNU   '
          WRITE(LU,*)'         POUR LE FICHIER DES DONNEES :',INDIC
        ELSE
          WRITE(LU,*)'LECDON : INDICATOR OF FORMAT FOR THE   '
          WRITE(LU,*)'         DATA FILE UNKNOWN :',INDIC
        ENDIF
        WRITE(LU,*)'***********************************************'
          
          CALL PLANTE(0)
      ENDIF

C
C-----------------------------------------------------------------------
C   LE COURANT EST CALCULE PAR INTERPOLATION SUR LES POINTS INTERIEURS
C                         AU DOMAINE
C-----------------------------------------------------------------------
C
      IF (COURAN) THEN
        CALL FASP(X,Y,U,NPOIN2,XRELV,YRELV,UR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,1.D-6)
        CALL FASP(X,Y,V,NPOIN2,XRELV,YRELV,VR,NP,NBOR,MESH%KP1BOR%I,
     *                                                    NPTFR,1.D-6)
      ENDIF
      IF (DONTEL) 
     *CALL FASP(X,Y,TRA03,NPOIN2,XRELV,YRELV,TRA01,NP,NBOR,
     *                                      MESH%KP1BOR%I,NPTFR,1.D-6)
C
C-----------------------------------------------------------------------
C
C     FORMATS
C
10    FORMAT (2I4,4F9.3,2I2)
20    FORMAT (10f6.2)
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
      CALL PLANTE(1)
C
      DEALLOCATE(W)
      RETURN
      END
