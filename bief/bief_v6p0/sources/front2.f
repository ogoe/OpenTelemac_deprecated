C                       *****************
                        SUBROUTINE FRONT2
C                       *****************
C
     *(NFRLIQ,NFRSOL,DEBLIQ,FINLIQ,DEBSOL,FINSOL,LIHBOR,LIUBOR,
     * X,Y,NBOR,KP1BOR,DEJAVU,NPOIN,NPTFR,KLOG,LISTIN,NUMLIQ,MAXFRO)
C
C***********************************************************************
C BIEF VERSION 5.6           27/02/04    J-M HERVOUET  30 87 80 18
C***********************************************************************
C
C  FONCTION  : REPERAGE, NUMEROTATION DES FRONTIERES LIQUIDES ET SOLIDES
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   NFRLIQ       |<-- | NOMBRE DE FRONTIERES LIQUIDES
C |   NFRSOL       |<-- | NOMBRE DE FRONTIERES SOLIDES
C |   DEBLIQ       |<-- | DEBUTS DES FRONTIERES LIQUIDES
C |   FINLIQ       |<-- | FINS DES FRONTIERES LIQUIDES
C |   DEBSOL       |<-- | DEBUTS DES FRONTIERES SOLIDES
C |   FINSOL       |<-- | FINS DES FRONTIERES SOLIDES
C |   LIHBOR       | -->| CONDITIONS AUX LIMITES SUR H
C |   X , Y        | -->| COORDONNEES DU MAILLAGE.
C |   NBOR         | -->| NUMEROS GLOBAUX DES POINTS DE BORD
C |   KP1BOR       | -->| NUMEROS DES EXTREMITES DES SEGMENTS DE BORD
C |                |    | DANS LA NUMEROTATION DES POINTS DE BORD
C |   DEJAVU       | -- | TABLEAU DE TRAVAIL
C |   NPOIN        | -->| NOMBRE DE POINTS DU MAILLAGE
C |   NPTFR        | -->| NOMBRE DE POINTS FRONTIERE
C |   KLOG         | -->| LIHBOR(K)=KLOG : FRONTIERE SOLIDE
C |   LISTIN       | -->| IMPRESSIONS SUR LISTING (OU NON)
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C  PRECAUTIONS D'EMPLOI : LES FRONTIERES SOLIDES SONT REPEREES PAR LE
C                         FAIT QUE LIHBOR(K) = KLOG POUR UN POINT DE
C                         BORD DE NUMERO K.
C                         UN SEGMENT COMPRIS ENTRE UN POINT LIQUIDE ET
C                         UN POINT SOLIDE EST SOLIDE.
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)  :: NPOIN,NPTFR,KLOG,MAXFRO
      INTEGER, INTENT(OUT) :: NFRLIQ,NFRSOL
      INTEGER, INTENT(OUT) :: DEBLIQ(MAXFRO),FINLIQ(MAXFRO)
      INTEGER, INTENT(OUT) :: DEBSOL(MAXFRO),FINSOL(MAXFRO)
      INTEGER , INTENT(IN) :: LIHBOR(NPTFR),LIUBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN) :: X(NPOIN) , Y(NPOIN)
      INTEGER, INTENT(IN) :: NBOR(NPTFR),KP1BOR(NPTFR)
      INTEGER, INTENT(OUT) :: DEJAVU(NPTFR)
      LOGICAL, INTENT(IN) :: LISTIN
      INTEGER, INTENT(OUT) :: NUMLIQ(NPTFR)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C 
      INTEGER K,KPREV,IDEP,SOL1,LIQ1,L1,L2,L3,NILE
C
      LOGICAL SOLF,LIQF,SOLD,LIQD
C
      DOUBLE PRECISION MINNS,MAXNS,EPS,YMIN,NS
C
      INTRINSIC ABS
C
C-----------------------------------------------------------------------
C
C  INITIALISATIONS
C
C  DEJAVU : MARQUE D'UN 1 LES POINTS DEJA TRAITES
C  NILE   : NOMBRE D'ILES
C
      DO 10 K=1,NPTFR
        DEJAVU(K) = 0
10    CONTINUE
C
      NILE = 0
      IDEP = 1
      NFRLIQ = 0
      NFRSOL = 0
C
C-----------------------------------------------------------------------
C
C  ON REVIENDRA A L'ETIQUETTE 20 S'IL Y A AU MOINS UNE ILE
C
20    CONTINUE
C
C  RECHERCHE DU POINT LE PLUS SUD-OUEST (IL PEUT Y EN AVOIR PLUSIEURS)
C
      MINNS = X(NBOR(IDEP)) + Y(NBOR(IDEP))
      MAXNS = MINNS
      YMIN  = Y(NBOR(IDEP))
C
      DO 30 K = 1 , NPTFR
      IF(DEJAVU(K).EQ.0) THEN
        NS = X(NBOR(K)) + Y(NBOR(K))
        IF(NS.LT.MINNS) THEN
         IDEP = K
         MINNS = NS
         YMIN = Y(NBOR(K))
        ENDIF
        IF(NS.GT.MAXNS) MAXNS = NS
      ENDIF
30    CONTINUE
C
      EPS = (MAXNS-MINNS) * 1.D-4
C
C  CHOIX DU POINT LE PLUS SUD PARMI LES CANDIDATS SUD-OUEST
C
      DO 40 K = 1 , NPTFR
      IF(DEJAVU(K).EQ.0) THEN
        NS = X(NBOR(K)) + Y(NBOR(K))
        IF(ABS(MINNS-NS).LT.EPS) THEN
          IF(Y(NBOR(K)).LT.YMIN) THEN
           IDEP = K
           YMIN = Y(NBOR(K))
          ENDIF
        ENDIF
      ENDIF
40    CONTINUE
C
C-----------------------------------------------------------------------
C
C  NUMEROTATION ET REPERAGE DES FRONTIERES DU CONTOUR COMMENCANT
C  AU POINT IDEP.
C
C  SOLD = .TRUE. : LA FRONTIERE AU DEPART DE IDEP EST SOLIDE
C  LIQD = .TRUE. : LA FRONTIERE AU DEPART DE IDEP EST LIQUIDE
C  SOLF = .TRUE. : LA FRONTIERE AU RETOUR A IDEP EST SOLIDE
C  LIQF = .TRUE. : LA FRONTIERE AU RETOUR A IDEP EST LIQUIDE
C  LIQ1 : NUMERO DE LA PREMIERE FRONTIERE LIQUIDE DU CONTOUR
C  SOL1 : NUMERO DE LA PREMIERE FRONTIERE SOLIDE DU CONTOUR
C
      K = IDEP
C
      SOL1 = 0
      LIQ1 = 0
      LIQF = .FALSE.
      SOLF = .FALSE.
C
C NATURE DU PREMIER SEGMENT
C
C     LOI DE DOMINANCE DU SOLIDE SUR LE LIQUIDE
      IF(LIHBOR(K).EQ.KLOG.OR.LIHBOR(KP1BOR(K)).EQ.KLOG) THEN
C       LE PREMIER SEGMENT EST SOLIDE
        NFRSOL = NFRSOL + 1
        SOL1 = NFRSOL
        SOLD = .TRUE.
        LIQD = .FALSE.
      ELSE
C       LE PREMIER SEGMENT EST LIQUIDE
        NFRLIQ = NFRLIQ + 1
        LIQ1 = NFRLIQ
        LIQD = .TRUE.
        SOLD = .FALSE.
      ENDIF
C
      DEJAVU(K) = 1
      KPREV = K
      K = KP1BOR(K)
C
50    CONTINUE
C
C RECHERCHE DES POINTS DE TRANSITION A PARTIR DU POINT SUIVANT IDEB
C
C ON CHERCHE AUSSI LES CAS DE POINTS ISOLES POUR DETECTER LES ERREURS
C DANS LES DONNEES.
C
      L1 = LIHBOR(KPREV)
      L2 = LIHBOR(K)
      L3 = LIHBOR(KP1BOR(K))
C
      IF(L1.EQ.KLOG.AND.L2.NE.KLOG.AND.L3.NE.KLOG) THEN
C     TRANSITION SOLIDE-LIQUIDE AU POINT K
        NFRLIQ = NFRLIQ + 1
        FINSOL(NFRSOL) = K
        DEBLIQ(NFRLIQ) = K
        LIQF = .TRUE.
        SOLF = .FALSE.
      ELSEIF(L1.NE.KLOG.AND.L2.NE.KLOG.AND.L3.EQ.KLOG) THEN
C     TRANSITION LIQUIDE-SOLIDE AU POINT K
        NFRSOL = NFRSOL + 1
        FINLIQ(NFRLIQ) = K
        DEBSOL(NFRSOL) = K
        LIQF = .FALSE.
        SOLF = .TRUE.
      ELSEIF(L1.NE.KLOG.AND.L2.NE.KLOG.AND.L3.NE.KLOG) THEN
C     RECHERCHE DES TRANSITIONS LIQUIDE-LIQUIDE AU POINT K
        IF(L2.NE.L3.OR.LIUBOR(K).NE.LIUBOR(KP1BOR(K))) THEN
          FINLIQ(NFRLIQ) = K
          NFRLIQ = NFRLIQ + 1
          DEBLIQ(NFRLIQ) = KP1BOR(K)
        ENDIF
      ELSEIF(L1.EQ.KLOG.AND.L2.NE.KLOG.AND.L3.EQ.KLOG) THEN
C     ERREUR DANS LES DONNEES
        IF(LNG.EQ.1) WRITE(LU,102) K
        IF(LNG.EQ.2) WRITE(LU,103) K
        CALL PLANTE(1)
        STOP
      ELSEIF(L1.NE.KLOG.AND.L2.EQ.KLOG.AND.L3.NE.KLOG) THEN
C     ERREUR DANS LES DONNEES
        IF(LNG.EQ.1) WRITE(LU,104) K
        IF(LNG.EQ.2) WRITE(LU,105) K
        CALL PLANTE(1)
        STOP
      ENDIF
C
      DEJAVU(K) = 1
      KPREV = K
      K = KP1BOR(K)
      IF(K.NE.IDEP) GO TO 50
C
C  CAS D'UN CHANGEMENT DE FRONTIERE AU POINT DE DEPART IDEP
C
      IF(SOLF) THEN
C       LA DERNIERE FRONTIERE ETAIT SOLIDE
        IF(SOLD) THEN
C         LA PREMIERE FRONTIERE ETAIT SOLIDE
          DEBSOL(SOL1) = DEBSOL(NFRSOL)
          NFRSOL = NFRSOL - 1
        ELSEIF(LIQD) THEN
C         LA PREMIERE FRONTIERE ETAIT LIQUIDE
          DEBLIQ(LIQ1) = IDEP
          FINSOL(NFRSOL) = IDEP
        ENDIF
C
      ELSEIF(LIQF) THEN
C       LA DERNIERE FRONTIERE DU CONTOUR ETAIT LIQUIDE
        IF(LIQD) THEN
C         LA PREMIERE FRONTIERE DU CONTOUR ETAIT LIQUIDE
          DEBLIQ(LIQ1) = DEBLIQ(NFRLIQ)
          NFRLIQ = NFRLIQ - 1
        ELSEIF(SOLD) THEN
C         LA PREMIERE FRONTIERE DU CONTOUR ETAIT SOLIDE
          DEBSOL(SOL1) = IDEP
          FINLIQ(NFRLIQ) = IDEP
        ENDIF
C
      ELSE
C     CAS OU TOUT LE CONTOUR EST DU MEME TYPE
        IF(SOL1.NE.0) THEN
          DEBSOL(SOL1) = IDEP
          FINSOL(SOL1) = IDEP
        ELSEIF(LIQ1.NE.0) THEN
          DEBLIQ(LIQ1) = IDEP
          FINLIQ(LIQ1) = IDEP
        ELSE
          IF(LISTIN.AND.LNG.EQ.1) THEN
           WRITE(LU,'(1X,A)') 'CAS IMPOSSIBLE DANS FRONT2'
          ENDIF
          IF(LISTIN.AND.LNG.EQ.2) THEN
           WRITE(LU,'(1X,A)') 'IMPOSSIBLE CASE IN FRONT2'
          ENDIF
          CALL PLANTE(1)
          STOP
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C
C  ON REGARDE S'IL RESTE DES CONTOURS :
C
      DO 60 K = 1 , NPTFR
        IF(DEJAVU(K).EQ.0) THEN
          IDEP = K
          NILE = NILE + 1
          GO TO 20
        ENDIF
60    CONTINUE
C
C-----------------------------------------------------------------------
C
      DO 79 K=1,NPTFR
        NUMLIQ(K)=0
79    CONTINUE
C
C  IMPRESSION DES RESULTATS ET CALCUL DE NUMLIQ
C
      IF(NILE.NE.0.AND.LISTIN.AND.LNG.EQ.1) WRITE(LU,69) NILE
      IF(NILE.NE.0.AND.LISTIN.AND.LNG.EQ.2) WRITE(LU,169) NILE
C
      IF(NFRLIQ.NE.0) THEN
        IF(LISTIN.AND.LNG.EQ.1) WRITE(LU,70) NFRLIQ
        IF(LISTIN.AND.LNG.EQ.2) WRITE(LU,170) NFRLIQ

        DO 80 K = 1, NFRLIQ
C
C  MARQUAGE DES NUMEROS DES FRONTIERES LIQUIDES
C
          L1=DEBLIQ(K)
          NUMLIQ(L1)=K
707       L1=KP1BOR(L1)
          NUMLIQ(L1)=K
          IF(L1.NE.FINLIQ(K)) GO TO 707
C
C  FIN DU MARQUAGE
C
          IF(LISTIN.AND.LNG.EQ.1) WRITE(LU,90)
     *                            K,DEBLIQ(K),NBOR(DEBLIQ(K)),
     *                            X(NBOR(DEBLIQ(K))),Y(NBOR(DEBLIQ(K))),
     *                            FINLIQ(K),NBOR(FINLIQ(K)),
     *                            X(NBOR(FINLIQ(K))),Y(NBOR(FINLIQ(K)))
          IF(LISTIN.AND.LNG.EQ.2) WRITE(LU,190)
     *                            K,DEBLIQ(K),NBOR(DEBLIQ(K)),
     *                            X(NBOR(DEBLIQ(K))),Y(NBOR(DEBLIQ(K))),
     *                            FINLIQ(K),NBOR(FINLIQ(K)),
     *                            X(NBOR(FINLIQ(K))),Y(NBOR(FINLIQ(K)))
80      CONTINUE
      ENDIF
C
      IF(NFRSOL.NE.0) THEN
        IF(LISTIN.AND.LNG.EQ.1) WRITE(LU,100) NFRSOL
        IF(LISTIN.AND.LNG.EQ.2) WRITE(LU,101) NFRSOL
        DO 110 K = 1, NFRSOL
          IF(LISTIN.AND.LNG.EQ.1) WRITE(LU,90)
     *                            K,DEBSOL(K),NBOR(DEBSOL(K)),
     *                            X(NBOR(DEBSOL(K))),Y(NBOR(DEBSOL(K))),
     *                            FINSOL(K),NBOR(FINSOL(K)),
     *                            X(NBOR(FINSOL(K))),Y(NBOR(FINSOL(K)))
          IF(LISTIN.AND.LNG.EQ.2) WRITE(LU,190)
     *                            K,DEBSOL(K),NBOR(DEBSOL(K)),
     *                            X(NBOR(DEBSOL(K))),Y(NBOR(DEBSOL(K))),
     *                            FINSOL(K),NBOR(FINSOL(K)),
     *                            X(NBOR(FINSOL(K))),Y(NBOR(FINSOL(K)))
110     CONTINUE
      ENDIF
C
C-----------------------------------------------------------------------
C
C  FORMATS
C
69    FORMAT(/,1X,'IL Y A ',1I3,' ILE(S) DANS LE DOMAINE')
169   FORMAT(/,1X,'THERE IS ',1I3,' ISLAND(S) IN THE DOMAIN')
70    FORMAT(/,1X,'IL Y A ',1I3,' FRONTIERE(S) LIQUIDE(S) :')
170   FORMAT(/,1X,'THERE IS ',1I3,' LIQUID BOUNDARIES:')
100   FORMAT(/,1X,'IL Y A ',1I3,' FRONTIERE(S) SOLIDE(S) :')
101   FORMAT(/,1X,'THERE IS ',1I3,' SOLID BOUNDARIES:')
102   FORMAT(/,1X,'FRONT2 : ERREUR AU POINT DE BORD ',1I5,
     *       /,1X,'         POINT LIQUIDE ENTRE DEUX POINTS SOLIDES')
103   FORMAT(/,1X,'FRONT2 : ERROR AT BOUNDARY POINT ',1I5,
     *       /,1X,'         LIQUID POINT BETWEEN TWO SOLID POINTS')
104   FORMAT(/,1X,'FRONT2 : ERREUR AU POINT DE BORD ',1I5,
     *       /,1X,'         POINT SOLIDE ENTRE DEUX POINTS LIQUIDES')
105   FORMAT(/,1X,'FRONT2 : ERROR AT BOUNDARY POINT ',1I5,
     *       /,1X,'         SOLID POINT BETWEEN TWO LIQUID POINTS')
90    FORMAT(/,1X,'FRONTIERE ',1I3,' : ',/,1X,
     *            ' DEBUT AU POINT DE BORD ',1I4,
     *            ' , DE NUMERO GLOBAL ',1I6,/,1X,
     *            ' ET DE COORDONNEES : ',G16.7,3X,G16.7,
     *       /,1X,' FIN AU POINT DE BORD ',1I4,
     *            ' , DE NUMERO GLOBAL ',1I6,/,1X,
     *            ' ET DE COORDONNEES : ',G16.7,3X,G16.7)
190   FORMAT(/,1X,'BOUNDARY ',1I3,' : ',/,1X,
     *            ' BEGINS AT BOUNDARY POINT: ',1I4,
     *            ' , WITH GLOBAL NUMBER: ',1I6,/,1X,
     *            ' AND COORDINATES: ',G16.7,3X,G16.7,
     *       /,1X,' ENDS AT BOUNDARY POINT: ',1I4,
     *            ' , WITH GLOBAL NUMBER: ',1I6,/,1X,
     *            ' AND COORDINATES: ',G16.7,3X,G16.7)
C
C-----------------------------------------------------------------------
C
      IF(NFRSOL.GT.MAXFRO.OR.NFRLIQ.GT.MAXFRO) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'FRONT2 : DEPASSEMENT DE TABLEAUX'
          WRITE(LU,*) '         AUGMENTER MAXFRO DANS LE CODE APPELANT' 
          WRITE(LU,*) '         A LA VALEUR ',MAX(NFRSOL,NFRLIQ)             
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'FRONT2: SIZE OF ARRAYS EXCEEDED'
          WRITE(LU,*) '        INCREASE MAXFRO IN THE CALLING PROGRAM' 
          WRITE(LU,*) '        UP TO THE VALUE ',MAX(NFRSOL,NFRLIQ)            
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
