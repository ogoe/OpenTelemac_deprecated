C                       *****************
                        SUBROUTINE TOPOGR
C                       *****************
C
     *(ZF,ZREF,ZFE,IKLE,IFABOR,NBOR,NELBOR,NULONE,
     * ITRA05,ITRA02,ITRA03,NELEM,NPTFR,NPOIN,MXPTVS)
C
C***********************************************************************
C  BIEF VERSION 5.1    17/08/94      J-M JANIN (LNH) 30 87 72 84
C
C***********************************************************************
C
C     FONCTION :
C
C       ANALYSE FINE DE LA TOPOGRAPHIE ET CONSTRUCTION DE ZFE
C
C       LE TABLEAU DES COTES DU FOND PAR ELEMENTS ZFE
C       GARANTIRA ENSUITE QUE L'ON N'AURA PAS DE DOMAINES
C       LIQUIDES CONNECTES PAR UN SEUL POINT
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    ZF          | -->| COTE DU FOND AUX NOEUDS                      |
C |    ZREF        | -- | CORRECTIF DE ZFE                             |
C |    ZFE         |<-- | COTE DU FOND PAR ELEMENT                     |
C |    IKLE        | -->| TABLE DE CONNECTIVITE                        |
C |    IFABOR      | -->| NUMERO DES ELEMENTS VOISINS                  |
C |    IFAN        | -- | NUMERO LOCAL DE LA FACE DANS L'ELEMENT VOISIN
C |    NBOR        | -->| NUMERO GLOBAL DES POINTS DE BORD             |
C |    NELBOR      | -->| NUMERO DES ELEMENTS DE BORD                  |
C |    NULONE      | -->| NUMERO LOCAL DES NOEUDS AU BORD              |
C |    ITRA01..05  | -- | TABLEAUX DE TRAVAIL D'ENTIERS                |
C |    NELEM       | -->| NOMBRE D'ELEMENTS                            |
C |    NPTFR       | -->| NOMBRE DE POINTS FRONTIERE                   |
C |    NPOIN       | -->| NOMBRE DE POINTS                             |
C |    MXPTVS      | -->| NOMBRE MAXIMUM DE POINTS VOISINS D'UN POINT
C |                |    | (A VERIFIER, LE NOMBRE D'ELEMENTS SUFFIT
C |                |    |  VRAISEMBLABLEMENT)
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C APPELE PAR:
C***********************************************************************
C                                                                      *
C TELMAC ET MITRID                                                     *
C                                                                      *
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER, INTENT(IN)    :: NELEM,NPTFR,NPOIN,MXPTVS
      INTEGER, INTENT(IN)    :: IKLE(NELEM,3),IFABOR(NELEM,3)
      INTEGER, INTENT(IN)    :: NBOR(NPTFR),NELBOR(NPTFR),NULONE(NPTFR)
      INTEGER, INTENT(INOUT) :: ITRA05(NPOIN),ITRA02(NPOIN)
      INTEGER, INTENT(INOUT) :: ITRA03(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: ZFE(NELEM),ZREF(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: ZF(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER IELEM,IPTFR,IPOIN,I,I1,I2,I3,N1,N2,ERR,IPREV(3),IMAX
      DOUBLE PRECISION EPSILO
      LOGICAL FLAG
C
C     POUR ALLOCATION DYNAMIQUE D''ENTIERS
C
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ITRA01,ITRA04,IFAN
C
      DATA IPREV  / 3 , 1 , 2 /
      DATA EPSILO / 1.D-6 /
C
C-----------------------------------------------------------------------
C
      ALLOCATE(IFAN(NELEM,3)          ,STAT=ERR)
      ALLOCATE(ITRA01(NPOIN,MXPTVS+1) ,STAT=ERR)
      ALLOCATE(ITRA04(NPOIN,6)        ,STAT=ERR)
      IF(ERR.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'TOPOGR : MAUVAISE ALLOCATION DE W'
        IF(LNG.EQ.2) WRITE(LU,*) 'TOPOGR: WRONG ALLOCATION OF W'
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C    REMPLISSAGE DE IFAN
C
C   IFAN(IELEM,IFACE) DONNE LE NUMERO LOCAL DANS L'ELEMENT DONNE
C   PAR IFABOR(IELEM,IFACE) DU NOEUD AYANT LE NUMERO LOCAL IFACE DANS
C   L'ELEMENT IELEM.
C
C    PREMIER JET POUR LE REMPLISSAGE DE ZFE
C
C    DEBUT DE REMPLISSAGE DE ITRA01 ET ITRA02
C
C   ITRA01 ET ITRA02 EFFECTUENT L'OPERATION INVERSE DE IKLE.
C   DANS CETTE BOUCLE ON REMPLIT ITRA01(IPOIN,1) QUI DONNE LE NUMERO LE
C   PLUS GRAND DES ELEMENTS CONTENANT IPOIN ET ITRA02(IPOIN) LE NUMERO
C   LOCAL DE IPOIN DANS CET ELEMENT.
C
C-----------------------------------------------------------------------
C
C>>>>
C  BOUCLES 5 ET 6 AJOUTEES PAR JMH LE 22/5/95
C  INITIALISATION POUR DETECTER LES TROUS DANS LE MAILLAGE
C  (POINTS NON RELIES A UN ELEMENT, CAS DES MAILLAGES CURVILIGNES)
C
      DO 6 I1    = 1 , MXPTVS+1
      DO 5 IPOIN = 1 , NPOIN
         ITRA01(IPOIN,I1) = 0
5     CONTINUE
6     CONTINUE
C
C<<<<
C
      DO 10 IELEM = 1,NELEM
C
         I1 = IKLE(IELEM,1)
         I2 = IKLE(IELEM,2)
         I3 = IKLE(IELEM,3)
C
         IFAN(IELEM,1) = 0
         N1 = IFABOR(IELEM,1)
         IF (N1.GT.0) THEN
            IF(IKLE(N1,1).EQ.I1) IFAN(IELEM,1) = 3
            IF(IKLE(N1,2).EQ.I1) IFAN(IELEM,1) = 1
            IF(IKLE(N1,3).EQ.I1) IFAN(IELEM,1) = 2
         ENDIF
C
         IFAN(IELEM,2) = 0
         N1 = IFABOR(IELEM,2)
         IF (N1.GT.0) THEN
            IF(IKLE(N1,1).EQ.I2) IFAN(IELEM,2) = 3
            IF(IKLE(N1,2).EQ.I2) IFAN(IELEM,2) = 1
            IF(IKLE(N1,3).EQ.I2) IFAN(IELEM,2) = 2
         ENDIF
C
         IFAN(IELEM,3) = 0
         N1 = IFABOR(IELEM,3)
         IF (N1.GT.0) THEN
            IF(IKLE(N1,1).EQ.I3) IFAN(IELEM,3) = 3
            IF(IKLE(N1,2).EQ.I3) IFAN(IELEM,3) = 1
            IF(IKLE(N1,3).EQ.I3) IFAN(IELEM,3) = 2
         ENDIF
C
         ZFE(IELEM) = MAX(ZF(I1),ZF(I2),ZF(I3))
         ITRA01(I1,1) = IELEM
         ITRA02(I1)   = 1
         ITRA01(I2,1) = IELEM
         ITRA02(I2)   = 2
         ITRA01(I3,1) = IELEM
         ITRA02(I3)   = 3
C
10    CONTINUE
C
C-----------------------------------------------------------------------
C
C    DEBUT DE REMPLISSAGE DE ITRA01 ET ITRA02 POUR LES POINTS FRONTIERE
C
C   POUR CES POINTS ITRA01(IPOIN,1) N'EST PAS UN ELEMENT QUELCONQUE,
C   C'EST L'ELEMENT QUI COMPORTE UNE FACE DE BORD ENTRE LE POINT IPOIN
C   ET LE POINT SUIVANT DANS LA NUMEROTATION LOCALE DE L'ELEMENT.
C
C-----------------------------------------------------------------------
C
      DO 20 IPTFR = 1,NPTFR
         ITRA01(NBOR(IPTFR),1) = NELBOR(IPTFR)
         ITRA02(NBOR(IPTFR))   = NULONE(IPTFR)
20    CONTINUE
C
C-----------------------------------------------------------------------
C
C    SUITE DU REMPLISSAGE DE ITRA01
C
C   ITRA01(IPOIN,I+1) EST LE NUMERO DE L'ELEMENT ADJACENT A L'ELEMENT
C   ITRA01(IPOIN,1) EN TOURNANT AUTOUR DU POINT IPOIN DANS LE SENS TRIGO
C   ITRA02 EST UN TABLEAU MONODIMENSIONNEL CAR ON N'AURA PAS BESOIN PAR
C   LA SUITE DE CONNAITRE LE NUMERO LOCAL DE IPOIN DANS UN QUELCONQUE
C   ELEMENT ITRA01(IPOIN,*).
C
C    REMPLISSAGE DE ITRA03
C
C   ITRA03(IPOIN) CORRESPOND AU NOMBRE D'ELEMENTS CONTENANT IPOIN
C   PRECEDE D'UN SIGNE NEGATIF SI IPOIN EST UN POINT DE BORD.
C
C    ATTENTION |||
C
C   ON AUTORISE A UN POINT D'ETRE CONTENU AU MAXIMUM DANS 10 ELEMENTS
C
C-----------------------------------------------------------------------
C
      DO 30 IPOIN = 1,NPOIN
         ITRA03(IPOIN) = 0
30    CONTINUE
C
      IMAX = 0
C
40    CONTINUE
      FLAG = .FALSE.
      IMAX = IMAX + 1
      IF (IMAX.GT.MXPTVS+1) THEN
        IF(LNG.EQ.1) WRITE(LU,23) MXPTVS
        IF(LNG.EQ.2) WRITE(LU,24) MXPTVS
23      FORMAT(1X,'TOPOGR : LE NOMBRE DE POINTS VOISINS MAXIMUM'/,1X,
     *            '         TROUVE DANS TOPOGR EST SUPERIEUR A ',/,1X,
     *            '         LA VALEUR ANNONCEE MXPTVS :',1I6)
24      FORMAT(1X,'TOPOGR : THE MAXIMUM NUMBER OF NEIGHBOURS TO'/,1X,
     *            '         A POINT IS GREATER THAN THE VALUE  ',/,1X,
     *            '         GIVEN BY MXPTVS :',1I6)
        CALL PLANTE(0)
        STOP
      ENDIF
C
      DO 50 IPOIN = 1,NPOIN
C
         IF (ITRA03(IPOIN).EQ.0) THEN
            N1 = ITRA01(IPOIN,IMAX)
            IF(N1.NE.0) THEN
              FLAG = .TRUE.
              N2 = IFABOR(N1,IPREV(ITRA02(IPOIN)))
C                          ICI IMAX N'EST JAMAIS A SON MAXIMUM
              ITRA01(IPOIN,IMAX+1) = N2
              ITRA02(IPOIN) = IFAN(N1,IPREV(ITRA02(IPOIN)))
              IF (N2.LE.0)               ITRA03(IPOIN) = -IMAX
              IF (N2.EQ.ITRA01(IPOIN,1)) ITRA03(IPOIN) =  IMAX
            ENDIF
         ENDIF
C
50    CONTINUE
C
      IF (FLAG) GOTO 40
C
60    CONTINUE
C
C-----------------------------------------------------------------------
C
C    RECHERCHE DES EXTREMA LOCAUX DE ZFE EN TOURNANT AUTOUR D'UN NOEUD
C
C   ITRA04(IPOIN,I) CORRESPOND AU IEME EXTREMA RENCONTRE SACHANT QUE
C   L'ELEMENT ASSOCIE EST DONNE PAR ITRA01(IPOIN,ITRA04(IPOIN,I)).
C   ITRA02(IPOIN) DONNE LE NOMBRE TOTAL DE PHASES DE CROISSANCE ET DE
C   DECROISSANCE TROUVEES PRECEDE D'UN SIGNE NEGATIF SI LA DERNIERE
C   TROUVEE EST DECROISSANTE
C
C-----------------------------------------------------------------------
C
      DO 70 IPOIN = 1,NPOIN
         ITRA02(IPOIN) = 0
         ITRA05(IPOIN) = 0
70    CONTINUE
C
      DO 80 I = 1,IMAX-1
C
         DO 90 IPOIN = 1,NPOIN
C
            IF (ITRA03(IPOIN).GE.I.OR.ITRA03(IPOIN).LT.-I) THEN
C
               N1 = ITRA01(IPOIN,I)
               N2 = ITRA01(IPOIN,I+1)
C
               IF (ZFE(N2).GT.ZFE(N1)+EPSILO) THEN
                  IF (ITRA02(IPOIN).LT.0) ITRA04(IPOIN,-ITRA02(IPOIN))=I
                  IF (ITRA02(IPOIN).LE.0) ITRA02(IPOIN)=-ITRA02(IPOIN)+1
               ELSEIF (ZFE(N2).LT.ZFE(N1)-EPSILO) THEN
                  IF (ITRA02(IPOIN).GT.0) ITRA04(IPOIN, ITRA02(IPOIN))=I
                  IF (ITRA02(IPOIN).GE.0) ITRA02(IPOIN)=-ITRA02(IPOIN)-1
               ENDIF
C
            ENDIF
C
90       CONTINUE
C
80    CONTINUE
C
      DO 95 IPOIN = 1,NPOIN
         IF((ITRA03(IPOIN).LT.0.AND.(ITRA02(IPOIN).LE.-4.OR.
     *       ITRA02(IPOIN).GE.5)).OR.ABS(ITRA02(IPOIN)).GE.6) THEN
           IF (LNG.EQ.1) THEN
            WRITE(LU,*) 'LE MAILLAGE AUTOUR DU POINT ',IPOIN,' EST TROP'
            WRITE(LU,*) 'GROSSIER PAR RAPPORT A LA BATHYMETRIE'
           ELSEIF(LNG.EQ.2) THEN
            WRITE(LU,*) 'THE MESH AROUND THE NODE ',IPOIN,' HAS TO'
            WRITE(LU,*) 'BE REFINED BECAUSE OF THE BATHYMETRY'
           ENDIF
           STOP
         ENDIF
95    CONTINUE
C
C-----------------------------------------------------------------------
C
C    CORRECTION DE ZFE EN FONCTION DE ITRA02
C
C-----------------------------------------------------------------------
C
      FLAG = .FALSE.
C
      DO 100 IPOIN = 1,NPOIN
C
         I1 = ITRA03(IPOIN)
C
         IF (I1.LT.0) THEN
C
C-----------------------------------------------------------------------
C
C    CORRECTION DE ZFE POUR LES POINTS DE BORD
C
C   SI ITRA02(IPOIN) VAUT 1, -1, 2 : PAS DE CORRECTION
C   SI ITRA02(IPOIN) VAUT -2, 3, -3, 4 : CORRECTION
C
C-----------------------------------------------------------------------
C
            IF (ITRA02(IPOIN).EQ.-2) THEN
C
               FLAG = .TRUE.
               IF (ZFE(ITRA01(IPOIN,-I1)).GT.ZFE(ITRA01(IPOIN,1))) THEN
                  ITRA02(IPOIN) = ITRA04(IPOIN,1) + 1
                  ITRA05(IPOIN) = -I1
               ELSE
                  ITRA02(IPOIN) = 1
                  ITRA05(IPOIN) = ITRA04(IPOIN,1) - 1
               ENDIF
               ZREF(IPOIN) = ZFE(ITRA01(IPOIN,ITRA04(IPOIN,1)))
C
            ELSEIF (ITRA02(IPOIN).EQ.3) THEN
C
               FLAG = .TRUE.
               IF (ZFE(ITRA01(IPOIN,ITRA04(IPOIN,2))).GT.
     *             ZFE(ITRA01(IPOIN,1))) THEN
                  ITRA02(IPOIN) = ITRA04(IPOIN,1) + 1
                  ITRA05(IPOIN) = -I1
               ELSE
                  ITRA02(IPOIN) = 1
                  ITRA05(IPOIN) = ITRA04(IPOIN,1) - 1
               ENDIF
               ZREF(IPOIN) = ZFE(ITRA01(IPOIN,ITRA04(IPOIN,1)))
C
            ELSEIF (ITRA02(IPOIN).EQ.-3) THEN
C
               FLAG = .TRUE.
               IF (ZFE(ITRA01(IPOIN,-I1)).GT.
     *             ZFE(ITRA01(IPOIN,ITRA04(IPOIN,1)))) THEN
                  ITRA02(IPOIN) = ITRA04(IPOIN,2) + 1
                  ITRA05(IPOIN) = -I1
               ELSE
                  ITRA02(IPOIN) = 1
                  ITRA05(IPOIN) = ITRA04(IPOIN,2) - 1
               ENDIF
               ZREF(IPOIN) = ZFE(ITRA01(IPOIN,ITRA04(IPOIN,2)))
C
            ELSEIF (ITRA02(IPOIN).EQ.4) THEN
C
               FLAG = .TRUE.
               IF (ZFE(ITRA01(IPOIN,ITRA04(IPOIN,3))).GT.
     *             ZFE(ITRA01(IPOIN,ITRA04(IPOIN,1)))) THEN
                  ITRA02(IPOIN) = ITRA04(IPOIN,2) + 1
                  ITRA05(IPOIN) = -I1
               ELSE
                  ITRA02(IPOIN) = 1
                  ITRA05(IPOIN) = ITRA04(IPOIN,2) - 1
               ENDIF
               ZREF(IPOIN) = ZFE(ITRA01(IPOIN,ITRA04(IPOIN,2)))
C
            ENDIF
C
         ELSE
C
C-----------------------------------------------------------------------
C
C    CORRECTION DE ZFE POUR LES POINTS INTERIEURS
C
C   SI ITRA02(IPOIN) VAUT 1, -1, 2, -2, 3, -3 : PAS DE CORRECTION
C   SI ITRA02(IPOIN) VAUT 4, -4, 5, -5 : CORRECTION
C
C-----------------------------------------------------------------------
C
            IF (ITRA02(IPOIN).EQ.4) THEN
C
               FLAG = .TRUE.
               IF (ZFE(ITRA01(IPOIN,ITRA04(IPOIN,3))).GT.
     *             ZFE(ITRA01(IPOIN,ITRA04(IPOIN,1)))) THEN
                  ITRA02(IPOIN) = ITRA04(IPOIN,2) + 1
                  ITRA05(IPOIN) = I1
               ELSE
                  ITRA02(IPOIN) = 2
                  ITRA05(IPOIN) = ITRA04(IPOIN,2) - 1
               ENDIF
               ZREF(IPOIN) = MIN(ZFE(ITRA01(IPOIN,ITRA04(IPOIN,2))),
     *                           ZFE(ITRA01(IPOIN,1)))
C
            ELSEIF (ITRA02(IPOIN).EQ.-4) THEN
C
               FLAG = .TRUE.
               IF (ZFE(ITRA01(IPOIN,ITRA04(IPOIN,2))).GT.
     *             ZFE(ITRA01(IPOIN,1))) THEN
                  ITRA02(IPOIN) = ITRA04(IPOIN,1) + 1
                  ITRA05(IPOIN) = ITRA04(IPOIN,3) - 1
               ELSE
                  ITRA02(IPOIN) = MOD(ITRA04(IPOIN,3),I1) + 1
                  ITRA05(IPOIN) = ITRA04(IPOIN,1) - 1
               ENDIF
               ZREF(IPOIN) = MIN(ZFE(ITRA01(IPOIN,ITRA04(IPOIN,1))),
     *                           ZFE(ITRA01(IPOIN,ITRA04(IPOIN,3))))
C
            ELSEIF (ITRA02(IPOIN).EQ.5) THEN
C
               FLAG = .TRUE.
               IF (ZFE(ITRA01(IPOIN,ITRA04(IPOIN,4))).GT.
     *             ZFE(ITRA01(IPOIN,ITRA04(IPOIN,2)))) THEN
                  ITRA02(IPOIN) = ITRA04(IPOIN,3) + 1
                  ITRA05(IPOIN) = ITRA04(IPOIN,1) - 1
               ELSE
                  ITRA02(IPOIN) = ITRA04(IPOIN,1) + 1
                  ITRA05(IPOIN) = ITRA04(IPOIN,3) - 1
               ENDIF
               ZREF(IPOIN) = MIN(ZFE(ITRA01(IPOIN,ITRA04(IPOIN,1))),
     *                           ZFE(ITRA01(IPOIN,ITRA04(IPOIN,3))))
C
            ELSEIF (ITRA02(IPOIN).EQ.-5) THEN
C
               FLAG = .TRUE.
               IF (ZFE(ITRA01(IPOIN,ITRA04(IPOIN,3))).GT.
     *             ZFE(ITRA01(IPOIN,ITRA04(IPOIN,1)))) THEN
                  ITRA02(IPOIN) = ITRA04(IPOIN,2) + 1
                  ITRA05(IPOIN) = ITRA04(IPOIN,4) - 1
               ELSE
                  ITRA02(IPOIN) = MOD(ITRA04(IPOIN,4),I1) + 1
                  ITRA05(IPOIN) = ITRA04(IPOIN,2) - 1
               ENDIF
               ZREF(IPOIN) = MIN(ZFE(ITRA01(IPOIN,ITRA04(IPOIN,2))),
     *                           ZFE(ITRA01(IPOIN,ITRA04(IPOIN,4))))
C
            ENDIF
C
         ENDIF
C
100   CONTINUE
C
      IF (FLAG) THEN
C
         DO 110 IPOIN = 1,NPOIN
C
            IF (ITRA05(IPOIN).NE.0) THEN
C
               IF (ITRA05(IPOIN).LT.ITRA02(IPOIN)) THEN
                  DO 120 I = ITRA02(IPOIN),ITRA03(IPOIN)
                     ZFE(ITRA01(IPOIN,I)) = MAX(ZFE(ITRA01(IPOIN,I)),
     *                                          ZREF(IPOIN))
120               CONTINUE
                  ITRA02(IPOIN) = 1
               ENDIF
C
               DO 130 I = ITRA02(IPOIN),ITRA05(IPOIN)
                  ZFE(ITRA01(IPOIN,I)) = MAX(ZFE(ITRA01(IPOIN,I)),
     *                                       ZREF(IPOIN))
130            CONTINUE
C
            ENDIF
C
110      CONTINUE
C
         GOTO 60
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      DEALLOCATE(IFAN)
      DEALLOCATE(ITRA01)
      DEALLOCATE(ITRA04)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
