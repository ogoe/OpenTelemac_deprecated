C                       *****************
                        SUBROUTINE MASKBD
C                       *****************
C
     *(MASKEL,ZFE,ZF,HN,HMIN,IKLE,IFABOR,ITRA01,NELEM,NPOIN)
C
C***********************************************************************
C  BIEF VERSION 5.1            11/08/94      J-M JANIN (LNH) 30 87 72 84
C
C***********************************************************************
C
C     FONCTION :
C
C    MASQUAGE DES ELEMENTS DECOUVERTS OU PARTIELLEMENT DECOUVERTS
C
C
C     ALGORITHME :
C
C    - UN ELEMENT EST MASQUE SI SA COTE DU FOND ZFE EST SUPERIEURE A LA
C      COTE DE SURFACE LIBRE DE SON CENTRE DE GRAVITE
C
C    - TOUT ELEMENT DONT LA COTE ZFE EST SUPERIEURE A LA COTE ZFE D'UN
C      VOISIN MASQUE EST MASQUE
C
C      LORSQU'ON TOURNE AUTOUR D'UN NOEUD LA FONCTION ZFE NE COMPORTE
C      QU'UN MINIMUM ET UN MAXIMUM (VOIR TOPOGR). LA SECONDE PHASE DE
C      L'ALGORITHME ASSURE DONC QU'IL N'Y A PAS 2 PARTIES DU DOMAINE
C      CONNECTEES UNIQUEMENT PAR UN SOMMET. CE TRAITEMENT EVITE
C      EGALEMENT DES MASQUAGES-DEMASQUAGES INTEMPESTIFS EN PARTICULIER
C      DANS LES ZONES PLATES DECOUVRANTES.
C
C     INCONVENIENTS :
C
C      CET ALGORITHME SUPPOSE LA SURFACE LIBRE QUASI HORIZONTALE.
C      IL EST BIEN ADAPTE POUR TRAITER DES EVOLUTIONS LIEES A LA MAREE
C      MAIS PAS POUR TRAITER UNE RUPTURE DE BARRAGE.
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   MASKEL       |<-- | TABLEAU DE MASQUAGE DES ELEMENTS             |
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE.        |
C |   ZFE          | -->| COTE DU FOND PAR ELEMENT.                    |
C |   ZF           | -->| COTE DU FOND PAR NOEUD.                      |
C |   IKLE         | -->| TABLE DE CONNECTIVITE.                       |
C |   IFABOR       | -->| NUMEROS DES ELEMENTS VOISINS.                |
C |   ITRA01       | -->| TABLEAU DE TRAVAIL D'ENTIERS.                |
C |   NELEM        | -->| NOMBRE D'ELEMENTS.                           |
C |   NPOIN        | -->| NOMBRE DE NOEUDS.                            |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C                                                                      *
C APPELE PAR: TELMAC ET MITRID
C                                                                      *
C***********************************************************************
C
      IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NELEM,NPOIN
      INTEGER, INTENT(IN)             :: IKLE(NELEM,3),IFABOR(NELEM,3)
      INTEGER, INTENT(INOUT)          :: ITRA01(NELEM)
      DOUBLE PRECISION, INTENT(IN)    :: ZFE(NELEM),ZF(NPOIN),HN(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: HMIN
      DOUBLE PRECISION, INTENT(INOUT) :: MASKEL(NELEM)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,I1,I2,I3,N      
C
      DOUBLE PRECISION EPSILO,ZSE
C
      LOGICAL FLAG
C
      DATA EPSILO / 1.D-6 /
C
C-----------------------------------------------------------------------
C
      FLAG = .FALSE.
C
      DO 10 IELEM = 1,NELEM
         I1 = IKLE(IELEM,1)
         I2 = IKLE(IELEM,2)
         I3 = IKLE(IELEM,3)
         ZSE = (ZF(I1)+HN(I1)+ZF(I2)+HN(I2)+ZF(I3)+HN(I3))/3.D0
         IF (ZFE(IELEM)+HMIN+EPSILO.GT.ZSE) THEN
            FLAG = .TRUE.
            MASKEL(IELEM) = 0.D0
         ENDIF
10    CONTINUE
C
20    CONTINUE
C
      IF (FLAG) THEN
C
         FLAG = .FALSE.
         DO 30 IELEM = 1,NELEM
C
            ITRA01(IELEM) = 0
            IF (MASKEL(IELEM).GT.0.5D0) THEN
C
               N=IFABOR(IELEM,1)
               IF (N.GT.0) THEN
                  IF (MASKEL(N).LT.0.5D0.AND.ZFE(IELEM).GT.
     *                ZFE(N)-EPSILO) ITRA01(IELEM) = 1
               ENDIF
               N=IFABOR(IELEM,2)
               IF (N.GT.0) THEN
                  IF (MASKEL(N).LT.0.5D0.AND.ZFE(IELEM).GT.
     *                ZFE(N)-EPSILO) ITRA01(IELEM) = 1
               ENDIF
               N=IFABOR(IELEM,3)
               IF (N.GT.0) THEN
                  IF (MASKEL(N).LT.0.5D0.AND.ZFE(IELEM).GT.
     *                ZFE(N)-EPSILO) ITRA01(IELEM) = 1
               ENDIF
C
            ENDIF
C
30       CONTINUE
C
         DO 40 IELEM = 1,NELEM
            IF (ITRA01(IELEM).EQ.1) THEN
               FLAG = .TRUE.
               MASKEL(IELEM) = 0.D0
            ENDIF
40       CONTINUE
C
         GOTO 20
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
