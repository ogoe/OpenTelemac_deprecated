!                       *****************
                        SUBROUTINE DECOUP
!                       *****************
!
     &(ISURC,X,Y,IKLE,NCOLOR,IFABOR, NELEM2,NPOIN2,COLOR)
!
!***********************************************************************
! PROGICIEL: STBTEL V5.2        19/04/91  J-C GALLAND  (LNH)
!                               19/02/93  J-M JANIN    (LNH)
!***********************************************************************
!
! FONCTION : DECOUPAGE DES TRIANGLES SURCONTRAINTS :
!            ILS SONT COUPES EN TROIS PAR AJOUT D'UN POINT A
!            LEUR BARYCENTRE
!
!
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________.
! |      NOM       |MODE|                   ROLE                       |
! |________________|____|______________________________________________|
! |   ISURC        | -->| NUMERO DE L'ELEMENT SURCONTRAINT A TRAITER
! |   X,Y          |<-->| COORDONNEES DU MAILLAGE .
! |   IKLE         |<-->| LISTE DES POINTS DE CHAQUE ELEMENT
! |   NCOLOR       |<-->| TABLEAU DES COULEURS DES POINTS
! |   IFABOR       | -->| TABLEAU DES VOISINS DES ELEMENTS
! |   NELEM2       |<-->| NOUVEAU NOMBRE D'ELEMENTS APRES DECOUP
! |   NPOIN2       |<-->| NOUVEAU NOMBRE DE POINTS APRES DECOUP
! |    COLOR       |<-->| STOCKAGE COULEURS DES NOEUDS SUR FICHIER GEO
! |________________|____|______________________________________________
! | COMMON:        |    |
! |  GEO:          |    |
! |    MESH        | -->| TYPE DES ELEMENTS DU MAILLAGE
! |    NDP         | -->| NOMBRE DE NOEUDS PAR ELEMENTS
! |    NPOIN       | -->| NOMBRE TOTAL DE NOEUDS DU MAILLAGE
! |    NELEM       | -->| NOMBRE TOTAL D'ELEMENTS DU MAILLAGE
! |    NPMAX       | -->| DIMENSION EFFECTIVE DES TABLEAUX X ET Y
! |                |    | (NPMAX = NPOIN + 0.1*NELEM)
! |    NELMAX      | -->| DIMENSION EFFECTIVE DES TABLEAUX CONCERNANT
! |                |    | LES ELEMENTS (NELMAX = NELEM + 0.2*NELEM)
! |                |    |
! |________________|____|______________________________________________|
! MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
!-----------------------------------------------------------------------
! APPELE PAR : SURCON
! APPEL DE : -
!***********************************************************************
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
      DOUBLE PRECISION X(*) , Y(*)
!
      INTEGER NELEM , NELEM2 , MESH , NDP , NPOIN , NPOIN2 , NPMAX
      INTEGER NELMAX , KELEM
      INTEGER ISURC , IFAC , ICOLOR , I , I1 , I2 , I3
      INTEGER IKLE(NELMAX,4) , NCOLOR(*)
      INTEGER IFABOR(NELMAX,*)
!
      LOGICAL COLOR
!
      COMMON/GEO/ MESH , NDP , NPOIN , NELEM , NPMAX , NELMAX
!
!=======================================================================
! CALCUL DES COORDONNEES DU NOUVEAU NOEUD 4 (DE NUMERO NPOIN2)
!=======================================================================
!
      NPOIN2 = NPOIN2 + 1
      I1 = IKLE(ISURC,1)
      I2 = IKLE(ISURC,2)
      I3 = IKLE(ISURC,3)
!
      X(NPOIN2) = (X(I1) + X(I2) + X(I3))/3.D0
      Y(NPOIN2) = (Y(I1) + Y(I2) + Y(I3))/3.D0
!
!=======================================================================
! DEFINITION DE LA COULEUR DU NOEUD CREE (C'EST CELLE DU NOEUD NON POINT
! DE BORD DE L'ELEMENT VOISIN)
!=======================================================================
!
      IF (COLOR) THEN
         DO 10 IFAC=1,3
            IF(IFABOR(ISURC,IFAC).GT.0) ICOLOR = IFABOR(ISURC,IFAC)
 10      CONTINUE
!
         DO 20 I=1,3
            IF(IKLE(ICOLOR,I).NE.I1.AND.IKLE(ICOLOR,I).NE.I2.AND.
     &         IKLE(ICOLOR,I).NE.I3)
     &         NCOLOR(NPOIN2) = NCOLOR(IKLE(ICOLOR,I))
 20      CONTINUE
       ENDIF
!
!=======================================================================
! CALCUL DES NOUVEAUX IKLE : L'ELEMENT (1,2,4) CONSERVE LE NUMERO ISURC
!                            L'ELEMENT (2,3,4) PREND LE NUMERO NELEM2+1
!                            L'ELEMENT (3,1,4) PREND LE NUMERO NELEM2+2
!=======================================================================
!
      IKLE(ISURC,3) = NPOIN2
!
      NELEM2 = NELEM2 + 1
      IKLE(NELEM2,1) = I2
      IKLE(NELEM2,2) = I3
      IKLE(NELEM2,3) = NPOIN2
!
      KELEM = IFABOR(ISURC,2)
      IFABOR(NELEM2,1) = KELEM
      IFABOR(NELEM2,2) = NELEM2+1
      IFABOR(NELEM2,3) = ISURC
      IF (KELEM.GT.0) THEN
         IF (IFABOR(KELEM,1).EQ.ISURC) IFABOR(KELEM,1) = NELEM2
         IF (IFABOR(KELEM,2).EQ.ISURC) IFABOR(KELEM,2) = NELEM2
         IF (IFABOR(KELEM,3).EQ.ISURC) IFABOR(KELEM,3) = NELEM2
      ENDIF
      IFABOR(ISURC,2) = NELEM2
!
      NELEM2 = NELEM2 + 1
      IKLE(NELEM2,1) = I3
      IKLE(NELEM2,2) = I1
      IKLE(NELEM2,3) = NPOIN2
!
      KELEM = IFABOR(ISURC,3)
      IFABOR(NELEM2,1) = IFABOR(ISURC,3)
      IFABOR(NELEM2,2) = ISURC
      IFABOR(NELEM2,3) = NELEM2-1
      IF (KELEM.GT.0) THEN
         IF (IFABOR(KELEM,1).EQ.ISURC) IFABOR(KELEM,1) = NELEM2
         IF (IFABOR(KELEM,2).EQ.ISURC) IFABOR(KELEM,2) = NELEM2
         IF (IFABOR(KELEM,3).EQ.ISURC) IFABOR(KELEM,3) = NELEM2
      ENDIF
      IFABOR(ISURC,3) = NELEM2
!
!=======================================================================
!
      RETURN
      END