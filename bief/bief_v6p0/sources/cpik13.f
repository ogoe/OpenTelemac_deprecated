C                       *****************
                        SUBROUTINE CPIK13
C                       *****************
C
     *(IKLE,IKLBOR,ELTSEG,NBOR,NELEM,NELMAX,NPOIN,NPTFR)
C
C***********************************************************************
C BIEF VERSION 5.9         20/03/08    J-M HERVOUET (LNH) 01 30 87 80 18
C***********************************************************************
C
C FONCTION : EXTENSION DU TABLEAU DES CONNECTIVITES
C            CAS DE L'EXTENSION A UN ELEMENT QUADRATIQUE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      IKLE      |<-->|  TABLEAU DES CONNECTIVITES                   |
C |      NELEM     | -->|  NOMBRE D'ELEMENTS
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS
C |      NPOIN     | -->|  NOMBRE DE SOMMETS DU MAILLAGE
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : INBIEF
C
C  SOUS-PROGRAMME APPELE : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NELEM,NELMAX,NPOIN,NPTFR
      INTEGER, INTENT(IN)    :: ELTSEG(NELMAX,3)
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,6),IKLBOR(NPTFR,*),NBOR(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,K
C
C-----------------------------------------------------------------------
C
C     CONNECTIVITY TABLE OF QUADRATIC GLOBAL POINTS
C
      DO IELEM = 1 , NELEM
C
C       NUMERO=NPOIN+NUMERO DU SEGMENT QUI CONTIENT LE POINT
C
        IKLE(IELEM,4) = NPOIN + ELTSEG(IELEM,1)
        IKLE(IELEM,5) = NPOIN + ELTSEG(IELEM,2)
        IKLE(IELEM,6) = NPOIN + ELTSEG(IELEM,3)
C
      ENDDO
C
C-----------------------------------------------------------------------
C
C     CONNECTIVITY TABLE OF QUADRATIC BOUNDARY POINTS
C     GLOBAL NUMBERS OF BOUNDARY QUADRATIC POINTS
C
      DO K=1,NPTFR
        IKLBOR(K,3)=K+NPTFR
C       SEGMENTS 1 TO NPTFR ARE THE BOUNDARY SEGMENTS
        NBOR(IKLBOR(K,3))=NPOIN+K
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
