C                       *****************
                        SUBROUTINE CPIK12
C                       *****************
C
     *(IKLE,NELEM,NELMAX,NPOIN)
C
C***********************************************************************
C BIEF VERSION 5.1           10/01/95    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FONCTION : EXTENSION DU TABLEAU DES CONNECTIVITES
C            CAS DE L'EXTENSION A UN ELEMENT QUASI-BULLE.
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
      INTEGER, INTENT(IN)    :: NELEM,NELMAX,NPOIN
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,4)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM
C
C-----------------------------------------------------------------------
C
      DO IELEM = 1 , NELEM
C
        IKLE(IELEM,4) = NPOIN + IELEM
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END     
 
