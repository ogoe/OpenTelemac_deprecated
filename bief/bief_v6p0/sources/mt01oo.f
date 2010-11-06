C                       *****************
                        SUBROUTINE MT01OO
C                       *****************
C
     *(A11,A12,A22,XMUL,LGSEG,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.1           06/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FONCTION :
C
C    CE SOUS-PROGRAMME CALCULE LES COEFFICIENTS DE LA MATRICE SUIVANTE:
C
C                              /
C                    A    =   /  (P *P )*J(X,Y) DX
C                     I J    /L    I  J
C
C     PAR MAILLE ELEMENTAIRE .
C
C     J(X,Y) : JACOBIEN DE LA TRANSFORMATION ISOPARAMETRIQUE
C
C     L'ELEMENT EST LE SEGMENT P1
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      A11,..,A22| -->|  ELEMENTS DE LA MATRICE.
C |      XMUL      | -->|  FACTEUR MULTIPLICATIF
C |      LGSEG     | -->|  LONGUEUR DES SEGMENTS.
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF, EX_MT01OO => MT01OO
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
C
      DOUBLE PRECISION, INTENT(INOUT) :: A11(NELMAX),A12(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: A22(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: LGSEG(NELMAX)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM
      DOUBLE PRECISION XSUR6,DET1
C
C-----------------------------------------------------------------------
C
      XSUR6  = XMUL/6.D0
C
C-----------------------------------------------------------------------
C
      DO 1 IELEM = 1 , NELEM
C
        DET1 = LGSEG(IELEM) * XSUR6
        A11(IELEM) = DET1 + DET1
        A12(IELEM) = DET1
        A22(IELEM) = DET1 + DET1
C
1     CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
