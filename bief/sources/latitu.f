C                       *****************
                        SUBROUTINE LATITU
C                       *****************
C
     *(COSLAT,SINLAT,LAMBD0,  Y,NPOIN)
C
C***********************************************************************
C BIEF VERSION 5.1           10/01/95    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C      FONCTION:    CALCUL DE TABLEAUX QUI DEPENDENT DE LA LATITUDE
C                   DU POINT CONSIDERE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    COSLAT      |<-- |  COS(LAMBDA)                                 |
C |    SINLAT      |<-- |  SIN(LAMBDA)                                 |
C |    LAMBD0      | -->|  LATITUDE DU POINT ORIGINE                   |
C |    Y           | -->|  COORDONNEES DES POINTS DU MAILLAGE          |
C |    NPOIN       | -->|  NOMBRE DE POINTS DANS LE MAILLAGE           |
C |    PRIVE       | -->|  TABLEAU PRIVE DE DIMENSION (NPOIN,NPRIV)    |
C |                |    |  NPRIV EST DONNE DANS LE PROGRAMME PRINCIPAL |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PREDON
C
C PROGRAMME APPELE : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(IN)    :: NPOIN
      DOUBLE PRECISION, INTENT(IN)    :: Y(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: COSLAT(NPOIN),SINLAT(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: LAMBD0
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C
      DOUBLE PRECISION LB2RAD,SURR,PISUR4,PISUR2,XLAMB
C
      INTRINSIC TAN,ATAN,SIN,COS,EXP
C
C-----------------------------------------------------------------------
C
C RAYON DE LA TERRE
C
      SURR = 1.D0 / 6400000.D0
C
C-----------------------------------------------------------------------
C
      PISUR4 = ATAN(1.D0)
      PISUR2 = PISUR4 + PISUR4
C
C  LAMBD0/2 EN RADIANS
C
      LB2RAD = LAMBD0 * PISUR4 / 90.D0
C
C  CONSTRUCTION DE 1/COS(LAMBDA),COS(LAMBDA),SIN(LAMBDA)
C
      DO 10 I = 1 , NPOIN
C
        XLAMB = 2.D0* ATAN(EXP(Y(I)*SURR)*TAN(LB2RAD+PISUR4))-PISUR2
        COSLAT(I) = COS(XLAMB)
        SINLAT(I) = SIN(XLAMB)
C
10    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
