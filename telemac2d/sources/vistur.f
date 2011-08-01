C                       *****************
                        SUBROUTINE VISTUR
C                       *****************
C
     *(VISC,AK,EP,NPOIN,CMU,PROPNU)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2    17/08/94    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C     FONCTION  : CALCUL DE LA VISCOSITE TURBULENTE EN FONCTION DE K ET
C                 EPSILON.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     VISC       |<-- | DIFFUSION TURBULENTE                         |
C |     AK         | -->| ENERGIE TURBULENTE                           |
C |     EP         | -->| DISSIPATION TURBULENTE                       |
C |     NPOIN      | -->| NOMBRE DE POINTS DANS LE MAILLAGE            |
C |     CMU        | -->| CONSTANTE DU MODELE K-EPSILON                |
C |     PROPNU     | -->| VISCOSITE LAMINAIRE                          |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)           :: NPOIN
      DOUBLE PRECISION, INTENT(IN)  :: CMU,PROPNU
      TYPE(BIEF_OBJ), INTENT(IN)    :: AK,EP
      TYPE(BIEF_OBJ), INTENT(INOUT) :: VISC
C   
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I  
C
C-----------------------------------------------------------------------
C
      DO I = 1 , NPOIN
C
        VISC%R(I) = PROPNU + CMU * AK%R(I)**2 / EP%R(I)
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
