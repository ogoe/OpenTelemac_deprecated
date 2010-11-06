C                       *****************
                        SUBROUTINE COSAKE
C                       *****************
C
     *(KARMAN,CMU,C1,C2,SIGMAK,SIGMAE,ESTAR,SCHMIT,KMIN,KMAX,EMIN,EMAX)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2    17/08/94    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C     FONCTION  : FIXE LES CONSTANTES DU MODELE K-EPSILON
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    KARMAN      |<-- | CONSTANTE DE KARMAN                          |
C |    CMU         |<-- | CONSTANTE DU MODELE K-EPSILON                |
C |    C1          |<-- | CONSTANTE DU MODELE K-EPSILON                |
C |    C2          |<-- | CONSTANTE DU MODELE K-EPSILON                |
C |    SIGMAK      |<-- | CONSTANTE DU MODELE K-EPSILON                |
C |    SIGMAE      |<-- | CONSTANTE DU MODELE K-EPSILON                |
C |    ESTAR       |<-- | CONSTANTE DU MODELE K-EPSILON                |
C |    SCHMITT     |<-- | NOMBRE DE SCHMITT                            |
C |    KMIN        |<-- | K MINIMUM EN CAS DE CLIPPING                 |
C |    KMAX        |<-- | K MAXIMUM EN CAS DE CLIPPING                 |
C |    EMIN        |<-- | EPSILON MINIMUM EN CAS DE CLIPPING           |
C |    EMAX        |<-- | EPSILON MAXIMUM EN CAS DE CLIPPING           |
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
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(OUT) :: KMIN,KMAX,EMIN,EMAX
      DOUBLE PRECISION, INTENT(OUT) :: KARMAN,CMU,C1,C2
      DOUBLE PRECISION, INTENT(OUT) :: SIGMAK,SIGMAE,ESTAR,SCHMIT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C CONSTANTE DE KARMAN
C
      KARMAN = 0.41D0
      CMU    = 0.09D0
      C1     = 1.44D0
      C2     = 1.92D0
      SIGMAK = 1.00D0
      SIGMAE = 1.30D0
      ESTAR  = 0.15D0
C
C NOMBRE DE SCHMIDT
C
      SCHMIT = 0.50D0
C
C VALEURS EXTREMES POUR LE CLIPPING DE K ET EPSILON
C
      KMIN = 1.D-8
      EMIN = 1.D-8
      KMAX = 1.D10
      EMAX = 1.D10
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
