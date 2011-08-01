C                       *****************
                        SUBROUTINE DIFFCL
C                       *****************
C
     *(LITBOR,TTILD,TBOR,NBOR,ICONV,NPOIN,NPTFR)
C
C***********************************************************************
C  BIEF VERSION 6.0     09/10/2009    J-M HERVOUET (LNHE) 01 30 87 80 18
C                                          
C  DEPLACE DE TELEMAC-2D POUR APPEL PAR SISYPHE
C***********************************************************************
C
C      FONCTION:
C      ========:
C
C      CE SOUS-PROGRAMME INITIALISE LA VALEUR DU TRACEUR POUR LES
C      CONDITIONS AUX LIMITES DE TYPE DIRICHLET DANS L'ETAPE DE
C      DIFFUSION.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   LITBOR       | -->| TYPES DE CONDITIONS AUX LIMITES DU TRACEUR.  |
C |   TTILD        | -->| TRACEUR.                                     |
C |   TBOR         |<-->| CONDITIONS AUX LIMITES SUR T.                |
C |   NBOR         | -->| ADRESSES DES POINTS FRONTIERES.
C |                |    | CONDITIONS AUX LIMITES (PHYSIQUE) .          |
C |   NPOIN        | -->| NOMBRE DE POINTS TOTAL .                     |
C |   NPTFR        | -->| NOMBRE DE POINTS FRONTIERES .                |
C |   KSORT        | -->| INDICATEUR DE POINT DE SORTIE FLUIDE .       |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE|NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMME APPELE : NEANT
C
C**********************************************************************
C
      USE DECLARATIONS_TELEMAC
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NPOIN,NPTFR,ICONV
      INTEGER, INTENT(IN)             :: NBOR(NPTFR)
      INTEGER, INTENT(IN)             :: LITBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: TTILD(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: TBOR(NPTFR)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K
C
C----------------------------------------------------------------------
C
C  INITIALISATION DE TBOR (VOIR LA CONSTRUCTION DE LIMTRA DANS DIFFIN)
C
      IF(ICONV.EQ.ADV_CAR) THEN
C
      DO 1 K=1,NPTFR
C
C  AUX SORTIES LIBRES AVEC LA METHODE DES CARACTERISTIQUES
C                     ON IMPOSE LE RESULTAT DE LA CONVECTION
C
        IF(LITBOR(K).EQ.KSORT) TBOR(K) = TTILD(NBOR(K))
C
1     CONTINUE
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
