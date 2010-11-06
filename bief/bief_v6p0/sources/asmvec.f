C                       *****************
                        SUBROUTINE ASMVEC
C                       *****************
C
     *(X, IKLE,NPOIN,NELEM,NELMAX,NDP,W,INIT,LV)
C
C***********************************************************************
C BIEF VERSION 5.1           17/08/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : ASSEMBLAGE MULTIPLICATIF DE VECTEUR.
C
C            ATTENTION, CE VECTEUR EST INITIALISE A 1 SI INIT = .TRUE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  X             |<-->| VECTEUR ASSEMBLE
C |  IKLE          | -->| CORRESPONDANCES NUMEROTATION LOCALE-GLOBALE
C |  NPOIN         | -->| DIMENSION DU TABLEAU X
C |  NELEM         | -->| NOMBRE D'ELEMENTS DANS LE MAILLAGE.
C |  NELMAX        | -->| PREMIERE DIMENSION DE IKLE ET W.
C |                |    | (CAS D'UN MAILLAGE ADAPTATIF)
C |  NDP           | -->| DEUXIEME DIMENSION DE IKLE.
C |  W             | -->| TABLEAUX DE TRAVAIL CONTENANT LE VECTEUR SOUS
C |                |    | FORME NON ASSEMBLEE
C |                |    | W EST DE DIMENSION NELMAX * NDP(IELM)
C |                |    | NDP EST LE NOMBRE DE POINTS DE L'ELEMENT
C |  INIT          | -->| LOGIQUE : SI VRAI : X EST INITIALISE A 0
C |   LV           | -->|  LONGUEUR DU VECTEUR POUR LA VECTORISATION
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE : ASMVE1,OV
C
C***********************************************************************
C
      USE BIEF, EX_ASMVEC => ASMVEC
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(IN)    :: NELMAX,NPOIN,NELEM,NDP,LV
      DOUBLE PRECISION, INTENT(INOUT) :: X(NPOIN)
      INTEGER         , INTENT(IN)    :: IKLE(NELMAX,NDP)
      DOUBLE PRECISION, INTENT(IN)    :: W(NELMAX,NDP)
      LOGICAL         , INTENT(IN)    :: INIT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IDP
C
C-----------------------------------------------------------------------
C   MISE A 1. EVENTUELLE DU VECTEUR X
C-----------------------------------------------------------------------
C
      IF(INIT) CALL OV( 'X=C     ' , X , X , X , 1.D0 , NPOIN )
C
C-----------------------------------------------------------------------
C   ASSEMBLAGE, CONTRIBUTION DES POINTS LOCAUX 1,... JUSQU'A NDP
C-----------------------------------------------------------------------
C
      DO IDP = 1 , NDP
C
        CALL ASMVE1(X, IKLE(1,IDP),W(1,IDP),NPOIN,NELEM,NELMAX,LV)
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
 
