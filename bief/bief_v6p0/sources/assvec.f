C                       *****************
                        SUBROUTINE ASSVEC
C                       *****************
C
     *(X, IKLE,NPOIN,NELEM,NELMAX,IELM,W,INIT,LV,MSK,MASKEL)
C
C***********************************************************************
C BIEF VERSION 5.9         29/02/08    J-M HERVOUET (LNH) 01 30 87 80 18
C                                       
C***********************************************************************
C
C FONCTION : ASSEMBLAGE DE VECTEUR.
C
C            ATTENTION, CE VECTEUR N'EST INITIALISE A 0 QUE SI
C            INIT = .TRUE.
C
C-----------------------------------------------------------------------
C
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
C |  IELM          | -->| TYPE D'ELEMENT (VOIR CI-DESSUS)
C |  W             | -->| TABLEAUX DE TRAVAIL CONTENANT LE VECTEUR SOUS
C |                |    | FORME NON ASSEMBLEE
C |                |    | W EST DE DIMENSION NELMAX * NDP(IELM)
C |                |    | NDP EST LE NOMBRE DE POINTS DE L'ELEMENT
C |  INIT          | -->| LOGIQUE : SI VRAI : X EST INITIALISE A 0
C |  LV            | -->|  LONGUEUR DU VECTEUR POUR LA VECTORISATION
C |  MSK           | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |  MASKEL        | -->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE : ASSVE1
C
C***********************************************************************
C
      USE BIEF, EX_ASSVEC => ASSVEC
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(INOUT) :: X(*)
      INTEGER         , INTENT(IN)    :: NELEM,NELMAX,NPOIN,IELM,LV
      INTEGER         , INTENT(IN)    :: IKLE(NELMAX,*)
      DOUBLE PRECISION, INTENT(IN)    :: W(NELMAX,*),MASKEL(NELMAX)
      LOGICAL         , INTENT(IN)    :: INIT,MSK
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IDP
C
C-----------------------------------------------------------------------
C   MISE A ZERO EVENTUELLE DU VECTEUR X
C-----------------------------------------------------------------------
C
      IF(INIT) CALL OV( 'X=C     ' , X , X , X , 0.D0 , NPOIN )
C
C-----------------------------------------------------------------------
C   ASSEMBLAGE, CONTRIBUTION DES POINTS LOCAUX 1,... JUSQU'A NDP
C-----------------------------------------------------------------------
C  
      DO IDP = 1 , NBPEL(IELM)
C
        CALL ASSVE1(X,IKLE(1,IDP),W(1,IDP),NELEM,NELMAX,LV,MSK,MASKEL)
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
