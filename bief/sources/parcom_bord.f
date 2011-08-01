C                       **********************
                        SUBROUTINE PARCOM_BORD
C                       **********************
C
     *( X , ICOM , MESH )
C
C***********************************************************************
C BIEF VERSION 5.9      24/10/2008    J-M HERVOUET (LNHE) 01 30 87 80 18
C COPYRIGHT 2008              AFTER REINHARD HINKELMANN (HANNOVER UNI.)
C***********************************************************************
C
C   FONCTION : COMPLEMENT D'UN VECTEUR AUX INTERFACES ENTRE
C              SOUS-DOMAINES. ICI VECTEUR DE BORD DE TYPE 1
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C |      X         |<-->| VECTEUR DE BORD.
C |      ICOM      | -->| COMMUNICATION MODE
C |                |    | = 1 : VALUE WITH MAXIMUM ABSOLUTE VALUE TAKEN 
C |                |    | = 2 : CONTRIBUTIONS ADDED
C |                |    | = 3 : MAXIMUM CONTRIBUTION RETAINED
C |                |    | = 4 : MINIMUM CONTRIBUTION RETAINED
C |      MESH      | -->| MAILLAGE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   :  PLANTE
C
C***********************************************************************
C
      USE BIEF, EX_PARCOM_BORD => PARCOM_BORD
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C-----------------------------------------------------------------------
C
      INTEGER, INTENT(IN) :: ICOM
C
C     STRUCTURES : VECTEURS OU BLOCS
C
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
      DOUBLE PRECISION, INTENT(INOUT) :: X(*)
C
C-----------------------------------------------------------------------
C
      INTEGER NPTFR,I,TELM,TDIM1,TDIM2,TDIMDISC
C
C***********************************************************************
C
      NPTFR=NBPTS(1)
C
      TELM     = MESH%T%ELM
      TDIM1    = MESH%T%DIM1
      TDIM2    = MESH%T%DIM2
      TDIMDISC = MESH%T%DIMDISC
C
      MESH%T%ELM     = 11
      MESH%T%DIM1    = NBPTS(11)
      MESH%T%DIM2    = 1
      MESH%T%DIMDISC = 0
C
      CALL OS('X=0     ',X=MESH%T)
C
      DO I=1,NPTFR
        MESH%T%R(MESH%NBOR%I(I))=X(I)
      ENDDO
C
      CALL PARCOM(MESH%T,ICOM,MESH)
C
      DO I=1,NPTFR
        X(I)=MESH%T%R(MESH%NBOR%I(I))
      ENDDO
C
      MESH%T%ELM     = TELM
      MESH%T%DIM1    = TDIM1
      MESH%T%DIM2    = TDIM2
      MESH%T%DIMDISC = TDIMDISC
C
C-----------------------------------------------------------------------
C
      RETURN
      END
