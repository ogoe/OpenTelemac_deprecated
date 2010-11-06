C                       ******************
                        SUBROUTINE PARCOM2
C                       ******************
C
     *( X1 , X2 , X3 , NPOIN , NPLAN , ICOM , IAN , MESH )
C
C***********************************************************************
C BIEF VERSION 5.9        13/08/08    J-M HERVOUET (LNHE) 01 30 87 80 18
C COPYRIGHT 2008              AFTER REINHARD HINKELMANN (HANNOVER UNI.)
C***********************************************************************
C
C   FONCTION : COMPLEMENT D'UN VECTEUR AUX INTERFACES ENTRE
C              SOUS-DOMAINES.
C
C              X PEUT ETRE UN BLOC DE VECTEURS, EN CE CAS, TOUS LES
C              VECTEURS DU BLOC SONT TRAITES.
C
C              ATTENTION ||||
C
C              SI LES VECTEURS ONT UNE DEUXIEME DIMENSION
C              ELLE EST POUR L'INSTANT IGNOREE
C
C   IMPORTANT NOTICE :
C
C   FROM VERSION 5.9 ON, IDENTICAL VALUES ARE ENSURED AT INTERFACE POINTS
C   SO THAT DIFFERENT PROCESSORS WILL ALWAYS MAKE THE SAME DECISION IN
C   TESTS ON REAL NUMBERS
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C |      X         |<-->| VECTEUR OU BLOC DE VECTEURS.
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
      USE BIEF, EX_PARCOM2 => PARCOM2
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C-----------------------------------------------------------------------
C
      INTEGER, INTENT(IN) :: ICOM,NPOIN,NPLAN,IAN
C
C     STRUCTURES : VECTEURS OU BLOCS
C
      TYPE(BIEF_MESH) , INTENT(INOUT) :: MESH
      DOUBLE PRECISION, INTENT(INOUT) :: X1(NPOIN,NPLAN)
      DOUBLE PRECISION, INTENT(INOUT) :: X2(NPOIN,NPLAN)
      DOUBLE PRECISION, INTENT(INOUT) :: X3(NPOIN,NPLAN)
C
C-----------------------------------------------------------------------
C
      CALL PARACO(X1,X2,X3,NPOIN,ICOM,IAN,NPLAN,
     *            MESH%NB_NEIGHB,MESH%NB_NEIGHB_PT%I,MESH%LIST_SEND%I,
     *            MESH%NH_COM%I,MESH%NH_COM%DIM1,MESH%BUF_SEND%R,
     *            MESH%BUF_RECV%R,MESH%BUF_SEND%DIM1)
C
C     VERSION 5.9 : ENSURING SAME VALUES AT INTERFACE POINTS
C                   SHARED BY SEVERAL SUB-DOMAINS                        
C
      IF(ICOM.EQ.2.AND.NCSIZE.GT.2) THEN
C
      CALL PARACO(X1,X2,X3,NPOIN,1,IAN,NPLAN,
     *            MESH%NB_NEIGHB,MESH%NB_NEIGHB_PT%I,MESH%LIST_SEND%I,
     *            MESH%NH_COM%I,MESH%NH_COM%DIM1,MESH%BUF_SEND%R,
     *            MESH%BUF_RECV%R,MESH%BUF_SEND%DIM1)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
