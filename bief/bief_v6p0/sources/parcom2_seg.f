C                       **********************
                        SUBROUTINE PARCOM2_SEG
C                       **********************
C
     *( X1 , X2 , X3 , NSEG , NPLAN , ICOM , IAN , MESH , OPT )
C
C***********************************************************************
C BIEF VERSION 6.0        19/10/2009  J-M HERVOUET (LNHE) 01 30 87 80 18
C COPYRIGHT EDF 2009              
C***********************************************************************
C
C   FONCTION : COMPLEMENT D'UN VECTEUR DE SEGMENTS AUX INTERFACES ENTRE
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
C              IN 3D, THE FINITE VOLUME SEGMENTS IN PRISMS ARE
C              CONSIDERED HERE, I.E. HORIZONTAL FIRST : NSEG*NPLAN
C              THEN, VERTICAL : NPOIN2*(NPLAN-1)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C |      X1,X2,X3  |<-->| 3 VECTORS CAN BE TREATED, SEE IAN
C |      NSEG      | -->| NUMBER OF 2D SEGMENTS
C |      NPLAN     | -->| NUMBER OF PLANES
C |      ICOM      | -->| COMMUNICATION MODE
C |                |    | = 1 : VALUE WITH MAXIMUM ABSOLUTE VALUE TAKEN 
C |                |    | = 2 : CONTRIBUTIONS ADDED
C |                |    | = 3 : MAXIMUM CONTRIBUTION RETAINED
C |                |    | = 4 : MINIMUM CONTRIBUTION RETAINED
C |      IAN       | -->| NUMBER OF VECTORS TO BE TREATED (X1, X2, X3)
C |      MESH      | -->| MAILLAGE 2D
C |      OPT       | -->| 1 : HORIZONTAL AND VERTICAL SEGMENTS ONLY
C |                |    | 2 : ALL SEGMENTS
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   :  PLANTE
C
C***********************************************************************
C
      USE BIEF, EX_PARCOM2_SEG => PARCOM2_SEG
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: ICOM,NSEG,NPLAN,IAN,OPT
C
C     STRUCTURES : VECTEURS OU BLOCS
C
      TYPE(BIEF_MESH) , INTENT(INOUT) :: MESH
C
C     IN 2D X1(NSEG)
C     IN 3D X1(NSEG*NPLAN+NPOIN2*(NPLAN-1))
C
      DOUBLE PRECISION, INTENT(INOUT) :: X1(*),X2(*),X3(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NPOIN2,IDEB,IFIN,IPLAN
C
C-----------------------------------------------------------------------
C
      IF(NPLAN.GT.1) THEN
C
C       1) HORIZONTAL FLUXES
C
        DO IPLAN=1,NPLAN
          IDEB=1+NSEG*(IPLAN-1)
          IFIN=  NSEG* IPLAN
          CALL PARACO(X1(IDEB:IFIN),X2(IDEB:IFIN),X3(IDEB:IFIN),
     *                NSEG,ICOM,IAN,1,
     *                MESH%NB_NEIGHB_SEG,MESH%NB_NEIGHB_PT_SEG%I,
     *                MESH%LIST_SEND_SEG%I,MESH%NH_COM_SEG%I,
     *                MESH%NH_COM_SEG%DIM1,MESH%BUF_SEND%R,
     *                MESH%BUF_RECV%R,MESH%BUF_SEND%DIM1)
        ENDDO
C
C       2) VERTICAL FLUXES
C
        NPOIN2=MESH%NPOIN
        DO IPLAN=1,NPLAN-1
          IDEB=NSEG*NPLAN + NPOIN2*(IPLAN-1) + 1
          IFIN=NSEG*NPLAN + NPOIN2* IPLAN
          CALL PARCOM2(X1(IDEB:IFIN),X2(IDEB:IFIN),X3(IDEB:IFIN),
     *                 NPOIN2,1,ICOM,IAN,MESH)
        ENDDO
C
C       3) CROSSED FLUXES
C
        IF(OPT.EQ.2) THEN
          DO IPLAN=1,NPLAN-1
            IDEB=NSEG*NPLAN + NPOIN2*(NPLAN-1) + 2*NSEG*(IPLAN-1) + 1
            IFIN=NSEG*NPLAN + NPOIN2*(NPLAN-1) + 2*NSEG* IPLAN
            CALL PARACO(X1(IDEB:IFIN),X2(IDEB:IFIN),X3(IDEB:IFIN),
     *                  NSEG,ICOM,IAN,2,
     *                  MESH%NB_NEIGHB_SEG,MESH%NB_NEIGHB_PT_SEG%I,
     *                  MESH%LIST_SEND_SEG%I,MESH%NH_COM_SEG%I,
     *                  MESH%NH_COM_SEG%DIM1,MESH%BUF_SEND%R,
     *                  MESH%BUF_RECV%R,MESH%BUF_SEND%DIM1)
          ENDDO
        ENDIF
C
      ELSE
C
        CALL PARACO(X1,X2,X3,
     *              NSEG,ICOM,IAN,1,
     *              MESH%NB_NEIGHB_SEG,MESH%NB_NEIGHB_PT_SEG%I,
     *              MESH%LIST_SEND_SEG%I,MESH%NH_COM_SEG%I,
     *              MESH%NH_COM_SEG%DIM1,MESH%BUF_SEND%R,
     *              MESH%BUF_RECV%R,MESH%BUF_SEND%DIM1)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
