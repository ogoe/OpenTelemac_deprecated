C                       *****************
                        SUBROUTINE PARCOM
C                       *****************
C
     *( X , ICOM , MESH )
C
C***********************************************************************
C BIEF VERSION 5.9        24/10/08    J-M HERVOUET (LNHE) 01 30 87 80 18
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
      USE BIEF, EX_PARCOM => PARCOM
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
      TYPE(BIEF_MESH), INTENT(INOUT)   :: MESH
      TYPE(BIEF_OBJ), INTENT(INOUT) :: X
C
C-----------------------------------------------------------------------
C
      TYPE(BIEF_OBJ), POINTER  :: X2,X3
      INTEGER NPOIN,NPLAN,IAN,NP11,NSEG
C
C***********************************************************************
C
C  INUTILE POUR UN SOUS-DOMAINE DECONNECTE DES AUTRES
C
      IF(NPTIR.EQ.0) RETURN
C
C-----------------------------------------------------------------------
C
      NPOIN = MESH%NPOIN
      NPLAN = 1
      IF(MESH%DIM.EQ.3) THEN
        NPOIN = NBPTS(11)
        NPLAN = MESH%NPOIN/NPOIN
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(X%TYPE.EQ.2) THEN
C
C     STRUCTURE DE VECTEUR
C
      IAN = 1
      CALL PARCOM2(X%R,X%R,X%R,NPOIN,NPLAN,ICOM,IAN,MESH)
C
      IF(X%ELM.EQ.13) THEN
        NP11=NBPTS(11)
        NSEG=MESH%NSEG
        CALL PARCOM2_SEG(X%R(NP11+1:NP11+NSEG),
     *                   X%R(NP11+1:NP11+NSEG),
     *                   X%R(NP11+1:NP11+NSEG),
     *                   NSEG,1,ICOM,IAN,MESH,1)
      ENDIF
C
      ELSEIF(X%TYPE.EQ.4) THEN
C
C     STRUCTURE DE BLOC
C
C     ATTENTION : NOMBRE LIMITE A 3 |||||||||||||||||||||||||
      IAN = X%N
      IF(IAN.EQ.1) THEN
        X2 => X%ADR(1)%P
        X3 => X%ADR(1)%P
      ELSEIF(IAN.EQ.2) THEN
        X2 => X%ADR(2)%P
        X3 => X%ADR(2)%P
      ELSEIF(IAN.EQ.3) THEN
        X2 => X%ADR(2)%P
        X3 => X%ADR(3)%P
      ELSE
        IF(LNG.EQ.1) WRITE(LU,*) 'PARCOM PREVU JUSQU''A 3 VECTEURS'
        IF(LNG.EQ.2) WRITE(LU,*) 'PARCOM: NO MORE THAN 3 VECTORS'
        CALL PLANTE(1)
        STOP
      ENDIF
C
      CALL PARCOM2(X%ADR(1)%P%R,X2%R,X3%R,NPOIN,NPLAN,ICOM,IAN,MESH)
C
C     PROVISIONNALY 1 BY 1, COULD BE OPTIMISED
C
      IF(X%ADR(1)%P%ELM.EQ.13) THEN
        NP11=NBPTS(11)
        NSEG=MESH%NSEG
        CALL PARCOM2_SEG(X%ADR(1)%P%R(NP11+1:NP11+NSEG),
     *                   X%ADR(1)%P%R(NP11+1:NP11+NSEG),
     *                   X%ADR(1)%P%R(NP11+1:NP11+NSEG),
C    *                   NSEG,1,ICOM,IAN,MESH)
     *                   NSEG,1,ICOM,1  ,MESH,1)
      ENDIF
      IF(IAN.GE.2.AND.X2%ELM.EQ.13) THEN
        NP11=NBPTS(11)
        NSEG=MESH%NSEG
        CALL PARCOM2_SEG(X2%R(NP11+1:NP11+NSEG),
     *                   X2%R(NP11+1:NP11+NSEG),
     *                   X2%R(NP11+1:NP11+NSEG),
C    *                   NSEG,1,ICOM,IAN,MESH)
     *                   NSEG,1,ICOM,1  ,MESH,1)
      ENDIF
      IF(IAN.EQ.3.AND.X3%ELM.EQ.13) THEN
        NP11=NBPTS(11)
        NSEG=MESH%NSEG
        CALL PARCOM2_SEG(X3%R(NP11+1:NP11+NSEG),
     *                   X3%R(NP11+1:NP11+NSEG),
     *                   X3%R(NP11+1:NP11+NSEG),
C    *                   NSEG,1,ICOM,IAN,MESH)
     *                   NSEG,1,ICOM,1  ,MESH,1)
      ENDIF
C
      ELSE
C
C     ERREUR SUR LA STRUCTURE
C
      IF (LNG.EQ.1) WRITE(LU,50) X%NAME,X%TYPE
      IF (LNG.EQ.1) WRITE(LU,53)
50    FORMAT(1X,'PARCOM (BIEF) : NOM DE X : ',A6,'  TYPE : ',1I6)
53    FORMAT(1X,'                CAS NON PREVU')
      IF (LNG.EQ.2) WRITE(LU,51) X%NAME,X%TYPE
      IF (LNG.EQ.2) WRITE(LU,54)
51    FORMAT(1X,'PARCOM (BIEF): NAME OF X: ',A6,'  TYPE : ',1I6)
54    FORMAT(1X,'               UNEXPECTED CASE')
      CALL PLANTE(1)
      STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
