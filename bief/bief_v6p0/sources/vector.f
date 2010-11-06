C                       *****************
                        SUBROUTINE VECTOR
C                       *****************
C
     *(VEC,OP,FORMUL,IELM1,XMUL,F,G,H,U,V,W,MESH,MSK,MASKEL)
C
C***********************************************************************
C BIEF VERSION 5.9        25/06/2008     JM HERVOUET (LNHE) 01 30 87 80 18
C                         22/09/05     REGINA NEBAUER
C
C 25/06/2008 JMH : PAS D'APPEL DE VECTOS SI NOMBRE D'ELEMENTS NUL
C
C***********************************************************************
C
C  FONCTION : CALCULS DE VECTEURS
C
C             LE VECTEUR EST IDENTIFIE PAR LA FORMULE CONTENUE DANS
C             LA CHAINE DE CARACTERES FORMUL
C
C             OP VAUT = OU +
C
C-----------------------------------------------------------------------
C
C  SIGNIFICATION DE IELM1
C
C  TYPE D'ELEMENT      NOMBRE DE POINTS
C
C   0 : SEGMENT P0             1
C   1 : SEGMENT P1             2
C
C  10 : TRIANGLE P0            1
C  11 : TRIANGLE P1            3
C  12 : TRIANGLE QUASI-BULLE   4
C
C  20 : QUADRILATERE Q0        1
C  21 : QUADRILATERE Q1        4
C
C  40 : PRISMES TELEMAC-3D P0  1
C  41 : PRISMES TELEMAC-3D P1  6
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  VEC           |<-->|  VECTEUR A REMPLIR OU MODIFIER
C |  OP            | -->|  '=' : ON FAIT VEC= LE VECTEUR
C |                |    |  '+' : ON FAIT VEC=VEC+ LE VECTEUR
C |  FORMUL        | -->|  FORMULE DECRIVANT LE VECTEUR
C |  IELM1         | -->|  TYPE D'ELEMENT DU VECTEUR.
C |  XMUL          | -->|  COEFFICIENT MULTIPLICATEUR DU RESULTAT
C |  F,G,H         | -->|  FONCTIONS INTERVENANT DANS LA FORMULE
C |  U,V,W         | -->|  COMPOSANTES D'UN VECTEUR U DANS LA FORMULE
C |  MESH          | -->|  STRUCTURE DE MAILLAGE : BLOC DES ENTIERS.
C |  MSK           | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |  MASKEL        | -->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : VECTOS
C
C**********************************************************************
C
      USE BIEF, EX_VECTOR => VECTOR
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ),    INTENT(INOUT) :: VEC
      DOUBLE PRECISION,  INTENT(IN)    :: XMUL
      INTEGER,           INTENT(IN)    :: IELM1
      LOGICAL,           INTENT(IN)    :: MSK
      CHARACTER(LEN=16), INTENT(IN)    :: FORMUL
      CHARACTER(LEN=1),  INTENT(IN)    :: OP
      TYPE(BIEF_OBJ),    INTENT(IN)    :: F,G,H,U,V,W,MASKEL
      TYPE(BIEF_MESH),   INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER  :: NPT   ! Nombre de points par element
      LOGICAL  :: LEGO  ! Assemblage ou pas
      INTEGER  :: IELM0 ! Discretisation P0
C
C-----------------------------------------------------------------------
C  CHANGEMENT EVENTUEL DE DISCRETISATION.
C-----------------------------------------------------------------------
C Au cas ou on a alloue un vecteur avec le statut 1, cad on ne veut 
C pas qu'il change de discretisation, il faut tester la coherence
C entre la discretisation du vecteur et celle proposee en argument.
C autrement, on peut avoir des mauvaises surprises!!!!
C
      IF(VEC%STATUS.EQ.1.AND.VEC%ELM.NE.IELM1) THEN
        IF(LNG.EQ.1) WRITE(LU,1001) VEC%NAME,VEC%ELM,IELM1
        IF(LNG.EQ.2) WRITE(LU,1002) VEC%NAME,VEC%ELM,IELM1
1001    FORMAT(1X,'VECTOR : CHANGEMENT DE DISCRETISATION IMPOSSIBLE',
     &  ' POUR VECTEUR ',A6,' : ',1I6,' <=> ',1I6)
1002    FORMAT(1X,'VECTOR: CHANGING DISCRETIZATION FORBIDDEN',
     &  ' FOR THE VECTOR ',A6,' : ',1I6,' <=> ',1I6)
        CALL PLANTE(1)
        STOP
      ELSEIF(VEC%STATUS.EQ.2.OR.VEC%STATUS.EQ.1) THEN
        NPT = NBPTS(IELM1)
        VEC%ELM = IELM1
        VEC%DIM1= NPT
      ELSE
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'LE VECTEUR ',VEC%NAME,' A UN STATUT EGAL A',
     *                VEC%STATUS,' IL NE PEUT ETRE UTILISE DANS VECTOR'
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'VECTOR ',VEC%NAME,' HAS A STATUS ',VEC%STATUS,
     *                ' IT CANNOT BE USED IN SUBROUTINE VECTOR'
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
C
C Assemblage : il n'est pas fait pour des vecteurs
C resultat de discretisation P0.
C On met LEGO a vrai si la discretisation du vecteur resultat
C est P1, a faux autrement. 
C
      IELM0 = 10*(IELM1/10)
      LEGO  = IELM0 .NE. IELM1 
C
C Changement : dans la suite on met la variable LEGO en argument.
C
C-----------------------------------------------------------------------
C  APPEL DU SOUS-PROGRAMME D'AIGUILLAGE ET D'ASSEMBLAGE
C-----------------------------------------------------------------------
C
      IF(DIMENS(IELM1).EQ.MESH%DIM) THEN
C       VECTEUR NORMAL : APPEL AVEC SURFAC, IKLE, NELEM, NELMAX
C                                   XEL,YEL,ZEL
        CALL VECTOS(VEC%R,OP,FORMUL,XMUL,
     *              F%R,G%R,H%R,U%R,V%R,W%R,
     *              F,G,H,U,V,W,MESH%W%R,LEGO,
     *              MESH%XEL%R  , MESH%YEL%R  , MESH%ZEL%R  ,
     *              MESH%SURFAC%R,MESH%IKLE%I,MESH%NBOR%I,
     *              MESH%XSGBOR%R, MESH%YSGBOR%R, MESH%ZSGBOR%R,
     *              NPT,MESH%NELEM,MESH%NELMAX,
     *              IELM1,MESH%LV,MSK,MASKEL%R,MESH)
      ELSE
C       VECTEUR DE BORD : APPEL AVEC LGSEG, IKLBOR, NELEB, NELEBX
C                                    X,Y,Z
        IF(MESH%NELEB.GT.0) THEN
          CALL VECTOS(VEC%R,OP,FORMUL,XMUL,
     *                F%R,G%R,H%R,U%R,V%R,W%R,
     *                F,G,H,U,V,W,MESH%W%R,LEGO,
     *                MESH%X%R,MESH%Y%R,MESH%Z%R  ,
     *                MESH%LGSEG%R,MESH%IKLBOR%I,MESH%NBOR%I,
     *                MESH%XSGBOR%R,MESH%YSGBOR%R,MESH%ZSGBOR%R,
     *                NPT,MESH%NELEB,MESH%NELEBX,
     *                IELM1,MESH%LV,MSK,MASKEL%R,MESH)
        ENDIF
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
