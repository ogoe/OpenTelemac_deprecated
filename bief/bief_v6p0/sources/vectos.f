C                       *****************
                        SUBROUTINE VECTOS
C                       *****************
C
     *(VEC,OP,FORMUL,
     * XMUL,F,G,H,U,V,W,SF,SG,SH,SU,SV,SW,
     * T,LEGO,
     * XEL,YEL,ZEL, SURFAC,IKLE,NBOR,
     * XNOR,YNOR,ZNOR,NPT,NELEM,NELMAX,IELM1,LV,MSK,MASKEL,MESH)
C
C***********************************************************************
C BIEF VERSION 6.0         16/07/07    JM HERVOUET (LNHE) 01 30 87 80 18
C                          22/09/05    REGINA NEBAUER
C
C***********************************************************************
C
C  FONCTION : CALCULS DE VECTEURS
C
C             LE VECTEUR EST IDENTIFIE PAR LA FORMULE CONTENUE DANS
C             LA CHAINE DE CARACTERES FORMUL
C
C-----------------------------------------------------------------------
C
C  SIGNIFICATION DE IELM1
C
C  TYPE D'ELEMENT      NOMBRE DE POINTS
C
C  10 : TRIANGLE P0            1
C  11 : TRIANGLE P1            3
C  12 : TRIANGLE QUASI-BULLE   4
C  13 : TRIANGLE P2            6
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
C |  OP            | -->|  OPERATION A EFFECTUER.01TT0
C |  FORMUL        | -->|  FORMULE DECRIVANT LE VECTEUR
C |  XMUL          | -->|  COEFFICIENT MULTIPLICATEUR DU RESULTAT
C |  F,G,H         | -->|  FONCTIONS INTERVENANT DANS LA FORMULE
C |  U,V,W         | -->|  COMPOSANTES D'UN VECTEUR U DANS LA FORMULE
C |  SF,SG,SH      | -->|  STRUCTURES DE F,G ET H
C |  SU,SV,SW      | -->|  STRUCTURES DE U,V ET W
C |  T             | -->|  TABLEAU DE TRAVAIL QUI
C |                |    |  CONTIENDRA LE VECTEUR NON ASSEMBLE
C |  LEGO          | -->|  LOGIQUE : POUR ASSEMBLER LE VECTEUR
C |  XEL,YEL,ZEL   | -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |  SURFAC        | -->|  SURFACE DES ELEMENTS.
C |  IKLE          | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
C |  NBOR          | -->|  NUMEROTATION GLOBALE DES POINTS DE BORD
C |  XNOR,YNOR,ZNOR| -->|  COMPOSANTES DES NORMALES.
C |  NPT           | -->|  NOMBRE DE POINTS DU VECTEUR.
C |  NELEM         | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |  NELMAX        | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS DU MAILLAGE ADAPTATIF)
C |  IELM1         | -->|  TYPE D'ELEMENT
C |  LV            | -->|  LONGUEUR DU VECTEUR POUR LA VECTORISATION
C |  MSK           | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |  MASKEL        | -->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : VECTOR
C PROGRAMMES APPELES : TOUS LES SOUS-PROGRAMMES DE VECTEURS ELEMENTAIRES
C                      ASSVEC
C
C**********************************************************************
C
      USE BIEF !, EX_VECTOS => VECTOS
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELMAX,NPT,NELEM,IELM1,LV
      INTEGER, INTENT(IN) :: IKLE(NELMAX,*),NBOR(*)
C
      DOUBLE PRECISION, INTENT(IN)    :: SURFAC(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: XEL(*),YEL(*),ZEL(*)
      DOUBLE PRECISION, INTENT(IN)    :: XNOR(*),YNOR(*),ZNOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: T(NELMAX,*),VEC(*)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL,MASKEL(NELMAX)
C
C     STRUCTURES DES FONCTIONS F,G,H,U,V,W ET DONNEES REELLES
C
      DOUBLE PRECISION, INTENT(IN) :: F(*),G(*),H(*),U(*),V(*),W(*)
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SF,SG,SH,SU,SV,SW
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
      LOGICAL, INTENT(IN) :: MSK,LEGO
C
      CHARACTER(LEN=16), INTENT(IN) :: FORMUL
      CHARACTER(LEN=1), INTENT(IN) ::  OP
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER ICOORD
      LOGICAL INIT,SPECAD
C
C-----------------------------------------------------------------------
C
C     TEST DU TYPE DE VECTEUR
C
C=======================================================================
C     VECTEUR MATRICE DE MASSE X VECTEUR F EGAL A 1
C=======================================================================
C
      IF(FORMUL(1:16).EQ.'MASBAS          '.OR.
     *   FORMUL(1:16).EQ.'MASBAS2         '     ) THEN
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1
C
        IF(IELM1.EQ.11) THEN
C
             CALL VC00AA(XMUL,SURFAC,NELEM,NELMAX,T(1,1),T(1,2),T(1,3))
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE QUASI-BULLE
C
        ELSEIF(IELM1.EQ.12) THEN
C
             CALL VC00BB(XMUL,SURFAC,NELEM,NELMAX,
     *                   T(1,1),T(1,2),T(1,3),T(1,4))
C
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P2
C
        ELSEIF(IELM1.EQ.13) THEN
C
             CALL VC00CC(XMUL,SURFAC,NELEM,NELMAX,
     *                   T(1,1),T(1,2),T(1,3),T(1,4),T(1,5),T(1,6))
C
C-----------------------------------------------------------------------
C
C       ELEMENT SEGMENT P1
C
C       ELSEIF(IELM1.EQ.1) THEN
C
C            CALL VC00OO(XMUL,SURFAC,NELEM,NELMAX,T(1,1),T(1,2))
C
C-----------------------------------------------------------------------
C
C       ELEMENT PRISME P1
C
        ELSEIF(IELM1.EQ.41) THEN
C
             IF(FORMUL(7:7).EQ.'2') THEN
C              COMPATIBLE WITH MASS-LUMPING IN 2D
               CALL VC00PP2(XMUL,ZEL,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   IKLE(1,5),IKLE(1,6),NELEM,NELMAX,T(1,1),T(1,2),
     *                   T(1,3),T(1,4),T(1,5),T(1,6))
             ELSE
C              NORMAL COMPUTATION OF INTEGRAL OF BASIS FUNCTIONS
               CALL VC00PP(XMUL,ZEL,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   IKLE(1,5),IKLE(1,6),NELEM,NELMAX,T(1,1),T(1,2),
     *                   T(1,3),T(1,4),T(1,5),T(1,6))
             ENDIF
C
C-----------------------------------------------------------------------
C
C       ELEMENT TETRAEDRE T1
C
        ELSEIF(IELM1.EQ.31.OR.IELM1.EQ.51) THEN
C
             CALL VC00TT(XMUL,XEL,YEL,ZEL,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   NELEM,NELMAX,T(1,1),T(1,2),T(1,3),T(1,4))
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1 (FACES LATERALES D'UN MAILLAGE DE PRISMES
C                            DECOUPES EN TETRAEDRES)
C
        ELSEIF(IELM1.EQ.61) THEN
C
             CALL VC00FT(XMUL,XEL,YEL,ZEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),NBOR,
     *                   NELEM,NELMAX,T(1,1),T(1,2),T(1,3))
C
C
C-----------------------------------------------------------------------
C
C       ELEMENT QUADRILATERE P1 
C
        ELSEIF(IELM1.EQ.71) THEN
C
C              POUR FACES RECTANGULAIRES VERTICALES DES PRISMES
               CALL VC00FF(XMUL,XEL,YEL,ZEL,
     *                     IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),NBOR,
     *                     NELEM,NELMAX,T(1,1),T(1,2),T(1,3),T(1,4))
C 
C-----------------------------------------------------------------------
C       AUTRE ELEMENT
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT
C-----------------------------------------------------------------------
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
C=======================================================================
C     VECTEUR MATRICE DE MASSE X VECTEUR F
C=======================================================================
C
      ELSEIF(FORMUL(1:6).EQ.'MASVEC') THEN
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1
C
        IF(IELM1.EQ.11) THEN
C                   
             CALL VC01AA(XMUL,SF,F,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *                   NELEM,NELMAX,T(1,1),T(1,2),T(1,3) )
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1 (FACES LATERALES D'UN MAILLAGE DE PRISMES
C                            DECOUPES EN TETRAEDRES)
C
        ELSEIF(IELM1.EQ.61.OR.IELM1.EQ.81) THEN
C
           IF(FORMUL(7:7).NE.'2') THEN
C
              CALL VC01FT(XMUL,SF,F,XEL,YEL,ZEL,
     *             IKLE(1,1),IKLE(1,2),IKLE(1,3),NBOR,
     *             NELEM,NELMAX,T(1,1),T(1,2),T(1,3)) 
           ELSE
C
              CALL VC01FT2(XMUL,SF,F,SG,G,XEL,YEL,ZEL,
     *             IKLE(1,1),IKLE(1,2),IKLE(1,3),NBOR,
     *             NELEM,NELMAX,T(1,1),T(1,2),T(1,3))
           ENDIF
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE QUASI-BULLE
C
        ELSEIF(IELM1.EQ.12) THEN
C
             CALL VC01BB(XMUL,SF,F,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   NELEM,NELMAX,T(1,1),T(1,2),T(1,3),T(1,4) )
C
C-----------------------------------------------------------------------
C
C       ELEMENT SEGMENT P1
C
        ELSEIF(IELM1.EQ.1) THEN
C
             CALL VC01OO(XMUL,SF,F,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),NBOR,NELEM,NELMAX,
     *                   T(1,1),T(1,2)  )
C
C-----------------------------------------------------------------------
C
C       ELEMENT QUADRILATERE P1 
C
        ELSEIF(IELM1.EQ.71) THEN
C
C              POUR FACES RECTANGULAIRES VERTICALES DES PRISMES
               CALL VC01FF(XMUL,SF,F,XEL,YEL,ZEL,
     *                     IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),NBOR,
     *                     NELEM,NELMAX,T(1,1),T(1,2),T(1,3),T(1,4))
C
C-----------------------------------------------------------------------
C
C       ELEMENT PRISME P1
C
        ELSEIF(IELM1.EQ.41) THEN
C
             CALL VC01PP(XMUL,SF,F,ZEL,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   IKLE(1,5),IKLE(1,6),NELEM,NELMAX,T(1,1),T(1,2),
     *                   T(1,3),T(1,4),T(1,5),T(1,6))
C
C-----------------------------------------------------------------------
C
C       ELEMENT TETRAEDRE T1
C
        ELSEIF(IELM1.EQ.31.OR.IELM1.EQ.51) THEN
C
             CALL VC01TT(XMUL,SF,F,XEL,YEL,ZEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   NELEM,NELMAX,T(1,1),T(1,2),T(1,3),T(1,4))
     
C
C-----------------------------------------------------------------------
C     
C       ELEMENT TETRAEDRE T0
C
        ELSEIF(IELM1.EQ.30) THEN
C
             CALL VC01TT0(XMUL,SF,F,XEL,YEL,ZEL,
     *                   IKLE(:,1),IKLE(:,2),IKLE(:,3),IKLE(:,4),
     *                   NELEM,NELMAX,VEC)    
     
          
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT
C-----------------------------------------------------------------------
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
C=======================================================================
C                 -    ---->   --->
C     VECTEUR :   F  . GRAD(U).GRAD(PSI)
C
C     F TENSEUR DE DIFFUSION DE COMPOSANTES DIAGONALES F,G,H
C
C     EQUIVALENT DE MATDIF * U
C                
C=======================================================================
C
      ELSEIF(FORMUL(1:16).EQ.'VECDIF          ') THEN
C
C       ELEMENT SEGMENT P1
C
C       IF(IELM1.EQ.41) THEN
C
C            CALL VC02PP(XMUL,SF,SG,SH,SU,F,G,H,U,XEL,YEL,ZEL,SURFAC,
C    *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
C    *                   IKLE(1,5),IKLE(1,6),NELEM,NELMAX,T(1,1),T(1,2),
C    *                   T(1,3),T(1,4),T(1,5),T(1,6)) 
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT
C-----------------------------------------------------------------------
C
C       ELSE
C
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
C
C       ENDIF
C
C=======================================================================
C     VECTEUR K GRAD(PSI) U.GRAD(F)
C=======================================================================
C
      ELSEIF(FORMUL(1:16).EQ.'SUPG            ') THEN
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1
C
        IF(IELM1.EQ.11) THEN
C
             CALL VC03AA(XMUL,SF,SG,SH,SU,SV,F,G,H,U,V,
     *                   XEL,YEL,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),NELEM,NELMAX,
     *                   T(1,1),T(1,2),T(1,3) )
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE QUASI-BULLE
C
        ELSEIF(IELM1.EQ.12) THEN
C
             CALL VC03BB(XMUL,SF,SG,SH,SU,SV,F,G,H,U,V,
     *                   XEL,YEL,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   NELEM,NELMAX,T(1,1),T(1,2),T(1,3),T(1,4) )
C
C-----------------------------------------------------------------------
C
C       ELEMENT SEGMENT P1
C
C       ELSEIF(IELM1.EQ.1) THEN
C
C            CALL VC03OO(XMUL,SF,SG,SH,SU,SV,SW,F,G,H,U,V,W,
C    *                   XEL,YEL,ZEL,SURFAC,
C    *                   IKLE(1,1),IKLE(1,2),NELEM,NELMAX,
C    *                   T(1,1),T(1,2)  )
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT
C-----------------------------------------------------------------------
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
C=======================================================================
C     VECTEUR U GRAD(PSI)
C=======================================================================
C
      ELSEIF(FORMUL(1:6).EQ.'VGRADP') THEN
C
        SPECAD = .FALSE.
        IF(FORMUL(8:8).EQ.'2') SPECAD = .TRUE.
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1
C
        IF(IELM1.EQ.11) THEN
C
             CALL VC04AA(XMUL,SU,SV,U,V,XEL,YEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *                   NELEM,NELMAX,T(1,1),T(1,2),T(1,3),SPECAD)
C
C-----------------------------------------------------------------------
C
C       ELEMENT PRISME P1
C
        ELSEIF(IELM1.EQ.41) THEN
C
             CALL VC04PP(XMUL,SU,SV,SW,U,V,W,SF,SG,SH,F,G,H,XEL,YEL,ZEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   IKLE(1,5),IKLE(1,6),NELEM,NELMAX,T(1,1),T(1,2),
     *                   T(1,3),T(1,4),T(1,5),T(1,6),SPECAD,FORMUL)
C
C-----------------------------------------------------------------------
C
C       ELEMENT TETRAEDRE T1
C
        ELSEIF(IELM1.EQ.31.OR.IELM1.EQ.51) THEN
C
             CALL VC04TT(XMUL,SU,SV,SW,U,V,W,XEL,YEL,ZEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   NELEM,NELMAX,T(1,1),T(1,2),
     *                   T(1,3),T(1,4),FORMUL)
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT
C-----------------------------------------------------------------------
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
C=======================================================================
C     VECTEUR U.N    (N VECTEUR NORMAL A L'ELEMENT)
C=======================================================================
C
      ELSEIF(FORMUL(1:16).EQ.'FLUBOR          ') THEN
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1 (FOND OU SURFACE D'UN MAILLAGE 3D ?)
C
        IF(IELM1.EQ.11) THEN
C
             CALL VC05AA(XMUL,SW,W,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),NELEM,NELMAX,
     *                   T(1,1),T(1,2),T(1,3) )
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1 POUR FACES VERTICALES DES PRISMES DECOUPES
C       EN TETRAEDRES
C
        ELSEIF(IELM1.EQ.61) THEN
C
             CALL VC05FT(XMUL,SU,SV,U,V,XEL,YEL,ZEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),NBOR,
     *                   NELEM,NELMAX,T(1,1),T(1,2),T(1,3))
C
C       ELEMENT QUADRILATERE P1 POUR FACES VERTICALES DES PRISMES
C
        ELSEIF(IELM1.EQ.71) THEN
C
             CALL VC05FF(XMUL,SU,SV,U,V,XEL,YEL,ZEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),NBOR,
     *                   NELEM,NELMAX,T(1,1),T(1,2),T(1,3),T(1,4))
C
C       ELEMENT SEGMENT LINEAIRE
C
        ELSEIF(IELM1.EQ.1) THEN
C
             CALL VC05OO(XMUL,SU,SV,U,V,XNOR,YNOR,SURFAC,
     *                   IKLE,NBOR,NELEM,NELMAX,T(1,1),T(1,2))
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT
C-----------------------------------------------------------------------
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
C=======================================================================
C     VECTEUR U GRAD(F)
C=======================================================================
C
      ELSEIF(FORMUL(1:13).EQ.'VGRADF       ') THEN
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1
C
        IF(IELM1.EQ.11) THEN
C
             CALL VC08AA(XMUL,SF,SU,SV,F,U,V,XEL,YEL,IKLE,
     *                   NELEM,NELMAX,T(1,1),T(1,2),T(1,3) , FORMUL )
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P2
C
        ELSEIF(IELM1.EQ.13) THEN
C
             CALL VC08CC(XMUL,SF,SU,SV,F,U,V,XEL,YEL,
     *            IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *            IKLE(1,5),IKLE(1,6),NELEM,NELMAX,
     *            T(1,1),T(1,2),T(1,3),T(1,4),T(1,5),T(1,6),FORMUL)

C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE QUASI-BULLE
C
        ELSEIF(IELM1.EQ.12) THEN
C
             CALL VC08BB(XMUL,SF,SU,SV,F,U,V,XEL,YEL,
     *            IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *            NELEM,NELMAX,T(1,1),T(1,2),T(1,3),T(1,4),FORMUL)
C
C-----------------------------------------------------------------------
C
C       ELEMENT SEGMENT P1
C
C       ELSEIF(IELM1.EQ.1) THEN
C
C            CALL VC08OO(XMUL,SF,SG,SH,SU,SV,SW,F,G,H,U,V,W,
C    *                   XEL,YEL,ZEL,SURFAC,
C    *                   IKLE(1,1),IKLE(1,2),NELEM,NELMAX,
C    *                   T(1,1),T(1,2)  )
C
C
C-----------------------------------------------------------------------
C
C       ELEMENT PRISME P1
C
        ELSEIF(IELM1.EQ.41) THEN
C
             CALL VC08PP(XMUL,SF,SU,SV,SW,F,U,V,W,XEL,YEL,ZEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   IKLE(1,5),IKLE(1,6),NELEM,NELMAX,T(1,1),T(1,2),
     *                   T(1,3),T(1,4),T(1,5),T(1,6))
C
C-----------------------------------------------------------------------
C
C       ELEMENT TETRAEDRE T1
C
        ELSEIF(IELM1.EQ.31.OR.IELM1.EQ.51) THEN
C
             CALL VC08TT(XMUL,SF,SU,SV,SW,F,U,V,W,XEL,YEL,ZEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   NELEM,NELMAX,T(1,1),T(1,2),
     *                   T(1,3),T(1,4))
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT
C-----------------------------------------------------------------------
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
! A.D. MODIF 25/11/04
!
C=======================================================================
C     VECTEUR U GRAD(F) 2 
C=======================================================================
C
      ELSEIF(FORMUL(1:13).EQ.'VGRADF2      ') THEN
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE QUASI-BULLE
C
C-----------------------------------------------------------------------
C
C       ELEMENT SEGMENT P1
C
C-----------------------------------------------------------------------
C
C       ELEMENT PRISME P1
C
        IF(IELM1.EQ.41) THEN
C
             CALL VC18PP(XMUL,SF,SU,SV,F,U,V,XEL,YEL,ZEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   IKLE(1,5),IKLE(1,6),NELEM,NELMAX,T(1,1),T(1,2),
     *                   T(1,3))
C
C-----------------------------------------------------------------------
C
C       ELEMENT TETRAEDRE T1
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT
C-----------------------------------------------------------------------
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
!  A.D. FIN MODIF 25/11/04
C=======================================================================
C     VECTEUR Q GRAD(F)
C=======================================================================
C
      ELSEIF(FORMUL(1:16).EQ.'QGRADF          ') THEN
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1
C
        IF(IELM1.EQ.11) THEN
C
             CALL VC09AA(XMUL,SF,SG,SU,SV,F,G,U,V,XEL,YEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),NELEM,NELMAX,
     *                   T(1,1),T(1,2),T(1,3) )
C
C-----------------------------------------------------------------------
C
C       ELEMENT SEGMENT P1
C
C       ELSEIF(IELM1.EQ.1) THEN
C
C            CALL VC09OO(XMUL,SF,SG,SH,SU,SV,SW,F,G,H,U,V,W,
C    *                   XEL,YEL,ZEL,SURFAC,
C    *                   IKLE(1,1),IKLE(1,2),NELEM,NELMAX,
C    *                   T(1,1),T(1,2)  )
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT
C-----------------------------------------------------------------------
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
C=======================================================================
C     VECTEUR F U.N    (N VECTEUR NORMAL A L'ELEMENT)
C=======================================================================
C
      ELSEIF(FORMUL(1:16).EQ.'FLUBDF          ') THEN
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1
C
C       IF(IELM1.EQ.11) THEN
C
C            CALL VC10AA(XMUL,SF,SG,SH,SU,SV,SW,F,G,H,U,V,W,
C    *                   XEL,YEL,ZEL,SURFAC,
C    *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),NELEM,NELMAX,
C    *                   T(1,1)   ,T(1,2)   ,T(1,3))
C
C-----------------------------------------------------------------------
C
C       ELEMENT SEGMENT P1
C
        IF(IELM1.EQ.1) THEN
C
             CALL VC10OO(XMUL,SF,SU,SV,F,U,V,XNOR,YNOR,SURFAC,
     *                   IKLE,NBOR,NELEM,NELMAX,T(1,1),T(1,2)  )
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT
C-----------------------------------------------------------------------
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
C=======================================================================
C     VECTEUR G GRADIENT(F)
C=======================================================================
C
      ELSEIF(FORMUL(1:15).EQ.'GGRADF         ') THEN
C
C     LE CARACTERE 16 EST LA COORDONNEE CHOISIE
C
      IF(FORMUL(16:16).EQ.'X') THEN
        ICOORD=1
      ELSEIF(FORMUL(16:16).EQ.'Y') THEN
        ICOORD=2
      ELSEIF(FORMUL(16:16).EQ.'Z') THEN
        ICOORD=3
      ENDIF
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1
C
        IF(IELM1.EQ.11) THEN
C
C TEST SI G EST P1 DISC
C
         IF (SG%DIM1.EQ.NELEM.AND.SG%DIM2.EQ.3.AND.
     *       SG%DIMDISC.EQ.11) THEN
C
            CALL VC11AA2(XMUL,SF,SG,SH,F,G,H,XEL,YEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),NELEM,NELMAX,
     *                   T(1,1),T(1,2),T(1,3) , ICOORD )
C
C CAS CLASSIQUE : G P1
C
         ELSE
C
             CALL VC11AA(XMUL,SF,SG,F,G,XEL,YEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),NELEM,NELMAX,
     *                   T(1,1),T(1,2),T(1,3) , ICOORD )
         ENDIF
C
C-----------------------------------------------------------------------
C
        ELSEIF(IELM1.EQ.12) THEN
C
             CALL VC11BB(XMUL,SF,SG,F,G,XEL,YEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   NELEM,NELMAX,
     *                   T(1,1),T(1,2),T(1,3),T(1,4) , ICOORD )
C
C-----------------------------------------------------------------------
C
C       ELEMENT SEGMENT P1
C
C       ELSEIF(IELM1.EQ.1) THEN
C
C            CALL VC11OO(XMUL,SF,SG,SH,SU,SV,SW,F,G,H,U,V,W,
C    *                   XEL,YEL,ZEL,XNOR,YNOR,ZNOR,SURFAC,
C    *                   IKLE(1,1),IKLE(1,2),NBOR,NELEM,NELMAX,
C    *                   T(1,1),T(1,2)  )
C
C-----------------------------------------------------------------------
C
C       ELEMENT PRISME P1
C
        ELSEIF(IELM1.EQ.41) THEN
C
             CALL VC11PP(XMUL,SF,SG,F,G,
     *                   XEL,YEL,ZEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   IKLE(1,5),IKLE(1,6),NELEM,NELMAX,T(1,1),T(1,2),
     *                   T(1,3),T(1,4),T(1,5),T(1,6),ICOORD)
C
C-----------------------------------------------------------------------
C
C       ELEMENT TETRAEDRE T1
C
        ELSEIF(IELM1.EQ.31.OR.IELM1.EQ.51) THEN
C
             CALL VC11TT(XMUL,SF,SG,F,G,
     *                   XEL,YEL,ZEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   NELEM,NELMAX,T(1,1),T(1,2),
     *                   T(1,3),T(1,4),ICOORD)
C     
        ELSEIF(IELM1.EQ.30) THEN 
C       
! Comment JP Renaud 21/09/2005
! Attention: this call create a P0 vector. Therefore
! there is no need to "assemble" it afterwards. Note 
! that VEC itself (and not T) is given to VC11TT0 and
! that LEGO is also set to .FALSE. to avoid
! the call ASSVEC at the end of the subroutine.
C
             CALL VC11TT0(XMUL,SF,SG,F,G,
     *                   XEL,YEL,ZEL,
     *                   IKLE(1:NELEM,1),IKLE(1:NELEM,2),
     *                   IKLE(1:NELEM,3),IKLE(1:NELEM,4),
     *                   NELEM,MESH%NPOIN,
     *                   VEC,ICOORD)
C     
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT
C-----------------------------------------------------------------------
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
C=======================================================================
C     VECTEUR GRADIENT(F)
C=======================================================================
C
      ELSEIF(FORMUL(1:5).EQ.'GRADF') THEN
C
C     LE CARACTERE 16 EST LA COORDONNEE CHOISIE
C
      IF(FORMUL(16:16).EQ.'X') THEN
        ICOORD=1
      ELSEIF(FORMUL(16:16).EQ.'Y') THEN
        ICOORD=2
      ELSEIF(FORMUL(16:16).EQ.'Z') THEN
        ICOORD=3
      ENDIF
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1
C
        IF(IELM1.EQ.11) THEN
C
             CALL VC13AA(XMUL,SF,F,XEL,YEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   NELEM,NELMAX,T(1,1),T(1,2),T(1,3),ICOORD)
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE QUASI-BULLE
C
        ELSEIF(IELM1.EQ.12) THEN
C
             CALL VC13BB(XMUL,SF,F,XEL,YEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),NELEM,
     *                   NELMAX,T(1,1),T(1,2),T(1,3),T(1,4),ICOORD)
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P2
C
        ELSEIF(IELM1.EQ.13) THEN
C
            CALL VC13CC(XMUL,SF,F,XEL,YEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *                   IKLE(1,4),IKLE(1,5),IKLE(1,6),
     *                   NELEM,NELMAX,
     *                   T(1,1),T(1,2),T(1,3), 
     *                   T(1,4),T(1,5),T(1,6),ICOORD)
C
C-----------------------------------------------------------------------
C
C       ELEMENT SEGMENT P1
C
C       ELSEIF(IELM1.EQ.1) THEN
C
C            CALL VC13OO(XMUL,SF,SG,SH,SU,SV,SW,F,G,H,U,V,W,
C    *                   XEL,YEL,ZEL,XNOR,YNOR,ZNOR,SURFAC,
C    *                   IKLE(1,1),IKLE(1,2),NBOR,NELEM,NELMAX,
C    *                   T(1,1),T(1,2)  )
C
C
C
C-----------------------------------------------------------------------
C
C       ELEMENT SEGMENT P2
C
C       ELSEIF(IELM1.EQ.2) THEN
C
C            CALL VC13OC(XMUL,SF,SG,SH,SU,SV,SW,F,G,H,U,V,W,
C    *                   XEL,YEL,ZEL,XNOR,YNOR,ZNOR,SURFAC,
C    *                   IKLE(1,1),IKLE(1,2),NBOR,NELEM,NELMAX,
C    *                   T(1,1),T(1,2),T(1,3)  )
C
C
C-----------------------------------------------------------------------
C
C       ELEMENT PRISME P1
C
        ELSEIF(IELM1.EQ.41) THEN
C
             IF(FORMUL(1:15).EQ.'GRADF(X,Y)     ') THEN
C            SIMPLIFIED FORMULA FOR EFFICIENCY AND ACCURACY
             CALL VC13PP2(XMUL,SF,F,XEL,YEL,ZEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   IKLE(1,5),IKLE(1,6),NELEM,NELMAX,T(1,1),T(1,2),
     *                   T(1,3),T(1,4),T(1,5),T(1,6),ICOORD)
             ELSEIF( FORMUL(8:15).EQ.'        ') THEN
C                    FORMUL(6:7) IS LEFT FOR OPTIONS
             CALL VC13PP(XMUL,SF,F,XEL,YEL,ZEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   IKLE(1,5),IKLE(1,6),NELEM,NELMAX,T(1,1),T(1,2),
     *                   T(1,3),T(1,4),T(1,5),T(1,6),ICOORD,FORMUL)
C           
             ELSE
               IF(LNG.EQ.1) WRITE(LU,1000) FORMUL
               IF(LNG.EQ.2) WRITE(LU,1001) FORMUL
               CALL PLANTE(1)
               STOP             
             ENDIF
C
C-----------------------------------------------------------------------
C
C       ELEMENT TETRAEDRE T1
C
        ELSEIF(IELM1.EQ.31.OR.IELM1.EQ.51) THEN
C
             IF(FORMUL(1:15).EQ.'GRADF          '.OR.
     *          FORMUL(1:15).EQ.'GRADF2         '  ) THEN
             CALL VC13TT(XMUL,SF,F,XEL,YEL,ZEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   NELEM,NELMAX,T(1,1),T(1,2),
     *                   T(1,3),T(1,4),ICOORD,FORMUL)
C      
             ELSE
               IF(LNG.EQ.1) WRITE(LU,1000) FORMUL
               IF(LNG.EQ.2) WRITE(LU,1001) FORMUL
               CALL PLANTE(1)
               STOP             
             ENDIF
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT
C-----------------------------------------------------------------------
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
C=======================================================================
C     VECTEUR PRODUCTION DE TURBULENCE
C=======================================================================
C
      ELSEIF(FORMUL(1:16).EQ.'PRODF           ') THEN
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1
C
        IF(IELM1.EQ.11) THEN
C
             CALL VC14AA(XMUL,SU,SV,U,V,XEL,YEL,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),NELEM,NELMAX,
     *                   T(1,1),T(1,2),T(1,3))
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE QUASI-BULLE
C
C       ELSEIF(IELM1.EQ.12) THEN
C
C            CALL VC14BB(XMUL,SF,SG,SH,SU,SV,SW,F,G,H,U,V,W,
C    *                   XEL,YEL,ZEL,SURFAC,
C    *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),NELEM,
C    *                   NELMAX,T(1,1),T(1,2),T(1,3),T(1,4))
C
C-----------------------------------------------------------------------
C
C       ELEMENT SEGMENT P1
C
C       ELSEIF(IELM1.EQ.1) THEN
C
C            CALL VC14OO(XMUL,SF,SG,SH,SU,SV,SW,F,G,H,U,V,W,
C    *                   XEL,YEL,ZEL,XNOR,YNOR,ZNOR,SURFAC,
C    *                   IKLE(1,1),IKLE(1,2),NBOR,NELEM,NELMAX,
C    *                   T(1,1),T(1,2))
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT
C-----------------------------------------------------------------------
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
C=======================================================================
C     VECTEUR DIV(HU)
C=======================================================================
C
      ELSEIF(FORMUL(1:16).EQ.'DIVQ            ') THEN
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1
C
        IF(IELM1.EQ.11) THEN
C
             CALL VC15AA(XMUL,SF,SU,SV,F,U,V,XEL,YEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   NELEM,NELMAX,T(1,1),T(1,2),T(1,3))
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE QUASI-BULLE
C
C       ELSEIF(IELM1.EQ.12) THEN
C
C            CALL VC15BB(XMUL,SF,SG,SH,SU,SV,SW,F,G,H,U,V,W,
C    *            XEL,YEL,ZEL,SURFAC,
C    *            IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
C    *            NELEM,NELMAX,T(1,1),T(1,2),T(1,3),T(1,4),FORMUL)
C
C-----------------------------------------------------------------------
C
C       ELEMENT SEGMENT P1
C
C       ELSEIF(IELM1.EQ.1) THEN
C
C            CALL VC15OO(XMUL,SF,SG,SH,SU,SV,SW,F,G,H,U,V,W,
C    *                   XEL,YEL,ZEL,SURFAC,
C    *                   IKLE(1,1),IKLE(1,2),NELEM,NELMAX,
C    *                   T(1,1),T(1,2)  )
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT
C-----------------------------------------------------------------------
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
C=======================================================================
C     VECTEUR K.GRAD(PSI) DIV(U)
C=======================================================================
C
      ELSEIF(FORMUL(1:16).EQ.'SUPGDIVU        ') THEN
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1
C
        IF(IELM1.EQ.11) THEN
C
             CALL VC16AA(XMUL,SF,SG,SU,SV,F,G,U,V,XEL,YEL,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *                   NELEM,NELMAX,
     *                   T(1,1),T(1,2),T(1,3) )
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE QUASI-BULLE
C
C       ELSEIF(IELM1.EQ.12) THEN
C
C            CALL VC16BB(XMUL,SF,SG,SH,SU,SV,SW,F,G,H,U,V,W,
C    *            XEL,YEL,ZEL,SURFAC,
C    *            IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
C    *            NELEM,NELMAX,T(1,1),T(1,2),T(1,3),T(1,4),FORMUL)
C
C-----------------------------------------------------------------------
C
C       ELEMENT SEGMENT P1
C
C       ELSEIF(IELM1.EQ.1) THEN
C
C            CALL VC16OO(XMUL,SF,SG,SH,SU,SV,SW,F,G,H,U,V,W,
C    *                   XEL,YEL,ZEL,SURFAC,
C    *                   IKLE(1,1),IKLE(1,2),NELEM,NELMAX,
C    *                   T(1,1),T(1,2)  )
C
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT
C-----------------------------------------------------------------------
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
C
C=======================================================================
C     VECTEUR H U.GRAD(PSI)
C=======================================================================
C
      ELSEIF(FORMUL(1:7).EQ.'HUGRADP') THEN
C
C-----------------------------------------------------------------------
C
C       ELEMENT TRIANGLE P1
C
        IF(IELM1.EQ.11) THEN
C
             CALL VC19AA(XMUL,SF,SG,SH,SU,SV,F,G,H,U,V,
     *                   XEL,YEL,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),NELEM,NELMAX,
     *                   T(1,1),T(1,2),T(1,3),FORMUL)
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT
C-----------------------------------------------------------------------
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
C=======================================================================
C     AUTRE VECTEUR
C=======================================================================
C
C
C=======================================================================
C     ERREUR : TYPE DE VECTEUR NON PROGRAMME
C=======================================================================
C
      ELSE
        IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
        IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C  ASSEMBLAGE EVENTUEL DU VECTEUR
C
      IF(LEGO) THEN
C
        IF(OP(1:1).EQ.'=') THEN
          INIT = .TRUE.
        ELSEIF(OP(1:1).EQ.'+') THEN
          INIT = .FALSE.
        ELSE
          IF (LNG.EQ.1) WRITE(LU,3000) OP
          IF (LNG.EQ.2) WRITE(LU,3001) OP
3000      FORMAT(1X,'VECTOS (BIEF) : OP NON PREVU :',A1)
3001      FORMAT(1X,'VECTOS (BIEF) : OP NOT RECOGNISED:',A1)
          CALL PLANTE(1)
          STOP
        ENDIF
C
        CALL ASSVEC(VEC, IKLE, NPT ,NELEM,NELMAX,IELM1,
     *              T,INIT,LV,MSK,MASKEL)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
1000  FORMAT(1X,'VECTOS (BIEF) : VECTEUR NON PREVU : ',A16)
1001  FORMAT(1X,'VECTOS (BIEF) : VECTOR NOT IMPLEMENTED:',A16)
2000  FORMAT(1X,'                POUR IELM1 = ',1I6)
2001  FORMAT(1X,'                FOR IELM1 = ',1I6)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
