C                       *******************
                        SUBROUTINE FILTER_H
C                       *******************
C
     *(VEC,T1,MESH,MSK,MASKEL,N,FLODEL,YAFLODEL,DT,W1,UNSV2D)
C
C***********************************************************************
C BIEF VERSION 5.9        20/05/08    J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C  FONCTION  : LISSAGE DES HAUTEURS NEGATIVES ET CALCUL DES FLUX
C              CORRESPONDANTS DANS L'EQUATION DE CONTINUITE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |    VEC         |<-->| VECTEUR A FILTRER
C |    BLDMAT      | -->| LOGIQUE : ON CONSTRUIT LA MATRICE OU PAS.
C |    T1          | -->| TABLEAU DE TRAVAIL.
C |    T2          |<-->| TABLEAU DE TRAVAIL. MATRICE A MASS-LUMPEE
C |                |    | EN SORTIE (VOIR AUSSI XMUL)
C |    A           |<-->| MATRICE (DONNEE OU CONSTRUITE SUIVANT BLDMAT)
C |    FORMUL      | -->| FORMULE DECRIVANT LA MATRICE
C |                |    | (MEMES CONVENTIONS QUE DANS MATRIX)
C |    F,G,H,U,V,W | -->| FONCTIONS INTERVENANT DANS LA MATRICE
C |    MESH,       | -->| BLOCS DU MAILLAGE.
C |    MSK,MASKEL  | -->| LOGIQUE ET TABLEAU POUR LE MASQUAGE
C |    N           | -->| NOMBRE DE FOIS OU ON FAIT L'OPERATION.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C
C  SOUS-PROGRAMMES APPELES :  LUMP
C
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)           :: N
      DOUBLE PRECISION, INTENT(IN)  :: DT
      LOGICAL, INTENT(IN)           :: MSK,YAFLODEL
      TYPE(BIEF_MESH), INTENT(INOUT):: MESH
      TYPE(BIEF_OBJ), INTENT(INOUT) :: VEC,T1,FLODEL,W1
      TYPE(BIEF_OBJ), INTENT(IN)    :: MASKEL,UNSV2D
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,NELEM
C
C-----------------------------------------------------------------------
C
      IF(YAFLODEL) CALL OV('X=C     ',W1%R,W1%R,W1%R,0.D0,3*MESH%NELMAX)
C
      DO 10 I=1,N
C
C-----------------------------------------------------------------------
C
C     CALCUL DES FLUX DUS AU LISSAGE (VOIR RELEASE NOTES 5.9)
C
      IF(YAFLODEL) THEN
        NELEM=MESH%NELEM
        CALL SMOOTHING_FLUX(-1.D0/DT,VEC,VEC%R,MESH%SURFAC%R,
     *                      MESH%IKLE%I(      1  :  NELEM),
     *                      MESH%IKLE%I(NELEM+1  :2*NELEM),
     *                      MESH%IKLE%I(2*NELEM+1:3*NELEM),
     *                      NELEM,MESH%NELMAX,W1%R)
      ENDIF
C
C-----------------------------------------------------------------------
C
C     CALCUL DU PRODUIT MATRICE DE MASSE x VEC
C
      CALL VECTOR(T1 ,'=','MASVEC          ',VEC%ELM,
     *            1.D0,VEC,VEC,VEC,VEC,VEC,VEC,MESH,MSK,MASKEL)
      IF(NCSIZE.GT.1) CALL PARCOM(T1,2,MESH)
C
C-----------------------------------------------------------------------
C
C     DIVISION PAR LA MASSE DES BASES : F = M * F / (AGGLOMEREE DE M)
C
      CALL OS('X=YZ    ',X=VEC,Y=T1,Z=UNSV2D)
C
C-----------------------------------------------------------------------
C
10    CONTINUE
C
C     PRISE EN COMPTE DES FLUX DUS AU LISSAGE DES HAUTEURS NEGATIVES
C
      IF(YAFLODEL) THEN
C                                                             
        CALL FLUX_EF_VF(FLODEL%R,W1%R,MESH%NSEG,MESH%NELEM,
     *                  MESH%ELTSEG%I,MESH%ORISEG%I,
     *                  MESH%IKLE%I,.TRUE.,0)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
