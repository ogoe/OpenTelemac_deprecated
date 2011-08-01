C                       ****************
                        SUBROUTINE PUOG1
C                       ****************
C
     *(X, A,B ,DITR,MESH,COPY)
C
C***********************************************************************
C BIEF VERSION 5.1           23/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR X = U B     (ELEMENT PAR ELEMENT)
C
C            ICI LA MATRICE L EST LE RESULTAT D'UNE DECOMPOSITION
C            EFFECTUEE PAR LE SOUS-PROGRAMME DECLDU.
C
C            CHAQUE MATRICE ELEMENTAIRE A ETE DECOMPOSEE SOUS LA FORME :
C
C            LE X DE X UE
C
C            LE : TRIANGULAIRE INFERIEURE AVEC DES 1 SUR LA DIAGONALE.
C            DE : DIAGONALE
C            UE : TRIANGULAIRE SUPERIEURE AVEC DES 1 SUR LA DIAGONALE.
C
C                                                T
C            SI LA MATRICE EST SYMETRIQUE : LE =  UE
C
C-----------------------------------------------------------------------
C  SIGNIFICATION DE IELM :
C
C  TYPE D'ELEMENT      NOMBRE DE POINTS          PROGRAMME ICI
C
C  11 : TRIANGLE P1            3                       OUI
C  12 : TRIANGLE QUASI-BULLE   4                       OUI
C  21 : QUADRILATERE Q1        4                       OUI
C  41 : PRISMES TELEMAC-3D     6                       OUI
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      X         |<-- |  SOLUTION DU SYSTEME AX = B
C |      A         |<-- |  MATRICE A SOUS FORME LDU
C |      B         |<-- |  SECOND MEMBRE DU SYSTEME A RESOUDRE.
C |      DITR      | -->|  CARACTERE  'D' : ON CALCULE AVEC A
C |                |    |             'T' : ON CALCULE AVEC A TRANSPOSEE
C |      MESH      | -->|  BLOC DES TABLEAUX D'ENTIERS DU MAILLAGE.
C |      COPY      | -->|  SI .TRUE. B EST RECOPIE SUR X.
C |                |    |  AU PREALABLE.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : TNOMER , PLANTE
C
C**********************************************************************
C
      USE BIEF, EX_PUOG1 => PUOG1
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      CHARACTER(LEN=1), INTENT(IN) :: DITR
C
      LOGICAL, INTENT(IN) :: COPY
C
C-----------------------------------------------------------------------
C
C  STRUCTURES DE VECTEURS
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: X,B
C
C-----------------------------------------------------------------------
C
C  STRUCTURE DE MAILLAGE
C
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C-----------------------------------------------------------------------
C
C  STRUCTURE DE MATRICE
C
      TYPE(BIEF_OBJ), INTENT(IN) :: A
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELM,NPOIN,NELEM,NELMAX
      CHARACTER(LEN=1) :: TYPX
C
C-----------------------------------------------------------------------
C
      TYPX  = A%TYPEXT
      NPOIN = A%D%DIM1
      IELM  = A%ELMLIN
      NELEM = MESH%NELEM
      NELMAX= MESH%NELMAX
      CALL CPSTVC(B,X)
C
C-----------------------------------------------------------------------
C
C 1) DESCENTE AVEC RECOPIE DE B DANS X
C
      CALL TNOMER(X%R,A%X%R,TYPX,
     *     B%R,MESH%IKLE%I,NELEM,NELMAX,NPOIN,IELM,DITR,COPY,MESH%LV)
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
