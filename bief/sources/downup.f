C                       *****************
                        SUBROUTINE DOWNUP
C                       *****************
C
     *(X, A,B ,DITR,MESH)
C
C***********************************************************************
C BIEF VERSION 5.1           24/04/97    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FONCTION : RESOLUTION DU SYSTEME A X = B
C
C            ICI LA MATRICE A EST LE RESULTAT D'UNE DECOMPOSITION
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
C            LES MATRICES "DE" SONT CONSIDEREES COMME DES DIAGONALES
C            DE TAILLE NPOIN X NPOIN QU'IL FAUT IMAGINER COMPLETEES
C            AVEC DES 1 POUR LES POINTS QUI N'APPARTIENNENT PAS A
C            L'ELEMENT CONSIDERE.
C
C            ON A EFFECTUE ENSUITE LE PRODUIT DE TOUTES CES DIAGONALES
C            CE QUI A DONNE LA DIAGONALE DB.
C
C !!!!!!!!!  ENFIN : DB A ETE INVERSEE CAR C'EST SOUS CETTE FORME
C                    QU'ELLE EST UTILISEE ICI.
C
C            LA MATRICE A EST ICI :
C
C            LE PRODUIT DE 1 A NELEM DE TOUTES LES MATRICES LE
C
C            MULTIPLIE PAR :
C
C            LA DIAGONALE DB
C
C            MULTIPLIE PAR :
C
C            LE PRODUIT DE NELEM A 1 DE TOUTES LES MATRICES UE.
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
C |      DA        |<-- |  DIAGONALE DE LA MATRICE A
C |      TYPDIA    |<-- |  TYPE DE DIAGONALE ( 'Q', 'I' , OU '0' )
C |      XA        |<-- |  TERMES EXTRADIAGONAUX DE LA MATRICE A
C |      B         |<-- |  SECOND MEMBRE DU SYSTEME A RESOUDRE.
C |      IKLE      | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |      NPOIN     | -->|  DIMENSION DES TABLEAUX
C |      IELM      | -->|  TYPE D'ELEMENT
C |      DITR      | -->|  CARACTERE  'D' : ON CALCULE AVEC A
C |                |    |             'T' : ON CALCULE AVEC A TRANSPOSEE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : DESCEN , REMONT , PLANTE
C
C**********************************************************************
C
      USE BIEF, EX_DOXNUP => DOWNUP
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: X
      TYPE(BIEF_OBJ), INTENT(IN)    :: A,B
      TYPE(BIEF_MESH), INTENT(IN)   :: MESH
      CHARACTER(LEN=1), INTENT(IN)  :: DITR
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER S,SA,I
C
C-----------------------------------------------------------------------
C
      IF(X%TYPE.EQ.4) THEN
        S = X%N
      ELSE
        S = 0
      ENDIF
C
C     ON PREVOIT LE CAS OU LE SYSTEME EST UN BLOC MAIS OU ON N'UTILISE
C     QU'UNE SEULE MATRICE DE PRECONDITIONNEMENT.
C
      IF(A%TYPE.EQ.3) THEN
        SA = 0
      ELSEIF(A%TYPE.EQ.4) THEN
        SA = A%N
      ELSE
        IF (LNG.EQ.1) WRITE(LU,300) A%TYPE
        IF (LNG.EQ.2) WRITE(LU,400) A%TYPE
300     FORMAT(1X,'DOWNUP (BIEF) :',1I6,' TYPE DE A NON PREVU.')
400     FORMAT(1X,'DOWNUP (BIEF) :',1I6,' UNEXPECTED TYPE FOR A.')
        CALL PLANTE(0)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(S.EQ.0.AND.SA.EQ.0) THEN
C
C     CAS OU A EST UNE MATRICE SIMPLE ET X UN VECTEUR SIMPLE.
C
        CALL DWNUP1(X, A,B ,DITR,MESH)
C
      ELSEIF(S.GT.0.AND.S.EQ.SA) THEN
C
C     CAS OU LE BLOC A NE CONTIENT QUE LES DIAGONALES
C
        DO 10 I=1,S
          CALL DWNUP1(X%ADR(I)%P,
     *                A%ADR(I)%P,
     *                B%ADR(I)%P,
     *                DITR,MESH)
10      CONTINUE
C
      ELSEIF(S.GT.0.AND.S**2.EQ.SA) THEN
C
C     CAS OU LE BLOC A CONTIENT AUTANT DE MATRICES QUE LE SYSTEME
C     COMPLET : ON NE PREND QUE LES DIAGONALES.
C
        DO 11 I=1,S
          CALL DWNUP1(X%ADR(I)%P,
     *                A%ADR(1+(S+1)*(I-1))%P,
     *                B%ADR(I)%P,
     *                DITR,MESH)
11      CONTINUE
C
C     CAS OU A EST UNE SEULE MATRICE ET X UN BLOC
C
      ELSEIF(S.GT.0.AND.SA.EQ.0) THEN
C
        DO 12 I=1,S
          CALL DWNUP1(X%ADR(I)%P,
     *                A,
     *                B%ADR(I)%P,
     *                DITR,MESH)
12      CONTINUE
C
      ELSE
        IF (LNG.EQ.1) WRITE(LU,301)
        IF (LNG.EQ.2) WRITE(LU,401)
301     FORMAT(1X,'DOWNUP (BIEF) : CAS NON PREVU')
401     FORMAT(1X,'DOWNUP (BIEF) : UNEXPECTED CASE')
        CALL PLANTE(0)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
 
