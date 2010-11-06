C                       ***************
                        SUBROUTINE GOUP
C                       ***************
C
     *(X, A,B ,DITR,MESH,COPY)
C
C***********************************************************************
C BIEF VERSION 5.1           23/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FONCTION : RESOLUTION DU SYSTEME U X = B (ELEMENT PAR ELEMENT)
C
C            ICI LA MATRICE U EST LE RESULTAT D'UNE DECOMPOSITION
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
C PROGRAMMES APPELES : DESCEN , REMONT , PLANTE
C
C**********************************************************************
C
      USE BIEF, EX_GOUP => GOUP
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
      LOGICAL, INTENT(IN) :: COPY
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
300     FORMAT(1X,'GOUP (BIEF) :',1I6,' TYPE DE A NON PREVU.')
400     FORMAT(1X,'GOUP (BIEF) :',1I6,' UNEXPECTED TYPE FOR A.')
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
        CALL GOUP1(X, A,B ,DITR,MESH,COPY)
C
      ELSEIF(S.GT.0.AND.S.EQ.SA) THEN
C
C     CAS OU LE BLOC A NE CONTIENT QUE LES DIAGONALES
C
        DO 10 I=1,S
          CALL GOUP1( X%ADR(I)%P,A%ADR(I)%P,B%ADR(I)%P,DITR,MESH,COPY)
10      CONTINUE
C
      ELSEIF(S.GT.0.AND.S**2.EQ.SA) THEN
C
C     CAS OU LE BLOC A CONTIENT AUTANT DE MATRICES QUE LE SYSTEME
C     COMPLET : ON NE PREND QUE LES DIAGONALES.
C
        DO 11 I=1,S
          CALL GOUP1( X%ADR(I)%P,
     *                A%ADR(1+(S+1)*(I-1))%P,
     *                B%ADR(I)%P,
     *                DITR,MESH,COPY)
11      CONTINUE
C
C     CAS OU A EST UNE SEULE MATRICE ET X UN BLOC
C
      ELSEIF(S.GT.0.AND.SA.EQ.0) THEN
C
        DO 12 I=1,S
          CALL GOUP1(X%ADR(I)%P,
     *                A,
     *                B%ADR(I)%P,
     *                DITR,MESH,COPY)
12      CONTINUE
C
      ELSE
        IF (LNG.EQ.1) WRITE(LU,301)
        IF (LNG.EQ.2) WRITE(LU,401)
301     FORMAT(1X,'GOUP (BIEF) : CAS NON PREVU')
401     FORMAT(1X,'GOUP (BIEF) : UNEXPECTED CASE')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
 
