C                       *****************
                        SUBROUTINE DESCEN
C                       *****************
C
     *(X, XA,TYPEXA,B,IKLE,NELEM,NELMAX,NPOIN,IELM,DITR,COPY,LV)
C
C***********************************************************************
C BIEF VERSION 5.3           18/08/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : RESOLUTION DU SYSTEME L X = B (ELEMENT PAR ELEMENT)
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
C
C-----------------------------------------------------------------------
C  SIGNIFICATION DE IELM :
C
C  TYPE D'ELEMENT      NOMBRE DE POINTS          PROGRAMME ICI
C
C  11 : TRIANGLE P1            3                       OUI
C  12 : TRIANGLE P2            6
C  13 : TRIANGLE P1-ISO P1     6
C  14 : TRIANGLE P2            7
C  21 : QUADRILATERE Q1        4                       OUI
C  22 : QUADRILATERE Q2        8
C  24 : QUADRILATERE Q2        9
C  31 : TETRAEDRE P1           4
C  32 : TETRAEDRE P2          10
C  41 : PRISMES TELEMAC-3D     6                       OUI
C  41 : PRISMES DECOUPES En TETRAEDRES     4           OUI
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      X         |<-- |  SOLUTION DU SYSTEME AX = B
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
C |      COPY      | -->|  SI .TRUE. B EST RECOPIE DANS X AU PREALABLE
C |      LV        | -->|  LONGUEUR DU VECTEUR POUR LA VECTORISATION
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : DES11 , DES21 , DES41 , PLANTE
C
C**********************************************************************
C
      USE BIEF, EX_DESCEN => DESCEN
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(IN)    :: IELM,NPOIN,NELEM,NELMAX,LV
      INTEGER         , INTENT(IN)    :: IKLE(NELMAX,*)
      DOUBLE PRECISION, INTENT(INOUT) :: X(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: XA(NELMAX,*),B(NPOIN)
      CHARACTER(LEN=1), INTENT(IN)    :: TYPEXA,DITR
      LOGICAL         , INTENT(IN)    :: COPY
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C 1) INITIALISATION : X = SECOND MEMBRE
C
      IF(COPY) CALL OV( 'X=Y     ' , X , B , B , 0.D0 , NPOIN )
C
C-----------------------------------------------------------------------
C
C 2) INVERSION DES MATRICES TRIANGULAIRES INFERIEURES (DESCENTE)
C
C     2.1) CAS TRANSPOSE
C
      IF(TYPEXA(1:1).EQ.'S' .OR.
     *  (TYPEXA(1:1).EQ.'Q'.AND.DITR(1:1).EQ.'T')) THEN
C
      IF(IELM.EQ.11) THEN
C
        CALL DES11(X,XA(1,1),XA(1,2),XA(1,3),
     *             IKLE(1,1),IKLE(1,2),IKLE(1,3),NELEM,NELMAX,NPOIN,LV)
C
      ELSEIF(IELM.EQ.21.OR.IELM.EQ.12.OR.IELM.EQ.31.OR.IELM.EQ.51) THEN
C
        CALL DES21(X,XA(1,1),XA(1,2),XA(1,3),XA(1,4),XA(1,5),XA(1,6),
     *             IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *             NELEM,NELMAX,NPOIN,LV)
C
      ELSEIF(IELM.EQ.41) THEN
C
        CALL DES41(X,XA(1,1),XA(1,2),XA(1,3),XA(1,4),XA(1,5),XA(1,6),
     *             XA(1,7),XA(1,8),XA(1,9),XA(1,10),XA(1,11),XA(1,12),
     *             XA(1,13),XA(1,14),XA(1,15),
     *             IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *             IKLE(1,4),IKLE(1,5),IKLE(1,6),NELEM,NELMAX,NPOIN,LV)
C
C  IELM NON PREVU : ERREUR
C
      ELSE
C
       IF (LNG.EQ.1) WRITE(LU,100) IELM
       IF (LNG.EQ.2) WRITE(LU,101) IELM
100    FORMAT(1X,'DESCEN (BIEF) : IELM = ',1I6,' ELEMENT NON PREVU')
101    FORMAT(1X,'DESCEN (BIEF) : IELM = ',1I6,' ELEMENT NOT AVAILABLE')
       CALL PLANTE(1)
       STOP
C
      ENDIF
C
C     2.2) CAS DIRECT
C
      ELSEIF(TYPEXA(1:1).EQ.'Q'.AND.DITR(1:1).EQ.'D') THEN
C
      IF(IELM.EQ.11) THEN
C
        CALL DES11(X,XA(1,4),XA(1,5),XA(1,6),
     *             IKLE(1,1),IKLE(1,2),IKLE(1,3),NELEM,NELMAX,NPOIN,LV)
C
      ELSEIF(IELM.EQ.21.OR.IELM.EQ.12.OR.IELM.EQ.31.OR.IELM.EQ.51) THEN
C
        CALL DES21(X,XA(1,7),XA(1,8),XA(1,9),XA(1,10),XA(1,11),XA(1,12),
     *             IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *             NELEM,NELMAX,NPOIN,LV)
C
      ELSEIF(IELM.EQ.41) THEN
C
        CALL DES41(X,
     *            XA(1,16),XA(1,17),XA(1,18),XA(1,19),XA(1,20),XA(1,21),
     *            XA(1,22),XA(1,23),XA(1,24),XA(1,25),XA(1,26),XA(1,27),
     *            XA(1,28),XA(1,29),XA(1,30),
     *            IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *            IKLE(1,4),IKLE(1,5),IKLE(1,6),NELEM,NELMAX,NPOIN,LV)
C
C  IELM NON PREVU : ERREUR
C
      ELSE
C
       IF (LNG.EQ.1) WRITE(LU,100) IELM
       IF (LNG.EQ.2) WRITE(LU,101) IELM
       CALL PLANTE(1)
       STOP
C
      ENDIF
C
C     2.3) CAS NON PREVU
C
      ELSE
         IF (LNG.EQ.1) WRITE(LU,200) TYPEXA(1:1),DITR(1:1)
         IF (LNG.EQ.2) WRITE(LU,201) TYPEXA(1:1),DITR(1:1)
200      FORMAT(1X,'DESCEN (BIEF) : TYPE DE MATRICE NON PREVU :',A1,/,
     *          1X,'AVEC DITR=',A1)
201      FORMAT(1X,'DESCEN (BIEF) : UNEXPECTED TYPE OF MATRIX :',A1,/,
     *          1X,'WITH DITR=',A1)
         CALL PLANTE(1)
         STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
