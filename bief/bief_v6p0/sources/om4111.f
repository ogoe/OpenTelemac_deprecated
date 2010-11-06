C                       *****************
                        SUBROUTINE OM4111
C                       *****************
C
     *(OP ,  DM,TYPDIM,XM,TYPEXM,   DN,TYPDIN,XN,TYPEXN,   C,
     * SIZDN,SZMDN,SIZXN,SZMXN,NETAGE, NELMAX3D)
C
C***********************************************************************
C BIEF VERSION 5.1           06/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FONCTION : OPERATIONS SUR DES MATRICES   M: TRIANGLE P1
C                                          N: MATRICE DE BORD
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES MATRICES M ET N ,LA MATRICE DIAGONALE D ET LA
C   CONSTANTE C.
C
C   LE RESULTAT EST LA MATRICE M.
C
C      OP = 'M=M+N   '  : ON AJOUTE N A M
C      OP = 'M=M+TN  '  : ON AJOUTE TRANSPOSEE DE N A M
C
C-----------------------------------------------------------------------
C
C  CONVENTION POUR LE STOCKAGE DES TERMES EXTRA-DIAGONAUX :
C
C      XM(     ,1)  ---->  M(1,2)
C      XM(     ,2)  ---->  M(1,3)
C      XM(     ,3)  ---->  M(2,3)
C      XM(     ,4)  ---->  M(2,1)
C      XM(     ,5)  ---->  M(3,1)
C      XM(     ,6)  ---->  M(3,2)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C |    OP          | -->| OPERATION A EFFECTUER
C |    DM,TYPDIM   |<-->| DIAGONALE ET TYPE DE DIAGONALE DE M
C |    XM,TYPEXM   | -->| TERMES EXTRA-DIAG. ET TYPE POUR M
C |    DN,TYPDIN   | -->| DIAGONALE ET TYPE DE DIAGONALE DE N
C |    XN,TYPEXN   | -->| TERMES EXTRA-DIAG. ET TYPE POUR N
C |    D           | -->| MATRICE DIAGONALE
C |    C           | -->| CONSTANTE DONNEE
C |    SIZDN       | -->| NOMBRE DE POINTS DU MAILLAGE 2D
C |    SZMDN       | -->| NOMBRE MAXIMUM DE POINTS DU MAILLAGE 2D
C |                |    | (CAS D'UN MAILLAGE ADAPTATIF)
C |    SIZXN       | -->| NOMBRE D'ELEMENTS DU MAILLAGE 2D
C |    SZMXN       | -->| NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE 2D
C |                |    | (CAS D'UN MAILLAGE ADAPTATIF)
C |    NETAGE      | -->| NOMBRE D'ETAGES DU MAILLAGE 3D
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   : OV
C
C***********************************************************************
C
      USE BIEF, EX_OM4111 => OM4111
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NETAGE,SIZDN,SZMXN,SZMDN,SIZXN
      integer, intent(in)             :: nelmax3d
      CHARACTER(LEN=8), INTENT(IN)    :: OP
ccc      DOUBLE PRECISION, INTENT(IN)    :: DN(*),XN(SZMXN,*)
      DOUBLE PRECISION, INTENT(IN)    :: DN(*),XN(NELMAX3D/NETAGE,*)
      DOUBLE PRECISION, INTENT(INOUT) :: DM(SZMDN,*)
ccc      DOUBLE PRECISION, INTENT(INOUT) :: XM(SZMXN,NETAGE,*)
      DOUBLE PRECISION, INTENT(INOUT) :: XM(NELMAX3D/NETAGE,NETAGE,*)
      CHARACTER(LEN=1), INTENT(INOUT) :: TYPDIM,TYPEXM,TYPDIN,TYPEXN
      DOUBLE PRECISION, INTENT(IN)    :: C
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K
C
      DOUBLE PRECISION Z(1)
C
C-----------------------------------------------------------------------
C
      IF(OP(1:8).EQ.'M=M+NF  ') THEN
C
        IF(TYPDIM.EQ.'Q'.AND.TYPDIN.EQ.'Q') THEN
          CALL OV( 'X=X+Y   ' , DM , DN , Z , C , SIZDN )
        ELSE
          IF (LNG.EQ.1) WRITE(LU,198) TYPDIM(1:1),OP(1:8),TYPDIN(1:1)
          IF (LNG.EQ.2) WRITE(LU,199) TYPDIM(1:1),OP(1:8),TYPDIN(1:1)
198       FORMAT(1X,'OM4111 (BIEF) : TYPDIM = ',A1,' NON PROGRAMME',
     *      /,1X,'POUR L''OPERATION : ',A8,' AVEC TYPDIN = ',A1)
199       FORMAT(1X,'OM4111 (BIEF) : TYPDIM = ',A1,' NOT IMPLEMENTED',
     *      /,1X,'FOR THE OPERATION : ',A8,' WITH TYPDIN = ',A1)
          CALL PLANTE(1)
          STOP
        ENDIF
C
        IF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'Q') THEN
C
C          CAS OU LES DEUX MATRICES SONT NON SYMETRIQUES
C
           DO 10 K = 1 , SIZXN
             XM(K,1, 1) = XM(K,1, 1) + XN(K,1)
             XM(K,1, 2) = XM(K,1, 2) + XN(K,2)
             XM(K,1, 6) = XM(K,1, 6) + XN(K,3)
             XM(K,1,16) = XM(K,1,16) + XN(K,4)
             XM(K,1,17) = XM(K,1,17) + XN(K,5)
             XM(K,1,21) = XM(K,1,21) + XN(K,6)
10         CONTINUE
C
        ELSEIF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU M EST QUELCONQUE ET N SYMETRIQUE
C
           DO 20 K = 1 , SIZXN
             XM(K,1, 1) = XM(K,1, 1) + XN(K,1)
             XM(K,1, 2) = XM(K,1, 2) + XN(K,2)
             XM(K,1, 6) = XM(K,1, 6) + XN(K,3)
             XM(K,1,16) = XM(K,1,16) + XN(K,1)
             XM(K,1,17) = XM(K,1,17) + XN(K,2)
             XM(K,1,21) = XM(K,1,21) + XN(K,3)
20         CONTINUE
C
        ELSEIF(TYPEXM(1:1).EQ.'S'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU LES DEUX MATRICES SONT SYMETRIQUES
C
           DO 30 K = 1 , SIZXN
             XM(K,1, 1) = XM(K,1, 1) + XN(K,1)
             XM(K,1, 2) = XM(K,1, 2) + XN(K,2)
             XM(K,1, 6) = XM(K,1, 6) + XN(K,3)
30         CONTINUE
C
        ELSE
           IF (LNG.EQ.1) WRITE(LU,98) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,99) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
98         FORMAT(1X,'OM4111 (BIEF) : TYPEXM = ',A1,' NE CONVIENT PAS',
     *       /,1X,'POUR L''OPERATION : ',A8,' AVEC TYPEXN = ',A1)
99         FORMAT(1X,'OM4111 (BIEF) : TYPEXM = ',A1,' DOES NOT GO',
     *       /,1X,'FOR THE OPERATION : ',A8,' WITH TYPEXN = ',A1)
           CALL PLANTE(1)
           STOP
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'M=M+TNF ') THEN
C
        CALL OV( 'X=X+Y   ' , DM , DN , Z , C , SIZDN )
C
        IF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'Q') THEN
C
C          CAS OU LES DEUX MATRICES SONT NON SYMETRIQUES
C
           DO 40 K = 1 , SIZXN
             XM(K,1, 1) = XM(K,1, 1) + XN(K,4)
             XM(K,1, 2) = XM(K,1, 2) + XN(K,5)
             XM(K,1, 6) = XM(K,1, 6) + XN(K,6)
             XM(K,1,16) = XM(K,1,16) + XN(K,1)
             XM(K,1,17) = XM(K,1,17) + XN(K,2)
             XM(K,1,21) = XM(K,1,21) + XN(K,3)
40         CONTINUE
C
        ELSEIF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU M EST QUELCONQUE ET N SYMETRIQUE
C
           DO 50 K = 1 , SIZXN
             XM(K,1, 1) = XM(K,1, 1) + XN(K,1)
             XM(K,1, 2) = XM(K,1, 2) + XN(K,2)
             XM(K,1, 6) = XM(K,1, 6) + XN(K,3)
             XM(K,1,16) = XM(K,1,16) + XN(K,1)
             XM(K,1,17) = XM(K,1,17) + XN(K,2)
             XM(K,1,21) = XM(K,1,21) + XN(K,3)
50         CONTINUE
C
        ELSEIF(TYPEXM(1:1).EQ.'S'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU LES DEUX MATRICES SONT SYMETRIQUES
C
           DO 60 K = 1 , SIZXN
             XM(K,1, 1) = XM(K,1, 1) + XN(K,1)
             XM(K,1, 2) = XM(K,1, 2) + XN(K,2)
             XM(K,1, 6) = XM(K,1, 6) + XN(K,3)
60         CONTINUE
C
        ELSE
           IF (LNG.EQ.1) WRITE(LU,98) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,99) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           CALL PLANTE(1)
           STOP
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'M=M+NS  ') THEN
C
        CALL OV( 'X=X+Y   ' , DM(1,NETAGE+1) , DN , Z , C , SIZDN )
C
        IF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'Q') THEN
C
C          CAS OU LES DEUX MATRICES SONT NON SYMETRIQUES
C
           DO 110 K = 1 , SIZXN
             XM(K,NETAGE,13) = XM(K,NETAGE,13) + XN(K,1)
             XM(K,NETAGE,14) = XM(K,NETAGE,14) + XN(K,2)
             XM(K,NETAGE,15) = XM(K,NETAGE,15) + XN(K,3)
             XM(K,NETAGE,28) = XM(K,NETAGE,28) + XN(K,4)
             XM(K,NETAGE,29) = XM(K,NETAGE,29) + XN(K,5)
             XM(K,NETAGE,30) = XM(K,NETAGE,30) + XN(K,6)
110        CONTINUE
C
        ELSEIF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU M EST QUELCONQUE ET N SYMETRIQUE
C
           DO 120 K = 1 , SIZXN
             XM(K,NETAGE,13) = XM(K,NETAGE,13) + XN(K,1)
             XM(K,NETAGE,14) = XM(K,NETAGE,14) + XN(K,2)
             XM(K,NETAGE,15) = XM(K,NETAGE,15) + XN(K,3)
             XM(K,NETAGE,28) = XM(K,NETAGE,28) + XN(K,1)
             XM(K,NETAGE,29) = XM(K,NETAGE,29) + XN(K,2)
             XM(K,NETAGE,30) = XM(K,NETAGE,30) + XN(K,3)
120        CONTINUE
C
        ELSEIF(TYPEXM(1:1).EQ.'S'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU LES DEUX MATRICES SONT SYMETRIQUES
C
           DO 130 K = 1 , SIZXN
             XM(K,NETAGE,13) = XM(K,NETAGE,13) + XN(K,1)
             XM(K,NETAGE,14) = XM(K,NETAGE,14) + XN(K,2)
             XM(K,NETAGE,15) = XM(K,NETAGE,15) + XN(K,3)
130        CONTINUE
C
        ELSE
           IF (LNG.EQ.1) WRITE(LU,98) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,99) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           CALL PLANTE(1)
           STOP
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'M=M+TNS ') THEN
C
        CALL OV( 'X=X+Y   ' , DM(1,NETAGE+1) , DN , Z , C , SIZDN )
C
        IF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'Q') THEN
C
C          CAS OU LES DEUX MATRICES SONT NON SYMETRIQUES
C
           DO 140 K = 1 , SIZXN
             XM(K,NETAGE,13) = XM(K,NETAGE,13) + XN(K,4)
             XM(K,NETAGE,14) = XM(K,NETAGE,14) + XN(K,5)
             XM(K,NETAGE,15) = XM(K,NETAGE,15) + XN(K,6)
             XM(K,NETAGE,28) = XM(K,NETAGE,28) + XN(K,1)
             XM(K,NETAGE,29) = XM(K,NETAGE,29) + XN(K,2)
             XM(K,NETAGE,30) = XM(K,NETAGE,30) + XN(K,3)
140        CONTINUE
C
        ELSEIF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU M EST QUELCONQUE ET N SYMETRIQUE
C
           DO 150 K = 1 , SIZXN
             XM(K,NETAGE,13) = XM(K,NETAGE,13) + XN(K,1)
             XM(K,NETAGE,14) = XM(K,NETAGE,14) + XN(K,2)
             XM(K,NETAGE,15) = XM(K,NETAGE,15) + XN(K,3)
             XM(K,NETAGE,28) = XM(K,NETAGE,28) + XN(K,1)
             XM(K,NETAGE,29) = XM(K,NETAGE,29) + XN(K,2)
             XM(K,NETAGE,30) = XM(K,NETAGE,30) + XN(K,3)
150        CONTINUE
C
        ELSEIF(TYPEXM(1:1).EQ.'S'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU LES DEUX MATRICES SONT SYMETRIQUES
C
           DO 160 K = 1 , SIZXN
             XM(K,NETAGE,13) = XM(K,NETAGE,13) + XN(K,1)
             XM(K,NETAGE,14) = XM(K,NETAGE,14) + XN(K,2)
             XM(K,NETAGE,15) = XM(K,NETAGE,15) + XN(K,3)
160        CONTINUE
C
        ELSE
           IF (LNG.EQ.1) WRITE(LU,98) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,99) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           CALL PLANTE(1)
           STOP
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,70) OP
        IF (LNG.EQ.2) WRITE(LU,71) OP
70      FORMAT(1X,'OM4111 (BIEF) : OPERATION INCONNUE : ',A8)
71      FORMAT(1X,'OM4111 (BIEF) : UNKNOWN OPERATION : ',A8)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
