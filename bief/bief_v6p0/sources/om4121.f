C                       *****************
                        SUBROUTINE OM4121
C                       *****************
C
     *(OP ,  DM,TYPDIM,XM,TYPEXM,   DN,TYPDIN,XN,TYPEXN,   C,
     * NULONE,NELBOR,NBOR,NELMAX,SIZDN,SIZXN,SZMXN)
C
C***********************************************************************
C BIEF VERSION 5.9      23/06/2008    J-M HERVOUET (LNHE) 01 30 87 80 18
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
C |    NULONE      | -->| NUMEROS LOCAUX DES NOEUDS DE BORD.
C |    NELBOR      | -->| NUMEROS DES ELEMENTS DE BORD.
C |    NBOR        | -->| NUMEROS GLOBAUX DES POINTS DE BORD.
C |    IKLE        | -->| CORRESPONDANCE NUMEROTATIONS LOCALE ET GLOBALE
C |    NELEM       | -->| NOMBRE D'ELEMENTS DU MAILLAGE
C |    NELMAX      | -->| NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    | (CAS D'UN MAILLAGE ADAPTATIF)
C |    SIZXN       | -->| NOMBRE D'ELEMENTS DE BORD.
C |    SZMXN       | -->| NOMBRE MAXIMUM D'ELEMENTS DE BORD.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   : OV
C
C***********************************************************************
C
      USE BIEF, EX_OM4121 => OM4121
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NELMAX,SIZDN,SIZXN,SZMXN
      CHARACTER(LEN=8), INTENT(IN)    :: OP
      INTEGER, INTENT(IN)             :: NULONE(*),NELBOR(*),NBOR(*)
      DOUBLE PRECISION, INTENT(IN)    :: DN(*),XN(SZMXN,*)
      DOUBLE PRECISION, INTENT(INOUT) :: DM(*),XM(NELMAX,*)
      CHARACTER(LEN=1), INTENT(INOUT) :: TYPDIM,TYPEXM,TYPDIN,TYPEXN
      DOUBLE PRECISION, INTENT(IN)    :: C
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K,IEL
C
      DOUBLE PRECISION Z(1)
C
      INTEGER CORNSY(3,12),CORSYM(3,6)
C
C-----------------------------------------------------------------------
C
      DATA CORNSY/ 1, 6,17, 4, 9,10, 3, 8,12, 8,12, 3, 7,11, 5,28,30,14,
     *            16,21, 2,19,24,25,18,23,27,23,27,18,22,26,20,13,15,29/
      DATA CORSYM/ 1, 6, 2, 4, 9,10, 3, 8,12, 8,12, 3, 7,11, 5,13,15,14/
C
C***********************************************************************
C
      IF(OP(1:8).EQ.'M=M+N   ') THEN
C
        IF(TYPDIM.EQ.'Q'.AND.TYPDIN.EQ.'Q') THEN
          CALL OVDB( 'X=X+Y   ' , DM , DN , Z , C , NBOR , SIZDN )
        ELSE
          IF (LNG.EQ.1) WRITE(LU,198) TYPDIM(1:1),OP(1:8),TYPDIN(1:1)
          IF (LNG.EQ.2) WRITE(LU,199) TYPDIM(1:1),OP(1:8),TYPDIN(1:1)
198       FORMAT(1X,'OM4121 (BIEF) : TYPDIM = ',A1,' NON PROGRAMME',
     *      /,1X,'POUR L''OPERATION : ',A8,' AVEC TYPDIN = ',A1)
199       FORMAT(1X,'OM4121 (BIEF) : TYPDIM = ',A1,' NOT IMPLEMENTED',
     *      /,1X,'FOR THE OPERATION : ',A8,' WITH TYPDIN = ',A1)
          CALL PLANTE(1)
          STOP
        ENDIF   
C
        IF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'Q') THEN
C
C          CAS OU LES DEUX MATRICES SONT NON SYMETRIQUES
C
           IF(NCSIZE.GT.1) THEN
           DO K = 1 , SIZXN
             IEL = NELBOR(K)
             IF(IEL.GT.0) THEN
             XM( IEL , CORNSY(NULONE(K), 1) ) =
     *       XM( IEL , CORNSY(NULONE(K), 1) ) + XN(K, 1)
             XM( IEL , CORNSY(NULONE(K), 2) ) =
     *       XM( IEL , CORNSY(NULONE(K), 2) ) + XN(K, 2)
             XM( IEL , CORNSY(NULONE(K), 3) ) =
     *       XM( IEL , CORNSY(NULONE(K), 3) ) + XN(K, 3)
             XM( IEL , CORNSY(NULONE(K), 4) ) =
     *       XM( IEL , CORNSY(NULONE(K), 4) ) + XN(K, 4)
             XM( IEL , CORNSY(NULONE(K), 5) ) =
     *       XM( IEL , CORNSY(NULONE(K), 5) ) + XN(K, 5)
             XM( IEL , CORNSY(NULONE(K), 6) ) =
     *       XM( IEL , CORNSY(NULONE(K), 6) ) + XN(K, 6)
             XM( IEL , CORNSY(NULONE(K), 7) ) =
     *       XM( IEL , CORNSY(NULONE(K), 7) ) + XN(K, 7)
             XM( IEL , CORNSY(NULONE(K), 8) ) =
     *       XM( IEL , CORNSY(NULONE(K), 8) ) + XN(K, 8)
             XM( IEL , CORNSY(NULONE(K), 9) ) =
     *       XM( IEL , CORNSY(NULONE(K), 9) ) + XN(K, 9)
             XM( IEL , CORNSY(NULONE(K),10) ) =
     *       XM( IEL , CORNSY(NULONE(K),10) ) + XN(K,10)
             XM( IEL , CORNSY(NULONE(K),11) ) =
     *       XM( IEL , CORNSY(NULONE(K),11) ) + XN(K,11)
             XM( IEL , CORNSY(NULONE(K),12) ) =
     *       XM( IEL , CORNSY(NULONE(K),12) ) + XN(K,12)
             ENDIF
           ENDDO
           ELSE
           DO K = 1 , SIZXN
             IEL = NELBOR(K)
             XM( IEL , CORNSY(NULONE(K), 1) ) =
     *       XM( IEL , CORNSY(NULONE(K), 1) ) + XN(K, 1)
             XM( IEL , CORNSY(NULONE(K), 2) ) =
     *       XM( IEL , CORNSY(NULONE(K), 2) ) + XN(K, 2)
             XM( IEL , CORNSY(NULONE(K), 3) ) =
     *       XM( IEL , CORNSY(NULONE(K), 3) ) + XN(K, 3)
             XM( IEL , CORNSY(NULONE(K), 4) ) =
     *       XM( IEL , CORNSY(NULONE(K), 4) ) + XN(K, 4)
             XM( IEL , CORNSY(NULONE(K), 5) ) =
     *       XM( IEL , CORNSY(NULONE(K), 5) ) + XN(K, 5)
             XM( IEL , CORNSY(NULONE(K), 6) ) =
     *       XM( IEL , CORNSY(NULONE(K), 6) ) + XN(K, 6)
             XM( IEL , CORNSY(NULONE(K), 7) ) =
     *       XM( IEL , CORNSY(NULONE(K), 7) ) + XN(K, 7)
             XM( IEL , CORNSY(NULONE(K), 8) ) =
     *       XM( IEL , CORNSY(NULONE(K), 8) ) + XN(K, 8)
             XM( IEL , CORNSY(NULONE(K), 9) ) =
     *       XM( IEL , CORNSY(NULONE(K), 9) ) + XN(K, 9)
             XM( IEL , CORNSY(NULONE(K),10) ) =
     *       XM( IEL , CORNSY(NULONE(K),10) ) + XN(K,10)
             XM( IEL , CORNSY(NULONE(K),11) ) =
     *       XM( IEL , CORNSY(NULONE(K),11) ) + XN(K,11)
             XM( IEL , CORNSY(NULONE(K),12) ) =
     *       XM( IEL , CORNSY(NULONE(K),12) ) + XN(K,12)
           ENDDO
           ENDIF
C
        ELSEIF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU M EST QUELCONQUE ET N SYMETRIQUE
C
           IF(NCSIZE.GT.1) THEN
           DO K = 1 , SIZXN
             IEL = NELBOR(K)
             IF(IEL.GT.0) THEN
             XM( IEL , CORNSY(NULONE(K), 1) ) =
     *       XM( IEL , CORNSY(NULONE(K), 1) ) + XN(K, 1)
             XM( IEL , CORNSY(NULONE(K), 2) ) =
     *       XM( IEL , CORNSY(NULONE(K), 2) ) + XN(K, 2)
             XM( IEL , CORNSY(NULONE(K), 3) ) =
     *       XM( IEL , CORNSY(NULONE(K), 3) ) + XN(K, 3)
             XM( IEL , CORNSY(NULONE(K), 4) ) =
     *       XM( IEL , CORNSY(NULONE(K), 4) ) + XN(K, 4)
             XM( IEL , CORNSY(NULONE(K), 5) ) =
     *       XM( IEL , CORNSY(NULONE(K), 5) ) + XN(K, 5)
             XM( IEL , CORNSY(NULONE(K), 6) ) =
     *       XM( IEL , CORNSY(NULONE(K), 6) ) + XN(K, 6)
             XM( IEL , CORNSY(NULONE(K), 7) ) =
     *       XM( IEL , CORNSY(NULONE(K), 7) ) + XN(K, 1)
             XM( IEL , CORNSY(NULONE(K), 8) ) =
     *       XM( IEL , CORNSY(NULONE(K), 8) ) + XN(K, 2)
             XM( IEL , CORNSY(NULONE(K), 9) ) =
     *       XM( IEL , CORNSY(NULONE(K), 9) ) + XN(K, 3)
             XM( IEL , CORNSY(NULONE(K),10) ) =
     *       XM( IEL , CORNSY(NULONE(K),10) ) + XN(K, 4)
             XM( IEL , CORNSY(NULONE(K),11) ) =
     *       XM( IEL , CORNSY(NULONE(K),11) ) + XN(K, 5)
             XM( IEL , CORNSY(NULONE(K),12) ) =
     *       XM( IEL , CORNSY(NULONE(K),12) ) + XN(K, 6)
             ENDIF
           ENDDO
           ELSE
           DO K = 1 , SIZXN
             IEL = NELBOR(K)
             XM( IEL , CORNSY(NULONE(K), 1) ) =
     *       XM( IEL , CORNSY(NULONE(K), 1) ) + XN(K, 1)
             XM( IEL , CORNSY(NULONE(K), 2) ) =
     *       XM( IEL , CORNSY(NULONE(K), 2) ) + XN(K, 2)
             XM( IEL , CORNSY(NULONE(K), 3) ) =
     *       XM( IEL , CORNSY(NULONE(K), 3) ) + XN(K, 3)
             XM( IEL , CORNSY(NULONE(K), 4) ) =
     *       XM( IEL , CORNSY(NULONE(K), 4) ) + XN(K, 4)
             XM( IEL , CORNSY(NULONE(K), 5) ) =
     *       XM( IEL , CORNSY(NULONE(K), 5) ) + XN(K, 5)
             XM( IEL , CORNSY(NULONE(K), 6) ) =
     *       XM( IEL , CORNSY(NULONE(K), 6) ) + XN(K, 6)
             XM( IEL , CORNSY(NULONE(K), 7) ) =
     *       XM( IEL , CORNSY(NULONE(K), 7) ) + XN(K, 1)
             XM( IEL , CORNSY(NULONE(K), 8) ) =
     *       XM( IEL , CORNSY(NULONE(K), 8) ) + XN(K, 2)
             XM( IEL , CORNSY(NULONE(K), 9) ) =
     *       XM( IEL , CORNSY(NULONE(K), 9) ) + XN(K, 3)
             XM( IEL , CORNSY(NULONE(K),10) ) =
     *       XM( IEL , CORNSY(NULONE(K),10) ) + XN(K, 4)
             XM( IEL , CORNSY(NULONE(K),11) ) =
     *       XM( IEL , CORNSY(NULONE(K),11) ) + XN(K, 5)
             XM( IEL , CORNSY(NULONE(K),12) ) =
     *       XM( IEL , CORNSY(NULONE(K),12) ) + XN(K, 6)
           ENDDO
           ENDIF
C
        ELSEIF(TYPEXM(1:1).EQ.'S'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU LES DEUX MATRICES SONT SYMETRIQUES
C
           IF(NCSIZE.GT.1) THEN
           DO K = 1 , SIZXN
             IEL = NELBOR(K)
             IF(IEL.GT.0) THEN
             XM( IEL , CORSYM(NULONE(K), 1) ) =
     *       XM( IEL , CORSYM(NULONE(K), 1) ) + XN(K, 1)
             XM( IEL , CORSYM(NULONE(K), 2) ) =
     *       XM( IEL , CORSYM(NULONE(K), 2) ) + XN(K, 2)
             XM( IEL , CORSYM(NULONE(K), 3) ) =
     *       XM( IEL , CORSYM(NULONE(K), 3) ) + XN(K, 3)
             XM( IEL , CORSYM(NULONE(K), 4) ) =
     *       XM( IEL , CORSYM(NULONE(K), 4) ) + XN(K, 4)
             XM( IEL , CORSYM(NULONE(K), 5) ) =
     *       XM( IEL , CORSYM(NULONE(K), 5) ) + XN(K, 5)
             XM( IEL , CORSYM(NULONE(K), 6) ) =
     *       XM( IEL , CORSYM(NULONE(K), 6) ) + XN(K, 6)
             ENDIF
           ENDDO
           ELSE
           DO K = 1 , SIZXN
             IEL = NELBOR(K)
             XM( IEL , CORSYM(NULONE(K), 1) ) =
     *       XM( IEL , CORSYM(NULONE(K), 1) ) + XN(K, 1)
             XM( IEL , CORSYM(NULONE(K), 2) ) =
     *       XM( IEL , CORSYM(NULONE(K), 2) ) + XN(K, 2)
             XM( IEL , CORSYM(NULONE(K), 3) ) =
     *       XM( IEL , CORSYM(NULONE(K), 3) ) + XN(K, 3)
             XM( IEL , CORSYM(NULONE(K), 4) ) =
     *       XM( IEL , CORSYM(NULONE(K), 4) ) + XN(K, 4)
             XM( IEL , CORSYM(NULONE(K), 5) ) =
     *       XM( IEL , CORSYM(NULONE(K), 5) ) + XN(K, 5)
             XM( IEL , CORSYM(NULONE(K), 6) ) =
     *       XM( IEL , CORSYM(NULONE(K), 6) ) + XN(K, 6)
           ENDDO
           ENDIF
C
        ELSE
           IF (LNG.EQ.1) WRITE(LU,98) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,99) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
98         FORMAT(1X,'OM4121 (BIEF) : TYPEXM = ',A1,' NE CONVIENT PAS',
     *       /,1X,'POUR L''OPERATION : ',A8,' AVEC TYPEXN = ',A1)
99         FORMAT(1X,'OM4121 (BIEF) : TYPEXM = ',A1,' DOES NOT GO',
     *       /,1X,'FOR THE OPERATION : ',A8,' WITH TYPEXN = ',A1)
           CALL PLANTE(1)
           STOP
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'M=M+TN  ') THEN
C
        CALL OVDB( 'X=X+Y   ' , DM , DN , Z , C , NBOR , SIZXN )
C
        IF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'Q') THEN
C
C          CAS OU LES DEUX MATRICES SONT NON SYMETRIQUES
C
           IF(NCSIZE.GT.1) THEN
           DO K = 1 , SIZXN
             IEL = NELBOR(K)
             IF(IEL.GT.0) THEN
             XM( IEL , CORNSY(NULONE(K), 1) ) =
     *       XM( IEL , CORNSY(NULONE(K), 1) ) + XN(K, 7)
             XM( IEL , CORNSY(NULONE(K), 2) ) =
     *       XM( IEL , CORNSY(NULONE(K), 2) ) + XN(K, 8)
             XM( IEL , CORNSY(NULONE(K), 3) ) =
     *       XM( IEL , CORNSY(NULONE(K), 3) ) + XN(K, 9)
             XM( IEL , CORNSY(NULONE(K), 4) ) =
     *       XM( IEL , CORNSY(NULONE(K), 4) ) + XN(K,10)
             XM( IEL , CORNSY(NULONE(K), 5) ) =
     *       XM( IEL , CORNSY(NULONE(K), 5) ) + XN(K,11)
             XM( IEL , CORNSY(NULONE(K), 6) ) =
     *       XM( IEL , CORNSY(NULONE(K), 6) ) + XN(K,12)
             XM( IEL , CORNSY(NULONE(K), 7) ) =
     *       XM( IEL , CORNSY(NULONE(K), 7) ) + XN(K, 1)
             XM( IEL , CORNSY(NULONE(K), 8) ) =
     *       XM( IEL , CORNSY(NULONE(K), 8) ) + XN(K, 2)
             XM( IEL , CORNSY(NULONE(K), 9) ) =
     *       XM( IEL , CORNSY(NULONE(K), 9) ) + XN(K, 3)
             XM( IEL , CORNSY(NULONE(K),10) ) =
     *       XM( IEL , CORNSY(NULONE(K),10) ) + XN(K, 4)
             XM( IEL , CORNSY(NULONE(K),11) ) =
     *       XM( IEL , CORNSY(NULONE(K),11) ) + XN(K, 5)
             XM( IEL , CORNSY(NULONE(K),12) ) =
     *       XM( IEL , CORNSY(NULONE(K),12) ) + XN(K, 6)
             ENDIF
           ENDDO
           ELSE
           DO K = 1 , SIZXN
             IEL = NELBOR(K)
             XM( IEL , CORNSY(NULONE(K), 1) ) =
     *       XM( IEL , CORNSY(NULONE(K), 1) ) + XN(K, 7)
             XM( IEL , CORNSY(NULONE(K), 2) ) =
     *       XM( IEL , CORNSY(NULONE(K), 2) ) + XN(K, 8)
             XM( IEL , CORNSY(NULONE(K), 3) ) =
     *       XM( IEL , CORNSY(NULONE(K), 3) ) + XN(K, 9)
             XM( IEL , CORNSY(NULONE(K), 4) ) =
     *       XM( IEL , CORNSY(NULONE(K), 4) ) + XN(K,10)
             XM( IEL , CORNSY(NULONE(K), 5) ) =
     *       XM( IEL , CORNSY(NULONE(K), 5) ) + XN(K,11)
             XM( IEL , CORNSY(NULONE(K), 6) ) =
     *       XM( IEL , CORNSY(NULONE(K), 6) ) + XN(K,12)
             XM( IEL , CORNSY(NULONE(K), 7) ) =
     *       XM( IEL , CORNSY(NULONE(K), 7) ) + XN(K, 1)
             XM( IEL , CORNSY(NULONE(K), 8) ) =
     *       XM( IEL , CORNSY(NULONE(K), 8) ) + XN(K, 2)
             XM( IEL , CORNSY(NULONE(K), 9) ) =
     *       XM( IEL , CORNSY(NULONE(K), 9) ) + XN(K, 3)
             XM( IEL , CORNSY(NULONE(K),10) ) =
     *       XM( IEL , CORNSY(NULONE(K),10) ) + XN(K, 4)
             XM( IEL , CORNSY(NULONE(K),11) ) =
     *       XM( IEL , CORNSY(NULONE(K),11) ) + XN(K, 5)
             XM( IEL , CORNSY(NULONE(K),12) ) =
     *       XM( IEL , CORNSY(NULONE(K),12) ) + XN(K, 6)
           ENDDO
           ENDIF
C
        ELSEIF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU M EST QUELCONQUE ET N SYMETRIQUE
C
           IF(NCSIZE.GT.1) THEN
           DO K = 1 , SIZXN
             IEL = NELBOR(K)
             IF(IEL.GT.0) THEN
             XM( IEL , CORNSY(NULONE(K), 1) ) =
     *       XM( IEL , CORNSY(NULONE(K), 1) ) + XN(K, 1)
             XM( IEL , CORNSY(NULONE(K), 2) ) =
     *       XM( IEL , CORNSY(NULONE(K), 2) ) + XN(K, 2)
             XM( IEL , CORNSY(NULONE(K), 3) ) =
     *       XM( IEL , CORNSY(NULONE(K), 3) ) + XN(K, 3)
             XM( IEL , CORNSY(NULONE(K), 4) ) =
     *       XM( IEL , CORNSY(NULONE(K), 4) ) + XN(K, 4)
             XM( IEL , CORNSY(NULONE(K), 5) ) =
     *       XM( IEL , CORNSY(NULONE(K), 5) ) + XN(K, 5)
             XM( IEL , CORNSY(NULONE(K), 6) ) =
     *       XM( IEL , CORNSY(NULONE(K), 6) ) + XN(K, 6)
             XM( IEL , CORNSY(NULONE(K), 7) ) =
     *       XM( IEL , CORNSY(NULONE(K), 7) ) + XN(K, 1)
             XM( IEL , CORNSY(NULONE(K), 8) ) =
     *       XM( IEL , CORNSY(NULONE(K), 8) ) + XN(K, 2)
             XM( IEL , CORNSY(NULONE(K), 9) ) =
     *       XM( IEL , CORNSY(NULONE(K), 9) ) + XN(K, 3)
             XM( IEL , CORNSY(NULONE(K),10) ) =
     *       XM( IEL , CORNSY(NULONE(K),10) ) + XN(K, 4)
             XM( IEL , CORNSY(NULONE(K),11) ) =
     *       XM( IEL , CORNSY(NULONE(K),11) ) + XN(K, 5)
             XM( IEL , CORNSY(NULONE(K),12) ) =
     *       XM( IEL , CORNSY(NULONE(K),12) ) + XN(K, 6)
             ENDIF
           ENDDO
           ELSE
           DO K = 1 , SIZXN
             IEL = NELBOR(K)
             XM( IEL , CORNSY(NULONE(K), 1) ) =
     *       XM( IEL , CORNSY(NULONE(K), 1) ) + XN(K, 1)
             XM( IEL , CORNSY(NULONE(K), 2) ) =
     *       XM( IEL , CORNSY(NULONE(K), 2) ) + XN(K, 2)
             XM( IEL , CORNSY(NULONE(K), 3) ) =
     *       XM( IEL , CORNSY(NULONE(K), 3) ) + XN(K, 3)
             XM( IEL , CORNSY(NULONE(K), 4) ) =
     *       XM( IEL , CORNSY(NULONE(K), 4) ) + XN(K, 4)
             XM( IEL , CORNSY(NULONE(K), 5) ) =
     *       XM( IEL , CORNSY(NULONE(K), 5) ) + XN(K, 5)
             XM( IEL , CORNSY(NULONE(K), 6) ) =
     *       XM( IEL , CORNSY(NULONE(K), 6) ) + XN(K, 6)
             XM( IEL , CORNSY(NULONE(K), 7) ) =
     *       XM( IEL , CORNSY(NULONE(K), 7) ) + XN(K, 1)
             XM( IEL , CORNSY(NULONE(K), 8) ) =
     *       XM( IEL , CORNSY(NULONE(K), 8) ) + XN(K, 2)
             XM( IEL , CORNSY(NULONE(K), 9) ) =
     *       XM( IEL , CORNSY(NULONE(K), 9) ) + XN(K, 3)
             XM( IEL , CORNSY(NULONE(K),10) ) =
     *       XM( IEL , CORNSY(NULONE(K),10) ) + XN(K, 4)
             XM( IEL , CORNSY(NULONE(K),11) ) =
     *       XM( IEL , CORNSY(NULONE(K),11) ) + XN(K, 5)
             XM( IEL , CORNSY(NULONE(K),12) ) =
     *       XM( IEL , CORNSY(NULONE(K),12) ) + XN(K, 6)
           ENDDO
           ENDIF
C
        ELSEIF(TYPEXM(1:1).EQ.'S'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU LES DEUX MATRICES SONT SYMETRIQUES
C
           IF(NCSIZE.GT.1) THEN
           DO K = 1 , SIZXN
             IEL = NELBOR(K)
             IF(IEL.GT.0) THEN
             XM( IEL , CORSYM(NULONE(K), 1) ) =
     *       XM( IEL , CORSYM(NULONE(K), 1) ) + XN(K, 1)
             XM( IEL , CORSYM(NULONE(K), 2) ) =
     *       XM( IEL , CORSYM(NULONE(K), 2) ) + XN(K, 2)
             XM( IEL , CORSYM(NULONE(K), 3) ) =
     *       XM( IEL , CORSYM(NULONE(K), 3) ) + XN(K, 3)
             XM( IEL , CORSYM(NULONE(K), 4) ) =
     *       XM( IEL , CORSYM(NULONE(K), 4) ) + XN(K, 4)
             XM( IEL , CORSYM(NULONE(K), 5) ) =
     *       XM( IEL , CORSYM(NULONE(K), 5) ) + XN(K, 5)
             XM( IEL , CORSYM(NULONE(K), 6) ) =
     *       XM( IEL , CORSYM(NULONE(K), 6) ) + XN(K, 6)
             ENDIF
           ENDDO
           ELSE
           DO K = 1 , SIZXN
             IEL = NELBOR(K)
             XM( IEL , CORSYM(NULONE(K), 1) ) =
     *       XM( IEL , CORSYM(NULONE(K), 1) ) + XN(K, 1)
             XM( IEL , CORSYM(NULONE(K), 2) ) =
     *       XM( IEL , CORSYM(NULONE(K), 2) ) + XN(K, 2)
             XM( IEL , CORSYM(NULONE(K), 3) ) =
     *       XM( IEL , CORSYM(NULONE(K), 3) ) + XN(K, 3)
             XM( IEL , CORSYM(NULONE(K), 4) ) =
     *       XM( IEL , CORSYM(NULONE(K), 4) ) + XN(K, 4)
             XM( IEL , CORSYM(NULONE(K), 5) ) =
     *       XM( IEL , CORSYM(NULONE(K), 5) ) + XN(K, 5)
             XM( IEL , CORSYM(NULONE(K), 6) ) =
     *       XM( IEL , CORSYM(NULONE(K), 6) ) + XN(K, 6)
           ENDDO
           ENDIF
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
70      FORMAT(1X,'OM4121 (BIEF) : OPERATION INCONNUE : ',A8)
71      FORMAT(1X,'OM4121 (BIEF) : UNKNOWN OPERATION : ',A8)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
