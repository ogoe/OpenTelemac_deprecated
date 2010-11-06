C                       *****************
                        SUBROUTINE OM1201
C                       *****************
C
     *(OP ,  DM,TYPDIM,XM,TYPEXM,   DN,TYPDIN,XN,TYPEXN,   C,
     * NULONE,NELBOR,NBOR,NELMAX,NDIAG,NPTFR,NPTFRX)
C
C***********************************************************************
C BIEF VERSION 5.9      23/06/2008    J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C FONCTION : OPERATIONS SUR DES MATRICES   M: TRIANGLE QUASI-BULLE
C                                          N: MATRICE DE BORD
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES MATRICES M ET N ,LA MATRICE DIAGONALE D ET LA
C   CONSTANTE C.
C
C   LE RESULTAT EST LA MATRICE M.
C
C      OP = 'M=M+N   '  : ON AJOUTE N A M
C
C-----------------------------------------------------------------------
C
C  CONVENTION POUR LE STOCKAGE DES TERMES EXTRA-DIAGONAUX :
C
C      XM(     , 1)  ---->  M(1,2)
C      XM(     , 2)  ---->  M(1,3)
C      XM(     , 3)  ---->  M(1,4)
C      XM(     , 4)  ---->  M(2,3)
C      XM(     , 5)  ---->  M(2,4)
C      XM(     , 6)  ---->  M(3,4)
C      XM(     , 7)  ---->  M(2,1)
C      XM(     , 8)  ---->  M(3,1)
C      XM(     , 9)  ---->  M(4,1)
C      XM(     ,10)  ---->  M(3,2)
C      XM(     ,11)  ---->  M(4,2)
C      XM(     ,12)  ---->  M(4,3)
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
C |    IKLE        | -->| CORRESPONDANCE NUMEROTATIONS LOCALE ET GLOBALE
C |    NELEM       | -->| NOMBRE D'ELEMENTS DU MAILLAGE
C |    NELMAX      | -->| NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    | (CAS D'UN MAILLAGE ADAPTATIF)
C |    NDIAG       | -->| NOMBRE DE VALEURS DE LA DIAGONALE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   : OV
C
C***********************************************************************
C
      USE BIEF, EX_OM1201 => OM1201
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NELMAX,NDIAG,NPTFR,NPTFRX
      CHARACTER(LEN=8), INTENT(IN)    :: OP
      INTEGER, INTENT(IN)             :: NULONE(*),NELBOR(*),NBOR(*)
      DOUBLE PRECISION, INTENT(IN)    :: DN(*),XN(*)
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
      INTEGER CORNSY(4,2),CORSYM(4)
C
C-----------------------------------------------------------------------
C     ATTENTION : CECI NE MARCHE QUE POUR QUASI-BULLE
      DATA CORNSY/ 1,4,8,0,  7,10,2,0/
      DATA CORSYM/ 1,4,2,0      /
C
C-----------------------------------------------------------------------
C
      IF(OP(1:8).EQ.'M=M+N   ') THEN
C
        IF(TYPDIM.EQ.'Q'.AND.TYPDIM.EQ.'Q'.AND.NDIAG.GE.NPTFR) THEN
          CALL OVDB( 'X=X+Y   ' , DM , DN , Z , C , NBOR , NPTFR )
        ELSE
          IF (LNG.EQ.1) WRITE(LU,198) TYPDIM(1:1),OP(1:8),TYPDIN(1:1)
          IF (LNG.EQ.2) WRITE(LU,199) TYPDIM(1:1),OP(1:8),TYPDIN(1:1)
198       FORMAT(1X,'OM1201 (BIEF) : TYPDIM = ',A1,' NON PROGRAMME',
     *      /,1X,'POUR L''OPERATION : ',A8,' AVEC TYPDIN = ',A1)
199       FORMAT(1X,'OM1201 (BIEF) : TYPDIM = ',A1,' NOT IMPLEMENTED',
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
             DO K = 1 , NPTFR
               IEL = NELBOR(K)
               IF(IEL.GT.0) THEN
                 XM( IEL , CORNSY(NULONE(K),1) ) =
     *           XM( IEL , CORNSY(NULONE(K),1) ) + XN(K)
                 XM( IEL , CORNSY(NULONE(K),2) ) =
     *           XM( IEL , CORNSY(NULONE(K),2) ) + XN(K+NPTFRX)
               ENDIF
             ENDDO
           ELSE
             DO K = 1 , NPTFR
               IEL = NELBOR(K)
               XM( IEL , CORNSY(NULONE(K),1) ) =
     *         XM( IEL , CORNSY(NULONE(K),1) ) + XN(K)
               XM( IEL , CORNSY(NULONE(K),2) ) =
     *         XM( IEL , CORNSY(NULONE(K),2) ) + XN(K+NPTFRX)
             ENDDO
           ENDIF
C
        ELSEIF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU M EST QUELCONQUE ET N SYMETRIQUE
C
           IF(NCSIZE.GT.1) THEN
             DO K = 1 , NPTFR
               IEL = NELBOR(K)
               IF(IEL.GT.0) THEN
                 XM( IEL , CORNSY(NULONE(K),1) ) =
     *           XM( IEL , CORNSY(NULONE(K),1) ) + XN(K)
                 XM( IEL , CORNSY(NULONE(K),2) ) =
     *           XM( IEL , CORNSY(NULONE(K),2) ) + XN(K)
               ENDIF
             ENDDO
           ELSE
             DO K = 1 , NPTFR
               IEL = NELBOR(K)
               XM( IEL , CORNSY(NULONE(K),1) ) =
     *         XM( IEL , CORNSY(NULONE(K),1) ) + XN(K)
               XM( IEL , CORNSY(NULONE(K),2) ) =
     *         XM( IEL , CORNSY(NULONE(K),2) ) + XN(K)
             ENDDO
           ENDIF
C
        ELSEIF(TYPEXM(1:1).EQ.'S'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU LES DEUX MATRICES SONT SYMETRIQUES
C
           IF(NCSIZE.GT.1) THEN
             DO K = 1 , NPTFR
               IEL = NELBOR(K)
               IF(IEL.GT.0) THEN
                 XM( IEL , CORSYM(NULONE(K)) ) =
     *           XM( IEL , CORSYM(NULONE(K)) ) + XN(K)
               ENDIF
             ENDDO
           ELSE
             DO K = 1 , NPTFR
               IEL = NELBOR(K)
               XM( IEL , CORSYM(NULONE(K)) ) =
     *         XM( IEL , CORSYM(NULONE(K)) ) + XN(K)
             ENDDO
           ENDIF
C
        ELSE
           IF (LNG.EQ.1) WRITE(LU,99) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,98) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
98         FORMAT(1X,'OM1201 (BIEF) : TYPEXM = ',A1,' NE CONVIENT PAS',
     *       /,1X,'POUR L''OPERATION : ',A8,' AVEC TYPEXN = ',A1)
99         FORMAT(1X,'OM1201 (BIEF) : TYPEXM = ',A1,' DOES NOT GO',
     *       /,1X,'FOR THE OPERATION : ',A8,' WITH TYPEXN = ',A1)
           CALL PLANTE(1)
           STOP
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'M=M+TN  ') THEN
C
        CALL OVDB( 'X=X+Y   ' , DM , DN , Z , C , NBOR , NPTFR )
C
        IF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'Q') THEN
C
C          CAS OU LES DEUX MATRICES SONT NON SYMETRIQUES
C
           IF(NCSIZE.GT.1) THEN
             DO K = 1 , NPTFR
               IEL = NELBOR(K)
               IF(IEL.GT.0) THEN
                 XM( IEL , CORNSY(NULONE(K),1) ) =
     *           XM( IEL , CORNSY(NULONE(K),1) ) + XN(K+NPTFRX)
                 XM( IEL , CORNSY(NULONE(K),2) ) =
     *           XM( IEL , CORNSY(NULONE(K),2) ) + XN(K)
               ENDIF
             ENDDO
           ELSE
             DO K = 1 , NPTFR
               IEL = NELBOR(K)
               XM( IEL , CORNSY(NULONE(K),1) ) =
     *         XM( IEL , CORNSY(NULONE(K),1) ) + XN(K+NPTFRX)
               XM( IEL , CORNSY(NULONE(K),2) ) =
     *         XM( IEL , CORNSY(NULONE(K),2) ) + XN(K)
             ENDDO
           ENDIF
C
        ELSEIF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU M EST QUELCONQUE ET N SYMETRIQUE
C
           IF(NCSIZE.GT.1) THEN
             DO K = 1 , NPTFR
               IEL = NELBOR(K)
               IF(IEL.GT.0) THEN
                 XM( IEL , CORNSY(NULONE(K),1) ) =
     *           XM( IEL , CORNSY(NULONE(K),1) ) + XN(K)
                 XM( IEL , CORNSY(NULONE(K),2) ) =
     *           XM( IEL , CORNSY(NULONE(K),2) ) + XN(K)
               ENDIF
             ENDDO
           ELSE
             DO K = 1 , NPTFR
               IEL = NELBOR(K)
               XM( IEL , CORNSY(NULONE(K),1) ) =
     *         XM( IEL , CORNSY(NULONE(K),1) ) + XN(K)
               XM( IEL , CORNSY(NULONE(K),2) ) =
     *         XM( IEL , CORNSY(NULONE(K),2) ) + XN(K)
             ENDDO
           ENDIF
C
        ELSEIF(TYPEXM(1:1).EQ.'S'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU LES DEUX MATRICES SONT SYMETRIQUES
C
           IF(NCSIZE.GT.1) THEN
             DO K = 1 , NPTFR
               IEL = NELBOR(K)
               IF(IEL.GT.0) THEN
                 XM( IEL , CORSYM(NULONE(K)) ) =
     *           XM( IEL , CORSYM(NULONE(K)) ) + XN(K)
               ENDIF
             ENDDO
           ELSE
             DO K = 1 , NPTFR
               IEL = NELBOR(K)
               XM( IEL , CORSYM(NULONE(K)) ) =
     *         XM( IEL , CORSYM(NULONE(K)) ) + XN(K)
             ENDDO
           ENDIF
C
        ELSE
           IF (LNG.EQ.1) WRITE(LU,99) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,98) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           CALL PLANTE(1)
           STOP
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,140) OP
        IF (LNG.EQ.2) WRITE(LU,141) OP
140     FORMAT(1X,'OM1201 (BIEF) : OPERATION INCONNUE : ',A8)
141     FORMAT(1X,'OM1201 (BIEF) : UNKNOWN OPERATION : ',A8)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
