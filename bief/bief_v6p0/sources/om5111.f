C                       *****************
                        SUBROUTINE OM5111
C                       *****************
C
     *(OP ,  DM,TYPDIM,XM,TYPEXM,   DN,TYPDIN,XN,TYPEXN,   C,
     * SIZDN,SZMDN,SIZXN,SZMXN,NETAGE, NELMAX3D)
C
C***********************************************************************
C BIEF VERSION 5.3           28/08/02    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FONCTION : OPERATIONS SUR DES MATRICES   M: TETRAEDRE
C                                          N: MATRICE DE TRIANGLES
C                                             AU FOND OU EN SURFACE
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
C  CONVENTION POUR LE STOCKAGE DES TERMES EXTRA-DIAGONAUX EN TRIANGLES
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
C     USE BIEF, EX_OM5111 => OM5111
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NETAGE,SIZDN,SZMXN,SZMDN,SIZXN
      INTEGER, INTENT(IN)             :: NELMAX3D
      CHARACTER(LEN=8), INTENT(IN)    :: OP
ccc   DOUBLE PRECISION, INTENT(IN)    :: DN(*),XN(SZMXN,*)
      DOUBLE PRECISION, INTENT(IN)    :: DN(*),XN(NELMAX3D/(3*NETAGE),*)
      DOUBLE PRECISION, INTENT(INOUT) :: DM(SZMDN,*)
ccc   DOUBLE PRECISION, INTENT(INOUT) :: XM(SZMXN,NETAGE,*)
      DOUBLE PRECISION,INTENT(INOUT)::XM(NELMAX3D/(3*NETAGE),3,NETAGE,*)
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
198       FORMAT(1X,'OM5111 (BIEF) : TYPDIM = ',A1,' NON PROGRAMME',
     *      /,1X,'POUR L''OPERATION : ',A8,' AVEC TYPDIN = ',A1)
199       FORMAT(1X,'OM5111 (BIEF) : TYPDIM = ',A1,' NOT IMPLEMENTED',
     *      /,1X,'FOR THE OPERATION : ',A8,' WITH TYPDIN = ',A1)
          CALL PLANTE(1)
          STOP
        ENDIF
C
        IF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'Q') THEN
C
C          CAS OU LES DEUX MATRICES SONT NON SYMETRIQUES
C
C          XM(K,1,1,  )  K : NUMERO DE TRIANGLE
C                        1 : TETRAEDRE T1, CELUI DONT UNE FACE EST AU FOND
C                        1 : PREMIER ETAGE, CELUI DU FOND
C
           DO K = 1 , SIZXN
             XM(K,1,1,01) = XM(K,1,1,01) + XN(K,1)
             XM(K,1,1,02) = XM(K,1,1,02) + XN(K,2)
             XM(K,1,1,04) = XM(K,1,1,04) + XN(K,3)
             XM(K,1,1,07) = XM(K,1,1,07) + XN(K,4)
             XM(K,1,1,08) = XM(K,1,1,08) + XN(K,5)
             XM(K,1,1,10) = XM(K,1,1,10) + XN(K,6)
           ENDDO
C
        ELSEIF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU M EST QUELCONQUE ET N SYMETRIQUE
C
           DO K = 1 , SIZXN
             XM(K,1,1,01) = XM(K,1,1,01) + XN(K,1)
             XM(K,1,1,02) = XM(K,1,1,02) + XN(K,2)
             XM(K,1,1,04) = XM(K,1,1,04) + XN(K,3)
             XM(K,1,1,07) = XM(K,1,1,07) + XN(K,1)
             XM(K,1,1,08) = XM(K,1,1,08) + XN(K,2)
             XM(K,1,1,10) = XM(K,1,1,10) + XN(K,3)
           ENDDO
C
        ELSEIF(TYPEXM(1:1).EQ.'S'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU LES DEUX MATRICES SONT SYMETRIQUES
C
           DO K = 1 , SIZXN
             XM(K,1,1,01) = XM(K,1,1,01) + XN(K,1)
             XM(K,1,1,02) = XM(K,1,1,02) + XN(K,2)
             XM(K,1,1,04) = XM(K,1,1,04) + XN(K,3)
           ENDDO
C
        ELSE
           IF (LNG.EQ.1) WRITE(LU,98) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,99) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
98         FORMAT(1X,'OM5111 (BIEF) : TYPEXM = ',A1,' NE CONVIENT PAS',
     *       /,1X,'POUR L''OPERATION : ',A8,' AVEC TYPEXN = ',A1)
99         FORMAT(1X,'OM5111 (BIEF) : TYPEXM = ',A1,' DOES NOT GO',
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
           DO K = 1 , SIZXN
             XM(K,1,1,01) = XM(K,1,1,01) + XN(K,4)
             XM(K,1,1,02) = XM(K,1,1,02) + XN(K,5)
             XM(K,1,1,04) = XM(K,1,1,04) + XN(K,6)
             XM(K,1,1,07) = XM(K,1,1,07) + XN(K,1)
             XM(K,1,1,08) = XM(K,1,1,08) + XN(K,2)
             XM(K,1,1,10) = XM(K,1,1,10) + XN(K,3)
           ENDDO
C
        ELSEIF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU M EST QUELCONQUE ET N SYMETRIQUE
C
           DO K = 1 , SIZXN
             XM(K,1,1,01) = XM(K,1,1,01) + XN(K,1)
             XM(K,1,1,02) = XM(K,1,1,02) + XN(K,2)
             XM(K,1,1,04) = XM(K,1,1,04) + XN(K,3)
             XM(K,1,1,07) = XM(K,1,1,07) + XN(K,1)
             XM(K,1,1,08) = XM(K,1,1,08) + XN(K,2)
             XM(K,1,1,10) = XM(K,1,1,10) + XN(K,3)
           ENDDO
C
        ELSEIF(TYPEXM(1:1).EQ.'S'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU LES DEUX MATRICES SONT SYMETRIQUES
C
           DO K = 1 , SIZXN
             XM(K,1,1,01) = XM(K,1,1,01) + XN(K,1)
             XM(K,1,1,02) = XM(K,1,1,02) + XN(K,2)
             XM(K,1,1,04) = XM(K,1,1,04) + XN(K,3)
           ENDDO
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
C          XM(K,1,1,  )  K      : NUMERO DE TRIANGLE
C                        2      : TETRAEDRE T2, CELUI DONT UNE FACE EST EN SURFACE
C                        NETAGE : DERNIER ETAGE, CELUI DE LA SURFACE
C
           DO K = 1 , SIZXN
             XM(K,2,NETAGE,02) = XM(K,2,NETAGE,02) + XN(K,1)
             XM(K,2,NETAGE,01) = XM(K,2,NETAGE,01) + XN(K,2)
             XM(K,2,NETAGE,10) = XM(K,2,NETAGE,10) + XN(K,3)
             XM(K,2,NETAGE,08) = XM(K,2,NETAGE,08) + XN(K,4)
             XM(K,2,NETAGE,07) = XM(K,2,NETAGE,07) + XN(K,5)
             XM(K,2,NETAGE,04) = XM(K,2,NETAGE,04) + XN(K,6)
           ENDDO
C
        ELSEIF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU M EST QUELCONQUE ET N SYMETRIQUE
C
           DO K = 1 , SIZXN
             XM(K,2,NETAGE,02) = XM(K,2,NETAGE,02) + XN(K,1)
             XM(K,2,NETAGE,01) = XM(K,2,NETAGE,01) + XN(K,2)
             XM(K,2,NETAGE,10) = XM(K,2,NETAGE,10) + XN(K,3)
             XM(K,2,NETAGE,08) = XM(K,2,NETAGE,08) + XN(K,1)
             XM(K,2,NETAGE,07) = XM(K,2,NETAGE,07) + XN(K,2)
             XM(K,2,NETAGE,04) = XM(K,2,NETAGE,04) + XN(K,3)
           ENDDO
C
        ELSEIF(TYPEXM(1:1).EQ.'S'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU LES DEUX MATRICES SONT SYMETRIQUES
C
           DO K = 1 , SIZXN
             XM(K,2,NETAGE,02) = XM(K,2,NETAGE,02) + XN(K,1)
             XM(K,2,NETAGE,01) = XM(K,2,NETAGE,01) + XN(K,2)
             XM(K,2,NETAGE,04) = XM(K,2,NETAGE,04) + XN(K,3)
           ENDDO
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
           DO K = 1 , SIZXN
             XM(K,2,NETAGE,02) = XM(K,2,NETAGE,02) + XN(K,4)
             XM(K,2,NETAGE,01) = XM(K,2,NETAGE,01) + XN(K,5)
             XM(K,2,NETAGE,10) = XM(K,2,NETAGE,10) + XN(K,6)
             XM(K,2,NETAGE,08) = XM(K,2,NETAGE,08) + XN(K,1)
             XM(K,2,NETAGE,07) = XM(K,2,NETAGE,07) + XN(K,2)
             XM(K,2,NETAGE,04) = XM(K,2,NETAGE,04) + XN(K,3)
           ENDDO
C
        ELSEIF(TYPEXM(1:1).EQ.'Q'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU M EST QUELCONQUE ET N SYMETRIQUE
C
           DO K = 1 , SIZXN
             XM(K,2,NETAGE,02) = XM(K,2,NETAGE,02) + XN(K,1)
             XM(K,2,NETAGE,01) = XM(K,2,NETAGE,01) + XN(K,2)
             XM(K,2,NETAGE,10) = XM(K,2,NETAGE,10) + XN(K,3)
             XM(K,2,NETAGE,08) = XM(K,2,NETAGE,08) + XN(K,1)
             XM(K,2,NETAGE,07) = XM(K,2,NETAGE,07) + XN(K,2)
             XM(K,2,NETAGE,04) = XM(K,2,NETAGE,04) + XN(K,3)
           ENDDO
C
        ELSEIF(TYPEXM(1:1).EQ.'S'.AND.TYPEXN(1:1).EQ.'S') THEN
C
C          CAS OU LES DEUX MATRICES SONT SYMETRIQUES
C
           DO K = 1 , SIZXN
             XM(K,2,NETAGE,02) = XM(K,2,NETAGE,02) + XN(K,1)
             XM(K,2,NETAGE,01) = XM(K,2,NETAGE,01) + XN(K,2)
             XM(K,2,NETAGE,04) = XM(K,2,NETAGE,04) + XN(K,3)
           ENDDO
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
70      FORMAT(1X,'OM5111 (BIEF) : OPERATION INCONNUE : ',A8)
71      FORMAT(1X,'OM5111 (BIEF) : UNKNOWN OPERATION : ',A8)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
