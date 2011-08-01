C                       ****************
                        SUBROUTINE OMSEG
C                       ****************
C
     *(OP ,  DM,TYPDIM,XM,TYPEXM,   DN,TYPDIN,XN,TYPEXN,   D,C,
     * NDIAG,NSEG1,NSEG2,GLOSEG,SIZGLO)
C
C***********************************************************************
C BIEF VERSION 5.6         29/12/05    J-M HERVOUET (LNH) 01 30 87 80 18
C                                       
C***********************************************************************
C
C FONCTION : OPERATIONS ON MATRICES WITH AN EDGE-BASED STORAGE
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES MATRICES M ET N ,LA MATRICE DIAGONALE D ET LA
C   CONSTANTE C.
C
C   LE RESULTAT EST LA MATRICE M.
C
C      OP = 'M=N     '  : COPIE DE N DANS M
C      OP = 'M=CN    '  : PRODUIT DE N PAR LA CONSTANTE C
C      OP = 'M=M+CN  '  : ON AJOUTE CN A M
C      OP = 'M=TN    '  : TRANSPOSEE DE N DANS M
C      OP = 'M=M+TN  '  : TRANSPOSEE DE N AJOUTEE A N
C      OP = 'M=M+CTN '  : ON AJOUTE C TRANSPOSEE(N) A M
C      OP = 'M=M+N   '  : ON AJOUTE CN A M
C      OP = 'M=MD    '  : PRODUIT DE M PAR D A DROITE
C      OP = 'M=DM    '  : PRODUIT DE M PAR D A GAUCHE
C      OP = 'M=M-ND  '  : PRODUIT DE N PAR D A DROITE RETRANCHE A M
C      OP = 'M=M-DN  '  : PRODUIT DE N PAR D A GAUCHE RETRANCHE A M
C      OP = 'M=DMD   '  : PRODUIT DE M A DROITE ET A GAUCHE PAR D
C      OP = 'M=0     '  : ANNULATION DE M
C      OP = 'M=X(M)  '  : PASSAGE A UNE FORME NON SYMETRIQUE
C                         (ANCIEN MATSNS)
C      OP = 'M=MSK(M)'  : MASQUAGE DES TERMES EXTRADIAGONAUX
C                         (ANCIEN MASKEX)
C                         LE MASQUE EST PRIS DANS D
C      OP = 'M=M+D   '  : DIAGONALE D AJOUTEE A M
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
      USE BIEF, EX_OMSEG => OMSEG
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NDIAG,NSEG1,NSEG2,SIZGLO
      INTEGER, INTENT(IN) :: GLOSEG(SIZGLO,2)
      CHARACTER(LEN=8), INTENT(IN)    :: OP
      DOUBLE PRECISION, INTENT(IN)    :: DN(*),D(*)
C     XM AND XN MAY BE ONLY OF SIZE NSEG1 IF SYMMETRIC MATRIX 
C     SIZE GIVEN HERE ONLY FOR CHECKING BOUNDS
      DOUBLE PRECISION, INTENT(INOUT) :: XM(NSEG1+NSEG2)
      DOUBLE PRECISION, INTENT(IN)    :: XN(NSEG1+NSEG2)
      CHARACTER(LEN=1), INTENT(INOUT) :: TYPDIM,TYPEXM,TYPDIN,TYPEXN
      DOUBLE PRECISION, INTENT(INOUT) :: DM(*)
      DOUBLE PRECISION, INTENT(IN)    :: C
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTRINSIC MIN
C
      INTEGER ISEG,DIMX
C
      DOUBLE PRECISION Y(1),Z(1)
C
C-----------------------------------------------------------------------
C
C     ARRAYS XM AND XN ARE BASICALLY OF SIZE XM(DIMX,1 OR 2)
C     BUT IN CAS OF RECTANGULAR MATRICES OTHER DATA ARE STORED BEYOND
C     XM(2*DIMX)
C
      DIMX=MIN(NSEG1,NSEG2)
C
      IF(OP(3:8).EQ.'N     ') THEN
C
        IF(TYPDIN(1:1).EQ.'Q') THEN
          CALL OV( 'X=Y     ' , DM , DN , Z , C , NDIAG )
        ELSEIF(TYPDIN(1:1).EQ.'I'.OR.TYPDIN(1:1).EQ.'0') THEN
C         RIEN A FAIRE, JUSTE TYPDIN A RECOPIER.
        ELSE
           IF (LNG.EQ.1) WRITE(LU,5) TYPDIN(1:1)
           IF (LNG.EQ.2) WRITE(LU,6) TYPDIN(1:1)
5          FORMAT(1X,'OMSEG (BIEF) : TYPDIN INCONNU :',A1)
6          FORMAT(1X,'OMSEG (BIEF) : TYPDIN UNKNOWN :',A1)
           CALL PLANTE(1)
           STOP
        ENDIF
        TYPDIM(1:1)=TYPDIN(1:1)
C
        IF(TYPEXN(1:1).EQ.'S') THEN
           CALL OV( 'X=Y     ' , XM , XN , Z , C , NSEG1 )
        ELSEIF(TYPEXN(1:1).EQ.'Q') THEN
           CALL OV( 'X=Y     ' , XM , XN , Z , C , NSEG1+NSEG2 )
        ELSEIF(TYPEXN(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,10) TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,11) TYPEXN(1:1)
10         FORMAT(1X,'OMSEG (BIEF) : TYPEXN INCONNU :',A1)
11         FORMAT(1X,'OMSEG (BIEF) : TYPEXN UNKNOWN :',A1)
           CALL PLANTE(1)
           STOP
        ENDIF
        TYPEXM(1:1)=TYPEXN(1:1)
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(3:8).EQ.'CN    ') THEN
C
        CALL OV( 'X=CY    ' , DM , DN , Z , C , NDIAG )
C
        IF(TYPEXN(1:1).EQ.'S') THEN
           CALL OV( 'X=CY    ' , XM , XN , Z , C , NSEG1 )
        ELSEIF(TYPEXN(1:1).EQ.'Q') THEN
           CALL OV( 'X=CY    ' , XM , XN , Z , C , NSEG1+NSEG2 )
        ELSEIF(TYPEXN(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,10) TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,11) TYPEXN(1:1)
           CALL PLANTE(1)
           STOP
        ENDIF
C
        TYPDIM(1:1)=TYPDIN(1:1)
        TYPEXM(1:1)=TYPEXN(1:1)
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(3:8).EQ.'M+CN  ' .OR.
     *      (OP(3:8).EQ.'M+CTN ').AND.TYPEXN(1:1).NE.'Q') THEN
C
        IF(TYPDIN(1:1).EQ.'I') THEN
          CALL OV( 'X=X+C   ' , DM , DN , Z , C , NDIAG )
        ELSEIF(TYPDIN(1:1).NE.'0') THEN
          CALL OV( 'X=X+CY  ' , DM , DN , Z , C , NDIAG )
        ENDIF
C
        IF(TYPEXN(1:1).EQ.'S') THEN
           CALL OV( 'X=X+CY  ' , XM , XN , Z , C , NSEG1 )
           IF(TYPEXM(1:1).EQ.'Q') THEN
             CALL OV( 'X=X+CY  ' , XM(DIMX+1:DIMX+NSEG1) ,XN,Z,C,NSEG1)
           ENDIF
        ELSEIF(TYPEXN(1:1).EQ.'Q') THEN
           IF(TYPEXM(1:1).NE.'Q') THEN
             IF (LNG.EQ.1) WRITE(LU,99) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
             IF (LNG.EQ.2) WRITE(LU,98) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
99          FORMAT(1X,'OMSEG (BIEF) : TYPEXM = ',A1,' NE CONVIENT PAS',
     *       /,1X,'POUR L''OPERATION : ',A8,' AVEC TYPEXN = ',A1)
98          FORMAT(1X,'OMSEG (BIEF) : TYPEXM = ',A1,' DOES NOT GO',
     *       /,1X,'FOR THE OPERATION : ',A8,' WITH TYPEXN = ',A1)
             CALL PLANTE(1)
             STOP
           ENDIF
           CALL OV( 'X=X+CY  ' , XM , XN , Z , C , NSEG1+NSEG2 )
        ELSEIF(TYPEXN(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,10) TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,11) TYPEXN(1:1)
           CALL PLANTE(1)
           STOP
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(3:8).EQ.'M+CTN ') THEN
C
C  LES CAS OU N EST SYMETRIQUE SONT TRAITES AVEC M=M+CN
C
        CALL OV( 'X=X+CY  ' , DM , DN , Z , C , NDIAG )
C
        IF(NSEG1.NE.NSEG2) THEN
          WRITE(LU,*) 'M+CTN : RECTANGULAR MATRIX NOT IMPLEMENTED'
          STOP
        ENDIF
        IF(TYPEXN(1:1).EQ.'Q') THEN
           IF(TYPEXM(1:1).NE.'Q') THEN
             IF (LNG.EQ.1) WRITE(LU,99) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
             IF (LNG.EQ.2) WRITE(LU,98) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
             CALL PLANTE(1)
             STOP
           ENDIF
           CALL OV( 'X=X+CY  ' , XM , XN(DIMX+1:DIMX+NSEG1) ,Z,C,NSEG1)
           CALL OV( 'X=X+CY  ' , XM(DIMX+1:DIMX+NSEG1) , XN ,Z,C,NSEG1)
        ELSE
           IF (LNG.EQ.1) WRITE(LU,10) TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,11) TYPEXN(1:1)
           CALL PLANTE(1)
           STOP
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(3:8).EQ.'TN    ') THEN
C
        CALL OV( 'X=Y     ' , DM , DN , Z , C , NDIAG )
C
        IF(TYPEXN(1:1).EQ.'S') THEN
           CALL OV( 'X=Y     ' , XM , XN , Z , C , NSEG1 )
        ELSEIF(TYPEXN(1:1).EQ.'Q') THEN
           IF(TYPEXM(1:1).NE.'S'.AND.NSEG1.NE.NSEG2) THEN
             WRITE(LU,*) 'TN : RECTANGULAR MATRIX NOT IMPLEMENTED'
             STOP
           ENDIF
           CALL OV( 'X=Y     ' , XM , XN(DIMX+1:DIMX+NSEG1) ,Z,C,NSEG1)
           CALL OV( 'X=Y     ' , XM(DIMX+1:DIMX+NSEG1) , XN ,Z,C,NSEG1)
        ELSEIF(TYPEXN(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,10) TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,11) TYPEXN(1:1)
           CALL PLANTE(1)
           STOP
        ENDIF
        TYPDIM(1:1)=TYPDIN(1:1)
        TYPEXM(1:1)=TYPEXN(1:1)
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(3:8).EQ.'M+N   '.OR.
     *      (OP(3:8).EQ.'M+TN  ').AND.TYPEXN(1:1).NE.'Q') THEN
C
        CALL OV( 'X=X+Y   ' , DM , DN , Z , C , NDIAG )
C
        IF(TYPEXN(1:1).EQ.'S') THEN
           CALL OV( 'X=X+Y   ' , XM , XN , Z , C , NSEG1 )
           IF(TYPEXM(1:1).EQ.'Q') THEN
           CALL OV( 'X=X+Y   ' , XM(DIMX+1:DIMX+NSEG1) , XN ,Z,C,NSEG1)
           ENDIF
        ELSEIF(TYPEXN(1:1).EQ.'Q') THEN
           IF(TYPEXM(1:1).NE.'Q') THEN
             IF (LNG.EQ.1) WRITE(LU,99) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
             IF (LNG.EQ.2) WRITE(LU,98) TYPEXM(1:1),OP(1:8),TYPEXN(1:1)
             CALL PLANTE(1)
             STOP
           ENDIF
           CALL OV( 'X=X+Y   ' , XM , XN , Z , C , NSEG1+NSEG2 )
        ELSEIF(TYPEXN(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,10) TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,11) TYPEXN(1:1)
           CALL PLANTE(1)
           STOP
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(3:8).EQ.'M+TN  ') THEN
C
C     LE CAS N SYMETRIQUE A ETE TRAITE PLUS HAUT
C
        CALL OV( 'X=X+Y   ' , DM , DN , Z , C , NDIAG )
C
        IF(NSEG1.NE.NSEG2) THEN
          WRITE(LU,*) 'M+TN : RECTANGULAR MATRIX NOT IMPLEMENTED'
          STOP
        ENDIF
        IF(TYPEXM(1:1).EQ.'Q') THEN
           CALL OV( 'X=X+Y   ' , XM , XN(DIMX+1:DIMX+NSEG1) ,Z,C,NSEG1)
           CALL OV( 'X=X+Y   ' , XM(DIMX+1:DIMX+NSEG1) , XN ,Z,C,NSEG1)
        ELSEIF(TYPEXN(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,10) TYPEXN(1:1)
           IF (LNG.EQ.2) WRITE(LU,11) TYPEXN(1:1)
           CALL PLANTE(1)
           STOP
        ENDIF
        TYPDIM(1:1)=TYPDIN(1:1)
        TYPEXM(1:1)=TYPEXN(1:1)
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(3:8).EQ.'MD    ') THEN
C
C   TERMES DIAGONAUX
C
         IF(TYPDIM(1:1).EQ.'Q') THEN
           CALL OV( 'X=XY    ' , DM , D , Z , C , NDIAG )
         ELSEIF(TYPDIM(1:1).EQ.'I') THEN
           CALL OV( 'X=Y     ' , DM , D , Z , C , NDIAG )
           TYPDIM(1:1)='Q'
         ELSEIF(TYPDIM(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,12) TYPDIM(1:1)
           IF (LNG.EQ.2) WRITE(LU,13) TYPDIM(1:1)
           CALL PLANTE(1)
           STOP
         ENDIF
C
C   TERMES EXTRADIAGONAUX
C
         IF(TYPEXM(1:1).EQ.'Q') THEN
C
         DO ISEG = 1 , MIN(NSEG1,NSEG2)
           XM(ISEG) = XM(ISEG) * D(GLOSEG(ISEG,2))
           XM(ISEG+DIMX) = XM(ISEG+DIMX) * D(GLOSEG(ISEG,1))
         ENDDO
         IF(NSEG1.GT.NSEG2) THEN
           DO ISEG = MIN(NSEG1,NSEG2)+1,MAX(NSEG1,NSEG2)
             XM(ISEG+DIMX)=XM(ISEG+DIMX)*D(GLOSEG(ISEG,1))
           ENDDO
         ELSEIF(NSEG2.GT.NSEG1) THEN
           DO ISEG = MIN(NSEG1,NSEG2)+1,MAX(NSEG1,NSEG2)
             XM(ISEG+DIMX)=XM(ISEG+DIMX)*D(GLOSEG(ISEG,2))
           ENDDO
         ENDIF
C
         ELSEIF(TYPEXM(1:1).EQ.'S') THEN
          IF (LNG.EQ.1) WRITE(LU,170)
          IF (LNG.EQ.2) WRITE(LU,171)
170       FORMAT(1X,'OMSEG (BIEF) : M=MD , M DOIT ETRE NON SYMETRIQUE')
171       FORMAT(1X,'OMSEG (BIEF) : M=MD , M MUST BE NON-SYMMETRIC')
          CALL PLANTE(1)
          STOP
         ELSEIF(TYPEXM(1:1).NE.'0') THEN
          IF (LNG.EQ.1) WRITE(LU,172) TYPEXM(1:1)
          IF (LNG.EQ.2) WRITE(LU,173) TYPEXM(1:1)
172       FORMAT(1X,'OMSEG (BIEF) : TYPEXM NON PREVU : ',A1)
173       FORMAT(1X,'OMSEG (BIEF) : TYPEXM NOT AVAILABLE : ',A1)
          CALL PLANTE(1)
          STOP
         ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(3:8).EQ.'DM    ') THEN
C
C   TERMES DIAGONAUX
C
         IF(TYPDIM(1:1).EQ.'Q') THEN
           CALL OV( 'X=XY    ' , DM , D , Z , C , NDIAG )
         ELSEIF(TYPDIM(1:1).EQ.'I') THEN
           CALL OV( 'X=Y     ' , DM , D , Z , C , NDIAG )
           TYPDIM(1:1)='Q'
         ELSEIF(TYPDIM(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,12) TYPDIM(1:1)
           IF (LNG.EQ.2) WRITE(LU,13) TYPDIM(1:1)
           CALL PLANTE(1)
           STOP
         ENDIF
C
C   TERMES EXTRADIAGONAUX
C
         IF(TYPEXM(1:1).EQ.'Q') THEN
C
         DO ISEG = 1 , MIN(NSEG1,NSEG2)
           XM(ISEG) = XM(ISEG) * D(GLOSEG(ISEG,1))
           XM(ISEG+DIMX) = XM(ISEG+DIMX) * D(GLOSEG(ISEG,2))
         ENDDO
         IF(NSEG1.GT.NSEG2) THEN
           DO ISEG = MIN(NSEG1,NSEG2)+1,MAX(NSEG1,NSEG2)
             XM(ISEG+DIMX)=XM(ISEG+DIMX)*D(GLOSEG(ISEG,2))
           ENDDO
         ELSEIF(NSEG2.GT.NSEG1) THEN
           DO ISEG = MIN(NSEG1,NSEG2)+1,MAX(NSEG1,NSEG2)
             XM(ISEG+DIMX)=XM(ISEG+DIMX)*D(GLOSEG(ISEG,1))
           ENDDO
         ENDIF
C
         ELSEIF(TYPEXM(1:1).EQ.'S') THEN
          IF (LNG.EQ.1) WRITE(LU,180)
          IF (LNG.EQ.2) WRITE(LU,181)
180       FORMAT(1X,'OMSEG (BIEF) : M=DM N''EST PAS SYMETRIQUE')
181       FORMAT(1X,'OMSEG (BIEF) : M=MD IS NOT SYMMETRIC')
          CALL PLANTE(1)
          STOP
         ELSEIF(TYPEXM(1:1).NE.'0') THEN
          IF (LNG.EQ.1) WRITE(LU,172) TYPEXM(1:1)
          IF (LNG.EQ.2) WRITE(LU,173) TYPEXM(1:1)
          CALL PLANTE(1)
          STOP
         ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(3:8).EQ.'M-DN  ') THEN
C
C   TERMES DIAGONAUX
C
         IF(TYPDIM(1:1).EQ.'Q') THEN
           CALL OV( 'X=X-YZ  ' , DM , DN , D , C , NDIAG )
         ELSEIF(TYPDIM(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,12) TYPDIM(1:1)
           IF (LNG.EQ.2) WRITE(LU,13) TYPDIM(1:1)
           CALL PLANTE(1)
           STOP
         ENDIF
C
C   TERMES EXTRADIAGONAUX
C
         IF(TYPEXM(1:1).EQ.'Q') THEN
           IF(TYPEXN(1:1).EQ.'Q') THEN
           DO ISEG = 1 , NSEG1
            XM(ISEG     )=XM(ISEG     )-XN(ISEG     )*D(GLOSEG(ISEG,1))
            XM(ISEG+DIMX)=XM(ISEG+DIMX)-XN(ISEG+DIMX)*D(GLOSEG(ISEG,2))
           ENDDO
           ELSEIF(TYPEXN(1:1).EQ.'S') THEN
           DO ISEG = 1 , NSEG1
            XM(ISEG     ) = XM(ISEG     ) - XN(ISEG) * D(GLOSEG(ISEG,1))
            XM(ISEG+DIMX) = XM(ISEG+DIMX) - XN(ISEG) * D(GLOSEG(ISEG,2))
           ENDDO
           ELSEIF(TYPEXN(1:1).NE.'0') THEN
            IF (LNG.EQ.1) WRITE(LU,10) TYPEXN(1:1)
            IF (LNG.EQ.2) WRITE(LU,11) TYPEXN(1:1)
            CALL PLANTE(1)
            STOP
           ENDIF
         ELSE
           IF (LNG.EQ.1) WRITE(LU,172) TYPEXM(1:1)
           IF (LNG.EQ.2) WRITE(LU,173) TYPEXM(1:1)
           CALL PLANTE(1)
           STOP
         ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(3:8).EQ.'M-ND  ') THEN
C
C   TERMES DIAGONAUX
C
         IF(TYPDIM(1:1).EQ.'Q') THEN
           CALL OV( 'X=X-YZ  ' , DM , DN , D , C , NDIAG )
         ELSEIF(TYPDIM(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,12) TYPDIM(1:1)
           IF (LNG.EQ.2) WRITE(LU,13) TYPDIM(1:1)
           CALL PLANTE(1)
           STOP
         ENDIF
C
C   TERMES EXTRADIAGONAUX
C
         IF(TYPEXM(1:1).EQ.'Q') THEN
           IF(TYPEXN(1:1).EQ.'Q') THEN
           DO ISEG = 1 , NSEG1
            XM(ISEG     )=XM(ISEG     )-XN(ISEG     )*D(GLOSEG(ISEG,2))
            XM(ISEG+DIMX)=XM(ISEG+DIMX)-XN(ISEG+DIMX)*D(GLOSEG(ISEG,1))
           ENDDO
           ELSEIF(TYPEXN(1:1).EQ.'S') THEN
           DO ISEG = 1 , NSEG1
            XM(ISEG     ) = XM(ISEG     ) - XN(ISEG) * D(GLOSEG(ISEG,2))
            XM(ISEG+DIMX) = XM(ISEG+DIMX) - XN(ISEG) * D(GLOSEG(ISEG,1))
           ENDDO
           ELSEIF(TYPEXN(1:1).NE.'0') THEN
            IF (LNG.EQ.1) WRITE(LU,10) TYPEXN(1:1)
            IF (LNG.EQ.2) WRITE(LU,11) TYPEXN(1:1)
            CALL PLANTE(1)
            STOP
           ENDIF
         ELSE
           IF (LNG.EQ.1) WRITE(LU,172) TYPEXM(1:1)
           IF (LNG.EQ.2) WRITE(LU,173) TYPEXM(1:1)
           CALL PLANTE(1)
           STOP
         ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(3:8).EQ.'DMD   ') THEN
C
C   DIAGONALE
C
         IF(TYPDIM(1:1).EQ.'Q') THEN
           CALL OV( 'X=XY    ' , DM , D , Z , C , NDIAG )
           CALL OV( 'X=XY    ' , DM , D , Z , C , NDIAG )
         ELSEIF(TYPDIM(1:1).EQ.'I') THEN
           CALL OV( 'X=YZ    ' , DM , D , D , C , NDIAG )
           TYPDIM(1:1)='Q'
         ELSEIF(TYPDIM(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,12) TYPDIM(1:1)
           IF (LNG.EQ.2) WRITE(LU,13) TYPDIM(1:1)
12         FORMAT(1X,'OMSEG (BIEF) : TYPDIM INCONNU :',A1)
13         FORMAT(1X,'OMSEG (BIEF) : TYPDIM UNKNOWN :',A1)
           CALL PLANTE(1)
           STOP
         ENDIF
C
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX
C
         IF(TYPEXM(1:1).EQ.'S') THEN
C
         DO ISEG = 1 , NSEG1
           XM(ISEG)=XM(ISEG)*D(GLOSEG(ISEG,1))*D(GLOSEG(ISEG,2))
         ENDDO
C
         ELSEIF(TYPEXM(1:1).EQ.'Q') THEN
C
         DO ISEG = 1 , NSEG1
           XM(ISEG     )=XM(ISEG     )
     *                      *D(GLOSEG(ISEG,1))*D(GLOSEG(ISEG,2))
           XM(ISEG+DIMX)=XM(ISEG+DIMX)
     *                      *D(GLOSEG(ISEG,1))*D(GLOSEG(ISEG,2))
         ENDDO
C
        ELSEIF(TYPEXM(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,20) TYPEXM(1:1)
           IF (LNG.EQ.2) WRITE(LU,21) TYPEXM(1:1)
20         FORMAT(1X,'OMSEG (BIEF) : TYPEXM INCONNU :',A1)
21         FORMAT(1X,'OMSEG (BIEF) : TYPEXM UNKNOWN :',A1)
           CALL PLANTE(1)
           STOP
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(3:8).EQ.'M+D   ') THEN
C
        CALL OV( 'X=X+Y   ' , DM , D , Z , 0.D0 , NDIAG )
C       ICI IL Y A UN DOUTE SUR TYPDIM
        TYPDIM(1:1)='Q'
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(3:8).EQ.'0     ') THEN
C
        CALL OV( 'X=C     ' , DM , Y , Z , 0.D0 , NDIAG )
C
        IF(TYPEXM(1:1).EQ.'S') THEN
           CALL OV( 'X=C     ' , XM , Y , Z , 0.D0 , NSEG1 )
        ELSEIF(TYPEXM(1:1).EQ.'Q') THEN
           CALL OV( 'X=C     ' , XM , Y , Z , 0.D0 , NSEG1+NSEG2 )
        ELSEIF(TYPEXM(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,710) TYPEXM(1:1)
           IF (LNG.EQ.2) WRITE(LU,711) TYPEXM(1:1)
710        FORMAT(1X,'OMSEG (BIEF) : TYPEXM INCONNU :',A1)
711        FORMAT(1X,'OMSEG (BIEF) : TYPEXM UNKNOWN :',A1)
           CALL PLANTE(1)
           STOP
        ENDIF
C       ON NE CHANGE PAS TYPDIM
C       TYPDIM(1:1)='0'
C       ON NE CHANGE PAS TYPEXM
C       TYPEXM(1:1)='0'
C-----------------------------------------------------------------------
C
      ELSEIF(OP(3:8).EQ.'X(M)  ') THEN
C
        IF(TYPEXM(1:1).EQ.'S') THEN
          CALL OV('X=Y     ',XM(DIMX+1:DIMX+NSEG1),
     *                       XM(     1:     NSEG1),Z,C,NSEG1)
        ELSEIF(TYPEXM(1:1).NE.'0') THEN
          IF (LNG.EQ.1) WRITE(LU,810) TYPEXM(1:1)
          IF (LNG.EQ.2) WRITE(LU,811) TYPEXM(1:1)
810       FORMAT(1X,'OMSEG (BIEF) : MATRICE DEJA NON SYMETRIQUE : ',A1)
811       FORMAT(1X,'OMSEG (BIEF): MATRIX ALREADY NON SYMMETRICAL: ',A1)
          CALL PLANTE(1)
          STOP
        ENDIF
        TYPEXM(1:1)='Q'
C
C-----------------------------------------------------------------------
C
C     ELSEIF(OP(3:8).EQ.'MSK(M)') THEN
C
C     IF(TYPEXM(1:1).EQ.'S') THEN
C       J = 3
C     ELSEIF(TYPEXM(1:1).EQ.'Q') THEN
C       J = 6
C     ELSEIF(TYPEXM(1:1).EQ.'0') THEN
C       J = 0
C     ELSE
C       IF(LNG.EQ.1) WRITE(LU,710) TYPEXM
C       IF(LNG.EQ.2) WRITE(LU,711) TYPEXM
C       J=0
C       CALL PLANTE(1)
C       STOP
C     ENDIF
C
C     IF(J.GT.0) THEN
C        DO I = 1,J
C           CALL OV ( 'X=XY    ' , XM(1,I) , D , Z , C , NELEM )
C        ENDDO
C     ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,40) OP
        IF (LNG.EQ.2) WRITE(LU,41) OP
40      FORMAT(1X,'OMSEG (BIEF) : OPERATION INCONNUE : ',A8)
41      FORMAT(1X,'OMSEG (BIEF) : UNKNOWN OPERATION : ',A8)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
