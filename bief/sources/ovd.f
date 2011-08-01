C                       **************
                        SUBROUTINE OVD
C                       **************
C
     * ( OP , X , Y , Z , C , NPOIN , IOPT , D , EPS )
C
C***********************************************************************
C BIEF VERSION 5.2           26/11/93    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C  FONCTION : OPERATIONS SUR LES VECTEURS AVEC DES DIVISIONS
C             LES DIVISIONS PAR ZERO PEUVENT ETRE TESTEES OU NON
C             EN CAS DE DIVISION PAR ZERO, ON PEUT S'ARRETER OU
C             METTRE UNE VALEUR D.
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES VECTEURS X,Y ET Z ET LA CONSTANTE C. LE RESULTAT
C   EST LE VECTEUR X.
C
C   OP = 'X=1/Y   '     :  INVERSE DE Y MIS DANS X
C   OP = 'X=Y/Z   '     :  DIVISION DE Y PAR Z MIS DANS X
C   OP = 'X=CY/Z  '     :  PRODUIT DE C ET Y DIVISE PAR Z ET MIS DANS X
C   OP = 'X=CXY/Z '     :  PRODUIT DE C, X ET Y DIVISE PAR Z ET MIS DS X
C   OP = 'X=X+CY/Z'     :  PRODUIT DE C ET Y DIVISE PAR Z AJOUTE A X
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |      OP        | -->| CHAINE DE CARACTERES INDIQUANT L'OPERATION
C |                |   >| A EFFECTUER.
C |      X         |<-- | VECTEUR RESULTAT
C |      Y         | -->| VECTEUR OPERANDE
C |      Z         | -->| VECTEUR OPERANDE
C |      C         | -->| CONSTANTE DONNEE
C |      NPOIN     | -->| DIMENSION DES VECTEURS
C |      IOPT      | -->| OPTION : 1 ON NE FAIT PAS DE TEST
C |                |    |          2 LES TERMES INFINIS SONT REMPLACES
C |                |    |            PAR LA CONSTANTE D.
C |                |    |          3 ARRET EN CAS DE DIVISION PAR ZERO
C |                |    |          4 LES TERMES INFINIS SONT TRONQUES
C |      EPS       | -->| CRITERE DE DIVISION PAR ZERO
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELANTS : BEAUCOUP
C PROGRAMMES APPELES   : NEANT
C
C PRECAUTIONS D'EMPLOI : LES OPERATIONS 1/Y ET Y/Z
C                        SUPPRIMENT LES DIVISIONS PAR 0.
C                        UN PASSAGE CORRECT N'EST DONC PAS
C                        UNE PREUVE QUE Y OU Z N'EST JAMAIS NUL.
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER,          INTENT(IN)    :: NPOIN,IOPT
      DOUBLE PRECISION, INTENT(INOUT) :: X(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: Y(NPOIN),Z(NPOIN),C,D,EPS
      CHARACTER(LEN=8), INTENT(IN)    :: OP
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C
C-----------------------------------------------------------------------
C
      IF(OP(1:8).EQ.'X=1/Y   ') THEN
C
        IF(IOPT.EQ.1) THEN
C
        DO 30 I=1,NPOIN
            X(I) = 1.D0/Y(I)
30      CONTINUE
C
        ELSEIF(IOPT.EQ.2) THEN
C
        DO 31 I=1,NPOIN
C
          IF (ABS(Y(I)).GT.EPS) THEN
            X(I) = 1.D0/Y(I)
          ELSE
            X(I) = D
          ENDIF
C
31      CONTINUE
C
        ELSEIF(IOPT.EQ.3) THEN
C
        DO 32 I=1,NPOIN
C
          IF (ABS(Y(I)).GT.EPS) THEN
            X(I) = 1.D0/Y(I)
          ELSE
            IF(LNG.EQ.1) WRITE(LU,1000) I,OP,EPS
            IF(LNG.EQ.2) WRITE(LU,2000) I,OP,EPS
            CALL PLANTE(1)
            STOP
          ENDIF
C
32      CONTINUE
C
        ELSEIF(IOPT.EQ.4) THEN
C
        DO 33 I=1,NPOIN
C
          IF (ABS(Y(I)).GT.EPS) THEN
            X(I) = 1.D0/Y(I)
          ELSEIF (Y(I).GE.0.D0) THEN
            X(I) =  1.D0/EPS
          ELSE
            X(I) = -1.D0/EPS
          ENDIF
C
33      CONTINUE
C
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=Y/Z   ') THEN
C
        IF(IOPT.EQ.1) THEN
C
        DO 40 I=1,NPOIN
            X(I) = Y(I) / Z(I)
40      CONTINUE
C
        ELSEIF(IOPT.EQ.2) THEN
C
        DO 41 I=1,NPOIN
C
          IF (ABS(Z(I)).GT.EPS) THEN
            X(I) = Y(I) / Z(I)
          ELSE
            X(I) = D
          ENDIF
C
41      CONTINUE
C
        ELSEIF(IOPT.EQ.3) THEN
C
        DO 42 I=1,NPOIN
C
          IF (ABS(Z(I)).GT.EPS) THEN
            X(I) = Y(I) / Z(I)
          ELSE
            IF(LNG.EQ.1) WRITE(LU,1000) I,OP,EPS
            IF(LNG.EQ.2) WRITE(LU,2000) I,OP,EPS
            CALL PLANTE(1)
            STOP
          ENDIF
C
42      CONTINUE
C>>>>
        ELSEIF(IOPT.EQ.4) THEN
C
        DO 43 I=1,NPOIN
C
          IF (ABS(Z(I)).GT.EPS) THEN
            X(I) = Y(I) / Z(I)
          ELSEIF (ABS(Y(I)).LT.EPS) THEN
            X(I) = D
          ELSEIF (Z(I).GE.0.D0) THEN
            X(I) =  Y(I) / EPS
          ELSE
            X(I) = -Y(I) / EPS
          ENDIF
C
43      CONTINUE
C<<<<
C
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=CY/Z  ') THEN
C
        IF(IOPT.EQ.1) THEN
C
        DO 50 I=1,NPOIN
            X(I) = C*Y(I) / Z(I)
50      CONTINUE
C
        ELSEIF(IOPT.EQ.2) THEN
C
        DO 51 I=1,NPOIN
C
          IF (ABS(Z(I)).GT.EPS) THEN
            X(I) = C*Y(I) / Z(I)
          ELSE
            X(I) = D
          ENDIF
C
51      CONTINUE
C
        ELSEIF(IOPT.EQ.3) THEN
C
        DO 52 I=1,NPOIN
C
          IF (ABS(Z(I)).GT.EPS) THEN
            X(I) = C*Y(I) / Z(I)
          ELSE
            IF(LNG.EQ.1) WRITE(LU,1000) I,OP,EPS
            IF(LNG.EQ.2) WRITE(LU,2000) I,OP,EPS
            CALL PLANTE(1)
            STOP
          ENDIF
C
52      CONTINUE
C>>>>
        ELSEIF(IOPT.EQ.4) THEN
C
        DO 53 I=1,NPOIN
C
          IF (ABS(Z(I)).GT.EPS) THEN
            X(I) = C*Y(I) / Z(I)
          ELSEIF (ABS(C*Y(I)).LT.EPS) THEN
            X(I) = D
          ELSEIF (Z(I).GE.0.D0) THEN
            X(I) =  C*Y(I) / EPS
          ELSE
            X(I) = -C*Y(I) / EPS
          ENDIF
C
53      CONTINUE
C<<<<
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=CXY/Z ') THEN
C
        IF(IOPT.EQ.1) THEN
C
        DO 60 I=1,NPOIN
            X(I) = C*X(I)*Y(I) / Z(I)
60      CONTINUE
C
        ELSEIF(IOPT.EQ.2) THEN
C
        DO 61 I=1,NPOIN
C
          IF (ABS(Z(I)).GT.EPS) THEN
            X(I) = C*X(I)*Y(I) / Z(I)
          ELSE
            X(I) = D
          ENDIF
C
61      CONTINUE
C
        ELSEIF(IOPT.EQ.3) THEN
C
        DO 62 I=1,NPOIN
C
          IF (ABS(Z(I)).GT.EPS) THEN
            X(I) = C*X(I)*Y(I) / Z(I)
          ELSE
            IF(LNG.EQ.1) WRITE(LU,1000) I,OP,EPS
            IF(LNG.EQ.2) WRITE(LU,2000) I,OP,EPS
            CALL PLANTE(1)
            STOP
          ENDIF
C
62      CONTINUE
C>>>>
        ELSEIF(IOPT.EQ.4) THEN
C
        DO 63 I=1,NPOIN
C
          IF (ABS(Z(I)).GT.EPS) THEN
            X(I) = C*X(I)*Y(I) / Z(I)
          ELSEIF (ABS(C*X(I)*Y(I)).LT.EPS) THEN
            X(I) = D
          ELSEIF (Z(I).GE.0.D0) THEN
            X(I) =  C*Y(I) / EPS
          ELSE
            X(I) = -C*Y(I) / EPS
          ENDIF
C
63      CONTINUE
C<<<<
C
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=X+CY/Z') THEN
C
        IF(IOPT.EQ.1) THEN
C
        DO 70 I=1,NPOIN
            X(I) = X(I) + C * Y(I) / Z(I)
70      CONTINUE
C
        ELSEIF(IOPT.EQ.2) THEN
C
        DO 71 I=1,NPOIN
C
          IF (ABS(Z(I)).GT.EPS) THEN
            X(I) = X(I) + C * Y(I) / Z(I)
          ELSE
            X(I) = D
          ENDIF
C
71      CONTINUE
C
        ELSEIF(IOPT.EQ.3) THEN
C
        DO 72 I=1,NPOIN
C
          IF (ABS(Z(I)).GT.EPS) THEN
            X(I) = X(I) + C * Y(I) / Z(I)
          ELSE
            IF(LNG.EQ.1) WRITE(LU,1000) I,OP,EPS
            IF(LNG.EQ.2) WRITE(LU,2000) I,OP,EPS
            CALL PLANTE(1)
            STOP
          ENDIF
C
72      CONTINUE
C>>>>
        ELSEIF(IOPT.EQ.4) THEN
C
        DO 73 I=1,NPOIN
C
          IF (ABS(Z(I)).GT.EPS) THEN
            X(I) = X(I) + C*Y(I) / Z(I)
          ELSEIF (ABS(C*Y(I)).LT.EPS) THEN
            X(I) = D
          ELSEIF (Z(I).GE.0.D0) THEN
            X(I) = X(I) + C*Y(I) / EPS
          ELSE
            X(I) = X(I) - C*Y(I) / EPS
          ENDIF
C
73      CONTINUE
C<<<<
C
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
C
         IF (LNG.EQ.1) WRITE(LU,3000) OP
         IF (LNG.EQ.2) WRITE(LU,4000) OP
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
1000     FORMAT(1X,'OVD (BIEF) : DIVISION PAR ZERO AU POINT ',1I6,
     *             ' LORS DE L''OPERATION ',A8,/,1X,
     *             'LE CRITERE EST ',G16.7)
2000     FORMAT(1X,'OVD (BIEF) : DIVIDE BY ZERO AT POINT ',1I6,
     *             ' FOR OPERATION ',A8,/,1X,
     *             'THE CRITERION IS ',G16.7)
3000     FORMAT(1X,'OVD (BIEF) : OPERATION INCONNUE : ',A8)
4000     FORMAT(1X,'OVD (BIEF) : UNKNOWN OPERATION: ',A8)
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
