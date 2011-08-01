C                       *************
                        SUBROUTINE OV
C                       *************
C
     * ( OP , X , Y , Z , C , NPOIN )
C
C***********************************************************************
C BIEF VERSION 6.0     17/03/2010     J-M HERVOUET (LNHE) 01 30 87 80 18
C                                         F  LEPEINTRE (LNH) 30 87 78 54
C                            IDEE EMPRUNTEE A ULYSSE. MERCI D. LAURENCE
C***********************************************************************
C
C  FONCTION : OPERATIONS SUR LES VECTEURS
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES VECTEURS X,Y ET Z ET LA CONSTANTE C. LE RESULTAT
C   EST LE VECTEUR X.
C
C   OP = 'X=C     '     :  X MIS A LA VALEUR C
C   OP = 'X=Y     '     :  Y COPIE DANS X
C   OP = 'X=+Y    '     :  IDEM
C   OP = 'X=-Y    '     : -Y COPIE DANS X
C   OP = 'X=1/Y   '     :  INVERSE DE Y MIS DANS X
C   OP = 'X=Y+Z   '     :  SOMME DE Y ET Z MISE DANS X
C   OP = 'X=Y-Z   '     :  DIFFERENCE DE Y ET Z MISE DANS X
C   OP = 'X=YZ    '     :  PRODUIT Y PAR  Z MIS DANS X
C   OP = 'X=-YZ   '     :  PRODUIT Y PAR  Z MIS DANS X
C   OP = 'X=XY    '     :  PRODUIT Y PAR  X MIS DANS X
C   OP = 'X=X+YZ  '     :  PRODUIT Y PAR  Z AJOUTE A X
C   OP = 'X=X-YZ  '     :  PRODUIT Y PAR  Z RETRANCHE A X
C   OP = 'X=CXY   '     :  PRODUIT DE C X ET Y MIS DANS X
C   OP = 'X=CYZ   '     :  PRODUIT DE C, Y ET Z MIS DANS X
C   OP = 'X=CXYZ  '     :  PRODUIT DE C, X, Y ET Z MIS DANS X
C   OP = 'X=X+CYZ '     :  PRODUIT DE C, Y ET Z AJOUTE A X
C   OP = 'X=Y/Z   '     :  DIVISION DE Y PAR Z MIS DANS X
C   OP = 'X=CY/Z  '     :  PRODUIT DE C ET Y DIVISE PAR Z ET MIS DANS X
C   OP = 'X=CXY/Z '     :  PRODUIT DE C, X ET Y DIVISE PAR Z ET MIS DS X
C   OP = 'X=X+CY/Z'     :  PRODUIT DE C ET Y DIVISE PAR Z AJOUTE A X
C   OP = 'X=X+Y   '     :  Y AJOUTE A X
C   OP = 'X=X-Y   '     :  Y RETRANCHE A X
C   OP = 'X=CX    '     :  X MULTIPLIE PAR C
C   OP = 'X=CY    '     :  CY MIS DANS X
C   OP = 'X=Y+CZ  '     :  CZ AJOUTE A Y ET MIS DANS X
C   OP = 'X=X+CY  '     :  CY AJOUTE A X
C   OP = 'X=SQR(Y)'     :  RACINE DE Y MIS DANS X
C   OP = 'X=ABS(Y)'     :  VALEUR ABSOLUE DE Y MIS DANS X
C   OP = 'X=N(Y,Z)'     :  X MODULE DU VECTEUR DE COMPOSANTES Y ET Z
C   OP = 'X=Y+C   '     :  C AJOUTE A Y MIS DANS X
C   OP = 'X=X+C   '     :  C AJOUTE A X
C   OP = 'X=Y**C  '     :  Y A LA PUISSANCE C MIS DANS X
C   OP = 'X=COS(Y)'     :  COSINUS DE Y MIS DANS X
C   OP = 'X=SIN(Y)'     :  SINUS DE Y MIS DANS X
C   OP = 'X=ATN(Y)'     :  ARC TANGENTE DE Y MIS DANS X
C   OP = 'X=A(Y,Z)'     :  ATAN2(Y,Z) MIS DANS X
C   OP = 'X=+(Y,C)'     :  X = MAX DE Y ET DE C
C   OP = 'X=-(Y,C)'     :  X = MIN DE Y ET DE C
C   OP = 'X=+(Y,Z)'     :  X = MAX DE Y ET DE Z
C   OP = 'X=-(Y,Z)'     :  X = MIN DE Y ET DE Z
C   OP = 'X=YIFZ<C'     :  X = Y SI Z < C , POUR CHAQUE POINT
C   OP = 'X=C(Y-Z)'     :  X = C*(Y-Z)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |      OP        | -->| CHAINE DE CARACTERES INDIQUANT L'OPERATION
C |                |    | A EFFECTUER.
C |      X         |<-- | VECTEUR RESULTAT
C |      Y         | -->| VECTEUR OPERANDE
C |      Z         | -->| VECTEUR OPERANDE
C |      C         | -->| CONSTANTE DONNEE
C |      NPOIN     | -->| DIMENSION DES VECTEURS
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELANTS : BEAUCOUP
C PROGRAMMES APPELES   : NEANT
C
C PRECAUTIONS D'EMPLOI : LES OPERATIONS AVEC DIVISION
C                        SUPPRIMENT LES DIVISIONS PAR 0.
C                        UN PASSAGE CORRECT N'EST DONC PAS
C                        UNE PREUVE QUE Y OU Z N'EST JAMAIS NUL.
C                        POUR CES OPERATIONS, UTILISER PLUTOT OVD.
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NPOIN
      DOUBLE PRECISION, INTENT(IN)    :: Y(NPOIN),Z(NPOIN),C
      DOUBLE PRECISION, INTENT(INOUT) :: X(NPOIN)
      CHARACTER(LEN=8), INTENT(IN)    :: OP
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C
      INTRINSIC SQRT,ABS,COS,SIN,ATAN,MAX,MIN
C
C-----------------------------------------------------------------------
C
      SELECT CASE(OP(3:8))
C
C-----------------------------------------------------------------------
C
        CASE('C     ')
C
        DO I=1,NPOIN
          X(I) = C
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('0     ')
C
        DO I=1,NPOIN
          X(I) = 0.D0
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('Y     ')
C
        DO I=1,NPOIN
          X(I) = Y(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('+Y    ')
C
        DO I=1,NPOIN
          X(I) = Y(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('-Y    ')
C
        DO I=1,NPOIN
          X(I) = - Y(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('1/Y   ')
C
        DO I=1,NPOIN
          X(I) = 1.D0/Y(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('Y+Z   ')
C
        DO I=1,NPOIN
          X(I) = Y(I) + Z(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('Y-Z   ')
C
        DO I=1,NPOIN
          X(I) = Y(I) - Z(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('YZ    ')
C
        DO I=1,NPOIN
          X(I) = Y(I) * Z(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('-YZ   ')
C
        DO I=1,NPOIN
          X(I) = - Y(I) * Z(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('XY    ')
C
        DO I=1,NPOIN
          X(I) = X(I) * Y(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('X+YZ  ')
C
        DO I=1,NPOIN
          X(I) = X(I) + Y(I) * Z(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('X-YZ  ')
C
        DO I=1,NPOIN
          X(I) = X(I) - Y(I) * Z(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('CXY   ')
C
         DO I=1,NPOIN
           X(I) = C * X(I) * Y(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('CYZ   ')
C
        DO I=1,NPOIN
          X(I) = C * Y(I) * Z(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('CXYZ  ')
C
        DO I=1,NPOIN
          X(I) = C * X(I) * Y(I) * Z(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('X+CYZ ')
C
        DO I=1,NPOIN
          X(I) = X(I) + C * Y(I) * Z(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('Y/Z   ')
C
        DO I=1,NPOIN
          X(I) = Y(I) / Z(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('CY/Z  ')
C
        DO I=1,NPOIN
          X(I) = C*Y(I) / Z(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('CXY/Z ')
C
        DO I=1,NPOIN
          X(I) = C*X(I)*Y(I) / Z(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('X+CY/Z')
C
        DO I=1,NPOIN
          X(I) = X(I) + C * Y(I) / Z(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('X+Y   ')
C
        DO I=1,NPOIN
          X(I) = X(I) + Y(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('X-Y   ')
C
        DO I=1,NPOIN
          X(I) = X(I) - Y(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('CX    ')
C
        DO I=1,NPOIN
          X(I) = C * X(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('CY    ')
C
        DO I=1,NPOIN
          X(I) = C * Y(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('Y+CZ  ')
C
        DO I=1,NPOIN
          X(I) = Y(I) + C * Z(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('X+CY  ')
C
        DO I=1,NPOIN
          X(I) = X(I) + C * Y(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('SQR(Y)')
C
        DO I=1,NPOIN
          X(I) = SQRT(Y(I))
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('ABS(Y)')
C
        DO I=1,NPOIN
          X(I) = ABS(Y(I))
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('N(Y,Z)')
C
        DO I=1,NPOIN
          X(I) = SQRT( Y(I)**2 + Z(I)**2 )
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('Y+C   ')
C
        DO I=1,NPOIN
          X(I) = Y(I) + C
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('X+C   ')
C
        DO I=1,NPOIN
          X(I) = X(I) + C
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('Y**C  ')
C
        DO I=1,NPOIN
          IF(Y(I).GE.0.D0) THEN
            X(I) = Y(I)**C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,100) 
            IF (LNG.EQ.2) WRITE(LU,101) 
100         FORMAT(1X,'OV (BIEF) : Y**C INTERDIT SI Y < 0')
101         FORMAT(1X,'OV (BIEF): Y**C FORBIDDEN IF Y < 0')
            CALL PLANTE(1)
            STOP
          ENDIF
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('COS(Y)')
C
        DO I=1,NPOIN
          X(I) = COS(Y(I))
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('SIN(Y)')
C
        DO I=1,NPOIN
          X(I) = SIN(Y(I))
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('ATN(Y)')
C
        DO I=1,NPOIN
          X(I) = ATAN(Y(I))
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('A(Y,Z)')
C
        DO I=1,NPOIN
          X(I) = ATAN2(Y(I),Z(I))
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('+(Y,C)')
C
        DO I=1,NPOIN
          X(I) = MAX(Y(I),C)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('-(Y,C)')
C
        DO I=1,NPOIN
          X(I) = MIN(Y(I),C)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('+(Y,Z)')
C
        DO I=1,NPOIN
          X(I) = MAX(Y(I),Z(I))
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('-(Y,Z)')
C
        DO I=1,NPOIN
          X(I) = MIN(Y(I),Z(I))
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('YIFZ<C')
C
        DO I=1,NPOIN
          IF ( Z(I).LT.C ) X(I) = Y(I)
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE('C(Y-Z)')
C
        DO I=1,NPOIN
          X(I) = C*(Y(I)-Z(I))
        ENDDO
C
C-----------------------------------------------------------------------
C
        CASE DEFAULT
C
          IF (LNG.EQ.1) WRITE(LU,1000) OP
          IF (LNG.EQ.2) WRITE(LU,1001) OP
1000      FORMAT(1X,'OV (BIEF) : OPERATION INCONNUE: ',A8)
1001      FORMAT(1X,'OV (BIEF) : UNKNOWN OPERATION: ',A8)
          CALL PLANTE(1)
          STOP
C
C-----------------------------------------------------------------------
C
      END SELECT
C
C-----------------------------------------------------------------------
C
      RETURN
      END
