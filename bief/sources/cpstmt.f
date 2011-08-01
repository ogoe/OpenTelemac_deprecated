C                       *****************
                        SUBROUTINE CPSTMT
C                       *****************
C
     *( X , Y , TRANS )
C
C***********************************************************************
C BIEF VERSION 6.0      03/02/2010    J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C  FONCTION  : COPIE D'UNE STRUCTURE DE MATRICE SUR UNE AUTRE
C              X COPIEE SUR Y
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   X            | -->| STRUCTURE COPIEE.
C |   Y            |<-- | STRUCTURE SUR LAQUELLE ON COPIE.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_CPSTMT => CPSTMT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(IN)    :: X
      TYPE(BIEF_OBJ), INTENT(INOUT) :: Y
      LOGICAL, INTENT(IN), OPTIONAL :: TRANS
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELM1,IELM2,IELN1,IELN2
      LOGICAL TR
C
C-----------------------------------------------------------------------
C  ONLY MATRICES ARE TREATED HERE :
C-----------------------------------------------------------------------
C
      IF(X%TYPE.NE.3.OR.Y%TYPE.NE.3) THEN
        IF(LNG.EQ.1) WRITE(LU,200) X%NAME,X%TYPE,Y%NAME,Y%TYPE
        IF(LNG.EQ.2) WRITE(LU,201) X%NAME,X%TYPE,Y%NAME,Y%TYPE
200     FORMAT(1X,'CPSTMT : CAS NON PREVU POUR X ET Y :',/,1X,
     *            'X=',A6,' TYPE :',1I6                 ,/,1X,
     *            'Y=',A6,' TYPE :',1I6)
201     FORMAT(1X,'CPSTMT : FORBIDDEN CASE FOR X AND Y:',/,1X,
     *            'X=',A6,' TYPE :',1I6                 ,/,1X,
     *            'Y=',A6,' TYPE :',1I6)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(PRESENT(TRANS)) THEN
        TR = TRANS
      ELSE
        TR = .FALSE.
      ENDIF
C
      IF(.NOT.TR) THEN
        IELM1 = X%ELMLIN
        IELM2 = X%ELMCOL
      ELSE
        IELM1 = X%ELMCOL
        IELM2 = X%ELMLIN
      ENDIF
C
C CONTROL OF MEMORY SIZE FOR DIAGONAL AND EXTRA-DIAGONAL TERMS :
C
      IF(X%D%DIM1.GT.Y%D%MAXDIM1.OR.
     *   X%X%DIM2*X%X%DIM1.GT.Y%X%MAXDIM1*Y%X%MAXDIM2) THEN
        IELN1 = Y%ELMLIN
        IELN2 = Y%ELMCOL
        IF(LNG.EQ.1) WRITE(LU,400) X%NAME,IELM1,IELM2,Y%NAME,IELN1,IELN2
        IF(LNG.EQ.2) WRITE(LU,401) X%NAME,IELM1,IELM2,Y%NAME,IELN1,IELN2
        IF(LNG.EQ.1) WRITE(LU,402) X%TYPDIA,X%D%DIM1,X%TYPEXT,
     *                 X%X%DIM2*X%X%DIM1,
     *                 Y%TYPDIA,Y%D%MAXDIM1,
     *                 Y%TYPEXT,Y%X%MAXDIM1*Y%X%MAXDIM2
        IF(LNG.EQ.2) WRITE(LU,403) X%TYPDIA,X%D%DIM1,X%TYPEXT,
     *                 X%X%DIM2*X%X%DIM1,
     *                 Y%TYPDIA,Y%D%MAXDIM1,
     *                 Y%TYPEXT,Y%X%MAXDIM1*Y%X%MAXDIM2
 400    FORMAT(1X,'CPSTMT : CAS IMPOSSIBLE POUR X ET Y :',/,1X,
     *            'X=',A6,/,1X,'ELEMENTS ',1I6,' ET ',1I6,/,1X,
     *            'Y=',A6,/,1X,'ELEMENTS ',1I6,' ET ',1I6,/,1X,
     *            'Y EST PLUS PETITE QUE X')
 402    FORMAT(1X,'X A UNE DIAGONALE DE TYPE ',A1,/,1X,
     *            'AVEC UNE TAILLE DE ',1I6,/,1X,
     *            'DES TERMES EXTRADIAGONAUX DE TYPE ',A1,/,1X,
     *            'AVEC UNE TAILLE DE ',1I6,/,1X,
     *            'Y A UNE DIAGONALE DE TYPE ',A1,/,1X,
     *            'ET DE TAILLE MAXIMUM ',1I6,/,1X,
     *            'DES TERMES EXTRADIAGONAUX DE TYPE ',A1,/,1X,
     *            'ET UNE TAILLE MAXIMUM DE ',1I6)
 401    FORMAT(1X,'CPSTMT : FORBIDDEN CASE FOR X AND Y:',/,1X,
     *            'X=',A6,/,1X,'ELEMENTS ',1I6,' AND ',1I6,/,1X,
     *            'Y=',A6,/,1X,'ELEMENTS ',1I6,' AND ',1I6,/,1X,
     *            'Y IS SMALLER THAN X')
 403    FORMAT(1X,'X HAS A DIAGONAL OF TYPE ',A1,/,1X,
     *            'WITH A SIZE OF ',1I6,/,1X,
     *            'AND OFF-DIAGONAL TERMS OF TYPE ',A1,/,1X,
     *            'WITH A SIZE OF ',1I6,/,1X,
     *            'Y HAS A DIAGONAL OF TYPE ',A1,/,1X,
     *            'AND A MAXIMUM SIZE OF ',1I6,/,1X,
     *            'AND OFF-DIAGONAL TERMS OF TYPE ',A1,/,1X,
     *            'AND A MAXIMUM SIZE OF ',1I6)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C COPY OF TYPES OF ELEMENTS
C
      Y%ELMLIN = IELM1
      Y%ELMCOL = IELM2
C
C 4) COPIE DES TYPES DE DIAGONALES ET DE TERMES EXTRADIAGONAUX
C
      CALL CPSTVC(X%D,Y%D)
      CALL CPSTVC(X%X,Y%X)
C
C 5) COPIE DES CARACTERISTIQUES DE LA MATRICE
C
      Y%TYPDIA = X%TYPDIA
      Y%TYPEXT = X%TYPEXT
C
C-----------------------------------------------------------------------
C
      RETURN
      END
