C                       *****************
                        SUBROUTINE CPSTVC
C                       *****************
C
     *( X , Y )
C
C***********************************************************************
C BIEF VERSION 5.1            01/03/95    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C  FONCTION  : COPIE D'UNE STRUCTURE DE VECTEUR SUR UNE AUTRE
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
      USE BIEF, EX_CPSTVC => CPSTVC
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(IN)    :: X
      TYPE(BIEF_OBJ), INTENT(INOUT) :: Y
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER SIZEX,SIZEY
C
C-----------------------------------------------------------------------
C  LE CAS DES VECTEURS SEULEMENT EST TRAITE ICI :
C-----------------------------------------------------------------------
C
      IF(Y%TYPE.NE.2.OR.X%TYPE.NE.2) THEN
        IF(LNG.EQ.1) WRITE(LU,200) X%NAME,X%TYPE,Y%NAME,Y%TYPE
        IF(LNG.EQ.2) WRITE(LU,201) X%NAME,X%TYPE,Y%NAME,Y%TYPE
 200    FORMAT(1X,'CPSTVC : CAS NON PREVU POUR X ET Y :',/,1X,
     *            'X=',A6,' TYPE :',1I6                 ,/,1X,
     *            'Y=',A6,' TYPE :',1I6)
 201    FORMAT(1X,'CPSTVC : FORBIDDEN CASE FOR X AND Y:',/,1X,
     *            'X=',A6,' TYPE :',1I6                 ,/,1X,
     *            'Y=',A6,' TYPE :',1I6)
        CALL PLANTE(1)
        STOP
      ENDIF
C
      SIZEX = X%DIM1*X%DIM2
      SIZEY = Y%MAXDIM1*Y%MAXDIM2
C
      IF(SIZEX.GT.SIZEY) THEN
        IF(LNG.EQ.1) WRITE(LU,300) X%NAME,SIZEX,Y%NAME,SIZEY
        IF(LNG.EQ.2) WRITE(LU,301) X%NAME,SIZEX,Y%NAME,SIZEY
 300    FORMAT(1X,'CPSTVC : CAS NON PREVU POUR X ET Y:',/,1X,
     *            'X=',A6,' TAILLE         :',1I6,/,1X,
     *            'Y=',A6,' TAILLE MAXIMUM :',1I6)
 301    FORMAT(1X,'CPSTVC : FORBIDDEN CASE FOR X AND Y:',/,1X,
     *            'X=',A6,' SIZE        :',1I6,/,1X,
     *            'Y=',A6,' MAXIMUM SIZE:',1I6)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C     DISCRETISATION
      IF(Y%ELM.NE.X%ELM.AND.Y%STATUS.EQ.1) THEN
        IF(LNG.EQ.1) WRITE(LU,400) X%NAME,Y%NAME
        IF(LNG.EQ.2) WRITE(LU,401) X%NAME,Y%NAME
400     FORMAT(1X,'CPSTVC : COPIE DE ',A6,' INTERDITE SUR ',A6)
401     FORMAT(1X,'CPSTVC : COPY OF ',A6,' FORBIDDEN ON ',A6)
        CALL PLANTE(1)
        STOP
      ELSE
        Y%ELM = X%ELM
      ENDIF
C     PREMIERE DIMENSION DU VECTEUR
      Y%DIM1 = X%DIM1
C     DEUXIEME DIMENSION DU VECTEUR
      Y%DIM2 = X%DIM2
C     CAS DES VECTEURS DISCONTINUS
      Y%DIMDISC = X%DIMDISC
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
