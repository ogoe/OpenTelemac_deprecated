C               *****************
                SUBROUTINE CALCQB
C               *****************
     *(Q1,Q2,Q3)
C
C***********************************************************************
C
C  ARTEMIS VERSION 5.1   04/06/99    D. AELBRECHT (LNH) 01 30 87 74 12
C
C***********************************************************************
C
C    FONCTION:  CALCULE LA RACINE DE L'EQUATION TRANSCENDANTE EN QB,
C               TAUX DE VAGUES DEFERLANTES OU AYANT DEFERLE, DANS LE
C               CAS D'UNE HOULE ALEATOIRE:
C                                  
C                   1 - QB     ( HE ) 2
C                   ------ = - ( -- )
C                   LOG QB     ( HM )
C
C               HE: HAUTEUR ENERGETIQUE DE HOULE
C               HM: HAUTEUR CRITIQUE DE DEFERLEMENT D'APRES MICHE
C
C               LE CALCUL DE CETTE RACINE SE FAIT PAR LA METHODE 
C               DE LA CORDE SUR LA FONCTION F(X) = 1-X + Cste*LOG(X)
C               QUI EST CONCAVE SUR ]0,1[ 
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     Q1         | -->| EXTREMITE GAUCHE INITIALE DE LA CORDE        |
C |     Q2         | -->| EXTREMITE DROITE INITIALE DE LA CORDE        |
C |     Q3         |<-- | SOLUTION APPROCHEE DE QB                     |
C |________________|____|______________________________________________|
C
C-----------------------------------------------------------------------
C APPELE PAR              :  BERKHO
C
C***********************************************************************
C         
      IMPLICIT NONE
C
      DOUBLE PRECISION Q1,Q2,Q3,FQ1,FQ2,FQ3,EPSIQB,RAP
C
      EPSIQB = 1.D-4
      RAP = Q2
C
      IF(Q2.GE.1.D0) THEN
         Q3 = 1.D0
      ELSE
         FQ3 = 1000.D0
C
C 10      FQ1 = (1.D0-Q1)+RAP*LOG(Q1)
 10      FQ1 = (1.D0-Q1)+RAP*LOG(ABS(Q1))
         FQ2 = (1.D0-Q2)+RAP*LOG(ABS(Q2))
C         FQ2 = (1.D0-Q2)+RAP*LOG(Q2)
         IF (FQ1.GE.0.D0) THEN
            Q3 = Q1
            FQ3 = EPSIQB/10.D0
         ELSE
            Q3 = Q1 - FQ1*(Q2-Q1)/(FQ2-FQ1)
            FQ3 = (1.D0-Q3)+RAP*LOG(ABS(Q3))
C            FQ3 = (1.D0-Q3)+RAP*LOG(Q3)
            IF ((FQ3*FQ1).GT.0.D0) THEN
               Q1 = Q3
            ELSE
               Q2 = Q3
            ENDIF
         ENDIF
         IF(ABS(FQ3).GE.EPSIQB) GOTO 10
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
