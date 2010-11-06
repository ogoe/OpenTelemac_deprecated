C                       ***************
                        SUBROUTINE HREF
C                       ***************
C
C***********************************************************************
C TELEMAC 2D VERSION 5.2      01/03/90    J-M HERVOUET
C***********************************************************************
C
C  FONCTION  : CALCUL DE LA HAUTEUR DE REFERNCE POUR LES EQUATIONS
C              DE BOUSSINESQ
C
C              PAR DEFAUT ON PREND LA HAUTEUR INITIALE
C     
C              CE SOUS-PROGRAMME PEUT ETRE MODIFIE
C
C              ON PEUT METTRE PAR EXEMPLE LA HAUTEUR DE LINEARISATION
C
C              SI ON VEUT RETOMBER SUR SAINT-VENANT, ON PEUT METTRE
C              H0 = 0
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|_______________________________________________
C |      H0        |<-- | HAUTEUR DE REFERENCE
C |      H         | -->| HAUTEUR INITIALE
C |      X,Y,(Z)   | -->| COORDONNEES DU MAILLAGE (Z N'EST PAS EMPLOYE).
C |      ZF        | -->| FOND A MODIFIER.
C |      HAULIN    | -->| PROFONDEUR DE LINEARISATION
C |      COTINI    | -->| COTE INITIALE
C |      MESH      | -->| MAILLAGE
C |      PRIVE     | -->| TABLEAU PRIVE POUR L'UTILISATEUR.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT :
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC2D
C
      IMPLICIT NONE
      INTEGER LNG,LU                                                 
      COMMON/INFO/LNG,LU        
C
C-----------------------------------------------------------------------
C
      CALL OS( 'X=Y     ' , H0 , H , H , 0.D0 )
C     NEXT LINE WILL DEGENERATE BOUSSINESQ INTO SAINT-VENANT
C     CALL OS( 'X=C     ' , H0 , H , H , 0.D0 )
C
C-----------------------------------------------------------------------
C
      RETURN
      END       
 
