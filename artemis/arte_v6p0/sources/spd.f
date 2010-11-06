C                       *****************************
                        DOUBLE PRECISION FUNCTION SPD
C                       *****************************
C
     *(TETA)
C
C***********************************************************************
C
C  ARTEMIS VERSION 2.0     01/06/93   F. LEPEINTRE (LNH) 01 30 87 78 54
C  ARTEMIS VERSION 5.1     04/06/99   D. AELBRECHT (LNH) 01 30 87 74 12
C
C***********************************************************************
C
C      FONCTION:    CALCULE LA DENSITE D'ENERGIE SUIVANT LA FORMULE
C                   DE GODA:RANDOM SEA AND DESIGN OF MARITIME STRUCTURES
C                           UNIVERSITY OF TOKYO PRESS - 1985
C
C
C SPD(TETA) = COS( (TETA)/2 )**(2*EXPO)
C
C
C OU TETA EST L'ANGLE DE PROPAGATION DE LA HOULE (LA DIRECTION
C    PRINCIPALE DE PROPAGATION EST TETA=0)
C    EXPO EST UN EXPOSANT DONNE PAR L'UTILISATEUR
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   TETA         | -->|  ANGLE DE PROPAGATION DE LA HOULE            |
C |                |    |  D'ENERGIE                                   |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : DIRALE
C
C***********************************************************************
C
c      USE INTERFACE_ARTEMIS                    
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      DOUBLE PRECISION TETA,PI,DEGRAD
C
      DOUBLE PRECISION EXPO
      COMMON /COEFHD/ EXPO
C
      INTRINSIC COS
C
C-----------------------------------------------------------------------
C
      PARAMETER( PI = 3.1415926535897932384626433D0 ,
     *           DEGRAD = PI/180.D0 )
C
C-----------------------------------------------------------------------
C
      SPD = COS ( TETA*DEGRAD / 2.D0 )**(2*EXPO)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
