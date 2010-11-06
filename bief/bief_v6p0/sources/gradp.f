C                       ****************
                        SUBROUTINE GRADP
C                       ****************
C
     *(NS,NT,NU,AIRT,X,Y,DPX,DPY)
C
C***********************************************************************
C BIEF VERSION 5.4                                                INRIA
C
C***********************************************************************
C
C     FONCTION  : COMPUTATION OF THE BASIS FUNCTIONS GRADIENTS
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    NS          | -->|  NOMBRE DE POINTS DU MAILLAGE                |
C |    NT          | -->|  NOMBRE DE TRIANGLES DU MAILLAGE             |
C |    NU          | -->|  NUMEROS DES NOEUDS PAR TRIANGLE             |
C |    AIRT        | -->|  AIRES DES TRIANGLES                         |
C |    X,Y         | -->|  COORDONNEES DES NOEUDS DU MAILLAGE          |
C |    DPX,DPY     |<-- |  GRADIENT DES FONCTIONS DE BASE P1           |
C |                |    |  PAR TRIANGLE                                |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C
C APPELE PAR : INBIEF
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                 
      INTEGER, INTENT(IN)           :: NS,NT,NU(NT,3)  
      DOUBLE PRECISION, INTENT(IN)  :: X(NS),Y(NS),AIRT(NT)
      DOUBLE PRECISION, INTENT(OUT) :: DPX(3,NT),DPY(3,NT)  
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER JT ,NUBO1,NUBO2,NUBO3
      DOUBLE PRECISION AIRJI,X1,X2,X3,Y1,Y2,Y3
C                                                          
C-----------------------------------------------------------------------
C
      DO JT=1,NT
C
         NUBO1 = NU(JT,1)
         NUBO2 = NU(JT,2)
         NUBO3 = NU(JT,3)
C
         AIRJI = 0.5D0/AIRT(JT)
C
C        COMPUTATION OF THE P1-GRADIENTS
C
         X1 = X(NUBO1)
         Y1 = Y(NUBO1)
         X2 = X(NUBO2)
         Y2 = Y(NUBO2)
         X3 = X(NUBO3)
         Y3 = Y(NUBO3)
C
         DPX(1,JT) = AIRJI*(Y2-Y3)
         DPX(2,JT) = AIRJI*(Y3-Y1)
         DPX(3,JT) = AIRJI*(Y1-Y2)
         DPY(1,JT) = AIRJI*(X3-X2)
         DPY(2,JT) = AIRJI*(X1-X3)
         DPY(3,JT) = AIRJI*(X2-X1)
C                  
      ENDDO
C                                                          
C-----------------------------------------------------------------------
C
      RETURN
      END

