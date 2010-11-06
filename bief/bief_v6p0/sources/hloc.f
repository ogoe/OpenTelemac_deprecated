C                       ***************
                        SUBROUTINE HLOC
C                       ***************
C
     *(NPOIN,NSEG,NPTFR,NUBO,NBOR,VNOCL,XNEBOR,YNEBOR,AIRS,DTHAUT)
C
C***********************************************************************
C BIEF VERSION 5.4                                               INRIA
C
C***********************************************************************
C
C  FONCTION : COMPUTATION OF LOCAL SPACE STEP [ |Ci|/Sum(Lij) ]
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |  NPOIN         | -->|  NOMBRE DE POINTS DU MAILLAGE                |
C |  NSEG          | -->|  NOMBRE D'ARETES DU MAILLAGE                 |
C |  NPTFR         | -->|  NOMBRE DE POINTS FRONTIERE                  |
C |  NUBO          | -->|  NUMEROS GLOBAUX DES EXTREMITES D'UNE ARETE  |
C |  NBOR          | -->|  NUMEROS GLOBAUX DES POINTS FRONTIERE        |
C |  VNOCL         | -->|  NORMALE A l'INTERFACE                       |
C !                !    !   (2 PREMIERES COMPOSANTES) ET               |
C !                !    !   LONGUEUR DE CE SEGMENT (3IEME COMPOSANTE)  |
C |  XNEBOR,YNEBOR | -->|  NORMALE AUX POINTS FRONTIERE                |
C !  AIRS          ! -->!  AIRES DES CELLULES DU MAILLAGE.             !
C |  DTHAUT        |<-- |  |Ci|/Sum(Lij) UTILISE POUR CFL              | 
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
      INTEGER, INTENT(IN)           :: NSEG,NPOIN,NPTFR,NUBO(2,*)
      INTEGER, INTENT(IN)           :: NBOR(*) 
      DOUBLE PRECISION, INTENT(IN)  :: VNOCL(3,*),XNEBOR(*),YNEBOR(*)
      DOUBLE PRECISION, INTENT(IN)  :: AIRS(NPOIN)
      DOUBLE PRECISION, INTENT(OUT) :: DTHAUT(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C 
      INTEGER I,K,NSG,NUBO1,NUBO2
      DOUBLE PRECISION VNX,VNY,VNL
C
C-----------------------------------------------------------------------
C                                                         
C     INITIALIZATION
C
      DO I=1,NPOIN
        DTHAUT(I) = 0.D0
      ENDDO
C
      DO NSG=1,NSEG 
C
         NUBO1     = NUBO(1,NSG)
         NUBO2     = NUBO(2,NSG)
C
         DTHAUT(NUBO1) = DTHAUT(NUBO1) + VNOCL(3,NSG)
         DTHAUT(NUBO2) = DTHAUT(NUBO2) + VNOCL(3,NSG) 
C
      ENDDO
C
      DO K=1,NPTFR
       I=NBOR(K)
       VNX=XNEBOR(K+NPTFR)
       VNY=YNEBOR(K+NPTFR)
       VNL=SQRT(VNX**2+VNY**2)
       DTHAUT(I) = DTHAUT(I) + VNL
C
      ENDDO
C
      DO I=1,NPOIN
         DTHAUT(I) = AIRS(I)/ DTHAUT(I)
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
