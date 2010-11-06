C                       *******************
                        SUBROUTINE GRADNODT
C                       *******************
C
     *(NS,NT,NU,AIRT,AIRS,H,T,DPX,DPY,DJX,DJY,
     * DX,DY,DIFT,CVIST,CE,DTT)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.4                              INRIA
C
C***********************************************************************
C
C   FONCTION  : CALCUL DES GRADIENTS PAR TRIANGLES ET PAR NOEUD 
C               ET DU TERME DE DIFFUSION  POUR  TRACEUR
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |  NS            | -->|  NOMBRE DE POINTS DU MAILLAGE                |
C |  NT            | -->|  NOMBRE D'ELEMENTS DU MAILLAGE               |
C |  NU            | -->|  NUMEROS DES NOEUDS PAR TRIANGLE             |
C |  AIRT          | -->|  AIRES DES TRIANGLES                         |
C |  AIRS          | -->|  AIRES DES CELLULES                          |
C |  H             | -->|  HAUTEURS D'EAU                              |
C |  T             | -->|  TRACEURS                                    |
C |  DPX, DPY      ! -->!  GRADIENTS DES FONCTIONS DE BASE             !
C |  DJX,DJY       |<-- |  GRADIENTS PAR TRIANGLES                     |
C |  DX,DY         |<-- |  GRADIENTS PAR NOEUDS                        |
C |  DIFT          | -->|  LOGIQUE INDIQUANT S'IL Y A DIFFUSION TRACEUR|
C |  CVIST         | -->|  COEFFICIENT DE DIFFUSION DU TRACEUR         |
C |  CE            |<-- |  TERME DE DIFFUSION                          |
C !  DTT           ! -->!  PAS DE TEMPS TRACEUR                        !
C !________________!____!______________________________________________!
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)   
C        -- (TABLEAU DE TRAVAIL)                                        
C-----------------------------------------------------------------------
C
C     - SOUS PROGRAMME(S) APPELANT : MAJTRAC                            
C 
C***********************************************************************
C
      IMPLICIT NONE
C     
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NS,NT
      INTEGER, INTENT(IN)             :: NU(NT,3)
      LOGICAL, INTENT(IN)             :: DIFT 
      DOUBLE PRECISION, INTENT(IN)    :: DPX(3,NT),DPY(3,NT)
      DOUBLE PRECISION, INTENT(IN)    :: AIRT(NT),AIRS(NS),H(NS),T(NS) 
      DOUBLE PRECISION, INTENT(INOUT) :: DJX(NT),DJY(NT),DX(NS),DY(NS)
      DOUBLE PRECISION, INTENT(INOUT) :: CE(NS)
      DOUBLE PRECISION, INTENT(IN)    :: DTT,CVIST
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C               
      INTEGER IS,JT,NUBO1,NUBO2,NUBO3            
      DOUBLE PRECISION AIRJ,UA1,UA2,UA3,AIS,HTT,AUX      
C
C-----------------------------------------------------------------------
C
C     INITIALIZING THE HERMITIAN NODAL GRADIENTS
C
      DO IS=1,NS
        DX(IS)  = 0.D0
        DY(IS)  = 0.D0
      ENDDO
C
C     LOOP ON GLOBAL LIST OF TRIANGLES
C
      DO JT=1,NT
C
         NUBO1 = NU(JT,1)
         NUBO2 = NU(JT,2)
         NUBO3 = NU(JT,3)
C
         AIRJ=   AIRT(JT)
         HTT = H(NUBO1)+H(NUBO2)+H(NUBO3)
         AUX =  CVIST*DTT*AIRJ*HTT/3.
C
C        COMPUTATION OF THE P1-GRADIENTS
C
           UA1=T(NUBO1)
           UA2=T(NUBO2)
           UA3=T(NUBO3)
C
C  GRADIENTS PAR TRIANGLES
C
            DJX(JT)      = UA1*DPX(1,JT) +
     &               UA2*DPX(2,JT) + UA3*DPX(3,JT)
            DJY(JT)      = UA1*DPY(1,JT) +
     &               UA2*DPY(2,JT) + UA3*DPY(3,JT)
C
C  GRADIENTS PAR NOEUDS
C
         DX(NUBO1)    = DX(NUBO1) + AIRJ*DJX(JT)
         DX(NUBO2)    = DX(NUBO2) + AIRJ*DJX(JT)
         DX(NUBO3)    = DX(NUBO3) + AIRJ*DJX(JT)
C
         DY(NUBO1)    = DY(NUBO1) + AIRJ*DJY(JT)
         DY(NUBO2)    = DY(NUBO2) + AIRJ*DJY(JT)
         DY(NUBO3)    = DY(NUBO3) + AIRJ*DJY(JT)
C
C   TERME DE DIFFUSION
C
       IF(DIFT.AND.CVIST.NE.0.) THEN
         CE(NUBO1)       = CE(NUBO1) -AUX*
     &  (DJX(JT)*DPX(1,JT)+DJY(JT)*DPY(1,JT))
         CE(NUBO2)       = CE(NUBO2) -AUX*
     &  (DJX(JT)*DPX(2,JT)+DJY(JT)*DPY(2,JT))
         CE(NUBO3)       = CE(NUBO3) -AUX*
     &  (DJX(JT)*DPX(3,JT)+DJY(JT)*DPY(3,JT))
C
        ENDIF
C
       ENDDO
C
C     COMPLETING THE COMPUTATION OF THE NODAL GRADIENTS 
C
      DO IS=1,NS
         AIS = 1.D0/(3.D0*AIRS(IS))
         DX(IS)       = DX(IS)*AIS
         DY(IS)       = DY(IS)*AIS
      ENDDO
C
C-----------------------------------------------------------------------
C 
      RETURN
      END
