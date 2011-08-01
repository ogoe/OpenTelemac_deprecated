C                       ******************
                        SUBROUTINE GRADNOD
C                       ******************
C
     *(NS,NT,NU,AIRT,AIRS,UA,DPX,DPY,DJX,DJY,DX,DY,IVIS,CVIS,CE,ZF)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.4                              INRIA
C
C***********************************************************************
C
C   FONCTION  : CALCUL DES GRADIENTS PAR TRIANGLES ET PAR NOEUD 
C               ET DES TERMES DE DIFFUSION
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
C |  UA            | -->|  UA(1,IS) = H,  UA(2,IS)=U  ,UA(3,IS)=V      |
C |  DPX, DPY      ! -->!  GRADIENTS DES FONCTIONS DE BASE             !
C |  DJX,DJY       |<-- |  GRADIENTS PAR TRIANGLES                     |
C |  DX,DY         |<-- |  GRADIENTS PAR NOEUDS                        |
C !  IVIS          ! -->!  OPTION DIFFUSION DES VITESSES               !
C !  CVIS          ! -->!  COEFFICIENT DE DIFFUSION DES VITESSES       !
C |  CE            |<-- |  TERMES DE DIFFUSION                         |
C |  ZF            | -->|  COTES DU FOND                               |
C !________________!____!______________________________________________!
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)   
C        -- (TABLEAU DE TRAVAIL)                                        
C-----------------------------------------------------------------------
C     - SOUS PROGRAMME(S) APPELANT : FLUHYD                             
C 
C***********************************************************************
C
      IMPLICIT NONE 
C  
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NS,NT,IVIS
      INTEGER, INTENT(IN)             :: NU(NT,3)
      DOUBLE PRECISION, INTENT(IN)    :: AIRT(NT),AIRS(NS),CVIS
      DOUBLE PRECISION, INTENT(INOUT) :: DJX(3,NT),DJY(3,NT)
      DOUBLE PRECISION, INTENT(INOUT) :: DX(3,NS),DY(3,NS)
      DOUBLE PRECISION, INTENT(INOUT) :: CE(3,NS)
      DOUBLE PRECISION, INTENT(IN)    :: UA(3,NS),ZF(NS)
      DOUBLE PRECISION, INTENT(IN)    :: DPX(3,NT),DPY(3,NT) 
C     
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                          
      INTEGER IS,JT,NUBO1,NUBO2,NUBO3,IVAR                  
      DOUBLE PRECISION AIRJ,UA1,UA2,UA3,AIS,HTT,AUX
C
C-----------------------------------------------------------------------
C
C     INITIALIZING THE HERMITIAN NODAL GRADIENTS
C
      DO IS=1,NS
        DO IVAR=1,3
          DX(IVAR,IS) = 0.D0
          DY(IVAR,IS) = 0.D0
        ENDDO
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
         HTT = UA(1,NUBO1)+UA(1,NUBO2)+UA(1,NUBO3)
         AUX = CVIS*AIRJ*HTT/3.D0
C
C        COMPUTATION OF THE P1-GRADIENTS
C
C
C   ON CALCULE LE GRADIENT DE H+Z
C
       IVAR=1
           UA1=UA(IVAR,NUBO1) + ZF(NUBO1)
           UA2=UA(IVAR,NUBO2) + ZF(NUBO2)
           UA3=UA(IVAR,NUBO3) + ZF(NUBO3)
C
            DJX(IVAR,JT)      = UA1*DPX(1,JT) +
     &               UA2*DPX(2,JT) + UA3*DPX(3,JT)
            DJY(IVAR,JT)      = UA1*DPY(1,JT) +
     &               UA2*DPY(2,JT) + UA3*DPY(3,JT)
C
         DX(IVAR,NUBO1)    = DX(IVAR,NUBO1) + AIRJ*DJX(IVAR,JT)
         DX(IVAR,NUBO2)    = DX(IVAR,NUBO2) + AIRJ*DJX(IVAR,JT)
         DX(IVAR,NUBO3)    = DX(IVAR,NUBO3) + AIRJ*DJX(IVAR,JT)
C
         DY(IVAR,NUBO1)    = DY(IVAR,NUBO1) + AIRJ*DJY(IVAR,JT)
         DY(IVAR,NUBO2)    = DY(IVAR,NUBO2) + AIRJ*DJY(IVAR,JT)
         DY(IVAR,NUBO3)    = DY(IVAR,NUBO3) + AIRJ*DJY(IVAR,JT)
C
C    CALCUL DES GRADIENTS DE VITESSES
C
         DO IVAR=2,3
C
           UA1=UA(IVAR,NUBO1)
           UA2=UA(IVAR,NUBO2)
           UA3=UA(IVAR,NUBO3)
C
            DJX(IVAR,JT)      = UA1*DPX(1,JT) +
     &               UA2*DPX(2,JT) + UA3*DPX(3,JT)
            DJY(IVAR,JT)      = UA1*DPY(1,JT) +
     &               UA2*DPY(2,JT) + UA3*DPY(3,JT)
C
         DX(IVAR,NUBO1)    = DX(IVAR,NUBO1) + AIRJ*DJX(IVAR,JT)
         DX(IVAR,NUBO2)    = DX(IVAR,NUBO2) + AIRJ*DJX(IVAR,JT)
         DX(IVAR,NUBO3)    = DX(IVAR,NUBO3) + AIRJ*DJX(IVAR,JT)
C
         DY(IVAR,NUBO1)    = DY(IVAR,NUBO1) + AIRJ*DJY(IVAR,JT)
         DY(IVAR,NUBO2)    = DY(IVAR,NUBO2) + AIRJ*DJY(IVAR,JT)
         DY(IVAR,NUBO3)    = DY(IVAR,NUBO3) + AIRJ*DJY(IVAR,JT)
         ENDDO
C
C  CALCUL DES TERMES DE DIFFUSION DES VITESSES
C
      IF(IVIS.EQ.0.OR.CVIS.EQ.0.) GOTO 10
         CE(2,NUBO1)       = CE(2,NUBO1) -AUX*
     &  (DJX(2,JT)*DPX(1,JT)+DJY(2,JT)*DPY(1,JT))
         CE(2,NUBO2)       = CE(2,NUBO2) -AUX*
     &  (DJX(2,JT)*DPX(2,JT)+DJY(2,JT)*DPY(2,JT))
         CE(2,NUBO3)       = CE(2,NUBO3) -AUX*
     &  (DJX(2,JT)*DPX(3,JT)+DJY(2,JT)*DPY(3,JT))
C
         CE(3,NUBO1)       = CE(3,NUBO1) -AUX*
     &  (DJX(3,JT)*DPX(1,JT)+DJY(3,JT)*DPY(1,JT))
         CE(3,NUBO2)       = CE(3,NUBO2) -AUX*
     &  (DJX(3,JT)*DPX(2,JT)+DJY(3,JT)*DPY(2,JT))
         CE(3,NUBO3)       = CE(3,NUBO3) -AUX*
     &  (DJX(3,JT)*DPX(3,JT)+DJY(3,JT)*DPY(3,JT))
 10       CONTINUE
         ENDDO
C
C     COMPLETING THE COMPUTATION OF THE NODAL GRADIENTS 
C
      DO IS=1,NS
C
         AIS = 1.D0/(3.D0*AIRS(IS))
C
         DO IVAR=1,3
           DX(IVAR,IS) = DX(IVAR,IS)*AIS
           DY(IVAR,IS) = DY(IVAR,IS)*AIS
         ENDDO
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
