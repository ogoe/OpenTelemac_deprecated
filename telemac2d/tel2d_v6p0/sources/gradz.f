C                       ****************                               
                        SUBROUTINE GRADZ
C                       ****************                               
C                                                                       
     *(NS,NT,NSEG,NU,NUBO,X,Y,AIRT,AIRS,CMI,JV,
     * ZF,DPX,DPY,DSZ,BETA,AIRST,DXIZ,DYIZ,DSP,DSM,CORR)
C                                                                       
C***********************************************************************
C  TELEMAC 2D VERSION 5.4                                         INRIA
C-----------------------------------------------------------------------
C
C                 SAINT-VENANT    CINETIQUE 
C  FONCTION : CALCUL DES VARIATIONS DE Z POUR ORDRE 2
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C |  NS            | -->|  NOMBRE DE POINTS DU MAILLAGE                |
C |  NT            | -->|  NOMBRE D'ELEMENTS DU MAILLAGE               |
C |  NSEG          | -->|  NOMBRE D'ARETES DU MAILLAGE                 |
C |  NU            | -->|  NUMEROS DES NOEUDS PAR TRIANGLE             |
C !  NUBO          ! -->!  NUMEROS DES DEUX SOMMETS D'UNE ARETE        !
C |  X,Y           | -->|  COORDONNEES DES NOEUDS DU MAILLAGE          |
C |  AIRT          | -->|  AIRES DES TRIANGLES                         |
C |  AIRS          | -->|  AIRES DES CELLULES                          |
C !  CMI           ! -->!  COORDONNEES DES POINTS MILIEUX D'INTERFACE  !
C !  JV            ! -->!  NUMERO DU TRIANGLE AUQUEL APPARTIENT LE     !
C !                !    !  POINT MILIEU D'UNE INTERFACE                !
C |  ZF            | -->|  COTES DU FOND                               |
C |  DPX,DPY       | -->|  GRADIENT DES FONCTIONS DE BASE P1           |
C |                |    |     PAR TRIANGLE                             |
C |  DSZ           !<-- !  VARIATION DE Z POUR ORDRE 2                 !
C |  BETA          ! -- !  COEFFICIENT EXTRAPOLATION                   !
C |  AIRST         ! -->! AIRES DES SOUS-TRIANGLES DANS CELLULES       !
C |  DXIZ,DYIZ,DSP | -- | TABLEAUX DE TRAVAIL                          |
C |  DSM,CORR      | -- !                                              !
C !________________!____!______________________________________________!
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)  | 
C        -- (TABLEAU DE TRAVAIL)                                       | 
C-----------------------------------------------------------------------
C     - SOUS PROGRAMME(S) APPELANT : RESOLU                            |
C***********************************************************************
C                                                                       
      IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NS,NT,NSEG
      INTEGER, INTENT(IN)             :: NU(NT,3),NUBO(2,NSEG),JV(*) 
      DOUBLE PRECISION, INTENT(INOUT) :: DSZ(2,*)
      DOUBLE PRECISION, INTENT(IN)    :: X(NS),Y(NS),AIRT(NT),AIRS(NS)
      DOUBLE PRECISION, INTENT(INOUT) :: DXIZ(NS),DYIZ(NS)
      DOUBLE PRECISION, INTENT(INOUT) :: DSP(NS),DSM(NS),CORR(NS),BETA
      DOUBLE PRECISION, INTENT(IN)    :: DPX(3,NT),DPY(3,NT)
      DOUBLE PRECISION, INTENT(IN)    :: CMI(2,*),AIRST(2,*),ZF(NS)
C  
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                          
      INTEGER IS,I1,I2,I3,JT,J,NSG,NUBO1,NUBO2,ILIM      
      DOUBLE PRECISION AIRJ,DXTZ,DYTZ,AIX,AIY,AJX,AJY
      DOUBLE PRECISION ZF1,ZF2,GRADI,GRADJ,GRIJ,GRJI,AMDS,DSH
C
      DOUBLE PRECISION EXLIM
      EXTERNAL         EXLIM
C
C-----------------------------------------------------------------------
C
      DO IS=1,NS
         DXIZ(IS) = 0.D0
         DYIZ(IS) = 0.D0
      ENDDO
C
      DO JT=1,NT
C
         I1 = NU(JT,1)
         I2 = NU(JT,2)
         I3 = NU(JT,3)
C
         AIRJ =   AIRT(JT)
         DXTZ =ZF(I1)*DPX(1,JT) +ZF(I2)*DPX(2,JT) + ZF(I3)*DPX(3,JT) 
         DYTZ =ZF(I1)*DPY(1,JT) +ZF(I2)*DPY(2,JT) + ZF(I3)*DPY(3,JT)
C
         DXIZ(I1) = DXIZ(I1) + AIRJ*DXTZ
         DXIZ(I2) = DXIZ(I2) + AIRJ*DXTZ
         DXIZ(I3) = DXIZ(I3) + AIRJ*DXTZ
C
         DYIZ(I1) = DYIZ(I1) + AIRJ*DYTZ
         DYIZ(I2) = DYIZ(I2) + AIRJ*DYTZ
         DYIZ(I3) = DYIZ(I3) + AIRJ*DYTZ
      ENDDO
C
      DO IS=1,NS
         DXIZ(IS) = DXIZ(IS)/(3.D0*AIRS(IS))
         DYIZ(IS) = DYIZ(IS)/(3.D0*AIRS(IS))
         DSP(IS)  = 0.D0
         DSM(IS)  = 0.D0
      ENDDO
C
C    RECONSTRUCTION PAR INTERFACE
C
      DO NSG=1,NSEG 
C
         J         = JV(NSG)
C
         NUBO1     = NUBO(1,NSG)
         NUBO2     = NUBO(2,NSG)
C
         ZF1   =    ZF(NUBO1)
         ZF2   =    ZF(NUBO2)
C
         AIX       = CMI(1,NSG)-X(NUBO1)
         AIY       = CMI(2,NSG)-Y(NUBO1) 
         AJX       = CMI(1,NSG)-X(NUBO2)
         AJY       = CMI(2,NSG)-Y(NUBO2) 
C
C        GRADIENT AUX NOEUDS
C
         GRADI  = AIX*DXIZ(NUBO1) + AIY*DYIZ(NUBO1) 
C
         GRADJ  = AJX*DXIZ(NUBO2) + AJY*DYIZ(NUBO2)
C
C
         I1 = NU(J,1)
         I2 = NU(J,2)
         I3 = NU(J,3)
C
C        GRADIENT PAR TRIANGLE
C
         DXTZ =ZF(I1)*DPX(1,J) +ZF(I2)*DPX(2,J) + ZF(I3)*DPX(3,J) 
         DYTZ =ZF(I1)*DPY(1,J) +ZF(I2)*DPY(2,J) + ZF(I3)*DPY(3,J)
C
         GRIJ  = AIX*DXTZ + AIY*DYTZ 
C
         GRJI  = AJX*DXTZ + AJY*DYTZ
C
C    EXTRAPOLATION ET LIMITEUR
C
       ILIM=1
       BETA=1.D0
         DSZ(1,NSG)  =  EXLIM(ILIM,BETA,GRADI,GRIJ )
         DSZ(2,NSG)  =  EXLIM (ILIM,BETA,GRADJ,GRJI )
C
         IF(DSZ(1,NSG).GE.0.D0) THEN
         DSP(NUBO1) = DSP(NUBO1) + AIRST(1,NSG)*DSZ(1,NSG)
         ELSE
         DSM(NUBO1) = DSM(NUBO1) - AIRST(1,NSG)*DSZ(1,NSG)
         ENDIF
         IF(DSZ(2,NSG).GE.0.) THEN
         DSP(NUBO2) = DSP(NUBO2) + AIRST(2,NSG)*DSZ(2,NSG)
         ELSE
         DSM(NUBO2) = DSM(NUBO2) - AIRST(2,NSG)*DSZ(2,NSG)
         ENDIF
C
       ENDDO
C
C  CALCUL DES CORRECTEURS POUR AVOIR LA CONSERVATION
C
      DO IS=1,NS
       CORR(IS) =  DSM(IS) - DSP(IS)
       AMDS =MAX(DSP(IS),DSM(IS))
        IF(AMDS.GT.0.D0) THEN
        CORR(IS) = CORR(IS)/AMDS
        ENDIF
      ENDDO 
C
      DO NSG=1,NSEG 
C
         NUBO1 = NUBO(1,NSG)
         NUBO2 = NUBO(2,NSG)
C
         DSH =  DSZ(1,NSG)  
         DSZ(1,NSG) =   DSH +
     & MIN(0.D0,CORR(NUBO1))*MAX(0.D0,DSH)+
     & MAX(0.D0,CORR(NUBO1))*MAX(0.D0,-DSH)
C      
         DSH     =  DSZ(2,NSG)  
         DSZ(2,NSG) =   DSH +
     & MIN(0.D0,CORR(NUBO2))*MAX(0.D0,DSH)+
     & MAX(0.D0,CORR(NUBO2))*MAX(0.D0,-DSH)
      ENDDO 
C
C-----------------------------------------------------------------------
C      
      RETURN
      END
