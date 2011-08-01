C                       ******************
                        SUBROUTINE FLUCINT
C                       ******************
C
     *(NS,NSEG,DIMT,NUBO,G,X,Y,UA,TN,ZF,VNOCL,CE,
     * NORDRE,CMI,JMI,DJX,DJY,DX,DY,DJXT,DJYT,DXT,DYT,EPSWL)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.3               INRIA
C
C***********************************************************************
C
C     FONCTION  : EDGEWISE CONVECTIVE FLUXES COMPUTATION USING THE
C                 KINETIC BOLTZMANN FLUX.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |                | -- |  
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NS,NSEG,DIMT,NORDRE
      INTEGER, INTENT(IN) :: NUBO(2,NSEG)
      DOUBLE PRECISION, INTENT(IN) :: G,X(NS),Y(NS),VNOCL(3,NSEG),ZF(NS)
      DOUBLE PRECISION, INTENT(INOUT) :: CE(3,NS)
      DOUBLE PRECISION, INTENT(IN)    :: UA(3,NS),TN(DIMT),CMI(2,*)
      DOUBLE PRECISION, INTENT(IN)    :: JMI(*),EPSWL
      DOUBLE PRECISION, INTENT(IN)    :: DJX(3,*),DJY(3,*)
      DOUBLE PRECISION, INTENT(IN)    :: DX(3,*),DY(3,*)
      DOUBLE PRECISION, INTENT(IN)    :: DJXT(*),DJYT(*),DXT(*),DYT(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NSG,NUBO1,NUBO2,J,IVAR
C      
      DOUBLE PRECISION E2,RA2,RA3,ALP,ZF1,ZF2,ZI,ZJ,XNN,YNN,RNN      
      DOUBLE PRECISION UAS11,UAS12,UAS21,UAS22,UAS31,UAS32
      DOUBLE PRECISION GRADI(4),GRADJ(4),GRADIJ(4),GRADJI(4),GRIJ,GRJI
      DOUBLE PRECISION GRI,GRJ,AUX1,AUX2,GRI2,GRIJ2,GRJ2,GRJI2
      DOUBLE PRECISION AIX,AIY,AJX,AJY,HI,HI0,HIJ,CIJ,ETAI,UAS110,UAS120
      DOUBLE PRECISION HJ,HJ0,HJI,CJI,ETAJ,EXT1,EXT2,FLU11,FLU12,FLU41
      DOUBLE PRECISION A01,A11,A02,A12,UAS41,UAS42
      DOUBLE PRECISION UAS410,UAS420,BETA,BETA1
C                                                          
C-----------------------------------------------------------------------
C
      RA2  = SQRT(2.D0)
      RA3  = SQRT(1.5D0*G)
      ALP  = 0.5D0/RA3
C
      BETA = 0.333D0
      BETA1 = 1.D0+BETA
      E2    = 1.E-6
C
C     LOOP ON GLOBAL LIST OF EDGES
C
      DO 500 NSG=1,NSEG 
         J         = JMI(NSG)
C
         NUBO1     = NUBO(1,NSG)
         NUBO2     = NUBO(2,NSG)
C
         ZF1   =    ZF(NUBO1)
         ZF2   =    ZF(NUBO2)
C
         XNN       = VNOCL(1,NSG)
         YNN       = VNOCL(2,NSG)
         RNN       = VNOCL(3,NSG) 
C
         UAS11     = UA(1,NUBO1)
         UAS12     = UA(1,NUBO2)
         UAS21     = UA(2,NUBO1)
         UAS22     = UA(2,NUBO2)
         UAS31     = UA(3,NUBO1)
         UAS32     = UA(3,NUBO2)
         UAS41     = TN(NUBO1)
         UAS42     = TN(NUBO2)
C
         HI0=UAS11*UAS11
         HJ0=UAS12*UAS12
C
C ROTATION
C
         UAS21  = XNN*UAS21+YNN*UAS31
C
         UAS22  = XNN*UAS22+YNN*UAS32
C
C
         IF(NORDRE.EQ.1)  GOTO 1234
C
         AIX       = CMI(1,NSG)-X(NUBO1)
         AIY       = CMI(2,NSG)-Y(NUBO1) 
         AJX       = CMI(1,NSG)-X(NUBO2)
         AJY       = CMI(2,NSG)-Y(NUBO2) 
C
       DO IVAR=1,3
C
         GRADI(IVAR)  = AIX*DX(IVAR,NUBO1) + AIY*DY(IVAR,NUBO1) 
C
         GRADJ(IVAR)  = AJX*DX(IVAR,NUBO2) + AJY*DY(IVAR,NUBO2)
C
         GRADIJ(IVAR)  = AIX*DJX(IVAR,J) + AIY*DJY(IVAR,J) 
C
         GRADJI(IVAR)  = AJX*DJX(IVAR,J) + AJY*DJY(IVAR,J)
C
       ENDDO
C
         GRADI(4)  = AIX*DXT(NUBO1) + AIY*DYT(NUBO1) 
C
         GRADJ(4)  = AJX*DXT(NUBO2) + AJY*DYT(NUBO2)
C
         GRADIJ(4)  = AIX*DJXT(J) + AIY*DJYT(J) 
C
         GRADJI(4)  = AJX*DJXT(J) + AJY*DJYT(J)
C
C ROTATION DES GRADIENTS
C
C
         GRADI(2)  = XNN*GRADI(2)+YNN*GRADI(3)
C
         GRADIJ(2)  = XNN*GRADIJ(2)+YNN*GRADIJ(3)
C
         GRADJ(2)  = XNN*GRADJ(2)+YNN*GRADJ(3)
C
         GRADJI(2)  = XNN*GRADJI(2)+YNN*GRADJI(3)
C
       DO IVAR=1,4
        IF(IVAR.EQ.3) GOTO 100
C
         GRIJ = GRADIJ(IVAR)
         GRJI = GRADJI(IVAR)
C 
         GRI = BETA1*GRADI(IVAR) - BETA*GRIJ
         GRJ = BETA1*GRADJ(IVAR) - BETA*GRJI
C
         AUX1 = 0.5*(1.0+  DSIGN(1.0D0, GRI*GRIJ))
         AUX2 = 0.5*(1.0 + DSIGN(1.0D0, GRJ*GRJI))
C
C    VAN ALBADA
C
         GRI2  = GRI*GRI    + E2
         GRIJ2 = GRIJ*GRIJ  + E2
         GRJ2  = GRJ*GRJ    + E2
         GRJI2 = GRJI*GRJI  + E2
C
         GRADI(IVAR)  = AUX1*
     &      (GRI2  *GRIJ  + GRIJ2  *GRI )/(GRI2 + GRIJ2)
C
         GRADJ(IVAR)  = AUX2*
     &      (GRJ2  *GRJI + GRJI2 *GRJ )/(GRJ2 + GRJI2)
C
 100     CONTINUE
       ENDDO
C
C
         UAS110 =   UAS11
         UAS11 = MIN(MAX(UAS110/RA2,UAS11 + GRADI(1)),RA2*UAS110)
         UAS21     = UAS21 + GRADI(2)
C
         UAS120 =   UAS12
         UAS12 = MIN(MAX(UAS120/RA2,UAS12 + GRADJ(1)),RA2*UAS120)
         UAS22     = UAS22 + GRADJ(2)
C
       IF(UAS11.LE.0.D0) THEN
         UAS11 =0.D0
         UAS21 =0.D0
       ENDIF
       IF(UAS12.LE.0.D0)  THEN
          UAS12 =0.D0
          UAS22 =0.D0
       ENDIF
C
C
1234   CONTINUE
C
         HI   = UAS11*UAS11
C
C   ETAI= COTE SURFACE LIBRE
C   SI HI < EPS  ETAI =MIN( ETAI, ETAJ) ==> FOND PLAT
C   ZI =COTE DU FOND ARTIFICIELLE
C
         ETAI = HI0+ZF1
         IF(HI0.LE.EPSWL) ETAI = MIN(ETAI,HJ0+ZF2)
         ZI= ETAI -HI0
C
C         HIJ= ETAI-MIN(ETAI,MAX(ZF1,ZF2))
         IF(ZF1.GE.ZF2) THEN
         HIJ=HI0
         ELSE
         HIJ=MAX(0.D0,ETAI-ZF2)
         ENDIF
C
         IF(HI.LE.0.D0) THEN
           CIJ=0.D0
           FLU11=0.D0
         ELSE
C
           IF(HIJ.LE.0.D0) THEN
             CIJ=0.D0
C
             IF(UAS21.GE.0.D0) THEN
                EXT1=-RA3
             ELSE
                EXT1=RA3
             ENDIF
C
           ELSE
             CIJ  =HIJ*SQRT(HIJ)/HI
             EXT1 = MIN(RA3,MAX(-RA3,-UAS21/CIJ))
           ENDIF
C
         A01  = ALP*(RA3-EXT1)
         A11  = ALP*(RA3**2-EXT1**2)/2.D0
C
         FLU11= HI*(UAS21*A01+CIJ*A11)
C
        ENDIF
C
C 
         HJ   = UAS12 *UAS12
C
         ETAJ = HJ0+ZF2
         IF(HJ0.LE.EPSWL) ETAJ = MIN(ETAJ,HI0+ZF1)
         ZJ=ETAJ-HJ0
C
C         HJI= ETAJ-MIN(ETAJ,MAX(ZF1,ZF2))
         IF(ZF2.GE.ZF1) THEN
         HJI=HJ0
         ELSE
         HJI=MAX(0.D0,ETAJ-ZF1)
         ENDIF
C
         IF(HJ.LE.0.D0) THEN
            CJI=0.D0
            FLU12=0.D0
         ELSE
C
            IF(HJI.LE.0.D0) THEN
              CJI=0.D0
C  
              IF(UAS22.GE.0.D0) THEN
                EXT2=-RA3
              ELSE
                EXT2=RA3
              ENDIF
C
            ELSE
              CJI  =HJI*SQRT(HJI)/HJ
              EXT2 = MIN(RA3,MAX(-RA3,-UAS22/CJI))
            ENDIF
C
         A02  = ALP*(RA3+EXT2)
         A12  = ALP*(EXT2**2-RA3**2)/2.D0
C
         FLU12= HJ*(UAS22*A02+CJI*A12)
        ENDIF
C
C
         FLU11=(FLU11 + FLU12)*RNN
C
         IF (FLU11.GE.0.D0) THEN 
       IF(NORDRE.GE.2) THEN
         UAS410 = UAS41
         UAS41 = MIN(MAX(0.5D0*UAS410,UAS41 + GRADI(4)),2.D0*UAS410)
        ENDIF
         FLU41 =  UAS41 * FLU11
         ELSE
       IF(NORDRE.GE.2) THEN
         UAS420 = UAS42
         UAS42 = MIN(MAX(0.5D0*UAS420,UAS42 + GRADJ(4)),2.D0*UAS420)
        ENDIF
         FLU41 =  UAS42 * FLU11
        ENDIF
C
         CE(1,NUBO1) = CE(1,NUBO1) - FLU41
C
         CE(1,NUBO2) = CE(1,NUBO2) + FLU41  
C
500   CONTINUE
C                                                         
C-----------------------------------------------------------------------
C
      RETURN
      END
