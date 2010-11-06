C                       **************                              
                        SUBROUTINE CDL
C                       **************                               
C                                                                       
     *(NS,NPTFR,NBOR,LIMPRO,XNEBOR,YNEBOR,KDIR,KNEU,G,HBOR,
     * UBOR,VBOR,UA,CE,FLUENT,FLUSORT,DTHAUT,DT,CFL,FLUHBTEMP,NTRAC)
C                                                                       
C***********************************************************************
C  TELEMAC 2D VERSION 5.8                                          INRIA
C-----------------------------------------------------------------------
C                 SAINT-VENANT    CINETIQUE 
C   
C     COMPUTATION OF THE CONVECTIVE FLUXES AT BOUNDARIES
C     
C
C     UA(1,IS) = H,  UA(2,IS)=U  ,UA(3,IS)=V
C
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C |  NS            | -->|  NOMBRE DE POINTS DU MAILLAGE                |
C |  NPTFR         | -->|  NOMBRE DE POINTS FRONTIERE                  |
C |  NBOR          | -->|  NUMEROS GLOBAUX DES POINTS DE BORD          | 
C |  LIMPRO        | -->|  TYPES DE CONDITIONS AUX LIMITES             |
C |  XNEBOR,YNEBOR | -->|  NORMALE AUX POINTS FRONTIERE                |
C |  KDIR          | -->|  CONVENTION POUR LES POINTS DIRICHLET        |
C |  KNEU          | -->|  CONVENTION POUR LES POINTS NEUMANN          |
C |  G             | -->|  CONSTANTE DE GRAVITE                        |
C |  HBOR          | -->|  VALEURS IMPOSEES DE H                       | 
C |  UBOR          | -->|  VALEURS IMPOSEES DE U                       | 
C |  VBOR          | -->|  VALEURS IMPOSEES DE V                       | 
C |  UA            | -->|  UA(1,IS) = H,  UA(2,IS)=U  ,UA(3,IS)=V      |
C |  CE            |<-->|  FLUX                                        |
C !  FLUENT,FLUSORT|<-- |  FLUX MASSE ENTREE ET SORTIE                 | 
C |  DTHAUT        ! -->!  UTILISE POUR CONDITION CFL                  !
C |  DT            |<-->|  PAS DE TEMPS                                |
C !  CFL           ! -->!  NOMBRE DE CFL                               !
C !  FLUHBTEMP     !<-- !  FLUX BORD POUR TRACEUR                      !
C |  TRAC          | -->|  LOGIQUE INDIQUANT LA PRESENCE D'UN TRACEUR  |
C !________________!____!______________________________________________!
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)   
C        -- (TABLEAU DE TRAVAIL)                                        
C-----------------------------------------------------------------------
C     - SOUS PROGRAMME(S) APPELANT : FLUHYD                              
C 
C***********************************************************************
C 
      USE BIEF
C  
      IMPLICIT NONE
C    
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NS,NPTFR,KDIR,KNEU,NTRAC
      INTEGER, INTENT(IN)             :: NBOR(NPTFR),LIMPRO(NPTFR,6)
      DOUBLE PRECISION, INTENT(IN)    :: XNEBOR(2*NPTFR),YNEBOR(2*NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: HBOR(NPTFR),UA(3,NS),DTHAUT(*)
      DOUBLE PRECISION, INTENT(IN)    :: UBOR(NPTFR),VBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: G,CFL
      DOUBLE PRECISION, INTENT(INOUT) :: DT
      DOUBLE PRECISION, INTENT(INOUT) :: CE(3,NS),FLUENT,FLUSORT
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: FLUHBTEMP
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IS,K,NIT,ITRAC      
C
      DOUBLE PRECISION RA3,RA32,RA33, ALP,ALP2,ALP3,SG,SQ2      
      DOUBLE PRECISION VNX,VNY,VNX1,VNY1,VNL,H,U,V,RUN
      DOUBLE PRECISION FLUH,FLUU,FLUV,AUX,FLUTMP,RH,HRH,UNN,VNN
      DOUBLE PRECISION FHPLUS,FUPLUS,FHMOINS,FUMOINS 
      DOUBLE PRECISION A,A1,A2,A3,ALPHA0,ALPHA1,ALPHA2,C,VP1,VP2 ,VP3
      DOUBLE PRECISION HG ,RHG,HRHG,UG,VG,DEST,RVG,CA1,AM
      DOUBLE PRECISION UIN,VIN,HUIN,HVIN,SIGMAX,DTL,UNORM 
C
      SQ2   = SQRT(2.D0)
      SG    = SQRT(G)
      RA3   = SQRT(1.5D0*G)
      RA32  = RA3**2
      RA33  = RA3*RA32
      ALP   = 0.5D0/RA3
      ALP2  = 0.5D0 *ALP
      ALP3  = ALP/3.D0
C
      FLUENT=0.D0
      FLUSORT=0.D0
C
      DO K=1,NPTFR
       IS=NBOR(K)
       VNX1=XNEBOR(K)
       VNY1=YNEBOR(K)
       VNX=XNEBOR(K+NPTFR)
       VNY=YNEBOR(K+NPTFR)
       VNL=SQRT(VNX**2+VNY**2)
C
       H   = UA(1,IS)
       RH  = SQRT(H)
       U   = UA(2,IS)
       V   = UA(3,IS)
C
C         PAROIS SOLIDES
C         **************
C
C      CONDITION DE GLISSEMENT
C      ***********************
C
       IF(LIMPRO(K,1).EQ.KNEU) THEN
C
         AUX=0.5D0*G*H**2
         FLUH = 0.D0
         FLUU = AUX*VNX
         FLUV = AUX*VNY
       ELSE
C
C        FRONTIERES LIQUIDES
C        *******************
C
C     CALCUL DE F+(H,U,V)
C
       HRH = RH * H
C
       IF(H.LE.0.D0) THEN
         U=0.D0
         V=0.D0
         UNN=0.D0
         VNN=0.D0
         FHPLUS = 0.D0
         FUPLUS = 0.D0
       ELSE
         UNN= +VNX1*U+VNY1*V
         VNN= -VNY1*U+VNX1*V
C
         A=MIN(RA3,MAX(-RA3,-UNN/RH))
         A2 =A * A
         A3 =A2 * A
         ALPHA0=ALP*(RA3-A)
         ALPHA1=ALP2*(RA32-A2)
         ALPHA2=ALP3*(RA33-A3)
C     
         FHPLUS = H*UNN*ALPHA0 + HRH*ALPHA1
         FUPLUS = UNN*(FHPLUS+HRH*ALPHA1) + H*H*ALPHA2
      ENDIF
C
C
C     CALCUL DE L'ETAT FICTIF (HG,UG,VG)
C     ----------------------------------
C
C     H DONNE 
C     ########
C
        IF(LIMPRO(K,1).EQ.KDIR) THEN
C
         C   = SG*RH
         VP1 = UNN
         VP2 = VP1  + C
         VP3 = VP1  - C
C
         HG     =HBOR(K)
         RHG    =SQRT (HG)
         HRHG   =RHG*HG
C
        IF (VP2*VP3.LE.0.D0.OR. VP1.LE.0.D0) THEN
C
         IF(HG.EQ.0.D0) THEN
           UG=0.D0
           VG=0.D0
           FHMOINS = 0.D0
           FUMOINS = 0.D0
           SIGMAX=1.D-2
         ELSE
C
C    REGIME FLUVIAL
C    --------------
C    
          IF (VP2*VP3.LE.0.D0) THEN
C
           UG=UNN+2.D0*SG*(RH-RHG)
           VG=VNN
C
C    REGIME TORRENTIEL 
C    -----------------
C
          ELSE
C
C   FLUX IMPOSE
C   -----------
          IF(LIMPRO(K,2).EQ.KDIR) THEN
C
            UIN = UBOR(K)
            VIN = VBOR(K)
            HUIN = H*UIN
            HVIN = H*VIN
C
            DEST=HUIN*VNX1+HVIN*VNY1
            RVG =-HUIN*VNY1+HVIN*VNX1
C 
            A1 = DEST-FHPLUS
            CA1= SQ2*A1/(SG*HG*RHG)
            CALL ZEROPHI(-1.D0,AM,NIT,CA1)
C
            UG= AM*SG*RHG
            VG=RVG/HG
C
          ELSE  
C
C   IL MANQUE UNE DONNEE,  ON SUPPOSE LE REPOS
C
            UG= 0.D0
            VG= 0.D0
C
          ENDIF
C
           ENDIF
C
         GOTO 220
      ENDIF
      GOTO 200
C
C  LA SORTIE EST EN FAIT TORRENTIELLE, ON NE PEUT MAINTENIR
C   LA CONDITION H DONNE
C
       ELSE
       GOTO 100
C
      ENDIF
C
C
C     VITESSES DONNEES 
C     ################
C
        ELSE IF(LIMPRO(K,2).EQ.KDIR) THEN
C
        UIN = UBOR(K)
        VIN = VBOR(K)
        HUIN = H*UIN
        HVIN = H*VIN
C
         DEST=HUIN*VNX1+HVIN*VNY1
         RVG =-HUIN*VNY1+HVIN*VNX1
C     ATTENTION CHANGEMENT DE SIGNE / RAPPORT INRIA
            A1 = -DEST+FHPLUS
            A2 = -UNN - 2.D0*SG*RH
C            
         IF (A1.LE.0.D0) THEN
C 
C     FH- =-A1 CANNOT BE SATISFIED
C
         FHMOINS = 0.D0
         FUMOINS = 0.D0
         VG=0.D0
        SIGMAX=1.E-2
           ELSE
           CA1= 1.D0/(G*SQ2*A1)**(1.D0/3.D0)
           CALL ZEROPSI(-0.5D0,AM,NIT,CA1,A2)
C
           RHG =A2/(SG*(AM-2.D0))
           HG= RHG * RHG
           HRHG= RHG * HG
C
         IF (HG.EQ.0.D0) THEN
         UG=0.D0
         VG=0.D0
         FHMOINS = 0.D0
         FUMOINS = 0.D0
               SIGMAX=1.D-2
         ELSE
            UG=-AM*A2/(AM-2.D0)
            VG=RVG/HG
        GOTO 220
      ENDIF
      ENDIF
        GOTO 200
C
C   PAS DE CONDITIONS 
C
        ELSE
C
        GOTO 100
C
        ENDIF
       GOTO 1000
C
C
C     CALCUL DE F-(HG,UG,VG)
C
 220   CONTINUE
C
         A=MIN(RA3,MAX(-RA3,-UG/RHG))
         A2 =A * A
         A3 =A2 * A
         ALPHA0=ALP*(A+RA3)
         ALPHA1=ALP2*(A2-RA32)
         ALPHA2=ALP3*(A3+RA33)
C
         FHMOINS = HG*UG*ALPHA0 + HRHG*ALPHA1
         FUMOINS = UG*(FHMOINS + HRHG*ALPHA1) 
     &  + HG*HG*ALPHA2
C
            SIGMAX= RHG
            UNORM=SQRT(UG *UG + VG*VG)
            SIGMAX=MAX( 1.D-2, RA3 *SIGMAX +UNORM )
C
C     CALCUL DES FLUX ET ROTATION INVERSE
C
 200   CONTINUE
         FLUH=(FHPLUS +FHMOINS)*VNL
         FLUU=(FUPLUS +FUMOINS)*VNL
C
         IF (FLUH.GE.0.D0) THEN 
         FLUV= VNN*FLUH 
         ELSE
         FLUV= VG*FLUH 
         ENDIF
C
      FLUTMP=FLUU
      FLUU = +VNX1*FLUTMP-VNY1*FLUV
      FLUV = +VNY1*FLUTMP+VNX1*FLUV
C
C       CORRECTION DU PAS DE PEMPS
C
       DTL = CFL*DTHAUT(IS)/SIGMAX
       DT  = MIN(DT, DTL)
C
       GOTO 1000
100    CONTINUE
       RUN     = H*UNN
C
       FLUH =  RUN* VNL
       FLUU =  (U *RUN + 0.5D0*G*H**2* VNX)*VNL
       FLUV =  (V *RUN + 0.5D0*G*H**2* VNY)*VNL
C
 1000  CONTINUE
       ENDIF
C
       IF(LIMPRO(K,1).EQ.KDIR)  FLUSORT = FLUSORT + FLUH
       IF(LIMPRO(K,2).EQ.KDIR)  FLUENT = FLUENT +FLUH
C
       CE(1,IS)  = CE(1,IS) - FLUH
       CE(2,IS)  = CE(2,IS) - FLUU
       CE(3,IS)  = CE(3,IS) - FLUV
C
       IF(NTRAC.GT.0) THEN
         DO ITRAC=1,NTRAC
           FLUHBTEMP%ADR(ITRAC)%P%R(K)=FLUH
         ENDDO
       ENDIF
C         
       ENDDO
C
C-----------------------------------------------------------------------
C
       RETURN
       END
