C                       ******************                            
                        SUBROUTINE MAJTRAC
C                       ******************                             
C                                                                       
     * (NS,NT,DIMT,DLIMT,NSEG,NPTFR,NUBO,
     * X,Y,AIRS,NU,AIRE,HT,HTN,TN,ZF,NBOR,
     * TBOR,FLUTENT,FLUTSOR,SMTR,NORDRE,CMI,JMI,
     * DJXT,DJYT,DXT,DYT,
     * DPX,DPY,DIFT,CVIST,BETA,DSZ,AIRST,HSTOK,
     * HCSTOK,FLUXT,FLUHBOR,MASSOU,DTT)
C                                                                       
C***********************************************************************
C  TELEMAC 2D VERSION 5.4                                         INRIA
C ----------------------------------------------------------------------
C
C          FONCTION  :  MISE A JOUR DU TRACEUR
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C |  NS            | -->|  NOMBRE DE NOEUDS DU MAILLAGE                |
C |  NT            | -->|  NOMBRE D'ELEMENTS                           |
C !  DIMT          ! -->!  DIMENSION DU TRACEUR                        !
C !  DLIMT         ! -->!  DIMENSION DU TRACEUR AU BORD                !
C |  NSEG          | -->|  NOMBRE D'ARETES                             |
C |  NPTFR         | -->|  NOMBRE DE POINTS FRONTIERE                  |
C !  NUBO          ! -->!  NUMEROS GLOBAUX DES EXTREMITES DES ARETES   !
C |  X,Y           | -->|  COORDONNEES DES NOEUDS DU MAILLAGE          |
C !  AIRS          ! -->!  AIRES DES CELLULES                          !
C |  NU            | -->|  NUMEROS DES NOEUDS PAR TRIANGLE             |
C !  HT            !<-- !  HT AU TEMPS N+1                             !
C !  HTN,TN        ! -->!  HT, T  AU TEMPS N                           !
C |  ZF            | -->|  COTES DU FOND                               |
C |  NBOR          | -->|  NUMEROS GLOBAUX DES POINTS DE BORD          |
C |  TBOR          | -->|  CONDITIONS AUX LIMITES SUR T                !
C ! FLUTENT,FLUTSOR|<-- |  FLUX TRACEUR ENTREE ET SORTIE               |
C |  SMTR          | -->!  TERMES SOURCES DU TRACEUR                   !
C |  NORDRE        | -->|  ORDRE DU SCHEMA                             |
C !  CMI           ! -->!  COORDONNEES DES POINTS MILIEUX D'INTERFACE  !
C !  JMI           ! -->!  NUMERO DU TRIANGLE AUQUEL APPARTIENT LE     !
C |                !    !  POINT MILIEU DE L'INTERFACE                 !
C |  DJXT,DJYT     | -- |  GRADIENTS PAR TRIANGLES                     |
C |  DXT,DYT       | -- |  GRADIENTS PAR NOEUDS                        |
C |  DPX, DPY      ! -->!  GRADIENTS DES FONCTIONS DE BASE             !
C |  DIFT          | -->|  LOGIQUE INDIQUANT S'IL Y A DIFFUSION TRACEUR|
C |  CVIST         | -->|  COEFFICIENT DE DIFFUSION DU TRACEUR         |
C |  BETA          ! -- !  COEFFICIENT EXTRAPOLATION POUR ORDRE 2      !
C |  DSZ           ! -->!  VARIATION DE Z POUR ORDRE 2                 !
C |  AIRST         ! -->!  AIRES DES SOUS-TRIANGLES DANS CELLULES      !
C |  HSTOK         | -->|  HAUTEURS D'EAU  STOCKEES                    |
C !  HCSTOK        ! -->!  H RECONSTRUIT ORDRE 2   CORRIGE  STOCKE     !
C |  FLUXT         | -->|  FLUX DE MASSE                               |
C |  FLUHBOR       | -->|  FLUX DE MASSE FRONTIERE                     |
C |  MASSOU        |<-- |  MASSE DE TRACEUR AJOUTEE PAR TERME SOURCE   |
C !  DTT           ! -->!  PAS DE TEMPS TRACEUR                        !
C !________________!____!______________________________________________!
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)   
C        -- (TABLEAU DE TRAVAIL)                                        
C-----------------------------------------------------------------------
C     - SOUS PROGRAMME(S) APPELANT : RESOLU                             
C     - SOUS PROGRAMME(S) APPELES  : GRADNODT 
C 
C***********************************************************************
C
      USE INTERFACE_TELEMAC2D, EX_MAJTRAC => MAJTRAC
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C     
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      LOGICAL, INTENT(IN) :: DIFT
      INTEGER, INTENT(IN) :: NSEG,NPTFR,NORDRE,DIMT,DLIMT,NS,NT
      INTEGER, INTENT(IN) :: NUBO(2,NSEG),NU(NT,3)
      INTEGER, INTENT(IN) :: NBOR(NPTFR),JMI(*)
      DOUBLE PRECISION, INTENT(INOUT) :: HT(DIMT),FLUTENT,FLUTSOR
      DOUBLE PRECISION, INTENT(INOUT) :: MASSOU
      DOUBLE PRECISION, INTENT(IN)    :: TBOR(DLIMT),DSZ(2,*)
      DOUBLE PRECISION, INTENT(IN)    :: X(NS),Y(NS),AIRS(NS),AIRE(NT)
      DOUBLE PRECISION, INTENT(IN)    :: HTN(DIMT),TN(DIMT),ZF(*)
      DOUBLE PRECISION, INTENT(IN)    :: SMTR(DIMT),DPX(3,NT),DPY(3,NT)
      DOUBLE PRECISION, INTENT(IN)    :: CMI(2,*),AIRST(2,*),CVIST
      DOUBLE PRECISION, INTENT(INOUT) :: DJXT(*),DJYT(*),DXT(*),DYT(*)
      DOUBLE PRECISION, INTENT(INOUT) :: BETA
      DOUBLE PRECISION, INTENT(IN)    :: HSTOK(*)
      DOUBLE PRECISION, INTENT(IN)    :: HCSTOK(2,*),FLUXT(*)
      DOUBLE PRECISION, INTENT(IN)    :: FLUHBOR(*),DTT
C       
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C         
      INTEGER IS,K,NSG,NUBO1,NUBO2,J,ILIM,ERR      
C                  
      DOUBLE PRECISION ZF1,ZF2,FLUH,FLUT,HI0,HJ0,AIX,AIY,AJX,AJY,AMDS
      DOUBLE PRECISION GRADI,GRADJ,GRADIJ,GRADJI,FLU11,FLU41,UAS41,UAS42
C
C     EX TABLEAUX AUTOMATIQUES !!!!!!!!
C
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: CET(:),DST(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: DSP(:),DSM(:),CORRT(:)
      LOGICAL DEJA
      DATA DEJA/.FALSE./
C                                                                       
C-----------------------------------------------------------------------
C
      IF(.NOT.DEJA) THEN
        ALLOCATE(CET(NS),STAT=ERR)
        IF(ERR.NE.0) GO TO 1001
        ALLOCATE(DST(2,NSEG),STAT=ERR)
        IF(ERR.NE.0) GO TO 1001
        ALLOCATE(DSP(NS),STAT=ERR)
        IF(ERR.NE.0) GO TO 1001
        ALLOCATE(DSM(NS),STAT=ERR)
        IF(ERR.NE.0) GO TO 1001
        ALLOCATE(CORRT(NS),STAT=ERR)
        IF(ERR.NE.0) GO TO 1001
        GO TO 1002
1001    CONTINUE
        IF(LNG.EQ.1) WRITE(LU,1000) ERR
        IF(LNG.EQ.2) WRITE(LU,2000) ERR
1000    FORMAT(1X,'MAJTRAC : ERREUR A L''ALLOCATION DE MEMOIRE : ',/,1X,
     *        'CODE D''ERREUR : ',1I6)
2000    FORMAT(1X,'MAJTRAC: ERROR DURING ALLOCATION OF MEMORY: ',/,1X,
     *        'ERROR CODE: ',1I6)
        CALL PLANTE(1)
        STOP
1002    CONTINUE
        DEJA=.TRUE.
      ENDIF
C                                                                       
C-----------------------------------------------------------------------
C                                                                               
C  INITIALISATION
C	
      DO IS=1,NS
        CET(IS) = 0.D0
      ENDDO
C
C   CALCUL DES GRADIENTS PAR TRIANGLE ET PAR NOEUD DU TRACEUR 
C                 ET DU TERME DE DIFFUSION
C
      IF(DIFT.OR.NORDRE.EQ.2) CALL GRADNODT(NS,NT,NU,AIRE,AIRS,
     *HSTOK,TN,DPX,DPY,DJXT,DJYT,DXT,DYT,DIFT,CVIST,CET,DTT)
C
      IF(NORDRE.EQ.2) THEN
C
C  RECONSTRUCTION 2ND ORDRE POUR TRACEUR
C  *************************************
C
      DO IS=1,NS
        DSP(IS)=0.D0
        DSM(IS)=0.D0
      ENDDO 
C      
      DO  NSG=1,NSEG 
C
         J         = JMI(NSG)
C
         NUBO1     = NUBO(1,NSG)
         NUBO2     = NUBO(2,NSG)
C
         ZF1   =    ZF(NUBO1)
         ZF2   =    ZF(NUBO2)
C
         HI0   =HSTOK(NUBO1)
         HJ0   =HSTOK(NUBO2)
C
C   POUR UNE ARETE EN RECOUVREMENT, ON RESTERA A L'ORDRE 1
C
         IF(ZF1.GE. (HJ0+ZF2) .OR. ZF2.GE. (HI0+ZF1)       
     &  .OR. 2.*ABS(DSZ(1,NSG)).GE.HI0 
     &  .OR. 2.*ABS(DSZ(1,NSG)).GE.HJ0
     &  .OR. 2.*ABS(DSZ(2,NSG)).GE.HI0
     &  .OR. 2.*ABS(DSZ(2,NSG)).GE.HJ0)  THEN
         DST(1,NSG) =0.D0
         DST(2,NSG) =0.D0
        ELSE
C
         AIX       = CMI(1,NSG)-X(NUBO1)
         AIY       = CMI(2,NSG)-Y(NUBO1)
         AJX       = CMI(1,NSG)-X(NUBO2)
         AJY       = CMI(2,NSG)-Y(NUBO2) 
C
         GRADI  = AIX*DXT(NUBO1) + AIY*DYT(NUBO1) 
         GRADJ  = AJX*DXT(NUBO2) + AJY*DYT(NUBO2)
         GRADIJ  = AIX*DJXT(J) + AIY*DJYT(J) 
         GRADJI  = AJX*DJXT(J) + AJY*DJYT(J)
C
C    EXTRAPOLATION DES GRADIENTS ET LIMITEUR DE PENTE 
C
         ILIM=2
         BETA=0.3333D0
         DST(1,NSG)  = EXLIM (ILIM,BETA,GRADI,GRADIJ)
         DST(2,NSG)  = EXLIM (ILIM,BETA,GRADJ,GRADJI)
C
       ENDIF
       ENDDO
C
      DO NSG=1,NSEG 
C
         NUBO1     = NUBO(1,NSG)
         NUBO2     = NUBO(2,NSG)
C
         IF(DST(1,NSG).GE.0.D0) THEN
         DSP(NUBO1) = DSP(NUBO1) + 
     *  AIRST(1,NSG)* HCSTOK(1,NSG)*DST(1,NSG)
         ELSE
         DSM(NUBO1) = DSM(NUBO1) - 
     *  AIRST(1,NSG)* HCSTOK(1,NSG)*DST(1,NSG)
         ENDIF
         IF(DST(2,NSG).GE.0.) THEN
         DSP(NUBO2) = DSP(NUBO2) + 
     *  AIRST(2,NSG)* HCSTOK(2,NSG)*DST(2,NSG)
         ELSE
         DSM(NUBO2) = DSM(NUBO2) - 
     *  AIRST(2,NSG)* HCSTOK(2,NSG)*DST(2,NSG)
         ENDIF
C
      ENDDO 
C
C  ON CALCULE LES CORRECTIONS POUR ASSURER LA CONSERVATION DE HT
C                 ***********                 ******************
C
      DO IS=1,NS
       CORRT(IS) =  DSM(IS) - DSP(IS)
       AMDS =MAX(DSP(IS),DSM(IS))
        IF(AMDS.GT.0.D0) THEN
        CORRT(IS) = CORRT(IS)/AMDS
        ENDIF
      ENDDO 
 12       CONTINUE
C
      ENDIF
C
C  CALCUL DES FLUX POUR LES INTERFACES INTERNES
C
      DO 500 NSG=1,NSEG 
C
         NUBO1     = NUBO(1,NSG)
         NUBO2     = NUBO(2,NSG)
C
         UAS41     = TN(NUBO1)
         UAS42     = TN(NUBO2)
C
         FLU11=FLUXT(NSG)
C
         IF (FLU11.GE.0.) THEN 
       IF(NORDRE.GE.2) THEN
          UAS41 = UAS41  + DST(1,NSG) +
     & MIN(0.D0,CORRT(NUBO1))*MAX(0.D0,DST(1,NSG))+
     & MAX(0.D0,CORRT(NUBO1))*MAX(0.D0,-DST(1,NSG))
        ENDIF
         FLU41 =  UAS41 * FLU11
         ELSE
       IF(NORDRE.GE.2) THEN
           UAS42 = UAS42 + DST(2,NSG) +
     & MIN(0.D0,CORRT(NUBO2))*MAX(0.D0,DST(2,NSG))+
     & MAX(0.D0,CORRT(NUBO2))*MAX(0.D0,-DST(2,NSG))
       ENDIF
         FLU41 =  UAS42 * FLU11
        ENDIF
C
         CET(NUBO1) = CET(NUBO1) - FLU41
C
         CET(NUBO2) = CET(NUBO2) + FLU41  
C
500   CONTINUE
C
C    FLUX FRONTIERE
C
      DO K=1,NPTFR
       IS=NBOR(K)
C
       FLUH =FLUHBOR(K)
C
       IF(FLUH.GE.0.D0) THEN
         FLUT= TN(IS)*FLUH  
         FLUTSOR = FLUTSOR +FLUT
       ELSE
         FLUT= TBOR(K)*FLUH  
         FLUTENT = FLUTENT +FLUT
       ENDIF
C
       CET(IS)  = CET(IS) - FLUT 
C               
      ENDDO
C
C   MISE A JOUR DE HT
C
      DO  IS =1,NS
C
        HT(IS)  = HTN(IS) +  (CET(IS)+SMTR(IS))/AIRS(IS)
        MASSOU = MASSOU + SMTR(IS)
C
        IF(HT(IS).LE.1.D-15) HT(IS)=0.D0
C
      ENDDO
C
C-----------------------------------------------------------------------
C                                                                       
      RETURN
      END
