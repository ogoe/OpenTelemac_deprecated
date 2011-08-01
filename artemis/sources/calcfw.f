C                       *****************                               
                        SUBROUTINE CALCFW
C                       *****************                               
     *(I,H,C,CG,K,HMU,            
     * NPOIN,OMEGA,GRAV,         
     * VISCO,DIAM90,DIAM50,MVSED,MVEAU,
     * FORMFR,REGIDO,RICOEF,
     * ENTREG,ENTRUG,FFW)
C
C***********************************************************************
C
C  ARTEMIS VERSION 5.1 02/06/99   D. AELBRECHT (LNH) 01 30 87 74 12 
C                                    D. PAUGAM (STAGE 1996)
C  LINKED TO BIEF VERS. 5.0          J-M HERVOUET (LNH) 01 30 87 80 18
C       
C***********************************************************************
C
C     FONCTION  : CALCULE LE COEFFICIENT DE FROTTEMENT DE FOND
C                 FW POUR DES FONDS SABLEUX 
C-------------------------------------------------------------------
C APPELE PAR:  BERKHO
C-------------------------------------------------------------------
C
      USE BIEF
      USE INTERFACE_ARTEMIS, EX_CALCFW => CALCFW
C
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                               
C                                                                       
      INTEGER NPOIN,I               
      INTEGER REGIME, REGIDO
      INTEGER FORMFR
C                                                                       
      DOUBLE PRECISION C(NPOIN),CG(NPOIN)       
      DOUBLE PRECISION K(NPOIN)
      DOUBLE PRECISION HMU(NPOIN) 
C                                                                       
      DOUBLE PRECISION H(NPOIN) 
      DOUBLE PRECISION GRAV,OMEGA      
C
      DOUBLE PRECISION VISCO, DIAM90, DIAM50
      DOUBLE PRECISION MVSED, MVEAU, RICOEF
      DOUBLE PRECISION KSRUGO, AEX, UEX, FFW
      DOUBLE PRECISION RAP1, RAP2
      DOUBLE PRECISION LL1 
      DOUBLE PRECISION RMAX1, RMIN1, RMAX2, RMIN2
      DOUBLE PRECISION SP, PSI, KS1, RH, RS, KS2 
      DOUBLE PRECISION DIAMAD, TETACR, TAUCR
C
      LOGICAL ENTREG,ENTRUG 
C                                                                       
C      
      INTRINSIC ABS,EXP                                                 
C                                                                       
C                                                                       
C  CALCUL DES VALEURS DE FW SUR LE DOMAINE
C--------------------------------------------------------
      RMAX1 = 0.D0
      RMIN1 = 1.D7
      RMAX2 = 0.D0
      RMIN2 = 1.D7
C
      KSRUGO = 3.D0 * DIAM90
C
C ON CALCULE AEX ET UEX EN CHAQUE POINT DU MAILLAGE
C AINSI LR REGIME HYDRAULIQUE SERA DETERMINE EN CHAQUE
C POINT DU MAILLAGE EGALEMENT
C
C
        AEX = HMU(I)/(2.D0*SINH(K(I)*H(I)))
        UEX = OMEGA * AEX
        RAP1 = (UEX*AEX)/VISCO
        RAP2 = AEX/KSRUGO
C
C   TESTS SUR LES VALEURS MAXIMALES DU NOMBRE DE REYNOLDS
C   ET DU RAPPORT EXCURSION SUR RUGOSITE
C
        IF (RAP1 .GT. RMAX1) THEN
            RMAX1 = RAP1
          ENDIF
C
        IF (RAP1 .LT. RMIN1) THEN
             RMIN1 = RAP1
          ENDIF
C
        IF (RAP2 .GT. RMAX2) THEN
              RMAX2 = RAP2
           ENDIF
C
        IF (RAP2 .LT. RMIN2) THEN
              RMIN2 = RAP2
           ENDIF
C
C-----------------------------------------------------------------
C DETERMINATION DU REGIME HYDRAULIQUE
C-----------------------------------------------------------------
C
      IF (ENTREG) THEN
         REGIME = REGIDO
      ELSE
C
C     ON INITIALISE LE REGIME A ZERO AVANT CHAQUE NOUVELLE ITERATION
C
         REGIME = 0
         LL1 = 0.D0
C
         IF (RAP1 .LE. 10250.D0) THEN
            LL1 = 0.0322D0*RAP1+3.33D0
            IF (RAP2 .GE. LL1) THEN
              REGIME = 1
C             WRITE(*,*) 'LE REGIME HYDRAULIQUE EST LAMINAIRE'
            ENDIF
         ENDIF
C
         IF (RAP1 .GE. 3.D4) THEN
            LL1 = 0.009792D0*RAP1+208.33D0
            IF (RAP2 .GE. LL1) THEN
              REGIME = 2
C             WRITE(*,*) 'LE REGIME HYDRAULIQUE EST TURBULENT LISSE'
            ENDIF
         ENDIF
C
         IF (RAP1 .GE. 5.D3 .AND. RAP1 .LE. 2.D4) THEN
            LL1 = 0.026D0*RAP1-12.D0
            IF (RAP2 .LE. LL1) THEN
              REGIME = 3
C             WRITE(*,*) 'LE REGIME HYDRAULIQUE EST TURBULENT RUGUEUX'
            ENDIF
         ENDIF
C
         IF (RAP1 .GT. 2.D4) THEN
            LL1 = 0.00099D0*RAP1+30.30D0
            IF (RAP2 .LE. LL1) THEN
              REGIME = 3
C             WRITE(*,*) 'LE REGIME HYDRAULIQUE EST TURBULENT RUGUEUX'
            ENDIF
         ENDIF
C
         IF (REGIME .EQ. 0) THEN
            REGIME = 3
C           WRITE(*,*) 'LE REGIME EST TURBULENT RUGEUX'
         ENDIF
C
      ENDIF
C
C
C--------------------------------------------------------------
C FORMULATION ET CALCUL DE FW : LE COEFFICIENT DE FROTTEMENT
C--------------------------------------------------------------
C
      IF (REGIME .EQ. 1) THEN
         FFW = 2.D0*(((AEX*UEX)/VISCO)**(-0.5D0))
      ENDIF
C
      IF (REGIME .EQ. 2) THEN
         FFW = 0.09*(((AEX*UEX)/VISCO)**(-0.2D0))
      ENDIF
C
      IF (REGIME .EQ. 3) THEN
C
         IF (ENTRUG) THEN
C
             KSRUGO = 3.D0 * DIAM90
C
         ELSE
C
            SP = MVSED/MVEAU
            PSI = (UEX**2.D0)/((SP-1.D0)*GRAV*DIAM50)
C
            IF (PSI .LT. 250.D0) THEN
                 KS1 = 3.D0*DIAM90
C
                 IF (PSI .LT. 10.D0) THEN
                      RH = 0.22D0*AEX
                      RS = 0.18D0
                 ELSE
                      RH = AEX*(2.4D-13)*((250.D0-PSI)**5.D0)
                      RS = (2.D-7)*((250.D0-PSI)**2.5D0)
                 ENDIF
C
                 KS2 = 20.D0*RICOEF*RH*RS
C
C  RICOEF = 1 SI IL Y A SEULEMENT DES RIDES
C  RICOEF = 0.7D0 SI LES RIDES SONT SUPERPOSES A DES VAGUES DE SABLE
C
            ELSE
                 KS1 = 3.D0*(0.04D0*PSI-9.D0)*DIAM90
                 KS2 = 0.D0
                 IF (KS1 .LT. 0.01D0) THEN
                      KS1 = 3.D0*DIAM90
                 ENDIF
            ENDIF
C
            KSRUGO = KS1+KS2
            KS1 = (AEX/KSRUGO)
C
         ENDIF
C
         FFW = EXP(-6.D0+5.2D0*((AEX/KSRUGO)**(-0.19D0)))
         IF (FFW .GT. 0.3D0) THEN
             FFW = 0.3D0
         ENDIF
      ENDIF
C
      IF (REGIME .EQ. 4) THEN
C
C     REGIME TRANSITOIRE NON PRIS EN COMPTE
C
         IF (ENTRUG) THEN
C
             KSRUGO = 3.D0 * DIAM90
C
         ELSE
C
            DIAMAD = (((MVSED/MVEAU)-1.D0)*GRAV)/(VISCO**2.D0)
            DIAMAD = DIAMAD**(1/3)
            DIAMAD = DIAMAD * DIAM50
C
            TETACR = 0.14D0*(DIAMAD**(-0.64D0))
C
            TAUCR = TETACR*(MVSED-MVEAU)*GRAV*DIAM50
C
            KSRUGO = (3.D0*DIAM90)+((3.3D0*VISCO)/
     *               ((TAUCR/MVEAU)**0.5D0))
C
         ENDIF
C
         FFW = EXP(-6.D0+5.2D0*((AEX/KSRUGO)**(-0.19D0)))
         IF (FFW .GT. 0.3D0) THEN
              FFW = 0.3D0
         ENDIF
C
      ENDIF      
C
      RETURN
      END
