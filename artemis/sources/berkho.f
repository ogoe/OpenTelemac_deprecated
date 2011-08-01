C                       *****************
                        SUBROUTINE BERKHO
C                       *****************
C
     *(LT)
C
C***********************************************************************
C
C  ARTEMIS VERSION 6.0   18/03/10    C. DENIS (SINETICS)
C
C  ARTEMIS VERSION 5.1   04/06/99    D. AELBRECHT (LNH) 01 30 87 74 12
C
C  LINKED TO BIEF VERS. 5.0         J-M HERVOUET (LNHE) 01 30 87 80 18
C
C  LE 02/04/2007 : INVERSION DE LA SECONDE EQUATION AVANT APPEL A SOLVE
C                  EN CAS DE SOLVEUR DIRECT
C
C
C***********************************************************************
C
C      FONCTION: RESOUT L'EQUATION DE BERKHOFF MODIFIEE PAR 
C      =========   L'INTRODUCTION DE TERMES DE DISSIPATION
C
C      DIV (C*CG*GRAD(PHI)) + C*CG*( K**2 + I*K*MU ) * PHI = 0
C                                           ------
C
C PHI EST UNE FONCTION COMPLEXE DE COMPOSANTES RELLE PHIR ET IMAGINAIRE
C PHII
C
C MU EST UN COEFFICIENT DE DISSIPATION A PRIORI INCONNU
C
C LES CONDITIONS AUX LIMITES COUPLENT LES EQUATIONS EN PHIR ET PHII
C ELLES S'ECRIVENT :
C
C D(PHI)/DN  - I*K*PHI = D(F)/DN - I*K*F  (N: NORMALE EXTERIEURE)
C POUR UNE FRONTIERE MARITIME OU L'ON IMPOSE UNE ONDE INCIDENTE DEFINIE
C PAR LE POTENTIEL F (F=0 POUR UNE SORTIE LIBRE)
C
C D(PHI)/DN - I*(1-R*EXP(I*ALFA))/(1+R*EXP(I*ALFA))*K*COS(TETA)*PHI = 0
C POUR UNE PAROI SOLIDE, AVEC R COEFFICIENT DE REFLEXION DE LA PAROI,
C TETA ANGLE D'INCIDENCE DE LA HOULE SUR LA PAROI ET ALFA DEPHASAGE
C INDUIT PAR LA PAROI.
C
C SOIT D'UNE FACON GENERALE :
C D(PHIR)/DN = APHIRB*PHII + BPHIRB*PHIR + CPHIRB
C D(PHII)/DN =-APHIRB*PHIR + BPHIRB*PHII + DPHIRB
C
C
C APRES FORMULATION VARIATIONNELLE, ON OBTIENT LE SYSTEME SUIVANT :
C
C         (  AM1          BM1     )  ( PHIR )   ( CV1 )
C         (                       )  (      ) = (     )
C         (                       )  (      )   (     )
C         (  -BM1         AM1     )  ( PHII )   ( CV2 )
C
C           /
C AM1 =    / C*CG * GRAD(PSII)*GRAD(PSIJ) DS
C         /S
C
C           /
C       -  / OMEGA**2 * CG/C * PSII*PSIJ  DS
C         /S
C
C           /
C       -  /  BPHIRB * PSII*PSIJ  DB
C         /B
C
C           /                         /
C BM1 =  - /  APHIR * PSII*PSIJ DB + /  C*CG* K * MU * PSII * PSIJ DS
C         /B                        /S
C
C          /
C CV1 =   /   CPHIR * PSII DB
C        /B
C
C          /
C CV2 =   /   CPHII * PSII DB
C        /B
C
C
C AVEC S DESIGNANT LE DOMAINE DE CALCUL ET B SA FRONTIERE
C      PSII ET PSIJ LES FONCTIONS DE BASE AUX POINTS I ET J DU MAILLAGE
C
C REMARQUE : COMME APHII=-APHIR, BM1 APPARAIT AUSSI DANS L'EQUATION EN
C --------   PHII.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     LT         | -->| INDICE DU CALCUL COURANT                     |
C |   IELM         | -->|  TYPE D'ELEMENT                              |
C |   IELMB        | -->|  TYPE D'ELEMENT DE BORD                      |
C |   OMEGA        | -->|  PULSATION DE LA HOULE                       |
C |   PER          | -->|  PERIODE DE LA HOULE                         |
C |--------------------------------------------------------------------|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C APPELE PAR              :  ARTEMIS
C SOUS-PROGRAMMES APPELES :  OS, OSDB, OM, MATRIX, VECTOR
C                            DIRICH, CNTPRE, LUMP, SOLVE
C                            CALCQB, FWSPEC, CALCFW, CALCUE
C**********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_ARTEMIS
      USE INTERFACE_ARTEMIS, EX_BERKHO => BERKHO
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C
      INTEGER I,LT
      INTEGER ITERMU 
      DOUBLE PRECISION MAX_ECRHMU
      DOUBLE PRECISION MAX_MODHMU 
      
      DOUBLE PRECISION HM,HMUE,HEFF,ECRHMU,MODHMU
      DOUBLE PRECISION Q1,Q2,Q3
C
      DOUBLE PRECISION CBID,FFW
      DOUBLE PRECISION PI,DEGRAD,RADDEG
C
C-----------------------------------------------------------------------
C
      PARAMETER( PI = 3.1415926535897932384626433D0 , DEGRAD=PI/180.D0 )
      PARAMETER( RADDEG = 180.D0 / PI )
C
C-----------------------------------------------------------------------
C
      INTRINSIC ABS,MIN,MAX,LOG
      DOUBLE PRECISION P_DMAX
      EXTERNAL P_DMAX
      DOUBLE PRECISION TDEB1,TFIN1
C
C----------------------------------------------------------------------
C
C INITIALISATION A 0 DES VARIABLES MU ET FW  
C         POUR LE PREMIER CALCUL
C
      ITERMU=0
      IF (LT.EQ.0) THEN
         CALL OS( 'X=C     ' , MU , SBID , SBID , 0.D0 )
         CALL OS( 'X=C     ' , FW , SBID , SBID , 0.D0 )
      ENDIF
C     
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C
C BOUCLE D'ITERATION SUR LA VARIABLE MU DE DISSIPATION
C
C
C     =========================================
C
C     CALCUL DES MATRICES & DES SECONDS MEMBRES
C
C     =========================================
C
C     ---------------------------
C     MATRICE DE DIFFUSION DE AM1
C     ---------------------------
C
98      CALL OS( 'X=YZ    ' , T1 , C , CG , CBID )
      CALL MATRIX(AM1,'M=N     ','MATDIF          ',IELM,IELM,          
     *            1.D0,S,S,S,T1,T1,S,MESH,MSK,MASKEL)
C
C-----------------------------------------------------------------------
C
C PANCHANG A REVOIR : 7 EST GMRES
C
C ON STOCKE LA MATRICE DE DIFFUSION QUI SERT POUR LE PRECONDITIONNEMENT
C SI LA METHODE EST CELLE DE PANCHANG ET AL (ISOLVE(1)=7)
C
C     IF (ISOLVE(1).EQ.7) THEN
C
C        CALL OM('M=CN    ',AM3,AM1,Z,1.D0/(RELAX*(2.D0-RELAX)),MESH)
C
C     ENDIF
C
C-----------------------------------------------------------------------
C
C     -----------------------
C     MATRICE DE MASSE DE AM1
C     -----------------------
C
      

      CALL OS( 'X=Y/Z   ' , T1 , CG , C , CBID )
      CALL MATRIX(AM2,'M=N     ','FMATMA          ', IELM , IELM ,
     *            OMEGA**2 , T1,S,S,S,S,S,MESH,MSK,MASKEL)
C
C     --------------------------------------------------
C     ON CALCULE MATRICE DE DIFFUSION - MATRICE DE MASSE
C     --------------------------------------------------
C
      CALL OM( 'M=M+CN  ' , AM1 , AM2 , C , -1.D0 , MESH )
C
C     --------------------------------
C     ON AJOUTE A AM1 LE TERME DE BORD 
C     --------------------------------
C
C     (T1 SERT ICI DE STRUCTURE BIDON)
C
C        ------------------------------
C        TERME DE BORD : ONDE INCIDENTE 
C     ------------------------------
C
      IF (NPTFR .GT. 0) THEN
         
         CALL MATRIX(MBOR,'M=N     ','FMATMA          ',IELMB,IELMB,       
     *        -1.D0,BPHI1B,S,S,S,S,S,MESH,.TRUE.,MASK1)
         CALL OM( 'M=M+N   ' , AM1 , MBOR , T1 , CBID , MESH )
      END IF
C     ------------------------------
C     TERME DE BORD : SORTIE LIBRE 
C     ------------------------------
      IF (NPTFR .GT. 0) THEN
         CALL MATRIX(MBOR,'M=N     ','FMATMA          ',IELMB,IELMB,       
     *        -1.D0,BPHI2B,S,S,S,S,S,MESH,.TRUE.,MASK2)
         CALL OM( 'M=M+N   ' , AM1 , MBOR , T1 , CBID , MESH )
      END IF
C     ------------------------------
C     TERME DE BORD : PAROI SOLIDE 
C     ------------------------------
      IF (NPTFR .GT. 0) THEN
         CALL MATRIX(MBOR,'M=N     ','FMATMA          ',IELMB,IELMB,       
     *        -1.D0,BPHI3B,S,S,S,S,S,MESH,.TRUE.,MASK3)
         
         CALL OM( 'M=M+N   ' , AM1 , MBOR , T1 , CBID , MESH )             
      END IF
C     ------------------------------
C     TERME DE BORD : ONDE IMPOSEE
C     ------------------------------
      IF (NPTFR .GT. 0) THEN
         CALL MATRIX(MBOR,'M=N     ','FMATMA          ',IELMB,IELMB,       
     *        -1.D0,BPHI4B,S,S,S,S,S,MESH,.TRUE.,MASK4)
         CALL OM( 'M=M+N   ' , AM1 , MBOR , T1 , CBID , MESH )              
      END IF
C     
C     ---------------------
C     SECONDS MEMBRES : CV1
C     ---------------------
C
         CALL OS( 'X=C     ' , CV1, SBID , SBID , 0.D0 )
C     ------------------------------
C     TERME DE BORD : ONDE INCIDENTE 
C     ------------------------------
         CALL OS( 'X=CY    ' , T1,TETAP,SBID,DEGRAD)
         CALL OS( 'X=COS(Y)' , T2,T1,SBID,0.D0)      
         CALL OS( 'X=YZ    ' , T3,CPHI1B,T2,0.D0)   
         
         CALL OS( 'X=C     ' , T1, SBID , SBID , 0.D0 )
         IF (NPTFR .GT. 0) THEN
            CALL VECTOR(T1,'=','MASVEC          ',IELMB,                      
     *           -1.D0,T3,S,S,S,S,S,MESH,.TRUE.,MASK1) 
         END IF
         CALL OSDB( 'X=X+Y   ' , CV1 , T1 , SBID , CBID , MESH )                 
         CALL OS( 'X=CY    ' , T1,TETAB,SBID,DEGRAD)      
         CALL OS( 'X=COS(Y)' , T2,T1,SBID,0.D0)      
         CALL OS( 'X=SIN(Y)' , T3,T1,SBID,0.D0)      
         CALL OS( 'X=C     ' , T1, SBID , SBID , 0.D0 )
         IF (NPTFR .GT. 0) THEN
            CALL VECTOR(T1,'=','FLUBDF          ',IELMB,                      
     *           1.D0,CPHI1B,S,S,T2,T3,S,MESH,.TRUE.,MASK1) 
         END IF
         CALL OSDB( 'X=X+Y   ' , CV1 , T1 , SBID , CBID , MESH )                 
C     ------------------------------
C     TERME DE BORD : SORTIE LIBRE 
C     ------------------------------
         IF (NPTFR .GT. 0) THEN
            CALL VECTOR(T1,'=','MASVEC          ',IELMB,                      
     *           1.D0,CPHI2B,S,S,S,S,S,MESH,.TRUE.,MASK2) 
C     END IF
            CALL OSDB( 'X=X+Y   ' , CV1 , T1 , SBID , CBID , MESH )                
         END IF
C     ------------------------------
C     TERME DE BORD : PAROI SOLIDE 
C     ------------------------------
         IF (NPTFR .GT. 0) THEN
            CALL VECTOR(T1,'=','MASVEC          ',IELMB,                      
     *           1.D0,CPHI3B,S,S,S,S,S,MESH,.TRUE.,MASK3) 
            
            CALL OSDB( 'X=X+Y   ' , CV1 , T1 , SBID , CBID , MESH )
         END IF
C     ------------------------------
C     TERME DE BORD : ONDE IMPOSEE
C     ------------------------------
         IF (NPTFR .GT. 0) THEN
            CALL VECTOR(T1,'=','MASVEC          ',IELMB,                      
     *           1.D0,CPHI4B,S,S,S,S,S,MESH,.TRUE.,MASK4) 
        
            CALL OSDB( 'X=X+Y   ' , CV1 , T1 , SBID , CBID , MESH )                 
         END IF
         
      
C     ---------------------
C     SECONDS MEMBRES : CV2
C     ---------------------
C
         CALL OS( 'X=C     ' , CV2, SBID , SBID , 0.D0 )
C     ------------------------------
C     TERME DE BORD : ONDE INCIDENTE 
C     ------------------------------
         CALL OS( 'X=CY    ' , T1,TETAP,SBID,DEGRAD)      
         CALL OS( 'X=COS(Y)' , T2,T1,SBID,0.D0)      
         CALL OS( 'X=YZ    ' , T3,DPHI1B,T2,0.D0)      
         CALL OS( 'X=C     ' , T1, SBID , SBID , 0.D0 )
         IF (NPTFR .GT. 0) THEN
            CALL VECTOR(T1,'=','MASVEC          ',IELMB,                      
     *           -1.D0,T3,S,S,S,S,S,MESH,.TRUE.,MASK1) 
         END IF
         CALL OSDB( 'X=X+Y   ' , CV2 , T1 , SBID , CBID , MESH )                 
         CALL OS( 'X=CY    ' , T1,TETAB,SBID,DEGRAD)      
         CALL OS( 'X=COS(Y)' , T2,T1,SBID,0.D0)      
         CALL OS( 'X=SIN(Y)' , T3,T1,SBID,0.D0)      
         CALL OS( 'X=C     ' , T1, SBID , SBID , 0.D0 )
         IF (NPTFR .GT. 0) THEN
            CALL VECTOR(T1,'=','FLUBDF          ',IELMB,                      
     *           1.D0,DPHI1B,S,S,T2,T3,S,MESH,.TRUE.,MASK1) 
         END IF
         CALL OSDB( 'X=X+Y   ' , CV2 , T1 , SBID , CBID , MESH )             
         
         
    
C        ------------------------------
C        TERME DE BORD : SORTIE LIBRE 
C        ------------------------------
         IF (NPTFR .GT. 0) THEN
            CALL VECTOR(T1,'=','MASVEC          ',IELMB,                      
     *           1.D0,DPHI2B,S,S,S,S,S,MESH,.TRUE.,MASK2) 
         END IF
         CALL OSDB( 'X=X+Y   ' , CV2 , T1 , SBID , CBID , MESH )                 
C     ------------------------------
C     TERME DE BORD : PAROI SOLIDE 
C     ------------------------------
         IF (NPTFR .GT. 0) THEN
            CALL VECTOR(T1,'=','MASVEC          ',IELMB,                      
     *           1.D0,DPHI3B,S,S,S,S,S,MESH,.TRUE.,MASK3)
            CALL OSDB( 'X=X+Y   ' , CV2 , T1 , SBID , CBID , MESH )
         END IF
C     ------------------------------
C     TERME DE BORD : ONDE IMPOSEE
C     ------------------------------
         IF (NPTFR .GT. 0) THEN
            CALL VECTOR(T1,'=','MASVEC          ',IELMB,                      
     *           1.D0,DPHI4B,S,S,S,S,S,MESH,.TRUE.,MASK4) 
            CALL OSDB( 'X=X+Y   ' , CV2 , T1 , SBID , CBID , MESH )                 
      
         END IF
         
C
C     ---------------------------------------------------------- 
C     CALCUL DE LA MATRICE BM1 POUR LES VALEURS DE MU SPECIFIEES 
C     POUR L'ITERATION 'ITERMU'
C     ----------------------------------------------------------
C
      CALL OS( 'X=YZ    ' , T1 , C  , CG , CBID )
      CALL OS( 'X=YZ    ' , T2 , K  , MU , CBID )
      CALL OS( 'X=YZ    ' , T1 , T1 , T2 , CBID )
      CALL MATRIX(BM1,'M=N     ','FMATMA          ', IELM , IELM ,
     *            1.D0 , T1,S,S,S,S,S,MESH,MSK,MASKEL)
C
C     -------------------------------------------
C     ON AJOUTE A LA MATRICE BM1 LE TERME DE BORD
C     -------------------------------------------
C     
c      IF (NPTFR .GT. 0) THEN
      IF (NPTFR .GT. 0) THEN
C        ------------------------------
C        TERME DE BORD : ONDE INCIDENTE 
C        ------------------------------
         CALL MATRIX(MBOR,'M=N     ','FMATMA          ',IELMB,IELMB,       
     *        -1.D0,APHI1B,S,S,S,S,S,MESH,.TRUE.,MASK1) 
         CALL OM( 'M=M+N   ' , BM1 , MBOR , T1 , CBID , MESH )
      END IF
C     ------------------------------
C        TERME DE BORD : SORTIE LIBRE 
C        ------------------------------
      IF (NPTFR .GT. 0) THEN
            CALL MATRIX(MBOR,'M=N     ','FMATMA          ',IELMB,IELMB,       
     *           -1.D0,APHI2B,S,S,S,S,S,MESH,.TRUE.,MASK2) 
         CALL OM( 'M=M+N   ' , BM1 , MBOR , T1 , CBID , MESH )
      END IF
C        ------------------------------
C        TERME DE BORD : PAROI SOLIDE 
C        ------------------------------
      IF (NPTFR .GT. 0) THEN
         CALL MATRIX(MBOR,'M=N     ','FMATMA          ',IELMB,IELMB,       
     *        -1.D0,APHI3B,S,S,S,S,S,MESH,.TRUE.,MASK3) 
         CALL OM( 'M=M+N   ' , BM1 , MBOR , T1 , CBID , MESH )
      END IF
C        ------------------------------
C        TERME DE BORD : ONDE IMPOSEE
C        ------------------------------
         IF (NPTFR .GT. 0) THEN
            CALL MATRIX(MBOR,'M=N     ','FMATMA          ',IELMB,IELMB,       
     *           -1.D0,APHI4B,S,S,S,S,S,MESH,.TRUE.,MASK4) 
         
            CALL OM( 'M=M+N   ' , BM1 , MBOR , T1 , CBID , MESH )
         END IF
C     --------- 
C     AM2 = AM1
C     ---------
C
      CALL OM( 'M=N     ' , AM2 , AM1 , SBID , CBID , MESH )
C
C     -------------------------- 
C     BM1 DEVIENT NON SYMETRIQUE
C     --------------------------
C 
      CALL OM( 'M=X(M)  ' , BM1 , BM1 , SBID , CBID , MESH )
C
C     ----------------------------
C     ESSAI DE MASS-LUMPING DE BM1
C     ----------------------------
C
C     MASLU = 1.D0
C     CALL LUMP(T1,BM1,MESH,XMESH,MASLU,MSK,MASKEL)                     
C     CALL OM( 'M=CN    ' , BM1 , BM1 , T1 , 1.D0-MASLU , MESH )        
C     CALL OM( 'M=M+D   ' , BM1 , BM1 , T1 , C          , MESH )        
C
C     ----------
C     BM2 = -BM1
C     ----------
C
      CALL OM( 'M=CN    ' , BM2 , BM1 , C , -1.D0 , MESH )
C
C     =======================================
C
C     PRISE EN COMPTE DES POINTS DE DIRICHLET
C
C     =======================================
C               
      IF (DEFERL .OR. FROTTE) THEN
         IF (LNG.EQ.1) WRITE(LU,220) ITERMU+1
         IF (LNG.EQ.2) WRITE(LU,221) ITERMU+1
 220     FORMAT(/,1X,'SOUS-ITERATION NUMERO :',1X,I3,/)
 221     FORMAT(/,1X,'SUB-ITERATION NUMBER :',1X,I3,/)
      ENDIF
      CALL DIRICH(UNK,MAT,RHS,PHIB,LIDIR%I,TB,MESH,KENT,MSK,MASKEL)
C
C     ===============================================================
C
C     SI UN ELEMENT DE DAM1 EST NEGATIF OU NUL, ON INHIBE UN EVENTUEL
C     PRECONDITIONNEMENT DIAGONAL
C
C     ===============================================================
C
      CALL CNTPRE(AM1%D%R,NPOIN,SLVART%PRECON,SLVART%PRECON)
c      IF (LNG.EQ.1) WRITE(LU,230) SLVART%PRECON
c      IF (LNG.EQ.2) WRITE(LU,231) SLVART%PRECON
c 230  FORMAT(/,1X,'PRECONDITIONNEMENT APRES CONTROLE :',1X,I3)
c 231  FORMAT(/,1X,'PRECONDITIONNING AFTER CONTROL :',1X,I3)
C
C     ==========================================================
C
C     PRECONDITIONNEMENT BLOC-DIAGONAL : LES MATRICES DEVIENNENT
C                                       NON SYMETRIQUES.
C
C     ==========================================================
C
      IF (3*(SLVART%PRECON/3).EQ.SLVART%PRECON) THEN
       CALL OM( 'M=X(M)  ' , AM1 , AM1 , SBID , CBID , MESH )
        CALL OM( 'M=X(M)  ' , AM2 , AM2 , SBID , CBID , MESH )
      ENDIF
C
C     ==============================
C
C     RESOLUTION DU SYSTEME LINEAIRE
C
C     ==============================
C
C     ----------------------------
C     INITIALISATION DES INCONNUES
C     ----------------------------
C
      IF(ITERMU.EQ.0.AND.LT.EQ.0) THEN
        CALL LUMP(T1,AM1,MESH,1.D0)                    
        CALL OS( 'X=Y/Z   ' , PHIR , CV1 , T1 , CBID )
        CALL LUMP(T1,AM2,MESH,1.D0)                    
        CALL OS( 'X=Y/Z   ' , PHII , CV2 , T1 , CBID )
      ENDIF
C
      IF (LNG.EQ.1) WRITE(LU,240)
      IF (LNG.EQ.2) WRITE(LU,241)
 240  FORMAT(/,1X,'RESOLUTION DU SYSTEME LINEAIRE (SOLVE)',/)
 241  FORMAT(/,1X,'LINEAR SYSTEM SOLVING (SOLVE)',/)
C
      IF(SLVART%SLV.EQ.8 .OR. SLVART%SLV.EQ.9 ) THEN
C
C       CHANGEMENT DE SIGNE DE LA SECONDE EQUATION
C
       CALL OS('X=-Y    ',X=MAT%ADR(3)%P%D,Y=MAT%ADR(3)%P%D)
       CALL OS('X=-Y    ',X=MAT%ADR(4)%P%D,Y=MAT%ADR(4)%P%D)
       CALL OS('X=-Y    ',X=MAT%ADR(3)%P%X,Y=MAT%ADR(3)%P%X)
       CALL OS('X=-Y    ',X=MAT%ADR(4)%P%X,Y=MAT%ADR(4)%P%X)
       CALL OS('X=-Y    ',X=RHS%ADR(2)%P,Y=RHS%ADR(2)%P)

      ENDIF
C
 


      CALL SOLVE(UNK,MAT,RHS,TB,SLVART,INFOGR,MESH,AM3)

C
C
C     ============================================================
C
C     SI DEFERLEMENT OU FROTTEMENT DE FOND, CALCUL DU COEFFICIENT
C               DE DISSIPATION TOTAL MU_DEFERL + MU_FROTTE
C                                      (MU2)       (T1)
C     ============================================================
C
     
      
      IF (DEFERL .OR. FROTTE) THEN
         ECRHMU = 0.D0
         MODHMU = 0.D0
C     
C     --------------------------------------------
C     INITIALISATION A 0 DES TABLEAUX MU2 ET T3
C        MU2: NOUVEAU COEF. DE DISSIPATION
C     T3 : VALEUR DE QB POUR LA PERIODE CONSIDEREE 
C     --------------------------------------------
C     
         CALL OS( 'X=C     ' , MU2 , SBID , SBID , 0.D0 )
         CALL OS( 'X=C     ' , T3  , SBID , SBID , 0.D0 )
C     
C        ----------------------------------------------------
C     CALCUL DE LA HAUTEUR DE HOULE HMU CORRESPONDANT A LA
C        SOLUTION DU SYSTEME PRECEDENT
C     
         
         CALL OS( 'X=N(Y,Z)', T1  , PHIR , PHII , CBID )
         CALL OS( 'X=CY    ', HMU , T1   , SBID , 2.D0*OMEGA/GRAV )
          
       
C
C        --------------
C        SI DEFERLEMENT
C        --------------
C 
         IF (DEFERL) THEN
C
C        ------------------------------------------------------
C        SI HOULE REGULIERE, ALORS ON TESTE SI HMU > HM (IL Y A
C        DEFERLEMENT) OU NON, ET ON CALCULE MU2 D'APRES LA
C        LOI DE DALLY OU DE BATTJES & JANSSEN
C        ------------------------------------------------------
C
            IF (.NOT. ALEMON .AND. .NOT. ALEMUL) THEN
               DO 20 I = 1,NPOIN
                  HM = 0.88D0/K%R(I)*TANH(GAMMAS*K%R(I)*H%R(I)/0.88D0)
C     
C     HAUTEUR ENERGETIQUE: HMUE = HMU / SQRT(2)
C     
                  HMUE = HMU%R(I)/1.4142D0
                  HEFF=MIN(HMUE,HM)
                  HEFF=MAX(HEFF,1.D-5)
                  Q1 = 1.D-10
                  Q2 = (HEFF/HM)**2.D0
C     AJOUT JMH A CAUSE DU LOG APRES
                  Q2 = MAX(Q2,1.D-9)
C     
C     ------------
C     CALCUL DE QB
C     ------------
C     
                  CALL CALCQB(Q1,Q2,Q3)
C     
     
C     ALGORITHME SPECIFIQUE A LA HOULE REGULIERE
C     POUR LE CALCUL DU TAUX DE DEFERLEMENT
C     
                  IF (ITERMU.EQ.0) THEN
                     IF (Q3.LT.0.19D0) THEN
                        T3%R(I) = 0.D0
                     ELSE
                        T3%R(I) = 1.D0
                     ENDIF
C     
C                 ON MET PROVISOIREMENT LA VALEUR DE T3 CALCULEE
C                 A ITERMU = 0 DANS LA VARIABLE QB
C
                     QB%R(I) = T3%R(I)
                  ELSE
                     IF (QB%R(I).EQ.1.D0) THEN
                        IF (Q3.LT.0.1D0) THEN
                           T3%R(I) = 0.D0
                        ELSE
                           T3%R(I) = 1.D0
                        ENDIF
                     ENDIF
                  ENDIF
 20            CONTINUE
             

C
C           --------------------------------  
C           FORMULATION DE DALLY ET AL. 1985
C           --------------------------------
C
               IF (IBREAK.EQ.2) THEN
                  DO 30 I = 1,NPOIN
                    HM = 0.88D0/K%R(I)*TANH(GAMMAS*K%R(I)*H%R(I)/0.88D0)
                     HEFF=MIN(HMU%R(I),HM)
                     HEFF=MAX(HEFF,1.D-5)
                     MU2%R(I)=T3%R(I)*KDALLY*
     *                    (1.D0-(GDALLY*H%R(I)/HEFF)**2.D0)/H%R(I)
 30               CONTINUE 
               ENDIF
C     
C     ------------------------------------- 
C           FORMULATION DE BATTJES & JANSSEN 1978
C     -------------------------------------
C     
               IF (IBREAK.EQ.1) THEN
                  DO 40 I = 1,NPOIN
                   HM = 0.88D0/K%R(I)*TANH(GAMMAS*K%R(I)*H%R(I)/0.88D0)
                     HEFF=MIN(HMU%R(I),HM)
                     MU2%R(I) = T3%R(I)*2.D0*HEFF/(H%R(I)*CG%R(I)*PER)
 40               CONTINUE 
               ENDIF
               
C     
C     -------------------------------------------------------------
C        SI HOULE ALEATOIRE, ON CALCULE D'ABORD QB=T3, FRACTION DE VAGUE
C     DEFERLANTES OU AYANT DEFERLE, PUIS MU2 D'APRES B&J 78
C     -------------------------------------------------------------
C     
            ELSE
               DO 50 I = 1,NPOIN
                  HM = 0.88D0/K%R(I)*TANH(GAMMAS*K%R(I)*H%R(I)/0.88D0)
C     
C     HAUTEUR ENERGETIQUE: HMUE = HMU / SQRT(2)
C
                  HMUE = HMU%R(I)/1.4142D0
               HEFF=MIN(HMUE,HM)
               HEFF=MAX(HEFF,1.D-5)
               Q1 = 1.D-10
               Q2 = (HEFF/HM)**2.D0
C               AJOUT JMH A CAUSE DU LOG APRES
               Q2 = MAX(Q2,1.D-9)
C
C              ------------
C              CALCUL DE QB
C              ------------
C
               CALL CALCQB(Q1,Q2,Q3)
               T3%R(I) = Q3
C
C              -------------------------
C              CALCUL DU COEFFICIENT MU2
C              -------------------------
C     
               HEFF = MIN((HMU%R(I)/1.4142D0),HM)
               MU2%R(I)=ALFABJ*OMEGA*T3%R(I)*((HM/HEFF)**2.D0)/
     *                (3.14159D0*CG%R(I))
 50         CONTINUE
              
          
         END IF
        
C   
C        ------------------
C        FIN SI DEFERLEMENT
C        ------------------ 
C
         ENDIF
C
C        --------------------------------
C        ON REINITIALISE T1 = 0 CAR APRES
C        T1 REPRESENTE MU_FROTTEMENT
C        --------------------------------
C
         CALL OS( 'X=C     ' , T1 , C , CG , 0.D0 )
C     
C        ---------------------
C        SI FROTTEMENT DE FOND
C        ---------------------
C
         IF (FROTTE) THEN 
C
C           ------------------------------------------------
C           SI ENTFW=VRAI, ON PEUT IMPOSER LA VALEUR DU
C           COEFFICIENT DE FROTTEMENT FW SUR TOUT LE DOMAINE
C           ------------------------------------------------
C
            IF (ENTFW) THEN
               CALL FWSPEC(FW%R,FWCOEF,MESH%X%R,MESH%Y%R,
     *                     NPOIN,PRIVE,ZF%R)
            ELSE
               DO 70 I = 1,NPOIN
                  CALL CALCFW
     *                   (I,H%R,C%R,CG%R,K%R,HMU%R,
     *                    NPOIN,OMEGA,GRAV,
     *                    VISCO,DIAM90,DIAM50,MVSED,MVEAU,
     *                    FORMFR,REGIDO,RICOEF,
     *                    ENTREG,ENTRUG,FFW)
                  FW%R(I) = FFW
 70            CONTINUE
            ENDIF 
C
C           -----------------------------------------
C           CALCUL DU COEFFICIENT DE DISSIPATION POUR
C           LE FROTTEMENT DE FOND 
C           -----------------------------------------
C
            IF (FORMFR .EQ. 1) THEN
C
C           ---------------------------------------------------
C           CALCUL D'UNE VITESSE EFFICACE
C           UE = 1.2*(0.5*((DPHIR/DX)**2 + (DPHIR/DY)**2
C                         +(DPHII/DX)**2 + (DPHII/DY)**2))**0.5
C           ICI, LA VALEUR DE UE EST DANS LA VARIABLE T4 
C           ---------------------------------------------------
C
               CALL CALCUE
C
C              ----------------------------------------
C              LE COEFFICIENT DE DISSIPATION MU POUR LE
C              FROTTEMENT EST MIS DANS LA VARIABLE T1
C              ----------------------------------------
C
               CALL OS( 'X=C     ' , T1 , SBID , SBID , 0.D0 )
C
               DO 80 I = 1,NPOIN 
                  T1%R(I) = (0.5D0*FW%R(I)*T4%R(I))/
     *                    (H%R(I)*((COSH(K%R(I)*H%R(I)))**2.D0))
                  T1%R(I) = T1%R(I)/CG%R(I)
 80            CONTINUE
            ENDIF
C
            IF (FORMFR .EQ. 2) THEN
               CALL OS( 'X=C     ' , T1 , SBID , SBID , 0.D0 )
               DO 90 I = 1,NPOIN
                  T1%R(I) = (2*FW%R(I)*HMU%R(I)*
     *                    ((OMEGA/SINH(K%R(I)*H%R(I)))**3.D0))
                  T1%R(I) = T1%R(I)/(3.D0*3.14159D0*GRAV)
                  T1%R(I) = T1%R(I)/CG%R(I)
 90            CONTINUE
            ENDIF
C
C        -------------------------
C        FIN SI FROTTEMENT DE FOND
C        -------------------------
C
         END IF
C
C        ------------------------------------------------------- 
C        RELAXATION SUR LA VARIABLE MU2 POUR TENTER D'EVITER LES
C        OSCILLATIONS SUR LA CONVERGENCE DU SOLVEUR
C        -------------------------------------------------------
C

         
         DO 60 I = 1,NPOIN
C           
C           --------------------------
C           MU = MU_DEFERL + MU_FROTTE
C           --------------------------
C
            MU2%R(I) = MU2%R(I) + T1%R(I)
C
C           ----------
C           RELAXATION
C           ----------
C
            MU2%R(I) = MU%R(I) + RELDIS * (MU2%R(I) - MU%R(I))
            IF (ITERMU.EQ.0) THEN
               HMUANC%R(I) = HMU%R(I)
               ECRHMU = 1000000.D0
               MODHMU = 1.D0
               MU%R(I) = MU2%R(I)
            ELSE
               ECRHMU = MAX(ECRHMU,ABS(HMU%R(I)-HMUANC%R(I)))
               MODHMU = MAX(MODHMU,ABS(HMU%R(I)))
               MU%R(I) = MU2%R(I)
               HMUANC%R(I) = HMU%R(I)
            ENDIF
 60      CONTINUE

   


C        POUR LA HOULE REGULIERE, ON AMORTIT LA RELAXATION
C        A CHAQUE SOUS-ITERATION POUR FAVORISER LA CONVERGENCE
C        DE L'ALGORITHME DE CALCUL DE LA DISSIPATION
C
         IF (.NOT. ALEMON .AND. .NOT. ALEMUL) THEN
            RELDIS = RELDIS * 0.85D0
         ENDIF
C
         IF (NCSIZE .NE. 0) then
            ECRHMU = P_DMAX(ECRHMU)
            MODHMU = P_DMAX(MODHMU)
         END IF
         IF (LNG.EQ.1) WRITE(LU,*) 'ECART ENTRE DEUX 
     c        SOUS-ITERATIONS (%)',
     c        100*ECRHMU/MODHMU
         IF (LNG.EQ.2) WRITE(LU,*) 'DIFF. BETWEEN TWO 
     c        SUB-ITERATIONS (%) ',
     c        100*ECRHMU/MODHMU
         ITERMU = ITERMU + 1

              
C
C        -----------------------------------------------------------
C        SI NOMBRE DE SOUS-ITER. SUR MU >= NOMBRE MAXI DE SOUS-ITER.
C        ON ARRETE LA BOUCLE SUR MU EN FIXANT L'ECART RELATIF 
C        ECRHMU/MODHMU A 10 % DE EPSDIS
C        -----------------------------------------------------------
C
         IF (ITERMU.GE.NITDIS) THEN
            IF (LNG.EQ.1) WRITE(LU,100) ITERMU
            IF (LNG.EQ.2) WRITE(LU,101) ITERMU
 100        FORMAT(/,1X,'BERKHO (ARTEMIS): NOMBRE DE SOUS-ITERATIONS',
     * 1X,'MAXIMUM ATTEINT :',1X,I3)
 101        FORMAT(/,1X,'BERKHO (ARTEMIS): YOU REACHED THE MAXIMUM',
     * 1X,'NUMBER OF SUB-ITERATIONS :)',1X,I3) 
            ECRHMU = EPSDIS*MODHMU/10.D0
         ENDIF
C     
C     ------------------------------------------------
C        TEST DE CONVERGENCE SUR LA BOUCLE DE DISSIPATION
C        ------------------------------------------------
C
         WRITE(LU,*) ' '
         WRITE(LU,*) '----------------------------------------------- '
         IF (ECRHMU.GT.EPSDIS*MODHMU) GOTO 98
C
         IF (.NOT. ALEMON .AND. .NOT. ALEMUL) THEN
            CALL OS( 'X=Y     ', QB,T3,SBID,CBID)
         ELSE
            CALL OS( 'X=X+Y   ', QB,T3,SBID,CBID)
         ENDIF
         
C
         IF (LNG.EQ.1) WRITE(LU,200) ITERMU
         IF (LNG.EQ.2) WRITE(LU,201) ITERMU
 200     FORMAT(/,1X,'NOMBRE DE SOUS-ITERATIONS POUR LA DISSIPATION:',
     *   1X,I3)  
 201     FORMAT(/,1X,'NUMBER OF SUB-ITERATIONS FOR DISSIPATION:',
     *   1X,I3)
C
C     ========================================
C
C     FIN SI DEFERLEMENT OU FROTTEMENT DE FOND
C
C     ========================================
C
      ENDIF
C
C FIN DE LA BOUCLE D'ITERATIONS SUR LE TERME MU DE DISSIPATION
C
C-----------------------------------------------------------------------
C
      RETURN
      END
