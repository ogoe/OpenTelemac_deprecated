C                       ****************
                        SUBROUTINE PHBOR
C                       ****************
C
C***********************************************************************
C  ARTEMIS V6P0          18/03/2010  C. DENIS (SINETICS)   
C
C  ARTEMIS VERSION 5.1   21/08/00    D. AELBRECHT (LNH) 01 30 87 74 12
C
C  LINKED TO BIEF VERS. 5.0          J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C      FONCTION:    TRADUIT LES CONDITIONS AUX LIMITES
C                   IMPOSEES PAR L'UTILISATEUR EN CALCULANT
C                   LES COEFFICIENTS APHIR, APHII, ...
C                   POUR CHAQUE SEGMENT DE FRONTIERE
C
C
C-----------------------------------------------------------------------
C
C APPELE PAR : ARTEMI
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_ARTEMIS
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      LOGICAL TRVDEB
C
      INTEGER I,IPREC,IG,IG0,IGP1,IFR,IOIDEB(5),IOIFIN(5),ITOTO,IFROI
C
      DOUBLE PRECISION PI,DEGRAD
      DOUBLE PRECISION AUXI1,AUXI2,PHASOI,AUXIC,AUXIS,RADDEG,BID
C
      DOUBLE PRECISION, ALLOCATABLE ::  APHI1BT(:)
      DOUBLE PRECISION, ALLOCATABLE ::  BPHI1BT(:) 
      DOUBLE PRECISION, ALLOCATABLE ::  CPHI1BT(:)
      DOUBLE PRECISION, ALLOCATABLE ::  DPHI1BT(:)
      DOUBLE PRECISION, ALLOCATABLE ::  APHI2BT(:)
      DOUBLE PRECISION, ALLOCATABLE ::  BPHI2BT(:)
      DOUBLE PRECISION, ALLOCATABLE ::  CPHI2BT(:)
      DOUBLE PRECISION, ALLOCATABLE ::  DPHI2BT(:)
      DOUBLE PRECISION, ALLOCATABLE ::  APHI3BT(:)
      DOUBLE PRECISION, ALLOCATABLE ::  BPHI3BT(:)
      DOUBLE PRECISION, ALLOCATABLE ::  CPHI3BT(:)
      DOUBLE PRECISION, ALLOCATABLE ::  DPHI3BT(:)
      DOUBLE PRECISION, ALLOCATABLE ::  APHI4BT(:)
      DOUBLE PRECISION, ALLOCATABLE ::  BPHI4BT(:)
      DOUBLE PRECISION, ALLOCATABLE ::  CPHI4BT(:)
      DOUBLE PRECISION, ALLOCATABLE ::  DPHI4BT(:)
      

      INTRINSIC COS,SIN
C
C-----------------------------------------------------------------------
C
      PARAMETER( PI = 3.1415926535897932384626433D0 , DEGRAD=PI/180.D0 )
      PARAMETER( RADDEG = 180.D0 / PI )
C
C-----------------------------------------------------------------------
C
C ON INITIALISE LIDIR A KSORT (UNE VALEUR DIFFERENTE DE KENT)
C AFIN DE NE PAS ACTIVER LA PRISE EN COMPTE DES POINTS IMPOSES DANS
C PRIDIH QUAND CE N'EST PAS DEMANDE.
C     
     

      IF (NCSIZE .GT. 1) THEN

      ALLOCATE(APHI1BT(NPTFR_TOT))
      ALLOCATE(BPHI1BT(NPTFR_TOT))
      ALLOCATE(CPHI1BT(NPTFR_TOT))
      ALLOCATE(DPHI1BT(NPTFR_TOT))
      ALLOCATE(APHI2BT(NPTFR_TOT))
      ALLOCATE(BPHI2BT(NPTFR_TOT))
      ALLOCATE(CPHI2BT(NPTFR_TOT))
      ALLOCATE(DPHI2BT(NPTFR_TOT))
      ALLOCATE(APHI3BT(NPTFR_TOT))
      ALLOCATE(BPHI3BT(NPTFR_TOT))
      ALLOCATE(CPHI3BT(NPTFR_TOT))
      ALLOCATE(DPHI3BT(NPTFR_TOT))
      ALLOCATE(APHI4BT(NPTFR_TOT))
      ALLOCATE(BPHI4BT(NPTFR_TOT))
      ALLOCATE(CPHI4BT(NPTFR_TOT))
      ALLOCATE(DPHI4BT(NPTFR_TOT))
      ALLOCATE(LIDIRT(2*NPTFR_TOT))

     
      DO I=1,MESH%NPTFR
         APHI1B%R(I) = 0.D0
         BPHI1B%R(I) = 0.D0
         CPHI1B%R(I) = 0.D0
         DPHI1B%R(I) = 0.D0
         APHI2B%R(I) = 0.D0
         BPHI2B%R(I) = 0.D0
         CPHI2B%R(I) = 0.D0
         DPHI2B%R(I) = 0.D0
         APHI3B%R(I) = 0.D0
         BPHI3B%R(I) = 0.D0
         CPHI3B%R(I) = 0.D0
         DPHI3B%R(I) = 0.D0
         APHI4B%R(I) = 0.D0
         BPHI4B%R(I) = 0.D0
         CPHI4B%R(I) = 0.D0
         DPHI4B%R(I) = 0.D0
      END DO
      


      
        DO I=1,NPTFR_TOT
           LIDIRT(I) = KSORT
C     ATTENTION ON SUPPOSE ICI QUE NPTFRX=NPTFR
           LIDIRT(I+NPTFR_TOT) = KSORT
           IF (LIHBORT(I).EQ.KENT) THEN
              LIHBORT(I) = KINC
           ENDIF
           APHI1BT(I) = 0.D0
           BPHI1BT(I) = 0.D0
           CPHI1BT(I) = 0.D0
           DPHI1BT(I) = 0.D0
           APHI2BT(I) = 0.D0
           BPHI2BT(I) = 0.D0
           CPHI2BT(I) = 0.D0
           DPHI2BT(I) = 0.D0
           APHI3BT(I) = 0.D0
           BPHI3BT(I) = 0.D0
           CPHI3BT(I) = 0.D0
           DPHI3BT(I) = 0.D0
           APHI4BT(I) = 0.D0
           BPHI4BT(I) = 0.D0
           CPHI4BT(I) = 0.D0
           DPHI4BT(I) = 0.D0
          
        END DO
      


C
C-----------------------------------------------------------------------
C
C
C     ************************************************
C     INITIALISATION DE LA PHASE POUR L'ONDE INCIDENTE
C     ************************************************
C
      PHASOI = 0.D0
C
C     ******************************************
C     TRAITEMENT PARTICULIER DE L'ONDE INCIDENTE
C     ******************************************
C
C     -------------------------------------------------
C     REPERAGE DES DEBUTS DES FRONTIERES ONDE INCIDENTE
C     -------------------------------------------------
C
      TRVDEB = .TRUE.
      IFROI = 0
C
      DO 10 I=1,NPTFR_TOT
         IF (LIHBORT(I).EQ.KINC) THEN
            ITOTO = KP1BOR_TOT(I+NPTFR_TOT)
               IF (LIHBORT(ITOTO).NE.KINC) THEN
                  IFROI = IFROI + 1
                  IOIDEB(IFROI) = I
               ENDIF
         ENDIF
 10   CONTINUE
     
C
      IF(LNG.EQ.1) WRITE(LU,11) IFROI
      IF(LNG.EQ.2) WRITE(LU,12) IFROI
11    FORMAT(1X,'PHBOR : IL Y A : ',1I3,' FRONTIERE(S) ',
     *       1X,'DE TYPE ONDE INCIDENTE ')
12    FORMAT(1X,'PHBOR : THERE ARE :',1I3,' BOUNDARIE(S) ',
     *       1X,'OF INCIDENT WAVE TYPE ')
C
C     --------------------------------------------------------------
C     CALCUL DES COEFFICIENTS POUR LES FRONTIERES D'ONDE INCIDENTE
C     A PARTIR DE IOIDEB (DEBUT ONDE INCIDENTE)
C     --------------------------------------------------------------
C    
   
      
      DO 15 IFR=1,IFROI
         I = IOIDEB(IFR)
         
C
 20   CONTINUE
C
C        ********************************
C        NUMERO GLOBAL DU POINT DE BORD I
C        ********************************
C
C         IG   = MESH%NBOR%I(I)
      IG   = NBOR_TOT(I)
C
C        ******************************************
C        NUMERO GLOBAL DU POINT DE BORD PRECEDENT I
C        ******************************************
C
      IG0  = NBOR_TOT(KP1BOR_TOT(I+NPTFR_TOT))
C
C        ****************************************
C        NUMERO GLOBAL DU POINT DE BORD SUIVANT I
C        ****************************************
C
        
      IGP1 = NBOR_TOT(KP1BOR_TOT(I))

             
         AUXIC      = COS(TETABT(I)*DEGRAD)
         AUXIS      = SIN(TETABT(I)*DEGRAD)
         AUXI1      = GRAV/OMEGA * HBT(I)/2.D0 * 
     *                CTT(IG) * CGT(IG) * KT(IG)
C
C
C           DEVELOPPEMENT POUR LA DIRECTION NON UNIFORME
C           PHASOI EST LA PHASE DE LA HOULE
C
C           ANCIENNE FORMULE : PHASE = K * X (ONDE PLANE) :
C      
C                     PHASOI = KM * ( XM*AUXIC + YM*AUXIS )
C
C           NOUVELLE FORMULATION (ONDE NON PLANE) :
C
C                                   M
C                                 / 
C                    PHASOI(M) = /        K(N) * dX(N) 
C                               / 
C                                Mdeb
C
C           Ou Mdeb est le point de debut du segment d'onde incidente,
C           a partir duquel on calcule la phase. L'integrale ci-dessus est
C           simplement calculee par approximation lineaire
C
C           On ajoute l'eventuel dephasage ALFAP pour assurer la coherence
C           des phases des differentes cretes de houle s'il y a plusieurs
C           frontieres d'ondes incidentes non continues
C
      PHASOI = PHASOI + KT(IG)*AUXIC*(XT(IG) - XT(IG0))
     *                + KT(IG)*AUXIS*(YT(IG) - YT(IG0))
C
      APHI1BT(I)  = - KT(IG) * CTT(IG) * CGT(IG)
     *             * COS(TETAPT(I)*DEGRAD)
C
      BPHI1BT(I)  = 0.D0
C
      CPHI1BT(I)  = AUXI1 * COS( PHASOI + ALFAPT(I)*DEGRAD )
C
      DPHI1BT(I)  = AUXI1 * SIN( PHASOI + ALFAPT(I)*DEGRAD )
C
      I = KP1BOR_TOT(I)
C
C     ON VA JUSQU'AU POINT SUIVANT LA FIN D'UNE FRONTIERE O.I.
C
      IF (LIHBORT(I).NE.KINC) THEN
         IOIFIN(IFR) = I
        
         IPREC      = KP1BOR_TOT(I+NPTFR_TOT)

         TETABT(I) = TETABT(IPREC)
         HBT(I)    = HBT(IPREC)
         AUXIC      = COS(TETABT(I)*DEGRAD)
         AUXIS      = SIN(TETABT(I)*DEGRAD)
         AUXI1      = GRAV/OMEGA * HBT(I)/2.D0 * 
     *                CTT(IG) * CGT(IG) * KT(IG)
         PHASOI = PHASOI + KT(IG)*AUXIC*(XT(IG) - XT(IG0))
     *                   + KT(IG)*AUXIS*(YT(IG) - YT(IG0))
C
         APHI1BT(I) = - KT(IG) * CTT(IG) * CGT(IG)
     *                 * COS(TETAPT(IPREC)*DEGRAD)
C
         BPHI1BT(I) = 0.D0
C
         CPHI1BT(I) = AUXI1*COS(PHASOI + ALFAPT(IPREC)*DEGRAD)
C
         DPHI1BT(I) = AUXI1*SIN(PHASOI + ALFAPT(IPREC)*DEGRAD)
C
         GOTO 15
C
      ELSE
         GOTO 20
      ENDIF
C
15    CONTINUE
C
C     ******************************************
C     FIN TRAITEMENT DE L'ONDE INCIDENTE
C     ******************************************
C
      DO 100 I=1,NPTFR_TOT
C
C        ********************************
C        NUMERO GLOBAL DU POINT DE BORD I
C        ********************************
C
         IG   = NBOR_TOT(I)
C
C        ******************************************
C        NUMERO GLOBAL DU POINT DE BORD PRECEDENT I
C        ******************************************
C
         IG0  = NBOR_TOT(KP1BOR_TOT(I+NPTFR_TOT))
C
C        ****************************************
C        NUMERO GLOBAL DU POINT DE BORD SUIVANT I
C        ****************************************
C
         IGP1 = NBOR_TOT(KP1BOR_TOT(I))
CCPHI1B%R
C        -------------------------------------------------
C        COEFFICIENTS POUR UN SEGMENT DE BORD SORTIE LIBRE
C        -------------------------------------------------
C
         IF (LIHBORT(I).EQ.KSORT) THEN
            APHI2BT(I)  = - KT(IG) * CTT(IG) * CGT(IG)
     *                   * COS(TETAPT(I)*DEGRAD)
C
             BPHI2BT(I)  = 0.D0
C
            CPHI2BT(I)  = 0.D0
C
            DPHI2BT(I)  = 0.D0
C
         ELSEIF (LIHBORT(KP1BOR_TOT(I)).EQ.KSORT) THEN
            APHI2BT(I)  = - KT(IG) * CTT(IG) * CGT(IG)
     *                   * COS(TETAPT(KP1BOR_TOT(I))*DEGRAD)
C
             BPHI2BT(I)  = 0.D0
C
            CPHI2BT(I)  = 0.D0
C
            DPHI2BT(I)  = 0.D0
C
         ELSEIF (LIHBORT(KP1BOR_TOT(I+NPTFR_TOT)).EQ.KSORT) THEN
            APHI2BT(I)  = - KT(IG) * CTT(IG) * CGT(IG)
     *              * COS(TETAPT(KP1BOR_TOT(I+NPTFR_TOT))*DEGRAD)
C
             BPHI2BT(I)  = 0.D0
C
            CPHI2BT(I)  = 0.D0
C
            DPHI2BT(I)  = 0.D0
C
         ELSE
            APHI2BT(I)  = 0.D0
C
             BPHI2BT(I)  = 0.D0
C
            CPHI2BT(I)  = 0.D0
C
            DPHI2BT(I)  = 0.D0
C
         ENDIF
C
C        -------------------------------------------
C        COEFFICIENTS POUR UN SEGMENT DE BORD SOLIDE
C        -------------------------------------------
C
         IF (LIHBORT(I).EQ.KLOG) THEN
          AUXI1 = KT(IG) * CTT(IG) * CGT(IG) * 
     *      COS(TETAPT(I)*DEGRAD) /
     *      ( 1.D0 + RPT(I)*RPT(I) + 
     *        2.D0*RPT(I)*COS(ALFAPT(I)*DEGRAD) )
C
          APHI3BT(I) = - (1.D0 - RPT(I) * RPT(I) ) * AUXI1
C
          BPHI3BT(I) = 2.D0*RPT(I)*SIN(ALFAPT(I)*DEGRAD) * AUXI1
C
          CPHI3BT(I)  = 0.D0
C
          DPHI3BT(I)  = 0.D0
C
         ELSEIF (LIHBORT(KP1BOR_TOT(I)).EQ.KLOG) THEN
          AUXI1 = KT(IG) * CTT(IG) * CGT(IG) *
     *      COS(TETAPT(KP1BOR_TOT(I))*DEGRAD) /
     *      (1.D0 + RPT(KP1BOR_TOT(I))*RPT(KP1BOR_TOT(I))
     *      +2.D0 * RPT(KP1BOR_TOT(I))*
     *       COS(ALFAPT(KP1BOR_TOT(I))*DEGRAD))
C
          APHI3BT(I) = - (1.D0-RPT(KP1BOR_TOT(I))*
     *      RPT(KP1BOR_TOT(I))) * AUXI1
C
          BPHI3BT(I) = 2.D0*RPT(KP1BOR_TOT(I))
     *                * SIN(ALFAPT(KP1BOR_TOT(I))*DEGRAD) * AUXI1
C
          CPHI3BT(I)  = 0.D0
C
          DPHI3BT(I)  = 0.D0
C
         ELSEIF (LIHBORT(KP1BOR_TOT(I+NPTFR_TOT)).EQ.KLOG) THEN
          AUXI1 = KT(IG) * CTT(IG) * CGT(IG) *
     *     COS(TETAPT(KP1BOR_TOT(I+NPTFR_TOT))*DEGRAD) /
     *     (1.D0 + RPT(KP1BOR_TOT(I+NPTFR_TOT))
     *      *RPT(KP1BOR_TOT(I+NPTFR_TOT))
     *      +2.D0 * RPT(KP1BOR_TOT(I+NPTFR_TOT))*
     *       COS(ALFAPT(KP1BOR_TOT(I+NPTFR_TOT))*DEGRAD))
C
          APHI3BT(I) = - (1.D0-RPT(KP1BOR_TOT(I+NPTFR_TOT))*
     *      RPT(KP1BOR_TOT(I+NPTFR_TOT))) * AUXI1
C
          BPHI3BT(I) = 2.D0*RPT(KP1BOR_TOT(I+NPTFR_TOT))
     *      * SIN(ALFAPT(KP1BOR_TOT(I+NPTFR_TOT))*DEGRAD) * AUXI1
C
          CPHI3BT(I)  = 0.D0
C
          DPHI3BT(I)  = 0.D0
C
         ELSE
          APHI3BT(I)  = 0.D0
C
          BPHI3BT(I)  = 0.D0
C
          CPHI3BT(I)  = 0.D0
C
          DPHI3BT(I)  = 0.D0
C
         ENDIF
C
C        -------------------------------------------------
C        COEFFICIENTS POUR UN SEGMENT DE BORD ONDE IMPOSEE
C        -------------------------------------------------
CDA      ----------------------------------- 
CDA      ON LAISSE CES LIGNES POUR MEMOIRE !
CDA      ----------------------------------- 
CDA
CDA         IF (LIHBOR(I).EQ.KENT) THEN
CDA         AUXIC      = COS(TETAB(I)*DEGRAD)
CDA         AUXIS      = SIN(TETAB(I)*DEGRAD)
CDA         AUXI1      = GRAV/OMEGA * HB(I)/2.D0 * C(IG) * CG(IG) *
CDA     *                K(IG) * ( AUXIC *XSGBOR(I) +
CDA     *                          AUXIS *YSGBOR(I) )
CDA         AUXI2      = K(IG) * ( X(IG)*AUXIC +
CDA     *                          Y(IG)*AUXIS )
CDA
CDA         APHI4B(I)  = 0.D0
CDA
CDA         BPHI4B(I)  = 0.D0
CDA
CDA         CPHI4B(I)  = AUXI1 * COS( AUXI2 )
CDA
CDA         DPHI4B(I)  = AUXI1 * SIN( AUXI2 )
CDA
CDA       VALEURS IMPOSEES AUX NOEUDS DU SEGMENT D'ENTREE
CDA         LIDIR(I)         = KENT
CDA
CDA         AUXI1 = GRAV/OMEGA * HB(I)/2.D0
CDA         AUXI2 = K(IG) * (X(IG)*AUXIC +
CDA     *                    Y(IG)*AUXIS )
CDA
CDA            PHIRB(I) =   AUXI1 * SIN( AUXI2 )
CDA            PHIIB(I) = - AUXI1 * COS( AUXI2 )
CDA         ENDIF
C
C
100   CONTINUE
C-----------------------------------------------------------------------
C
        CALL GLOBAL_TO_LOCAL_BOUND(APHI1BT,APHI1B,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(BPHI1BT,BPHI1B,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(CPHI1BT,CPHI1B,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(DPHI1BT,DPHI1B,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(APHI2BT,APHI2B,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(BPHI2BT,BPHI2B,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(CPHI2BT,CPHI2B,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(DPHI2BT,DPHI2B,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(APHI3BT,APHI3B,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(BPHI3BT,BPHI3B,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(CPHI3BT,CPHI3B,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(DPHI3BT,DPHI3B,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(APHI4BT,APHI4B,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(BPHI4BT,BPHI4B,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(CPHI4BT,CPHI4B,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(DPHI4BT,DPHI4B,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUNDI(LIDIRT,LIDIR,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(TETAPT,TETAP,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(TETABT,TETAB,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(RPT,RP,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(HBT,HB,MESH%NPTFR,NPTFR_TOT)
        CALL GLOBAL_TO_LOCAL_BOUND(ALFAPT,ALFAP,MESH%NPTFR,NPTFR_TOT)
c        CALL GLOBAL_TO_LOCAL_BOUND2(CTT,C,MESH%NPOIN,NPOIN_TOT)
c        CALL GLOBAL_TO_LOCAL_BOUND2(CGT,CG,MESH%NPOIN,NPOIN_TOT)
c        CALL GLOBAL_TO_LOCAL_BOUND2(KT,K,MESH%NPOIN,NPOIN_TOT)


        DEALLOCATE(APHI1BT)
        DEALLOCATE(BPHI1BT)
        DEALLOCATE(CPHI1BT)
        DEALLOCATE(DPHI1BT)
        DEALLOCATE(APHI2BT)
        DEALLOCATE(BPHI2BT)
        DEALLOCATE(CPHI2BT)
        DEALLOCATE(DPHI2BT)
        DEALLOCATE(APHI3BT)
        DEALLOCATE(BPHI3BT)
        DEALLOCATE(CPHI3BT)
        DEALLOCATE(DPHI3BT)
        DEALLOCATE(APHI4BT)
        DEALLOCATE(BPHI4BT)
        DEALLOCATE(CPHI4BT)
        DEALLOCATE(DPHI4BT)

        DEALLOCATE(LIDIRT)
c        DEALLOCATE(XT)
c        DEALLOCATE(YT)
c        DEALLOCATE(CTT)
c        DEALLOCATE(KT)
c        DEALLOCATE(CGT)
      DO 110 IFR = 1,IFROI
         I          = IOIFIN(IFR)
         IPREC      = KP1BOR_TOT(I+NPTFR_TOT)
         TETAPT(I) = TETAPT(IPREC)
110   CONTINUE
C 
      ELSE



  
        DO I=1,NPTFR
                LIHBOR%I(I)=LIHBORT(I)
                RP%R(I)=RPT(I)
                HB%R(I)=HBT(I)
                ALFAP%R(I)=ALFAPT(I)
                TETAB%R(I)=TETABT(I)
                TETAP%R(I)=TETAPT(I)
         END DO       



         DO 501 I=1,NPTFR
        LIDIR%I(I) = KSORT
C       ATTENTION ON SUPPOSE ICI QUE NPTFRX=NPTFR
        LIDIR%I(I+NPTFR) = KSORT
        IF (LIHBOR%I(I).EQ.KENT) THEN
           LIHBOR%I(I) = KINC
        ENDIF
        APHI1B%R(I) = 0.D0
        BPHI1B%R(I) = 0.D0
        CPHI1B%R(I) = 0.D0
        DPHI1B%R(I) = 0.D0
        APHI2B%R(I) = 0.D0
        BPHI2B%R(I) = 0.D0
        CPHI2B%R(I) = 0.D0
        DPHI2B%R(I) = 0.D0
        APHI3B%R(I) = 0.D0
        BPHI3B%R(I) = 0.D0
        CPHI3B%R(I) = 0.D0
        DPHI3B%R(I) = 0.D0
        APHI4B%R(I) = 0.D0
        BPHI4B%R(I) = 0.D0
        CPHI4B%R(I) = 0.D0
        DPHI4B%R(I) = 0.D0
 501    CONTINUE
C
C-----------------------------------------------------------------------
C
C
C     ************************************************
C     INITIALISATION DE LA PHASE POUR L'ONDE INCIDENTE
C     ************************************************
C
      PHASOI = 0.D0
C
C     ******************************************
C     TRAITEMENT PARTICULIER DE L'ONDE INCIDENTE
C     ******************************************
C
C     -------------------------------------------------
C     REPERAGE DES DEBUTS DES FRONTIERES ONDE INCIDENTE
C     -------------------------------------------------
C
      TRVDEB = .TRUE.
      IFROI = 0
C
      DO 101 I=1,NPTFR
         IF (LIHBOR%I(I).EQ.KINC) THEN
               ITOTO = MESH%KP1BOR%I(I+NPTFR)
               IF (LIHBOR%I(ITOTO).NE.KINC) THEN
                  IFROI = IFROI + 1
                  IOIDEB(IFROI) = I
               ENDIF
         ENDIF
 101  CONTINUE

     
C
      IF(LNG.EQ.1) WRITE(LU,111) IFROI
      IF(LNG.EQ.2) WRITE(LU,121) IFROI
 111  FORMAT(1X,'PHBOR : IL Y A : ',1I3,' FRONTIERE(S) ',
     *       1X,'DE TYPE ONDE INCIDENTE ')
 121  FORMAT(1X,'PHBOR : THERE ARE :',1I3,' BOUNDARIE(S) ',
     *       1X,'OF INCIDENT WAVE TYPE ')
C
C     --------------------------------------------------------------
C     CALCUL DES COEFFICIENTS POUR LES FRONTIERES D'ONDE INCIDENTE
C     A PARTIR DE IOIDEB (DEBUT ONDE INCIDENTE)
C     --------------------------------------------------------------
C

   
      DO 151 IFR=1,IFROI
         I = IOIDEB(IFR)
C
 201     CONTINUE
C
C        ********************************
C        NUMERO GLOBAL DU POINT DE BORD I
C        ********************************
C
         IG   = MESH%NBOR%I(I)
C
C        ******************************************
C        NUMERO GLOBAL DU POINT DE BORD PRECEDENT I
C        ******************************************
C
         IG0  = MESH%NBOR%I(MESH%KP1BOR%I(I+NPTFR))
C
C        ****************************************
C        NUMERO GLOBAL DU POINT DE BORD SUIVANT I
C        ****************************************
C
         IGP1 = MESH%NBOR%I(MESH%KP1BOR%I(I))
C
     

         AUXIC      = COS(TETAB%R(I)*DEGRAD)
         AUXIS      = SIN(TETAB%R(I)*DEGRAD)
         AUXI1      = GRAV/OMEGA * HB%R(I)/2.D0 * 
     *                C%R(IG) * CG%R(IG) * K%R(IG)
C
C
C           DEVELOPPEMENT POUR LA DIRECTION NON UNIFORME
C           PHASOI EST LA PHASE DE LA HOULE
C
C           ANCIENNE FORMULE : PHASE = K * X (ONDE PLANE) :
C      
C                     PHASOI = KM * ( XM*AUXIC + YM*AUXIS )
C
C           NOUVELLE FORMULATION (ONDE NON PLANE) :
C
C                                   M
C                                 / 
C                    PHASOI(M) = /        K(N) * dX(N) 
C                               / 
C                                Mdeb
C
C           Ou Mdeb est le point de debut du segment d'onde incidente,
C           a partir duquel on calcule la phase. L'integrale ci-dessus est
C           simplement calculee par approximation lineaire
C
C           On ajoute l'eventuel dephasage ALFAP pour assurer la coherence
C           des phases des differentes cretes de houle s'il y a plusieurs
C           frontieres d'ondes incidentes non continues
C
      PHASOI = PHASOI + K%R(IG)*AUXIC*(X(IG) - X(IG0))
     *                + K%R(IG)*AUXIS*(Y(IG) - Y(IG0))
C
      APHI1B%R(I)  = - K%R(IG) * C%R(IG) * CG%R(IG)
     *             * COS(TETAP%R(I)*DEGRAD)
C
      BPHI1B%R(I)  = 0.D0
C
      CPHI1B%R(I)  = AUXI1 * COS( PHASOI + ALFAP%R(I)*DEGRAD )
C
      DPHI1B%R(I)  = AUXI1 * SIN( PHASOI + ALFAP%R(I)*DEGRAD )
C
      I = MESH%KP1BOR%I(I)
C
C     ON VA JUSQU'AU POINT SUIVANT LA FIN D'UNE FRONTIERE O.I.
C
      IF (LIHBOR%I(I).NE.KINC) THEN
         IOIFIN(IFR) = I
         IPREC      = MESH%KP1BOR%I(I+NPTFR)
         TETAB%R(I) = TETAB%R(IPREC)
         HB%R(I)    = HB%R(IPREC)
         AUXIC      = COS(TETAB%R(I)*DEGRAD)
         AUXIS      = SIN(TETAB%R(I)*DEGRAD)
         AUXI1      = GRAV/OMEGA * HB%R(I)/2.D0 * 
     *                C%R(IG) * CG%R(IG) * K%R(IG)
         PHASOI = PHASOI + K%R(IG)*AUXIC*(X(IG) - X(IG0))
     *                   + K%R(IG)*AUXIS*(Y(IG) - Y(IG0))
C
         APHI1B%R(I) = - K%R(IG) * C%R(IG) * CG%R(IG)
     *                 * COS(TETAP%R(IPREC)*DEGRAD)
C
         BPHI1B%R(I) = 0.D0
C
         CPHI1B%R(I) = AUXI1*COS(PHASOI + ALFAP%R(IPREC)*DEGRAD)
C
         DPHI1B%R(I) = AUXI1*SIN(PHASOI + ALFAP%R(IPREC)*DEGRAD)
C
         GOTO 151
C
      ELSE
         GOTO 201
      ENDIF
C
 151  CONTINUE
C
C     ******************************************
C     FIN TRAITEMENT DE L'ONDE INCIDENTE
C     ******************************************
C
      DO 1001 I=1,NPTFR
C
C        ********************************
C        NUMERO GLOBAL DU POINT DE BORD I
C        ********************************
C
         IG   = MESH%NBOR%I(I)
C
C        ******************************************
C        NUMERO GLOBAL DU POINT DE BORD PRECEDENT I
C        ******************************************
C
         IG0  = MESH%NBOR%I(MESH%KP1BOR%I(I+NPTFR))
C
C        ****************************************
C        NUMERO GLOBAL DU POINT DE BORD SUIVANT I
C        ****************************************
C
         IGP1 = MESH%NBOR%I(MESH%KP1BOR%I(I))
C
C        -------------------------------------------------
C        COEFFICIENTS POUR UN SEGMENT DE BORD SORTIE LIBRE
C        -------------------------------------------------
C
         IF (LIHBOR%I(I).EQ.KSORT) THEN
            APHI2B%R(I)  = - K%R(IG) * C%R(IG) * CG%R(IG)
     *                   * COS(TETAP%R(I)*DEGRAD)
C
            BPHI2B%R(I)  = 0.D0
C
            CPHI2B%R(I)  = 0.D0
C
            DPHI2B%R(I)  = 0.D0
C
         ELSEIF (LIHBOR%I(MESH%KP1BOR%I(I)).EQ.KSORT) THEN
            APHI2B%R(I)  = - K%R(IG) * C%R(IG) * CG%R(IG)
     *                   * COS(TETAP%R(MESH%KP1BOR%I(I))*DEGRAD)
C
            BPHI2B%R(I)  = 0.D0
C
            CPHI2B%R(I)  = 0.D0
C
            DPHI2B%R(I)  = 0.D0
C
         ELSEIF (LIHBOR%I(MESH%KP1BOR%I(I+NPTFR)).EQ.KSORT) THEN
            APHI2B%R(I)  = - K%R(IG) * C%R(IG) * CG%R(IG)
     *              * COS(TETAP%R(MESH%KP1BOR%I(I+NPTFR))*DEGRAD)
C
            BPHI2B%R(I)  = 0.D0
C
            CPHI2B%R(I)  = 0.D0
C
            DPHI2B%R(I)  = 0.D0
C
         ELSE
            APHI2B%R(I)  = 0.D0
C
            BPHI2B%R(I)  = 0.D0
C
            CPHI2B%R(I)  = 0.D0
C
            DPHI2B%R(I)  = 0.D0
C
         ENDIF
C
C        -------------------------------------------
C        COEFFICIENTS POUR UN SEGMENT DE BORD SOLIDE
C        -------------------------------------------
C
         IF (LIHBOR%I(I).EQ.KLOG) THEN
          AUXI1 = K%R(IG) * C%R(IG) * CG%R(IG) * 
     *      COS(TETAP%R(I)*DEGRAD) /
     *      ( 1.D0 + RP%R(I)*RP%R(I) + 
     *        2.D0*RP%R(I)*COS(ALFAP%R(I)*DEGRAD) )
C
          APHI3B%R(I) = - (1.D0 - RP%R(I) * RP%R(I) ) * AUXI1
C
          BPHI3B%R(I) = 2.D0*RP%R(I)*SIN(ALFAP%R(I)*DEGRAD) * AUXI1
C
          CPHI3B%R(I)  = 0.D0
C
          DPHI3B%R(I)  = 0.D0
C
         ELSEIF (LIHBOR%I(MESH%KP1BOR%I(I)).EQ.KLOG) THEN
          AUXI1 = K%R(IG) * C%R(IG) * CG%R(IG) *
     *      COS(TETAP%R(MESH%KP1BOR%I(I))*DEGRAD) /
     *      (1.D0 + RP%R(MESH%KP1BOR%I(I))*RP%R(MESH%KP1BOR%I(I))
     *      +2.D0 * RP%R(MESH%KP1BOR%I(I))*
     *       COS(ALFAP%R(MESH%KP1BOR%I(I))*DEGRAD))
C
          APHI3B%R(I) = - (1.D0-RP%R(MESH%KP1BOR%I(I))*
     *      RP%R(MESH%KP1BOR%I(I))) * AUXI1
C
          BPHI3B%R(I) = 2.D0*RP%R(MESH%KP1BOR%I(I))
     *                * SIN(ALFAP%R(MESH%KP1BOR%I(I))*DEGRAD) * AUXI1
C
          CPHI3B%R(I)  = 0.D0
C
          DPHI3B%R(I)  = 0.D0
C
         ELSEIF (LIHBOR%I(MESH%KP1BOR%I(I+NPTFR)).EQ.KLOG) THEN
          AUXI1 = K%R(IG) * C%R(IG) * CG%R(IG) *
     *     COS(TETAP%R(MESH%KP1BOR%I(I+NPTFR))*DEGRAD) /
     *     (1.D0 + RP%R(MESH%KP1BOR%I(I+NPTFR))
     *      *RP%R(MESH%KP1BOR%I(I+NPTFR))
     *      +2.D0 * RP%R(MESH%KP1BOR%I(I+NPTFR))*
     *       COS(ALFAP%R(MESH%KP1BOR%I(I+NPTFR))*DEGRAD))
C
          APHI3B%R(I) = - (1.D0-RP%R(MESH%KP1BOR%I(I+NPTFR))*
     *      RP%R(MESH%KP1BOR%I(I+NPTFR))) * AUXI1
C
          BPHI3B%R(I) = 2.D0*RP%R(MESH%KP1BOR%I(I+NPTFR))
     *      * SIN(ALFAP%R(MESH%KP1BOR%I(I+NPTFR))*DEGRAD) * AUXI1
C
          CPHI3B%R(I)  = 0.D0
C
          DPHI3B%R(I)  = 0.D0
C
         ELSE
          APHI3B%R(I)  = 0.D0
C
          BPHI3B%R(I)  = 0.D0
C
          CPHI3B%R(I)  = 0.D0
C
          DPHI3B%R(I)  = 0.D0
C
         ENDIF
C
C        -------------------------------------------------
C        COEFFICIENTS POUR UN SEGMENT DE BORD ONDE IMPOSEE
C        -------------------------------------------------
CDA      ----------------------------------- 
CDA      ON LAISSE CES LIGNES POUR MEMOIRE !
CDA      ----------------------------------- 
CDA
CDA         IF (LIHBOR(I).EQ.KENT) THEN
CDA         AUXIC      = COS(TETAB(I)*DEGRAD)
CDA         AUXIS      = SIN(TETAB(I)*DEGRAD)
CDA         AUXI1      = GRAV/OMEGA * HB(I)/2.D0 * C(IG) * CG(IG) *
CDA     *                K(IG) * ( AUXIC *XSGBOR(I) +
CDA     *                          AUXIS *YSGBOR(I) )
CDA         AUXI2      = K(IG) * ( X(IG)*AUXIC +
CDA     *                          Y(IG)*AUXIS )
CDA
CDA         APHI4B(I)  = 0.D0
CDA
CDA         BPHI4B(I)  = 0.D0
CDA
CDA         CPHI4B(I)  = AUXI1 * COS( AUXI2 )
CDA
CDA         DPHI4B(I)  = AUXI1 * SIN( AUXI2 )
CDA
CDA       VALEURS IMPOSEES AUX NOEUDS DU SEGMENT D'ENTREE
CDA         LIDIR(I)         = KENT
CDA
CDA         AUXI1 = GRAV/OMEGA * HB(I)/2.D0
CDA         AUXI2 = K(IG) * (X(IG)*AUXIC +
CDA     *                    Y(IG)*AUXIS )
CDA
CDA            PHIRB(I) =   AUXI1 * SIN( AUXI2 )
CDA            PHIIB(I) = - AUXI1 * COS( AUXI2 )
CDA         ENDIF
C
C
 1001 CONTINUE
C-----------------------------------------------------------------------

    
      


      DO 1101 IFR = 1,IFROI
         I          = IOIFIN(IFR)
         IPREC      = MESH%KP1BOR%I(I+NPTFR)
         TETAP%R(I) = TETAP%R(IPREC)
 1101 CONTINUE
C 
c      DO I=1,NPTFR
c         ALFAP%R(I)=ALFAPT(I)
c         TETAP%R(I)=TETAPT(I)
c         TETAB%R(I)=TETABT(I)
c         RP%R(I)=RPT(I)
c         HB%R(I)=HBT(I)
c      END DO



      END IF
      
      
      write(lu,*) 'END PHBOR' 
      RETURN 
      CONTAINS
      SUBROUTINE GLOBAL_TO_LOCAL_BOUND(TAB1,OBJ,NPTFR,NPTFR_TOT)
      DOUBLE PRECISION, INTENT(INOUT)  :: TAB1(:)
      TYPE (BIEF_OBJ), INTENT(INOUT) :: OBJ
      INTEGER, INTENT(IN) :: NPTFR
      INTEGER, INTENT(IN) :: NPTFR_TOT
      INTEGER :: I,J
      OBJ%R=0.0
      DO I=1,NPTFR_TOT
         DO J=1,MESH%NPTFR
            IF (NBOR_TOT(I) .EQ. MESH%KNOLG%I(MESH%NBOR%I(J))) THEN
               OBJ%R(J)=TAB1(I)
            END IF
         END DO
      END DO
      OBJ%DIM1=NPTFR
   

      END SUBROUTINE

      SUBROUTINE GLOBAL_TO_LOCAL_BOUNDi(TAB1,OBJ,NPTFR,NPTFR_TOT)
      INTEGER, INTENT(INOUT)  :: TAB1(:)
      TYPE (BIEF_OBJ), INTENT(INOUT) :: OBJ
      INTEGER, INTENT(IN) :: NPTFR
      INTEGER, INTENT(IN) :: NPTFR_TOT
      INTEGER :: I,J
      OBJ%I=0.0
      DO I=1,NPTFR_TOT
         DO J=1,MESH%NPTFR
            IF (NBOR_TOT(I) .EQ. MESH%KNOLG%I(MESH%NBOR%I(J))) THEN
               OBJ%I(J)=TAB1(I)
            END IF
         END DO
      END DO
      DO  I=NPTFR_TOT+1,2*NPTFR_TOT
         DO J=1,MESH%NPTFR
            IF (NBOR_TOT(I/2) .EQ. MESH%KNOLG%I(MESH%NBOR%I(J))) THEN
               OBJ%I(2*J)=TAB1(I)
            END IF
         END DO
      END DO
      END SUBROUTINE

       SUBROUTINE GLOBAL_TO_LOCAL_BOUND2(TAB1,OBJ,NPOIN,NPOIN_TOT)
      DOUBLE PRECISION, INTENT(INOUT)  :: TAB1(:)
      TYPE (BIEF_OBJ), INTENT(INOUT) :: OBJ
      INTEGER, INTENT(IN) :: NPOIN_TOT
      INTEGER, INTENT(IN) :: NPOIN
      INTEGER :: I,J
      OBJ%R=0.0
      DO I=1,NPOIN_TOT
         DO J=1,NPOIN
            IF (I .EQ. MESH%KNOLG%I(J)) THEN
               OBJ%R(J)=TAB1(I)
            END IF
         END DO
      END DO
      OBJ%DIM1=NPOIN
   
      


      END SUBROUTINE













      
      END
