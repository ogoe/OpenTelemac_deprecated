C                       *****************
                        SUBROUTINE QTRIA2
C                       *****************
C
     *( F     , XK    , FREQ  , DFREQ , DEPTH , TETA  , SINTET, COSTET ,
     *  KSPB  , BDISPB, BDSSPB, RAISF , GRAVIT, NF    , NPLAN , NPOIN2 ,
     *  NBD   , INDI  , TSTOT , TSDER )
C
C**********************************************************************
C  TOMAWAC - V5P0                           (EDF/DER/LNH)  -  11/06/98
C**********************************************************************
C
C  FONCTION : TERME SOURCE LIE AUX INTERACTIONS NON-LINEAIRES ENTRE
C  ********** TRIPLETS DE FREQUENCES.
C             MODELE ISSU DES EQUATIONS DE BOUSSINESQ.
C
C  ARGUMENTS :
C  ***********
C  +-------------+----+--------------------------------------------+
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C  +-------------+----+--------------------------------------------+
C  ! F(-,-,-)    !<-->! SPECTRE DIRECTIONNEL DE VARIANCE           !
C  ! XK(-,-)     ! -->! TABLEAU DES NOMBRES D'ONDE                 !
C  ! FREQ(-)     ! -->! TABLEAU DES FREQUENCES DE DISCRETISATION   !
C  ! DFREQ(-)    ! -->! TABLEAU DES PAS DE FREQUENCE               !
C  ! DEPTH(-)    ! -->! TABLEAU DES PROFONDEURS (METRES)           !
C  ! TETA(-)     ! -->! VECTEUR DES DIRECTIONS DE DISCRETISATION   !
C  ! SINTET(-)   ! -->! VECTEUR DES   SINUS DES DIRECTIONS         !
C  ! COSTET(-)   ! -->! VECTEUR DES COSINUS DES DIRECTIONS         !
C  ! RAISF       ! -->! RAISON FREQUENTIELLE POUR DISCRETISATION   !
C  ! GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR               !
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !
C  ! TSTOT(-,-,-)!<-- ! CONTRIBUTION TERME SOURCE - PARTIE TOTALE  !
C  ! TSDER(-,-,-)!<-- ! CONTRIBUTION TERME SOURCE - PARTIE DERIVEE !
C  +-------------+----+--------------------------------------------+
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C  +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  SEMIMP
C  ********    - PROGRAMME(S) APPELE(S) :  KERBOU
C
C  REMARQUES :
C  ***********
C
C**********************************************************************
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NF, NPOIN2, NPLAN
      DOUBLE PRECISION  F(NPOIN2,NPLAN,NF), XK(NPOIN2,NF)
      DOUBLE PRECISION  RAISF , DFREQ(NF) , FREQ(NF), DEPTH(NPOIN2)
      DOUBLE PRECISION  TETA(NPLAN) , SINTET(NPLAN) , COSTET(NPLAN)
      DOUBLE PRECISION  TSTOT(NPOIN2,NPLAN,NF), TSDER(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION  GRAVIT , BMS , BMSP
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  IFF, JFF, IPL, JPL , IPO
      INTEGER  IFR, IPP, IPM
      DOUBLE PRECISION  DTETA, FR1 , AP2  , XK1   , XK3 , DEP
      DOUBLE PRECISION  TETA2, XK2 , K2NL , NRJ2
      DOUBLE PRECISION  FREQ0, FREQ1, FREQ2, FREQ3, LRAISF, RAISM1
      DOUBLE PRECISION  VR1 , VR2 , VR3 , TK1 , TK2 , TK3
      DOUBLE PRECISION  BK1 , BK2 , BK3 , D12 , D21 , DEP2
      DOUBLE PRECISION  FILT , BISP , BDISPB , BDSSPB
      DOUBLE PRECISION  PI  , DEUPI , DEUPI2 , KSPB
      DOUBLE PRECISION  VAR1 , XC1 , XC2 , XC3
      PARAMETER(PI=3.141592654D0, DEUPI=2.D0*PI )
C
      INTEGER           NBD, INDI(NBD), IP1, IP3
C
C
C.....FONCTIONS EXTERNES.
C     """""""""""""""""""
      DOUBLE PRECISION KERBOU
      EXTERNAL         KERBOU
C
C     """"""""""""""""""""""""""""""""""""""""""""""""""""""""""
C
      DTETA = TETA(2)-TETA(1)
      BMS  = 1.D0/15.D0
      BMSP = BMS + 1.D0/3.D0
      DEUPI2 = DEUPI*DEUPI
      FREQ0  = FREQ(1)
      LRAISF = DLOG(RAISF)
      RAISM1 = RAISF-1.D0
C
      DO 200 IFF = 1,NF
       FREQ3 = FREQ(IFF)
       DO 201 JFF = 1,IFF-1
        FREQ1 = FREQ(JFF)
        FREQ2 = FREQ3-FREQ1
        IF(FREQ2.LE.FREQ0) THEN
          GOTO 201
        END IF
        FR1    = 1.D0 + DLOG(FREQ2/FREQ0)/LRAISF
        IFR    = IDINT(FR1)
        FR1    = FR1 - DBLE(IFR)
        FR1    = (RAISF**FR1-1.D0)/RAISM1
        DO 202 IP3 = 1,NBD
         IPL = INDI(IP3)
         DO 203 IP1 = 1,NBD
          JPL = INDI(IP1)
          DO 205 IPO = 1,NPOIN2
C         CALCUL DU NOMBRE K2.
C         --------------------
          DEP = DEPTH(IPO)
          XK1 = XK(IPO,JFF)
          XK3 = XK(IPO,IFF)
          K2NL   = DSQRT((XK3*COSTET(IPL)-XK1*COSTET(JPL))**2
     *             +(XK3*SINTET(IPL)-XK1*SINTET(JPL))**2)
          XK2    = (1.D0-FR1)*XK(IPO,IFR) + FR1*XK(IPO,IFR+1)
C
          TETA2=DATAN2(XK3*SINTET(IPL)-XK1*SINTET(JPL)
     *                ,XK3*COSTET(IPL)-XK1*COSTET(JPL))
          IF(TETA2.LT.0.D0) TETA2 = DEUPI + TETA2
C
          IF(TETA2.LT.BDISPB .OR. TETA2.GT.BDSSPB) THEN
C         On ne prend pas en compte les interactions entre les composantes
C         dont les directions n'appartiennent pas au secteur angulaire
C         definit par l'utilisateur via les variables BDISPB et BDSSPB
             GOTO 205
          END IF
C
          AP2    = (TETA2-TETA(1))/DTETA
          IPM    = IDINT(AP2)
          AP2    = AP2 - DBLE(IPM)
          IPM    = IPM + 1
          IPP    = IPM + 1
          IF(IPP.EQ.NPLAN+1) IPP = 1
C
C.........CALCUL DES COEFFICIENTS DE COUPLAGE.
C         """"""""""""""""""""""""""""""""""""
C         R(p-m,m)
          VR1 = KERBOU
     *(   XK1 , XK2   , FREQ1 , FREQ2 , DEP , TETA(JPL) , TETA2     )
C         R(m-p,p)
          VR2 = KERBOU
     *(   -XK1, XK3   , -FREQ1, FREQ3 , DEP , TETA(JPL) , TETA(IPL) )
C         R(-m,p)
          VR3 = KERBOU
     *(   -XK2  , XK3 , -FREQ2, FREQ3 , DEP , TETA2     , TETA(IPL) )
C
          FILT = KSPB/((XK2-K2NL)**2+KSPB*KSPB)
          FILT = -0.5D0*FILT/XK2
C
          DEP2 = DEP*DEP
          VAR1 = 2.D0*BMS*DEP2
          XC1  = VAR1*XK3*XK3
          XC2  = VAR1*XK2*XK2
          XC3  = VAR1*XK1*XK1
          VAR1 = BMSP*DEUPI2*DEP2
          TK1 = (GRAVIT*DEP*(1.D0+XC1)-VAR1*FREQ3*FREQ3)
          TK2 = (GRAVIT*DEP*(1.D0+XC2)-VAR1*FREQ2*FREQ2)
          TK3 = (GRAVIT*DEP*(1.D0+XC3)-VAR1*FREQ1*FREQ1)
C
          BK1 = DEUPI*FREQ3*(1.D0+3.D0*XC1)
          BK2 = DEUPI*FREQ2*(1.D0+3.D0*XC2)
          BK3 = DEUPI*FREQ1*(1.D0+3.D0*XC3)
C
C
C.........PRISE EN COMPTE DU TERME SOURCE.
C         """"""""""""""""""""""""""""""""
C
          NRJ2  = (1.D0-AP2)*((1.D0-FR1)*F(IPO,IPM,IFR)+FR1*
     *            F(IPO,IPM,IFR+1)) + AP2*((1.D0-FR1)*F(IPO,IPP,IFR)
     *            +FR1*F(IPO,IPP,IFR+1))
C
          BISP  = FILT*
     *            ((VR2/TK2)*F(IPO,IPL,IFF)*F(IPO,JPL,JFF)
     *            +(VR3/TK3)*F(IPO,IPL,IFF)*NRJ2
     *            -(VR1/TK1)*F(IPO,JPL,JFF)*NRJ2)
C
C          D12   = FILT*((VR2/TK2)*F(IPO,JPL,JFF)
C     *                 +(VR3/TK3)*NRJ2)
C          D21   = FILT*((VR2/TK2)*F(IPO,IPL,IFF)
C     *                 -(VR1/TK1)*NRJ2)
          VR1   = DFREQ(JFF)*DTETA*VR1/BK1
          VR3   = 2.D0*DFREQ(IFF)*DTETA*VR3/BK3
C
          TSTOT(IPO,IPL,IFF) = TSTOT(IPO,IPL,IFF) + VR1*BISP
C          TSDER(IPO,IPL,IFF) = TSDER(IPO,IPL,IFF) + VR1*D12
          TSTOT(IPO,JPL,JFF) = TSTOT(IPO,JPL,JFF) - VR3*BISP
C          TSDER(IPO,JPL,JFF) = TSDER(IPO,JPL,JFF) - VR3*D21
 205      CONTINUE
 203     CONTINUE
 202    CONTINUE
 201   CONTINUE
 200  CONTINUE
C
      RETURN
      END
