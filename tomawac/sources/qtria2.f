!                    *****************
                     SUBROUTINE QTRIA2
!                    *****************
!
     &( F     , XK    , FREQ  , DFREQ , DEPTH , TETA  , SINTET, COSTET ,
     &  KSPB  , BDISPB, BDSSPB, RAISF , GRAVIT, NF    , NPLAN , NPOIN2 ,
     &  NBD   , INDI  , TSTOT , TSDER )
!
!***********************************************************************
! TOMAWAC   V6P1                                   27/06/2011
!***********************************************************************
!
!brief    COMPUTES THE CONTRIBUTION OF THE NON-LINEAR
!+                INTERACTIONS SOURCE TERM (FREQUENCY TRIADS).
!+<BR>           (INSPIRED FROM THE BOUSSINESQ EQUATIONS)
!
!history  EDF/DER/LNH
!+        11/06/98
!+        V5P0
!+
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!history  G.MATTAROLO (EDF - LNHE)
!+        27/06/2011
!+        V6P1
!+   Translation of French names of the variables in argument
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| BDISPB         |-->| LOWER DIRECTIONAL BOUND. OF SPB TRIAD MODEL
!| BDSSPB         |-->| UPPER DIRECTIONAL BOUND. OF SPB TRIAD MODEL
!| COSTET         |-->| CPSINE OF TETA ANGLE
!| DEPTH          |-->| WATER DEPTH
!| DFREQ          |-->| FREQUENCY STEPS BETWEEN DISCRETIZED FREQUENCIES
!| F              |-->| DIRECTIONAL SPECTRUM
!| FREQ           |-->| DISCRETIZED FREQUENCIES
!| GRAVIT         |-->| GRAVITY ACCELERATION
!| INDI           |-->| WORK TABLE
!| KSPB           |-->| COEFFICIENT K OF SPB TRIAD INTERACTION MODEL
!| NBD            |-->| NUMBER OF TRIAD CONFIGURATIONS
!| NF             |-->| NUMBER OF FREQUENCIES
!| NPLAN          |-->| NUMBER OF DIRECTIONS
!| NPOIN2         |-->| NUMBER OF POINTS IN 2D MESH
!| RAISF          |-->| RAISON FREQUENTIELLE POUR DISCRETISATION
!| SINTET         |-->| SINE OF TETA ANGLE
!| TETA           |-->| DISCRETIZED DIRECTIONS
!| TSDER          |<->| DERIVED PART OF THE SOURCE TERM CONTRIBUTION
!| TSTOT          |<->| TOTAL PART OF THE SOURCE TERM CONTRIBUTION
!| XK             |-->| DISCRETIZED WAVE NUMBER
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IMPLICIT NONE
!
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
!
!.....VARIABLES IN ARGUMENT
!     """"""""""""""""""""
      INTEGER  NF, NPOIN2, NPLAN
      DOUBLE PRECISION  F(NPOIN2,NPLAN,NF), XK(NPOIN2,NF)
      DOUBLE PRECISION  RAISF , DFREQ(NF) , FREQ(NF), DEPTH(NPOIN2)
      DOUBLE PRECISION  TETA(NPLAN) , SINTET(NPLAN) , COSTET(NPLAN)
      DOUBLE PRECISION  TSTOT(NPOIN2,NPLAN,NF), TSDER(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION  GRAVIT , BMS , BMSP
!
!.....LOCAL VARIABLES
!     """""""""""""""""
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
!
      INTEGER           NBD, INDI(NBD), IP1, IP3
!
!
!.....EXTERNAL FUNCTIONS
!     """""""""""""""""""
      DOUBLE PRECISION KERBOU
      EXTERNAL         KERBOU
!
!     """"""""""""""""""""""""""""""""""""""""""""""""""""""""""
!
      DTETA = TETA(2)-TETA(1)
      BMS  = 1.D0/15.D0
      BMSP = BMS + 1.D0/3.D0
      DEUPI2 = DEUPI*DEUPI
      FREQ0  = FREQ(1)
      LRAISF = DLOG(RAISF)
      RAISM1 = RAISF-1.D0
!
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
!         COMPUTES K2
!         --------------------
          DEP = DEPTH(IPO)
          XK1 = XK(IPO,JFF)
          XK3 = XK(IPO,IFF)
          K2NL   = DSQRT((XK3*COSTET(IPL)-XK1*COSTET(JPL))**2
     &             +(XK3*SINTET(IPL)-XK1*SINTET(JPL))**2)
          XK2    = (1.D0-FR1)*XK(IPO,IFR) + FR1*XK(IPO,IFR+1)
!
          TETA2=DATAN2(XK3*SINTET(IPL)-XK1*SINTET(JPL)
     &                ,XK3*COSTET(IPL)-XK1*COSTET(JPL))
          IF(TETA2.LT.0.D0) TETA2 = DEUPI + TETA2
!
          IF(TETA2.LT.BDISPB .OR. TETA2.GT.BDSSPB) THEN
!         INTERACTIONS BETWEEN COMPONENTS WHICH DIRECTIONS ARE NOT
!         WITHIN THE ANGULAR SECTOR DEFINED BY THE USER (VARIABLES
!         BDISPB AND BDSSPB) ARE NOT TAKEN INTO ACCOUNT
             GOTO 205
          END IF
!
          AP2    = (TETA2-TETA(1))/DTETA
          IPM    = IDINT(AP2)
          AP2    = AP2 - DBLE(IPM)
          IPM    = IPM + 1
          IPP    = IPM + 1
          IF(IPP.EQ.NPLAN+1) IPP = 1
!
!.........COMPUTES COUPLING COEFFICIENTS
!         """"""""""""""""""""""""""""""""""""
!         R(P-M,M)
          VR1 = KERBOU
     &(   XK1 , XK2   , FREQ1 , FREQ2 , DEP , TETA(JPL) , TETA2     )
!         R(M-P,P)
          VR2 = KERBOU
     &(   -XK1, XK3   , -FREQ1, FREQ3 , DEP , TETA(JPL) , TETA(IPL) )
!         R(-M,P)
          VR3 = KERBOU
     &(   -XK2  , XK3 , -FREQ2, FREQ3 , DEP , TETA2     , TETA(IPL) )
!
          FILT = KSPB/((XK2-K2NL)**2+KSPB*KSPB)
          FILT = -0.5D0*FILT/XK2
!
          DEP2 = DEP*DEP
          VAR1 = 2.D0*BMS*DEP2
          XC1  = VAR1*XK3*XK3
          XC2  = VAR1*XK2*XK2
          XC3  = VAR1*XK1*XK1
          VAR1 = BMSP*DEUPI2*DEP2
          TK1 = (GRAVIT*DEP*(1.D0+XC1)-VAR1*FREQ3*FREQ3)
          TK2 = (GRAVIT*DEP*(1.D0+XC2)-VAR1*FREQ2*FREQ2)
          TK3 = (GRAVIT*DEP*(1.D0+XC3)-VAR1*FREQ1*FREQ1)
!
          BK1 = DEUPI*FREQ3*(1.D0+3.D0*XC1)
          BK2 = DEUPI*FREQ2*(1.D0+3.D0*XC2)
          BK3 = DEUPI*FREQ1*(1.D0+3.D0*XC3)
!
!
!.........TAKES THE SOURCE TERM INTO ACCOUNT
!         """"""""""""""""""""""""""""""""
!
          NRJ2  = (1.D0-AP2)*((1.D0-FR1)*F(IPO,IPM,IFR)+FR1*
     &            F(IPO,IPM,IFR+1)) + AP2*((1.D0-FR1)*F(IPO,IPP,IFR)
     &            +FR1*F(IPO,IPP,IFR+1))
!
          BISP  = FILT*
     &            ((VR2/TK2)*F(IPO,IPL,IFF)*F(IPO,JPL,JFF)
     &            +(VR3/TK3)*F(IPO,IPL,IFF)*NRJ2
     &            -(VR1/TK1)*F(IPO,JPL,JFF)*NRJ2)
!
!          D12   = FILT*((VR2/TK2)*F(IPO,JPL,JFF)
!     *                 +(VR3/TK3)*NRJ2)
!          D21   = FILT*((VR2/TK2)*F(IPO,IPL,IFF)
!     *                 -(VR1/TK1)*NRJ2)
          VR1   = DFREQ(JFF)*DTETA*VR1/BK1
          VR3   = 2.D0*DFREQ(IFF)*DTETA*VR3/BK3
!
          TSTOT(IPO,IPL,IFF) = TSTOT(IPO,IPL,IFF) + VR1*BISP
!          TSDER(IPO,IPL,IFF) = TSDER(IPO,IPL,IFF) + VR1*D12
          TSTOT(IPO,JPL,JFF) = TSTOT(IPO,JPL,JFF) - VR3*BISP
!          TSDER(IPO,JPL,JFF) = TSDER(IPO,JPL,JFF) - VR3*D21
 205      CONTINUE
 203     CONTINUE
 202    CONTINUE
 201   CONTINUE
 200  CONTINUE
!
      RETURN
      END
