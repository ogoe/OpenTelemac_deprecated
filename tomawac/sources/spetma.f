!                    *****************
                     SUBROUTINE SPETMA
!                    *****************
!
     &( SPEC  , FREQ  , NF    , AL    , FP     , GAMMA , SIGMAA, SIGMAB,
     &  DEUPI , GRAVIT, E2FMIN, FPMIN , DEPTH  )
!
!***********************************************************************
! TOMAWAC   V6P1                                   28/06/2011
!***********************************************************************
!
!brief    COMPUTES A TMA FREQUENCY SPECTRUM BASED
!+                ON A SERIES OF FREQUENCIES.
!
!history  M. BENOIT (EDF/DER/LNH)
!+        10/01/96
!+        V1P1
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
!+        28/06/2011
!+        V6P1
!+   Translation of French names of the variables in argument
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| AL             |-->| PHILLIPS CONSTANT (ALPHA)
!| DEPTH          |-->| WATER DEPTH
!| DEUPI          |-->| 2.PI
!| E2FMIN         |-->| SPECTRUM ENERGY THRESHOLD
!| FP             |-->| JONSWAP SPECTRUM PEAK FREQUENCY 
!| FPMIN          |-->| MINIMUM PEAK FREQUENCY VALUE
!| FREQ           |-->| DISCRETIZED FREQUENCIES
!| GAMMA          |-->| JONSWAP SPECTRUM PEAK FACTOR
!| GRAVIT         |-->| GRAVITY ACCELERATION
!| NF             |-->| NUMBER OF FREQUENCIES
!| SIGMAA         |-->| VALUE OF SIGMA FOR JONSWAP SPECTRUM (F<FP)
!| SIGMAB         |-->| VALUE OF SIGMA FOR JONSWAP SPECTRUM (F>FP)
!| SPEC           |<--| TMA VARIANCE DENSITY FREQUENCY SPECTRUM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IMPLICIT NONE
!
!.....VARIABLES IN ARGUMENT
!     """"""""""""""""""""
      INTEGER  NF
      DOUBLE PRECISION GRAVIT, SIGMAA, SIGMAB, GAMMA , DEUPI , FPMIN
      DOUBLE PRECISION FP    , E2FMIN, AL    , DEPTH
      DOUBLE PRECISION SPEC(NF)      , FREQ(NF)
!
!.....LOCAL VARIABLES
!     """""""""""""""""
      INTEGER  JF
      DOUBLE PRECISION COEF  , ARG1   , ARG2  , ARG3  , SIG   , FF
      DOUBLE PRECISION ARG4  , OMEGH
!
!
      IF (FP.GT.FPMIN) THEN
        COEF=AL*GRAVIT**2/DEUPI**4
        DO 100 JF=1,NF
          FF=FREQ(JF)
          IF (FF.LT.FP) THEN
            SIG=SIGMAA
          ELSE
            SIG=SIGMAB
          ENDIF
          ARG1=0.5D0*((FF-FP)/(SIG*FP))**2
          IF (ARG1.LT.99.D0) THEN
            ARG1=GAMMA**EXP(-ARG1)
          ELSE
            ARG1=1.D0
          ENDIF
          ARG2=1.25D0*(FP/FF)**4
          IF (ARG2.LT.99.D0) THEN
            ARG2=EXP(-ARG2)
          ELSE
            ARG2=0.D0
          ENDIF
          ARG3=COEF/FF**5
          OMEGH=DEUPI*FF*SQRT(DEPTH/GRAVIT)
          IF (OMEGH.LT.1.D0) THEN
            ARG4=0.5D0*OMEGH*OMEGH
          ELSEIF (OMEGH.LT.2.D0) THEN
            ARG4=1.D0-0.5D0*(2.D0-OMEGH)**2
          ELSE
            ARG4=1.D0
          ENDIF
          SPEC(JF)=ARG1*ARG2*ARG3*ARG4
          IF (SPEC(JF).LT.E2FMIN) SPEC(JF)=0.D0
  100   CONTINUE
      ELSE
        DO 150 JF=1,NF
          SPEC(JF)=0.D0
  150   CONTINUE
      ENDIF
!
      RETURN
      END
