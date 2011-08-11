!                    *****************
                     SUBROUTINE PERALE
!                    *****************
!
     &(PALE,GAMMA,PERPIC,NPALE,TRA01,NPOIN,PRIVE,NPRIV,PMIN,PMAX)
!
!***********************************************************************
! ARTEMIS   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    DISCRETISES AN ENERGY SPECTRUM IN NPALE BANDS
!+                OF EQUAL ENERGY. THE RESULT IS A LIST OF
!+                PERIODS CORRESPONDING TO EACH BAND.
!
!history  J-M HERVOUET (LNH)
!+
!+
!+   LINKED TO BIEF 5.0
!
!history  D. AELBRECHT (LNH)
!+        02/06/1999
!+        V5P1
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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| GAMMA          |-->| COEFFICIENT IN THE SPECTRUM FORMULA
!| NPALE          |-->| NUMBER OF DISCRETISATION PERIOD
!| NPOIN          |-->| NUMBER OF POINTS
!| NPRIV          |<->| NUMBER OF PRIVATE TABLES
!| PALE           |<--| PERIODS FOR SPECTRUM DISCRETISATION
!| PERPIC         |-->| PEAK PERIOD FOR SPECTRUM
!| PMAX           |-->| MAXIMUM FREQUENCY FOR SPECTRUM
!| PMIN           |-->| MINIMUM FREQUENCY FOR SPECTRUM
!| PRIVE          |<->| USER PRIVATE TABLE
!| TRA01          |<->| WORK STRUCTURE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE INTERFACE_ARTEMIS, EX_PERALE=> PERALE
      USE BIEF
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
      INTEGER NPALE,NPOIN,NPRIV,NPAS,I,K
!
      DOUBLE PRECISION PALE(NPALE),TRA01(NPOIN)
      DOUBLE PRECISION PERPIC,GAMMA ,SUMB  ,SUMICI   ,DF    ,VAR
      DOUBLE PRECISION PMIN  ,PMAX  ,FMIN  ,FMAX
      DOUBLE PRECISION B     ,B1    ,B2
!
      TYPE(BIEF_OBJ) :: PRIVE
!
      DOUBLE PRECISION FP,GAM,DELTA
      COMMON /COEFHE/ FP,GAM,DELTA
!
!      DOUBLE PRECISION SPE
!      EXTERNAL SPE
!
      INTRINSIC LOG,FLOAT
!
!-----------------------------------------------------------------------
!
! PEAK FREQUENCY
      FP   = 1.D0 / PERPIC
      FMIN = 1.D0 / PMAX
      FMAX = 1.D0 / PMIN
      IF (FMAX.GE.99.D0) THEN
         FMAX = 2.5D0 * FP
      ENDIF
!
! GAMMA IS IN THE COMMON STATEMENT OF FUNCTION SPE (CANNOT BE
! CALLED GAMMA BECAUSE IT IS AN ARGUMENT OF THIS SUBROUTINE)
      GAM  = GAMMA
!
!-----------------------------------------------------------------------
!
      IF (GAMMA.GT.0.99D0 .AND. GAMMA.LT.1.01D0) THEN
!
!
!        PIERSON-MOSKOWITZ SPECTRUM
!        ----------------------------
!
         B1 = EXP(-1.25D0 * (FP/FMAX)**4)
         B2 = EXP(-1.25D0 * (FP/FMIN)**4)
         B  = B1 - B2
         DO 10 I=1,NPALE
            PALE(NPALE-I+1) = PERPIC *
     &    (-0.8D0*LOG( B2 + B*FLOAT(2*I-1)/FLOAT(2*NPALE) ))**(0.25D0)
10       CONTINUE
!
      ELSE
!
!
!        JONSWAP SPECTRUM
!        ------------------
!
!        THE FREQUENCIES LIMITING THE SPECTRUM TO THE LEFT AND RIGHT
!        ARE GIVEN BY KEYWORDS IN THE ARGUMENTS
!
      IF (FMAX.LE.FP) THEN
         FMAX = 2.5D0 * FP
         WRITE(LU,110) FMAX
 110     FORMAT(/,1X,'(PERALE) : FMAX < FP ??? =>',1X,
     &          'CORRECTION : FMAX =',1X,F5.3,' HZ',/)
      ENDIF
!
!        NUMBER OF INTEGRATION INTERVALS FOR THE TRAPEZOIDS METHOD
!
         NPAS = 2000*NPALE
!
!        WIDTH OF AN INTEGRATION INTERVAL
!
         DF = (FMAX-FMIN)/FLOAT(NPAS)
!
!        COEFFICIENT FOR THE FUNCTION OF THE SPECTRUM (COMPUTED HERE
!        SO THAT IT'S NOT RECOMPUTED WHEN CALLS SPE)
!
         DELTA = 0.0624D0 * FP**4 /
     &           ( 0.230D0 + 0.0336D0*GAM - 0.185D0 / (1.9D0+GAM) )
!
!        SURFACE OF THE SPECTRUM (TRAPEZOIDS METHOD)
!
         SUMB = (SPE(FMIN) + SPE(FMAX))/2.D0
         DO 20 I = 2,NPAS-1
            SUMB = SUMB + SPE(FMIN+FLOAT(I)*DF)
20       CONTINUE
!
!        DIVIDES THE SPECTRUM INTO 2*NPALE BANDS OF EQUAL ENERGY
!
         SUMB = SUMB/FLOAT(2*NPALE)
!
!        IDENTIFIES THE FREQUENCIES EVERY (2*I-1)*SUMB (I=1,NPALE)
!
         SUMICI = SPE(FMIN)/2.D0
         I   = 1
         DO 30 K=1,NPAS
            VAR = SPE(FMIN+DF*FLOAT(K))
            SUMICI = SUMICI + VAR/2.D0
            IF (SUMICI.GT.SUMB*FLOAT(2*I-1)) THEN
               PALE(NPALE-I+1) = 1.D0 / ( FMIN + DF*(FLOAT(K)-0.5D0) )
               I = I + 1
               IF (I.GT.NPALE) RETURN
            ENDIF
            SUMICI = SUMICI + VAR/2.D0
30       CONTINUE
!
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
