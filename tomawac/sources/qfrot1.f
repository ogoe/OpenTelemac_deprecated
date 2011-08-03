!                    *****************
                     SUBROUTINE QFROT1
!                    *****************
!
     &( TSTOT , TSDER , F     , XK    , DEPTH , CFROT1, GRAVIT, NF    ,
     &  NPLAN , NPOIN2, BETA  )
!
!***********************************************************************
! TOMAWAC   V6P1                                   23/06/2011
!***********************************************************************
!
!brief    COMPUTES THE CONTRIBUTION OF THE BOTTOM FRICTION
!+                SOURCE TERM BASED ON HASSELMANN ET AL.'S FORMULATION
!+                (1973), MODIFIED BY BOUWS ET KOMEN (1983).
!
!note     THIS SOURCE TERM IS LINEAR IN F(FREQ,TETA), AND THE LINEAR
!+          COEFFICIENT DOES NOT VARY WITH TIME.
!note   CFROT1 (USED IN WAM CYCLE 4) EQUALS 0.038 M2.S-3.
!
!reference  HASSELMANN ET AL. (1973) :
!+                     "MEASUREMENTS OF WIND-WAVE GROWTH AND SWELL
!+                      DECAY DURING THE JOINT NORTH SEA WAVE PROJECT
!+                     (JONSWAP)". DEUTSCHEN HYDROGRAPHISVHEN ZEITSCHRIFT, REIHE A(8), NUM 12.
!reference BOUWS E., KOMEN G.J. (1983) :
!+                     "ON THE BALANCE BETWEEN GROWTH AND DISSIPATION
!+                      IN AN EXTREME DEPTH-LIMITED WIND-SEA IN THE
!+                      SOUTHERN NORTH-SEA". JPO, VOL 13.
!
!history  P. THELLIER; M. BENOIT (EDF/DER/LNH)
!+        03/04/95
!+        V1P0
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
!+        23/06/2011
!+        V6P1
!+   Translation of French names of the variables in argument
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| BETA           |<--| WORK TABLE
!| CFROT1         |-->| BOTTOM FRICTION COEFFICIENT
!| DEPTH          |-->| WATER DEPTH
!| F              |-->| DIRECTIONAL SPECTRUM
!| GRAVIT         |-->| ACCELERATION DE LA PESANTEUR
!| NF             |-->| NOMBRE DE FREQUENCES DE DISCRETISATION
!| NPLAN          |-->| NOMBRE DE DIRECTIONS DE DISCRETISATION
!| NPOIN2         |-->| NOMBRE DE POINTS DU MAILLAGE SPATIAL
!| TSDER          |<->| DERIVED PART OF THE SOURCE TERM CONTRIBUTION
!| TSTOT          |<->| TOTAL PART OF THE SOURCE TERM CONTRIBUTION
!| XK             |-->| DISCRETIZED WAVE NUMBER
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IMPLICIT NONE
!
!.....VARIABLES IN ARGUMENT
!     """"""""""""""""""""
      INTEGER  NF    , NPLAN , NPOIN2
      DOUBLE PRECISION CFROT1, GRAVIT
      DOUBLE PRECISION DEPTH(NPOIN2), BETA(NPOIN2)
      DOUBLE PRECISION TSTOT(NPOIN2,NPLAN,NF), TSDER(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION     F(NPOIN2,NPLAN,NF),    XK(NPOIN2,NF)
!
!.....LOCAL VARIABLES
!     """""""""""""""""
      INTEGER  JP    , JF    , IP
      DOUBLE PRECISION COEF , DEUKD
!
!
      COEF=-2.D0*CFROT1/GRAVIT
!
!.....LOOP OVER DISCRETISED FREQUENCIES
!     """"""""""""""""""""""""""""""""""""""""""""
      DO JF=1,NF
!
!.......COMPUTES THE LINEAR COEFFICIENT BETA : QFROT1 = BETA * F
!       """""""""""""""""""""""""""""""""""""""""""""""""""""""
        DO IP=1,NPOIN2
          DEUKD = MIN(2.D0*DEPTH(IP)*XK(IP,JF),7.D2)
          BETA(IP) = COEF*XK(IP,JF)/SINH(DEUKD)
        ENDDO
!
!.......TAKES THE SOURCE TERM INTO ACCOUNT
!       """"""""""""""""""""""""""""""""
        DO JP=1,NPLAN
          DO IP=1,NPOIN2
            TSTOT(IP,JP,JF) = TSTOT(IP,JP,JF)+BETA(IP)*F(IP,JP,JF)
            TSDER(IP,JP,JF) = TSDER(IP,JP,JF)+BETA(IP)
          ENDDO
        ENDDO
      ENDDO
!
      RETURN
      END
