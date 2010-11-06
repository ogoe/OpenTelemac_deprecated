C                       *****************
                        SUBROUTINE ECRETE
C                       *****************
     *( F     , DEPTH , NPOIN2, NPLAN , NF    , PROMIN)
C***********************************************************************
C TOMAWAC  V5.4      19/01/2004     M. BENOIT (EDF LNHE) 01 30 87 83 51
C***********************************************************************
C
C      FONCTION:
C      =========
C      MET LE SPECTRE DE VARIANCE A ZERO EN TOUS LES POINTS OU LA
C      PROFONDEUR D'EAU EST INFERIEURE A PROMIN.
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    F(-,-,-)    !<-->!  DENSITE SPECTRALE                           !
C !    DEPTH(-)    ! -->!  VECTEUR DES PROFONDEURS D'EAU               !
C !    NPOIN2      ! -->!  NOMBRE DE POINTS 2D                         !
C !    NPLAN       ! -->!  NOMBRE DE DIRECTIONS                        !
C !    NF          ! -->!  NOMBRE DE FREQUENCES                        !
C !    PROMIN      ! -->!  PROFONDEUR MINIMALE                         !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE PAR : WAC
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
      INTEGER          NPOIN2 , NPLAN, NF
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF), DEPTH(NPOIN2)
      DOUBLE PRECISION PROMIN
      INTEGER          IP    , JP    , JF
C
      DO IP=1,NPOIN2
        IF (DEPTH(IP).LT.PROMIN) THEN
          DO JF=1,NF
            DO JP=1,NPLAN
              F(IP,JP,JF)=0.0D0
            ENDDO
          ENDDO
        ENDIF
      ENDDO
C
      RETURN
      END
