C                       **********************
                        SUBROUTINE CONFIG_CODE
C                       **********************
C
     &(ICODE)
C
C***********************************************************************
C BIEF VERSION 5.5               J-M HERVOUET (LNH) 01 30 87 80 18
C***********************************************************************
C
C   FONCTIONS: 1) RESETS CORRESPONDING LOGICAL UNITS AND FILE NAMES 
C   ========== WHEN THERE ARE SEVERAL PROGRAMS COUPLED
C
C              2) RESETS THE COMMON NDS CORRESPONDING TO THE CODE NUMBER
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    ICODE       |    | NUMERO DU CODE EN CAS DE COUPLAGE
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : HOMERE
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      USE BIEF, EX_CONFIG_CODE => CONFIG_CODE
      USE DECLARATIONS_TELEMAC
C
      IMPLICIT NONE
      INTEGER     LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER , INTENT(IN)    :: ICODE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER  NDS(0:81,7)
      INTEGER NNDS(0:81,7,3)
C
      INTEGER I,J
C
      COMMON/NODES/NDS
      COMMON/NNODES/NNDS
C
C-----------------------------------------------------------------------
C
C  2) RE-SETTING COMMON NDS          
C
      DO I=0,81
        DO J=1,7
          NDS(I,J)=NNDS(I,J,ICODE)
        ENDDO
      ENDDO
C
C-----------------------------------------------------------------------
C
C  3) SETTING NAME OF CURRENT CODE          
C
      NAMECODE = NNAMECODE(ICODE)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
