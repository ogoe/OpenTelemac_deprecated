C                       *******************
                        SUBROUTINE SAVE_NDS
C                       *******************
C
     *(ICODE)
C
C***********************************************************************
C BIEF VERSION 5.5          05/08/04   J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C
C FONCTION : 1) SAVING COMMON NDS (INITIALISED BY A CALL TO ININDS)
C               TO BE ABLE TO RETRIEVE IT BY A CALL CONFIG_CODE
C
C               THIS COMMON IS USED BY FUNCTIONS IN BIEF SUCH AS NBPTS
C
C               ININDS IS CALLED BY ALMESH
C               ALMESH IS CALLED BY POINT_NAMEOFCODE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  ICODE         | -->| CODE NUMBER
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER , INTENT(IN) :: ICODE
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
      DO I=0,81
        DO J=1,7
          NNDS(I,J,ICODE)=NDS(I,J)
        ENDDO
      ENDDO
C
C-----------------------------------------------------------------------
C      
      RETURN
      END
