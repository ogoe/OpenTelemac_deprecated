C                       ***********************
                        SUBROUTINE ATTEND(ISEC)
C                       ***********************
C
C***********************************************************************
C  BIEF VERSION 5.2    08/02/2001    NATHALY BARBRY (UNIVERSITE DE CAEN)
C                                        J-M HERVOUET (LNHE) 30 87 80 18
C***********************************************************************
C
C  FONCTION :  ATTEND PENDANT UNE DUREE ISEC EN SECONDES
C
C----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |                |    |                                              |
C |   D            | -->| DUREE                                        |
C |________________|____|______________________________________________|
C 
C-----------------------------------------------------------------------
C
C APPELE PAR : WAITFOR
C
C**********************************************************************
C
C     PLEASE ADD OR REMOVE COMMENTS ACCORDING TO YOUR COMPILER
C
C-----------------------------------------------------------------------
C     VERSION F90 STANDARD
C-----------------------------------------------------------------------
C 
      INTEGER, INTENT(IN) :: ISEC     
      INTEGER T1,T2
      INTEGER  TIME_IN_SECONDS
      EXTERNAL TIME_IN_SECONDS
      T1 = TIME_IN_SECONDS()
      T2 = T1
C                            .AND. : WHEN CLOCK RESET TO ZERO
      DO WHILE (T2.LT.T1+ISEC.AND.T2.GE.T1)
        T2 = TIME_IN_SECONDS()
      END DO
C
C-----------------------------------------------------------------------
C     VERSION F95 NAG
C-----------------------------------------------------------------------
C
C     USE F90_UNIX_PROC
C     INTEGER, INTENT(IN) :: ISEC
C     CALL SLEEP(ISEC)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
