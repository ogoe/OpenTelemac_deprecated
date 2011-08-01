C                       *************************
                        SUBROUTINE P_MPI_ADDRESS3
C                       *************************
C
     *(LOCATION,ADDRESS,IER)
C
C***********************************************************************
C  PARA VERSION 5.9    19/08/2008   J.-M. HERVOUET (LNHE) 01 30 87 80 18
C  
C***********************************************************************
C
C      FONCTIONS: APPEL DE LA FONCTION MPI_ADDRESS (ICI PREMIER ARGUMENT 
C      ==========                              TABLEAU DOUBLE PRECISION)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________|
C |                | -->|
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR :
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      DOUBLE PRECISION LOCATION(10)
      INTEGER ADDRESS,IER
C
C-----------------------------------------------------------------------
C
      CALL MPI_ADDRESS(LOCATION,ADDRESS,IER)
C
      IF(IER.NE.0) THEN
        WRITE(LU,*) 'P_MPI_ADDRESS3:'
        WRITE(LU,*) 'MPI ERROR ',IER
        STOP
      ENDIF
C
C----------------------------------------------------------------------
C
      RETURN
      END
