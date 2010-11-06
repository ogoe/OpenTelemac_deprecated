C                       ************************
                        SUBROUTINE P_MPI_ADDRESS
C                       ************************
C
     *(LOCATION,ADDRESS,IER)
C
C***********************************************************************
C  PARA VERSION 5.9    23/06/2008   J.-M. HERVOUET (LNHE) 01 30 87 80 18
C  
C***********************************************************************
C
C      FONCTIONS: APPEL DE LA FONCTION MPI_ADDRESS
C      ==========
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
      INTEGER LOCATION,ADDRESS,IER
C
C-----------------------------------------------------------------------
C
      IF(LNG.EQ.1) WRITE(LU,*) 'APPEL DE P_MPI_ADDRESS VERSION VIDE'
      IF(LNG.EQ.2) WRITE(LU,*) 'CALL OF P_MPI_ADDRESS VOID VERSION'
C
C----------------------------------------------------------------------
C
      RETURN
      END
