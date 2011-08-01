C                       *****************
                        SUBROUTINE P_MAIL
C                       *****************
C
     *(CHAINE,NCAR)
C
C***********************************************************************
C  PARA       VERSION 5.9         23/06/2008        J-M HERVOUET (LNHE)
C  ADAPTED FOR MPI              /10/99       RAINER JOHANNI (SGI MUNICH)
C  VERSION 5.0 MODIFIED       28/12/99    J.A. JANKOWSKI (BAW KARLSRUHE)
C***********************************************************************
C
C      FONCTIONS: PASSAGE D'UNE CHAINE DE CARACTERES DE LONGUEUR NCAR.
C      ========== DEPUIS LA STATION MAITRESSE A TOUTES LES AUTRES.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________|
C |  MYPART        | -->| CONTRIBUTION DU PROCESSEUR APPELANT.
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
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INCLUDE 'mpif.h'
C
      INTEGER NCAR
      INTEGER IER
C
      CHARACTER*250 CHAINE
C
C-----------------------------------------------------------------------
C
      CALL MPI_BCAST(CHAINE,NCAR,MPI_CHARACTER,0,MPI_COMM_WORLD,IER)
C
      IF (IER.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'P_MAIL: PROBLEME DANS MPI_BCAST'
        IF(LNG.EQ.2) WRITE(LU,*) 'P_MAIL: PROBLEM IN MPI_BCAST'
        WRITE (LU,*) 'MPI ERROR ',IER
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END            
