C                       *****************
                        SUBROUTINE ERRPVM
C                       *****************
C
     *(ERROR_NUMBER)
C
C***********************************************************************
C  PARA       VERSION 5.9        23/06/2008        J-M HERVOUET (LNHE)
c  adapted for MPI              /10/99       Rainer Johanni (SGI Munich)
c  version 5.0 modified       28/12/99    J.A. Jankowski (BAW Karlsruhe)
C***********************************************************************
C
C      FONCTIONS: IMPRESSION DES MESSAGES D'ERREUR DE PVM.
C      ==========
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________|
C |  ERROR_NUMBER  | -->| RETOUR D'UN APPEL A MPI
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
      INTEGER ERROR_NUMBER
C
C-----------------------------------------------------------------------
C
      WRITE(LU,*) 'MPI ERROR NUMBER: ',ERROR_NUMBER
C
      RETURN
      END
