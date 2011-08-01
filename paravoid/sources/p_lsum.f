C                           *****************
                            SUBROUTINE P_LSUM
C                           *****************
C
     *(IARG1,LARG2)
C
C***********************************************************************
C  PARAVOID    VERSION 5.8         01/07/2006      O.BOITEAU (SINETICS)
C***********************************************************************
C
C      FONCTIONS: REDUCTION D'UN VECTEUR DE LOGICAL AVEC DIFFUSION DU
C      ==========   RESULTAT SUR TOUS LES PROCESSEURS.
C                 VERSION BIDON
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________|
C |  LARG2         |<-->| CONTRIBUTION DU PROCESSEUR APPELANT.         |
C |  IARG1         | -->| TAILLE DU VECTEUR                            |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : READ_CONNECTIVITY (ESTEL3D)
C
C SOUS-PROGRAMMES APPELES : MPI_ALLREDUCE + MPI_LOR (OU LOGIQUE)
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: IARG1
      LOGICAL, DIMENSION(IARG1), INTENT(OUT) :: LARG2
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      IF(LNG.EQ.1) WRITE(LU,*) 'APPEL DE P_LSUM VERSION VIDE'
      IF(LNG.EQ.2) WRITE(LU,*) 'CALL OF P_LSUM IN ITS VOID VERSION'
C
C     JMH LE 14/01/2008, POUR EVITER UN WARNING DU AU INTENT(OUT)
      LARG2 = .FALSE.
C
      RETURN
      END
