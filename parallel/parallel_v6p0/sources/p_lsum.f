C                           *****************
                            SUBROUTINE P_LSUM
C                           *****************
C
     *(IARG1,LARG2)
C
C***********************************************************************
C  PARALLEL    VERSION 5.9         01/07/2006      O.BOITEAU (SINETICS)
C***********************************************************************
C
C      FONCTIONS: REDUCTION D'UN VECTEUR DE LOGICAL AVEC DIFFUSION DU
C      ==========   RESULTAT SUR TOUS LES PROCESSEURS.
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
      INCLUDE 'mpif.h'
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: IARG1
      LOGICAL, DIMENSION(IARG1), INTENT(INOUT) :: LARG2
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      LOGICAL, DIMENSION(:), ALLOCATABLE :: LAUX
      INTEGER IER,I
C
C-----------------------------------------------------------------------
C
      ALLOCATE(LAUX(IARG1),STAT=IER)
      IF (IER.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*)'P_LSUM: ERREUR DANS ALLOCATION MEMOIRE'
        IF(LNG.EQ.2) WRITE(LU,*)'P_LSUM: ERROR IN MEMORY ALLOCATION'
        CALL PLANTE(1)
        STOP 
      ENDIF
C
      DO I=1,IARG1
        LAUX(I)=LARG2(I)
      ENDDO
C
      CALL MPI_ALLREDUCE(LAUX,LARG2,IARG1,MPI_LOGICAL,            
     &                   MPI_LOR,MPI_COMM_WORLD,IER)
C
      IF(IER.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'P_LSUM: ERREUR DANS MPI_ALLREDUCE'
        IF(LNG.EQ.2) WRITE(LU,*) 'P_LSUM: ERROR IN MPI_ALLREDUCE'
        WRITE(LU,*) 'MPI ERROR: ',IER
        CALL PLANTE(1)
        STOP 
      ENDIF
C
      DEALLOCATE(LAUX)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
