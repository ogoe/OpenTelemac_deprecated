C                       *****************
                        SUBROUTINE SMTRAC
C                       *****************
C
     *(NPOIN,DIMT,AT,DT,SMTR,SMH,NREJET,ISCE,TSCE2,MAXSCE,MAXTRA,ITRAC)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.8           INRIA
C
C***********************************************************************
C
C     FONCTION  : CALCUL DU SECOND MEMBRE DU TRACEUR
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |  NPOIN         | -->|  NOMBRE DE POINTS DU MAILLAGE                |
C !  DIMT          ! -->!  DIMENSION DU TRACEUR                        |
C |  AT            | -->|  TEMPS                                       |
C |  DT            | -->|  PAS DE TEMPS HYDRO                          |
C |  SMTR          |<-->!  TERMES SOURCES DU TRACEUR                   !
C |  SMH           | -->|  TERMES SOURCES DE L'EQUATION DE CONTINUITE  |
C |  NREJET        | -->|  NOMBRE DE SOURCES/PUITS                     |
C |  ISCE          | -->|  NUMEROS GLOBAUX DES POINTS SOURCES          |
C |  TSCE          | -->|  VALEURS DU TRACEUR AUX SOURCES              |
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C     - SOUS PROGRAMME(S) APPELANT : RESOLU                             
C 
C***********************************************************************
C 
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NPOIN,NREJET,ISCE(*),DIMT,ITRAC
      INTEGER, INTENT(IN) :: MAXSCE,MAXTRA
      DOUBLE PRECISION, INTENT(IN)    :: AT,DT,SMH(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: TSCE2(MAXSCE,MAXTRA)
      DOUBLE PRECISION, INTENT(INOUT) :: SMTR(DIMT)
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,IS
C                                                          
C-----------------------------------------------------------------------
C
      IF(NREJET.NE.0) THEN
        DO I=1,NREJET
          IS =ISCE(I)
          SMTR(IS) = SMTR(IS) + DT*SMH(IS) * TSCE2(I,ITRAC)
        ENDDO
      ENDIF
C                                                          
C-----------------------------------------------------------------------
C                                                                        
      RETURN                                                            
      END
