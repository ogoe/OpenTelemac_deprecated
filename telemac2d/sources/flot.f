C                       ***************
                        SUBROUTINE FLOT
C                       ***************
C
     *(XFLOT,YFLOT,NFLOT,NITFLO,FLOPRD,X,Y,NPOIN,DEBFLO,FINFLO,NIT)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2    17/08/94    J-M JANIN    (LNH) 30 87 72 84
C
C***********************************************************************
C
C   FONCTION : L'UTILISATEUR DOIT DONNER ICI :
C
C   1) LE PAS DE TEMPS DE LARGAGE DE CHAQUE FLOTTEUR
C
C   2) LE PAS DE TEMPS DE FIN DE CALCUL DE DERIVE DE CHAQUE FLOTTEUR
C
C   3) LA POSITION DES FLOTTEURS AU MOMENT DU LARGAGE
C
C-----------------------------------------------------------------------
C
C   FUNCTION : THE USER MUST GIVE HERE
C
C   1) WHEN THE FLOATING BODY IS RELEASED (IN TERMS OF TIME STEP)
C
C   2) WHEN THE COMPUTATION IS STOPPED FOR THIS FLOATING BODY.
C
C   3) THE INITIAL POSITION OF THE FLOATING BODY
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    XFLOT,YFLOT |<-- | POSITIONS OF FLOATING BODIES
C |    NFLOT       | -->| NOMBRE DE FLOTTEURS.
C |    NITFLO      | -->| MAXIMUM NUMBER OF RECORDS OF SUCCESIVE
C |                |    | POSITIONS OF FLOATING BODIES.
C |    FLOPRD      | -->| NUMBER OF TIME-STEPS BETWEEN 2 RECORDS OF
C |                |    | SUCCESSIVE POSITIONS OF FLOATING BODIES.
C |    X,Y         | -->| COORDINATES OF POINTS IN THE MESH
C |    NPOIN       | -->| NUMBER OF POINTS IN THE MESH
C |    DEBFLO      |<-- | TIME STEP OF INITIAL RELEASE
C |    FINFLO      |<-- | TIME STEP FOR END OF FOLLOW UP
C |    NIT         | -->| NUMBEr OF TIME STEPS
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
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
      INTEGER, INTENT(IN)             :: NPOIN,NIT,NFLOT,NITFLO,FLOPRD
      INTEGER, INTENT(INOUT)          :: DEBFLO(NFLOT),FINFLO(NFLOT)
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: XFLOT(NITFLO,NFLOT)
      DOUBLE PRECISION, INTENT(INOUT) :: YFLOT(NITFLO,NFLOT)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IFLOT    
C
C-----------------------------------------------------------------------
C
C  1) STEP FOR THE BEGINNING OF RELEASE (DEBFLO)
C  2) STEP FOR THE END OF RELEASE (FINFLO)
C
      DO 10 IFLOT=1,NFLOT
         DEBFLO(IFLOT) = 1
         FINFLO(IFLOT) = NIT
10    CONTINUE
C
C-----------------------------------------------------------------------
C
C  3) COORDINATES OF FLOATING BODIES AT THE BEGINNING
C
C     INITIAL POSITION OF FLOATING BODIES WHEN RELEASED.
C     DEFAULT VALUE NUMBER "IFLOT" IS RELEASED AT POINT "IFLOT"
C
C-----------------------------------------------------------------------
C
C     XFLOT(1,1)=-300000.D0
C     YFLOT(1,1)= 300000.D0
C
C     XFLOT(1,2)= 0.D0
C     YFLOT(1,2)= 300000.D0
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
