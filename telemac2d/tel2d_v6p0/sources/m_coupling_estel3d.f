      MODULE M_COUPLING_ESTEL3D
C
C-----------------------------------------------------------------------
C This file is part of TELEMAC-2D v5p7
C-----------------------------------------------------------------------
C Written by JP Renaud

C Set of structures and subroutines that lets TELEMAC-2D and ESTEL-3D
C interact by reading/writing two arrays which live in this module:
C
C       - ESTEL-3D write flux values
C       - TELEMAC-2D read the flux values
C       - TELEMAC-2D writes the depth values
C       - ESTEL-3D read the depth values
C       - etc...
C
C As the arrays are private, they are accessed by public methods:
C DEPTH_FILL saves the depths from TELEMAC-2D and DEPTH_GET allows
C ESTEL-3D to recover the information.
C-----------------------------------------------------------------------
C
      IMPLICIT NONE

      PRIVATE

      PUBLIC :: INFILTRATION_INIT
      PUBLIC :: INFILTRATION_FINISH
      PUBLIC :: INFILTRATION_FILL
      PUBLIC :: INFILTRATION_GET
      PUBLIC :: DEPTH_FILL
      PUBLIC :: DEPTH_GET

      INTEGER :: NPOIN2D
      LOGICAL :: DOINFILTRATION

      DOUBLE PRECISION, ALLOCATABLE :: FLUX_FROM_ESTEL3D(:)
      DOUBLE PRECISION, ALLOCATABLE :: DEPTH_FROM_T2D(:)

C-----------------------------------------------------------------------
      CONTAINS
C-----------------------------------------------------------------------

      SUBROUTINE INFILTRATION_INIT(NPOIN,ACTIVATE)
C-----------------------------------------------------------------------
C Allocates the coupling arrays if required and fill them up with zeros.
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: NPOIN
      LOGICAL, INTENT(IN) :: ACTIVATE
C-----------------------------------------------------------------------
C
      LOGICAL DEJA
      DATA DEJA/.FALSE./
C
C-----------------------------------------------------------------------
C
      NPOIN2D = 1

      IF(ACTIVATE) THEN
        NPOIN2D        = NPOIN
        DOINFILTRATION = .TRUE.
      ENDIF

      IF(.NOT.DEJA) THEN
        ALLOCATE( FLUX_FROM_ESTEL3D( NPOIN2D ) )
        ALLOCATE( DEPTH_FROM_T2D( NPOIN2D ) )
        DEJA=.TRUE.
      ENDIF

      FLUX_FROM_ESTEL3D(:) = 0.D0
      DEPTH_FROM_T2D(:)    = 0.D0
C-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE INFILTRATION_INIT

      SUBROUTINE INFILTRATION_FINISH()
C-----------------------------------------------------------------------
C De-allocates the coupling arrays
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------

      DEALLOCATE( FLUX_FROM_ESTEL3D )
      DEALLOCATE( DEPTH_FROM_T2D )

C-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE INFILTRATION_FINISH

      SUBROUTINE INFILTRATION_FILL(ARRAY1,ARRAY2,COEFF)
C-----------------------------------------------------------------------
C Fill the array FLUX_FROM_ESTEL3D with the values from the arguments.
C
C This subroutine is called from ESTEL-3D.
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      DOUBLE PRECISION, INTENT(IN) :: ARRAY1(NPOIN2D)
      DOUBLE PRECISION, INTENT(IN) :: ARRAY2(NPOIN2D)
      DOUBLE PRECISION, INTENT(IN) :: COEFF
C-----------------------------------------------------------------------
              FLUX_FROM_ESTEL3D(:) = COEFF       * ARRAY1(:)
     &                             + (1 - COEFF) * ARRAY2(:)
C-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE INFILTRATION_FILL

      SUBROUTINE DEPTH_FILL(ARRAY_FROM_T2D)
C-----------------------------------------------------------------------
C Fill the array DEPTH_FROM_T2D with the values from the argument.
C
C This subroutine is called from TELEMAC-2D.
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      DOUBLE PRECISION, INTENT(INOUT) :: ARRAY_FROM_T2D(NPOIN2D)
C-----------------------------------------------------------------------
      DEPTH_FROM_T2D(:) = ARRAY_FROM_T2D(:)
C-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE DEPTH_FILL

      SUBROUTINE DEPTH_GET(ARRAY_FROM_ESTEL3D)
C-----------------------------------------------------------------------
C  It basically reads the array DEPTH_FROM_T2D so that ESTEL-3D can
C use it for its boundary conditions.
C
C This subroutine is called from ESTEL-3D.
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      DOUBLE PRECISION, INTENT(INOUT) :: ARRAY_FROM_ESTEL3D(NPOIN2D)
C-----------------------------------------------------------------------
      ARRAY_FROM_ESTEL3D(:) = DEPTH_FROM_T2D(:)
C-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE DEPTH_GET

      SUBROUTINE INFILTRATION_GET(SMH,UNSV2D,YASMH)
C-----------------------------------------------------------------------
C Adds the infiltration term to the source term SMH and switch YASMH to
C TRUE. Note that a mass-vector is required as argument because the flux
C calculated within ESTEL-3D is multiplied by a mass vector and the
C division is easier to do from within TELEMAC-2D for mesh reasons.
C
C This subroutine is called from TELEMAC-2D.
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      DOUBLE PRECISION, INTENT(INOUT) :: SMH(NPOIN2D)
      DOUBLE PRECISION, INTENT(IN)    :: UNSV2D(NPOIN2D)
      LOGICAL, INTENT(INOUT)          :: YASMH
C-----------------------------------------------------------------------
      IF(DOINFILTRATION) THEN
        YASMH  = .TRUE.
        SMH(:) = SMH(:) + FLUX_FROM_ESTEL3D(:)*UNSV2D(:)
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE INFILTRATION_GET

C-----------------------------------------------------------------------
      END MODULE M_COUPLING_ESTEL3D
