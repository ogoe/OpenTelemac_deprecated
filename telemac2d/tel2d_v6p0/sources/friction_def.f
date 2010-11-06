C     *******************
      MODULE FRICTION_DEF
C     ******************* 
C
C***********************************************************************
C TELEMAC 2D VERSION 5.6
C***********************************************************************
C
      TYPE POINTER_TO_FRICTION
C         SEQUENCE
         TYPE(FRICTION_OBJ), POINTER :: P
      END TYPE POINTER_TO_FRICTION
C
      TYPE FRICTION_OBJ
C         SEQUENCE
         INTEGER          :: gnumb(2) ! Global number of the zone
         INTEGER          :: rtype(2) ! Type of law used

         ! Use REAL because CHESTR IS saved as simple precision in selafin data
         ! --------------------------------------------------------------------
         DOUBLE PRECISION :: rcoef(2) ! Friction parameter
         DOUBLE PRECISION :: nDef(2)  ! Default Manning (for C-W law)
         DOUBLE PRECISION :: dp       ! Diameter of roughness element
         DOUBLE PRECISION :: sp       ! Spacing of roughness element
         TYPE(POINTER_TO_FRICTION), POINTER, DIMENSION(:) :: ADR
      END TYPE FRICTION_OBJ
C
      END MODULE FRICTION_DEF

