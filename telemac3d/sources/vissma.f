!                    *****************
                     SUBROUTINE VISSMA
!                    *****************
!
     &(VISCVI,VISCTA,DNUTAH,DNUVIH,DNUVIV,DNUTAV,
     & U,V,W,TRAV1,TRAV2,TRAV3,TRAV4,TRAV5,TRAV6,
     & SVIDE,MESH3,IELM3,NTRAC,MSK,MASKEL,ITURBV)
!
!***********************************************************************
! TELEMAC3D   V6P1                                  21/08/2010
!***********************************************************************
!
!brief    INITIALISES VISCOSITIES
!+                FOR THE SMAGORINSKI MODEL.
!
!history  JACEK A. JANKOWSKI PINXIT
!+        **/03/99
!+
!+   FORTRAN95 VERSION
!
!history  OLIVER GOETHEL    UNI-HAN
!+        **/02/04
!+        V5P7
!+   SMAGORINSKY 3D
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| DNUTAH         |-->| COEFFICIENT FOR HORIZONTAL DIFFUSION OF TRACER
!| DNUTAV         |-->| COEFFICIENT FOR VERTICAL DIFFUSION OF TRACER
!| DNUVIH         |-->| COEFFICIENT FOR HORIZONTAL DIFFUSION OF VELOCITIES
!| DNUVIV         |-->| COEFFICIENT FOR VERTICAL DIFFUSION OF VELOCITIES
!| IELM3          |-->| TYPE OF ELEMENT
!| ITURBV         |-->| VERTICAL TURBULENCE MODEL
!| MASKEL         |-->| MASKING OF ELEMENTS
!|                |   | =1. : NORMAL   =0. : MASKED ELEMENT
!| MESH3          |---| 3D MESH
!| MSK            |-->| IF YES, THERE IS MASKED ELEMENTS.
!| NTRAC          |-->| NUMBER OF ACTIVE TRACERS
!| SVIDE          |<->| VOID STRUCTURE
!| TRAV1          |<->| WORK ARRAY
!| TRAV2          |<->| WORK ARRAY
!| TRAV3          |<->| WORK ARRAY
!| TRAV4          |<->| WORK ARRAY
!| TRAV5          |<->| WORK ARRAY
!| TRAV6          |<->| WORK ARRAY
!| U              |-->| COMPONENT OF VELOCITY
!| V              |-->| COMPONENT OF VELOCITY
!| VISCTA         |<->| TURBULENT VISCOSITY COEFFICIENTS FOR TRACERS
!| VISCVI         |<->| TURBULENT VISCOSITY COEFFICIENTS FOR VELOCITIES
!| W              |-->| COMPONENT OF VELOCITY
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE INTERFACE_TELEMAC3D, EX_VISSMA => VISSMA
      USE DECLARATIONS_TELEMAC
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)            :: NTRAC,ITURBV
      INTEGER, INTENT(IN)            :: IELM3
      LOGICAL, INTENT(IN)            :: MSK
      TYPE (BIEF_OBJ), INTENT(IN)    :: U, V, W
      TYPE (BIEF_OBJ), INTENT(INOUT) :: VISCVI, VISCTA
      TYPE (BIEF_OBJ), INTENT(INOUT) :: TRAV1, TRAV2, TRAV3, TRAV4
      TYPE (BIEF_OBJ), INTENT(INOUT) :: TRAV5      !  NUSMAG
      TYPE (BIEF_OBJ), INTENT(INOUT) :: TRAV6
      TYPE (BIEF_OBJ), INTENT(IN)    :: MASKEL
      TYPE (BIEF_OBJ), INTENT(INOUT) :: SVIDE
      TYPE (BIEF_MESH)               :: MESH3
      DOUBLE PRECISION, INTENT(IN)   :: DNUVIH,DNUTAH,DNUVIV,DNUTAV
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER ITRAC
!
!***********************************************************************
!
! INITIALISES VISCOSITIES
!
! VISCVI%ADR(1)%P IS THE X HORIZONTAL VISCOSITY
! VISCVI%ADR(2)%P IS THE Y HORIZONTAL VISCOSITY
! VISCVI%ADR(3)%P IS THE Z (VERTICAL) VISCOSITY
!
! FOR THE TRACERS:
!
! VISCTA%ADR(ITRAC)%P%ADR(1)%P IS THE X HORIZONTAL DIFFUSIVITY FOR THE
!     TRACER NUMBER ITRAC, ETC...
!
!***********************************************************************
!
      IF(ITURBV.NE.4) THEN
!
      CALL SMAGO(U,V,TRAV1,TRAV2,TRAV3,TRAV4,TRAV5,MESH3,IELM3,
     &           MSK,MASKEL)
!
! VISCOSITY COMPUTED BY SMAGORINSKI : IN TRAV5
!
      CALL OS('X=Y+C   ',X=VISCVI%ADR(1)%P,Y=TRAV5,C=DNUVIH)
      CALL OS('X=Y+C   ',X=VISCVI%ADR(2)%P,Y=TRAV5,C=DNUVIH)
!
      IF (NTRAC.NE.0) THEN
!
        DO ITRAC=1,NTRAC
          CALL OS('X=Y+C   ',X=VISCTA%ADR(ITRAC)%P%ADR(1)%P,
     &                       Y=TRAV5,C=DNUTAH)
          CALL OS('X=Y+C   ',X=VISCTA%ADR(ITRAC)%P%ADR(2)%P,
     &                       Y=TRAV5,C=DNUTAH)
        ENDDO
!
      ENDIF
!
      ELSE !ITURBV=4 -> 3D SMAGORINSKI
!
      CALL SMAGO3D(U,V,W,TRAV1,TRAV2,TRAV3,TRAV4,TRAV5,TRAV6,
     &             SVIDE,MESH3,IELM3,MSK,MASKEL)
!
      CALL OS('X=Y+C   ',X=VISCVI%ADR(1)%P,Y=TRAV5,C=DNUVIH)
      CALL OS('X=Y+C   ',X=VISCVI%ADR(2)%P,Y=TRAV5,C=DNUVIH)
      CALL OS('X=Y+C   ',X=VISCVI%ADR(3)%P,Y=TRAV5,C=DNUVIV)
!
      IF(NTRAC.NE.0) THEN
!
        DO ITRAC=1,NTRAC
          CALL OS ('X=Y+C   ',X=VISCTA%ADR(ITRAC)%P%ADR(1)%P,
     &                        Y=TRAV5,C=DNUTAH)
          CALL OS ('X=Y+C   ',X=VISCTA%ADR(ITRAC)%P%ADR(2)%P,
     &                        Y=TRAV5,C=DNUTAH)
          CALL OS ('X=Y+C   ',X=VISCTA%ADR(ITRAC)%P%ADR(3)%P,
     &                        Y=TRAV5,C=DNUTAV)
        ENDDO
!
      ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
