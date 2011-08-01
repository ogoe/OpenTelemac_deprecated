!                       *******************
                        SUBROUTINE MASBAS2D
!                       *******************
!
     *(VOLU2D,V2DPAR,UNSV2D,IELM,MESH,MSK,MASKEL,T1,S)
!
!***********************************************************************
! TELEMAC 2D VERSION 5.7    14/06/2006     J.-M. HERVOUET 01 30 87 80 18
!                          
!***********************************************************************
!
!      FONCTION:
!      =========
!
!      COMPUTING VARIOUS VOLUMES OF 2D BASIS AND THE INVERSE
!
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________.
! !  NOM           !MODE!                  ROLE                        !
! !________________!____!______________________________________________!
! !  VOLU2D        !<-->! INTEGRAL OF TEST FUNCTIONS, WITH MASKING 
! !  V2DPAR        !<-->! AS VOLU2D IF NOT PARALLEL
! !                !    ! IN PARALLEL COMPLETED WITH OTHER SUBDOMAINS
! !  UNSV2D        !<-->! INVERSE OF INTEGRAL OF TEST FUNCTIONS 
! !                !    ! WITHOUT MASKING
! !  IELM          ! -->! TYPE OF ELEMENT (11 FOR LINEAR)
! !  MESH          !<-->! MESH
! !  MSK           ! -->! IF YES, THERE IS MASKING, MASKEL IS TO BE USED
! !  MASKEL        ! -->! ARRAY OF MASKS, PER ELEMENT
! !  T1            !<-->! BIEF_OBJ STRUCTURE FOR LOCAL WORK
! !  S             ! -->! EMPTY BIEF_OBJ STRUCTURE
! !________________!____!______________________________________________!
! MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
!
!-----------------------------------------------------------------------
!
      USE BIEF
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)            :: IELM
      LOGICAL, INTENT(IN)            :: MSK
      TYPE(BIEF_OBJ) , INTENT(INOUT) :: VOLU2D,V2DPAR,UNSV2D,T1
      TYPE(BIEF_OBJ) , INTENT(IN)    :: MASKEL,S
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     VOLU2D : VOLUME WITH POSSIBLE MASKING
C
      CALL VECTOR(VOLU2D,'=','MASBAS          ',IELM,1.D0, 
     *            S,S,S,S,S,S,MESH,MSK,MASKEL)
C
C     V2DPAR : LIKE VOLU2D BUT IN PARALLEL VALUES COMPLETED AT 
C              INTERFACES BETWEEN SUBDOMAINS
C
      CALL OS('X=Y     ',X=V2DPAR,Y=VOLU2D)
      IF(NCSIZE.GT.1) CALL PARCOM(V2DPAR,2,MESH)
C
C     INVERSE OF VOLUMES (DONE WITHOUT MASKING), THERE SHOULD BE
C     NO DIVISION BY ZERO, UNLESS ELEMENT WITH NO AREA
C
      IF(MSK) THEN
        CALL VECTOR(T1,'=','MASBAS          ',IELM,1.D0, 
     *              S,S,S,S,S,S,MESH,.FALSE.,MASKEL)
        IF(NCSIZE.GT.1) CALL PARCOM(T1,2,MESH)
        CALL OS('X=1/Y   ',X=UNSV2D,Y=T1)
      ELSE
        CALL OS('X=1/Y   ',X=UNSV2D,Y=V2DPAR)
      ENDIF      
!
!-----------------------------------------------------------------------
!
      RETURN
      END
