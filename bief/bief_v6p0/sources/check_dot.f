!                       ********************
                        SUBROUTINE CHECK_DOT
!                       ********************
!
     *(X,T,TEXTE,MESH)
!
!***********************************************************************
! TELEMAC 3D VERSION 5.9    18/11/08  J-M HERVOUET (LNHE) 01 30 87 80 18
! 
!***********************************************************************
!
!      FONCTION:
!      =========
!
!      IN PARALLEL MODE, PRINTS THE EUCLIDIAN NORM OF A VECTOR WHICH
!      HAS NOT BEEN ASSEMBLED WITH PARCOM. E.G. A RIGHT HAND SIDE
!      BEFORE CALLING SOLVE.
!
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________.
! !  NOM           !MODE!                  ROLE                        !
! !________________!____!______________________________________________!
! !  X             ! -->! THE BIEF_OBJ STRUCTURE WITH THE VECTOR
! !  T             !<-->! A WORK BIEF_OBJ STRUCTURE 
! !  TEXTE         ! -->! A TEXT TO BE PRINTED 
! !  MESH          ! -->! MESH STRUCTURE            
! !________________!____!______________________________________________!
! MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
!
!-----------------------------------------------------------------------
!
! SOUS-PROGRAMME APPELE PAR : 
! SOUS-PROGRAMME APPELE : PARCOM
!
!***********************************************************************
!
      USE BIEF
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER(LEN=*) :: TEXTE
      TYPE(BIEF_OBJ)   :: X,T
      TYPE(BIEF_MESH)  :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CALL OS('X=Y     ',X=T,Y=X)
      IF(NCSIZE.GT.1) CALL PARCOM(T,2,MESH)
      WRITE(LU,*) TEXTE,'=',P_DOTS(T,T,MESH)
!
!-----------------------------------------------------------------------
!
      RETURN
      END
