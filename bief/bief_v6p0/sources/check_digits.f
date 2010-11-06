!                       ***********************
                        SUBROUTINE CHECK_DIGITS
!                       ***********************
!
     *(F,T1,MESH)
!
!***********************************************************************
! TELEMAC 3D VERSION 5.9    02/06/08  J-M HERVOUET (LNHE) 01 30 87 80 18
! 
!***********************************************************************
!
!      FONCTION:
!      =========
!
!      IN PARALLEL MODE, CHECKING THAT PROCESSORS SHARING AN INTERFACE
!      POINT HAVE EXACTLY THE SAME VALUE OF ARRAY F.
!
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________.
! !  NOM           !MODE!                  ROLE                        !
! !________________!____!______________________________________________!
! !  F             ! -->! TABLEAU A VERIFIER
! !  T1            !<-->! TABLEAU DE TRAVAIL 
! !  MESH          ! -->! STRUCTURE DE MAILLAGE            
! !________________!____!______________________________________________!
! MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
!
!-----------------------------------------------------------------------
!
! SOUS-PROGRAMME APPELE PAR : LIMTEL , MITRID , TELSOU
! SOUS-PROGRAMME APPELE : OV
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
      TYPE(BIEF_OBJ),  INTENT(IN   ) :: F
      TYPE(BIEF_OBJ),  INTENT(INOUT) :: T1
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,ISTOP
C
      INTEGER  P_IMAX
      EXTERNAL P_IMAX
C
C-----------------------------------------------------------------------
C
      CALL OS('X=Y     ',X=T1,Y=F)
      CALL PARCOM(T1,3,MESH)
      ISTOP=0
      DO I=1,T1%DIM1
        IF(T1%R(I).NE.F%R(I)) THEN
          IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'CHECK_DIGITS : DIFFERENCE DANS ',F%NAME
          WRITE(LU,*) '               AU POINT LOCAL ',I
          WRITE(LU,*) '               =  POINT GLOBAL ',MESH%KNOLG%I(I)
          WRITE(LU,*) '               VALEUR ',F%R(I) 
          WRITE(LU,*) '               MINIMUM ',T1%R(I)
          WRITE(LU,*) '            DIFFERENCE ',F%R(I)-T1%R(I)
          ENDIF
          IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'CHECK_DIGITS : DIFFERENCE IN ',F%NAME
          WRITE(LU,*) '               AT LOCAL POINT ',I
          WRITE(LU,*) '               =  GLOBAL POINT ',MESH%KNOLG%I(I)
          WRITE(LU,*) '               VALUE ',F%R(I) 
          WRITE(LU,*) '               MINIMUM ',T1%R(I)
          WRITE(LU,*) '            DIFFERENCE ',F%R(I)-T1%R(I)
          ENDIF
          ISTOP=I
        ENDIF
      ENDDO
!
      IF(NCSIZE.GT.1) ISTOP=P_IMAX(ISTOP)
      IF(ISTOP.GT.0) THEN
        WRITE(LU,*) 'CHECK_DIGITS : ERREUR SUR VECTEUR ',F%NAME
        CALL PLANTE(1)
        STOP
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
