!                       **************************
                        SUBROUTINE UTIMP_TELEMAC2D
!                       **************************
!
     &(LTL,ATL,GRADEBL,GRAPRDL,LISDEBL,LISPRDL) 
!
!***********************************************************************
! TELEMAC 2D VERSION 5.4    AUGUST 2003       JACEK A. JANKOWSKI PINXIT
! BAW KARLSRUHE                                  jacek.jankowski@baw.de
!***********************************************************************
!
!      FONCTION:
!      =========
!
!    IMPRIME D'EVENTUELLES DONNEES DEMANDEES PAR L'UTILISATEUR
!
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________.
! !  NOM           !MODE!                  ROLE                        !
! !________________!____!______________________________________________!
! !  LTL           ! -->! NUMERO DU PAS DE TEMPS                       !
! !  ATL           ! -->! TEMPS DU PAS DE TEMPS                        !
! !  GRADEBL       ! -->! 1ER PAS DE TEMPS POUR LES SORTIES GRAPHIQUES !
! !  GRAPRDL       ! -->! PERIODE DE SORTIE SUR LE FICHIER DE RESULTAT !
! !  LISDEBL       ! -->! 1ER PAS DE TEMPS POUR LES SORTIES LISTING    !
! !  LISPRDL       ! -->! PERIODE DE SORTIE LISTING                    !
! !________________!____!______________________________________________!
! MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
!
!-----------------------------------------------------------------------
!
! SOUS-PROGRAMME APPELE PAR : TELEMAC2D
!
! please note this subroutine is called in the same places as the 
! main output subroutine of telemac2d named DESIMP, i.e. twice: 
!
! (1) once a run, when ltl==0, independently if
!     OUTPUT OF INITIAL CONDITIONS : YES is set or not
! (2) each time step just after DESIMP-output 
!
!***********************************************************************
!
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(IN) :: ATL
      INTEGER, INTENT(IN) :: LTL,GRADEBL,GRAPRDL,LISDEBL,LISPRDL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
!
!***********************************************************************
! USER OUTPUT  
!
!
!
!-----------------------------------------------------------------------
!
      RETURN
      END SUBROUTINE UTIMP_TELEMAC2D
