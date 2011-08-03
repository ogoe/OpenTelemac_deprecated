!                       *****************
                        SUBROUTINE INITRI
!                       *****************
!
     &( NPOIN1,TYPELE,NGEO,NFO1)
!
!***********************************************************************
! PROGICIEL : STBTEL V5.2      07.04.1993  P. LANG     (LHF)
!***********************************************************************
!
!     FONCTION  : INITIALISATION DES INFORMATIONS DANS LE CAS DE TRIGRID
!
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________
! !      NOM       |MODE!                   ROLE
! !________________|____!______________________________________________
! ! NPOIN1         !<-- ! NOMBRE DE POINTS DU MAILLAGE
! ! TYPELE         !<-- ! TYPE D'ELEMENT DU MAILLAGE (ICI TRIAGLES)
! !________________!____!______________________________________________
! ! COMMON:        !    !
! !  GEO:          !    !
! !    MESH        !<-- ! TYPE DES ELEMENTS DU MAILLAGE
! !    NDP         !<-- ! NOMBRE DE NOEUDS PAR ELEMENTS ( ICI FORCEMENT
! !    NPOIN       !<-- ! NOMBRE TOTAL DE NOEUDS DU MAILLAGE
! !    NELEM       !<-- ! NOMBRE TOTAL D'ELEMENTS DU MAILLAGE
! !    NPMAX       ! -->! DIMENSION EFFECTIVE DES TABLEAUX X ET Y
! !                !    ! (NPMAX = NPOIN + 100)
! !    NELMAX      ! -->! DIMENSION EFFECTIVE DES TABLEAUX CONCERNANT
! !                !    ! LES ELEMENTS (NELMAX = NELEM + 200)
! !  FICH:         !    !
! !    NRES        ! -->! NUMERO DU CANAL DU FICHIER DE SERAFIN
! !    NGEO       ! -->! NUMERO DU CANAL DU FICHIER MAILLEUR
! !    NLIM      ! -->! NUMERO DU CANAL DU FICHIER DYNAM DE TELEMAC
! �    NFO1      � -->! NUMERO DU CANAL DU FICHIER TRIANGLE DE TRIGRID
! !________________!____!______________________________________________
! MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
!-----------------------------------------------------------------------
! APPELE PAR : HOMERE
! APPEL DE :
!***********************************************************************
!
      IMPLICIT NONE
!
      INTEGER       MESH, NDP, NPOIN, NELEM, NPMAX, NELMAX
      INTEGER       NGEO
      INTEGER       NPOIN1, NFO1
      CHARACTER*(*) TYPELE
      CHARACTER*1   ZDUMMY
!
! COMMON
!
      COMMON/GEO/MESH, NDP, NPOIN, NELEM, NPMAX, NELMAX
!
      REWIND (NGEO)
      REWIND (NFO1)
      READ (NGEO,*) NPOIN1
      NPOIN = NPOIN1
      NELEM = 0
 1    CONTINUE
        READ (NFO1, '(A1)', END=9000) ZDUMMY
        NELEM = NELEM + 1
      GOTO 1
!
 9000 CONTINUE
      TYPELE = 'TRIANGLES  '
      NDP = 3
      MESH = 3
!
      RETURN
      END