!                       *****************
                        SUBROUTINE ELEB3D                              
!                       *****************
!
     *(IKLE3 , NBOR , KP1BOR , NELBOR, IKLBOR, NULONE,                              
     * NELEM2, NPOIN2, NPLAN, NETAGE, NPTFR )
!
!***********************************************************************
C BIEF VERSION 5.9      23/06/2008    J-M HERVOUET (LNHE) 01 30 87 80 18
!                                                    
!***********************************************************************
!
!      FONCTION:
!      =========
!
!    CONSTRUCTION DU MAILLAGE 3D : A L'ENTREE, LES TABLEAUX DU MAILLAGE
!    3D REMPLIS PAR UN APPEL PREALABLE DE ELEBD. A LA SORTIE, LES
!    TABLEAUX COMPLETES EN 3D.
!
!----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________.
! !  NOM           !MODE!                  ROLE                        !
! !________________!____!______________________________________________!
! !  IKLE3         !<-- ! CORRESPONDANCE LOCALE - GLOBALE EN 3D        !
! !  SURFAC        !<-->! SURFACE DES TRIANGLES ETENDUE EN 3D          !
! !  NBOR3         !<-- ! CORRESPONDANCE ENTRE LA NUMEROTATION DE BORD !
! !                !    ! ET LA NUMEROTATION GLOBALE (3D)              !
! !  NBOR          ! -->! CORRESPONDANCE ENTRE LA NUMEROTATION DE BORD !
! !                !    ! ET LA NUMEROTATION GLOBALE (2D)              !
! !  KP1BOR        ! -->! PT FRONT. SUIVANT LE PT FRONT. CONSIDERE     !
! !  NELBOR        ! -->! NUMERO GLOBAUX DES ELEMENTS DE BORD          !
! !  IKLBOR        !<-- ! TABLE DE CONNECTIVITE ELEMENTS DE BORD       !
! !  NELBO3        !<-- ! ASSOCIE A CHAQUE FACE DE BORD L'ELEMENT 3D   !
! !                !    ! AUQUEL ELLE APPARTIENT                       !
! !  NULONE        !<-- ! ASSOCIE LA NUMEROTATION LOCALE DE BORD A LA  !
! !                !    ! NUMEROTATION LOCALE 3D                       !
! !  IKLE2         ! -->! CORRESPONDANCE LOCALE - GLOGALE EN 2D        !
! !  NELEM2        ! -->! NOMBRE D'ELEMENTS EN 2D                      !
! !  NPOIN2        ! -->! NOMBRE DE POINTS 2D                          !
! !  NPLAN         ! -->! NOMBRE DE PLANS HORIZONTAUX                  !
! !  NETAGE        ! -->! NPLAN - 1                                    !
! !  NPTFR         ! -->! NOMBRE DE POINTS DE BORD 2D                  !
! !________________!____!______________________________________________!
! MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
!
!-----------------------------------------------------------------------
!
! SOUS-PROGRAMME APPELE PAR : MITRID
! SOUS-PROGRAMMES APPELES : OV , PLANTE
!
!***********************************************************************
!
      USE BIEF, EX_ELEB3D => ELEB3D
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NELEM2, NPOIN2, NPLAN, NETAGE, NPTFR
      INTEGER, INTENT(INOUT) :: IKLE3(NELEM2,NETAGE,6)
      INTEGER, INTENT(INOUT) :: IKLBOR(NPTFR,NETAGE,4)
      INTEGER, INTENT(INOUT) :: NULONE(NPTFR,NETAGE,4)
      INTEGER, INTENT(INOUT) :: NELBOR(NPTFR*NETAGE), NBOR(NPTFR*NPLAN)
      INTEGER, INTENT(INOUT) :: KP1BOR(NPTFR)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IPOIN
      INTEGER IETAGE,IPTFR
!
!***********************************************************************
!
! TABLES DE CONNECTIVITE POUR LES FACES DE BORDS --> IKLBOR , NBOR3 ,
! CORRESPONDANCE NUMEROTATION LOCALE BORD NUMEROTATION LOCALE 3D
! --> NULONE
! ET CALCUL DE NELBO3
!
!     BORDS LATERAUX
!
      DO IETAGE = 1,NETAGE
        DO IPTFR = 1,NPTFR
          IKLBOR(IPTFR,IETAGE,1) =        IPTFR  + (IETAGE-1)*NPTFR
          IKLBOR(IPTFR,IETAGE,2) = KP1BOR(IPTFR) + (IETAGE-1)*NPTFR
          IKLBOR(IPTFR,IETAGE,3) = IKLBOR(IPTFR,IETAGE,2) + NPTFR
          IKLBOR(IPTFR,IETAGE,4) = IKLBOR(IPTFR,IETAGE,1) + NPTFR
          IPOIN = NBOR(IPTFR)
          NBOR(IPTFR +(IETAGE-1)*NPTFR)=IPOIN+(IETAGE-1)*NPOIN2
          IELEM = NELBOR(IPTFR)
          IF(IELEM.GT.0) THEN
            NELBOR(IPTFR+(IETAGE-1)*NPTFR)=IELEM+(IETAGE-1)*NELEM2
            IF(IPOIN.EQ.IKLE3(IELEM,1,1)) THEN
              NULONE(IPTFR,IETAGE,1) = 1
              NULONE(IPTFR,IETAGE,2) = 2
              NULONE(IPTFR,IETAGE,3) = 5
              NULONE(IPTFR,IETAGE,4) = 4
            ELSEIF(IPOIN.EQ.IKLE3(IELEM,1,2)) THEN
              NULONE(IPTFR,IETAGE,1) = 2
              NULONE(IPTFR,IETAGE,2) = 3
              NULONE(IPTFR,IETAGE,3) = 6
              NULONE(IPTFR,IETAGE,4) = 5
            ELSEIF(IPOIN.EQ.IKLE3(IELEM,1,3)) THEN
              NULONE(IPTFR,IETAGE,1) = 3
              NULONE(IPTFR,IETAGE,2) = 1
              NULONE(IPTFR,IETAGE,3) = 4
              NULONE(IPTFR,IETAGE,4) = 6
            ELSE
              IF(LNG.EQ.1) WRITE(LU,101) IPOIN
              IF(LNG.EQ.2) WRITE(LU,102) IPOIN
              CALL PLANTE(1)
              STOP
            ENDIF
          ELSEIF(NCSIZE.GT.1) THEN
            NULONE(IPTFR,IETAGE,1) = 0
            NULONE(IPTFR,IETAGE,2) = 0
            NULONE(IPTFR,IETAGE,3) = 0
            NULONE(IPTFR,IETAGE,4) = 0 
            NELBOR(IPTFR+(IETAGE-1)*NPTFR)=0
          ELSE
            IF(LNG.EQ.1) WRITE(LU,101) IPOIN
            IF(LNG.EQ.2) WRITE(LU,102) IPOIN
            CALL PLANTE(1)
            STOP 
          ENDIF
        ENDDO
      ENDDO
!
!     COMPLETING NBOR IN VIEW OF 2D VALUES
!
      DO IPTFR = 1,NPTFR
         NBOR(IPTFR +(NPLAN-1)*NPTFR) = NBOR(IPTFR) + NETAGE*NPOIN2
      END DO
!
!-----------------------------------------------------------------------
!
101   FORMAT(' ELEB3D : PROBLEME A LA CONSTRUCTION DE NULONE, IPOIN =',
     &  I6)
102   FORMAT(' ELEB3D: PROBLEM WHEN BUILDING NULONE, IPOIN =',I6)
!
!-----------------------------------------------------------------------
!
      RETURN
      END
