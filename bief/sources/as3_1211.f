C                       *******************
                        SUBROUTINE AS3_1211
C                       *******************
C
     *(XM,NSEG11,NSEG12,XMT,NELMAX,NELEM,ELTSEG1,ELTSEG2,ELTSEG3,
     *                                   ELTSEG4,ELTSEG5,ELTSEG6,
     *                                   ORISEG1,ORISEG2,ORISEG3)
C
C***********************************************************************
C BIEF VERSION 5.6         29/12/05   J-M HERVOUET (LNH) 01 30 87 80 18
C***********************************************************************
C
C FONCTION : ASSEMBLING EXTRA-DIAGONAL TERMS OF MATRICES
C            IN THE CASE OF EDGE-BASED STORAGE
C
C            CASE OF QUASIBUBBLE-LINEAR ELEMENT
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !  XMAS          !<-- ! TERMES EXTRA-DIAGONAUX ASSEMBLES XA12,23,31
C !  XM2           ! -->! TERMES EXTRA-DIAGONAUX XA21,32,31
C !  TR            ! -->! TABLEAU DE TRAVAIL DE TAILLE > NPTFR
C !  NELMAX        ! -->! PREMIERE DIMENSION DE IKLE ET W.
C !                !    ! (CAS D'UN MAILLAGE ADAPTATIF)
C !  NELEM         ! -->! NOMBRE D'ELEMENTS DANS LE MAILLAGE.
C !  NPTFR         ! -->! NOMBRE DE POINTS FRONTIERES.
C !________________!____!_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELANT : ASSVEC
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(IN)    :: NELMAX,NELEM,NSEG11,NSEG12
      INTEGER         , INTENT(IN)    :: ELTSEG1(NELMAX)
      INTEGER         , INTENT(IN)    :: ELTSEG2(NELMAX)
      INTEGER         , INTENT(IN)    :: ELTSEG3(NELMAX)
      INTEGER         , INTENT(IN)    :: ELTSEG4(NELMAX)
      INTEGER         , INTENT(IN)    :: ELTSEG5(NELMAX)
      INTEGER         , INTENT(IN)    :: ELTSEG6(NELMAX)
      INTEGER         , INTENT(IN)    :: ORISEG1(NELMAX)
      INTEGER         , INTENT(IN)    :: ORISEG2(NELMAX)
      INTEGER         , INTENT(IN)    :: ORISEG3(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: XMT(NELMAX,*)
      DOUBLE PRECISION, INTENT(INOUT) :: XM(NSEG11+NSEG12)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER ISEG,IELEM
C
C-----------------------------------------------------------------------
C
C  INITIALISATION
C
      DO ISEG = 1 , 2*NSEG11
        XM(ISEG) = 0.D0
      ENDDO
C
C  ASSEMBLAGE PARTIE P1 (ENTRE 1 ET 2*NSEG11)
C
      DO IELEM = 1,NELEM
C         TERME 12
          XM(ELTSEG1(IELEM)+NSEG11*(ORISEG1(IELEM)-1)) 
     *  = XM(ELTSEG1(IELEM)+NSEG11*(ORISEG1(IELEM)-1)) + XMT(IELEM,01)
C         TERME 23
          XM(ELTSEG2(IELEM)+NSEG11*(ORISEG2(IELEM)-1))
     *  = XM(ELTSEG2(IELEM)+NSEG11*(ORISEG2(IELEM)-1)) + XMT(IELEM,04)
C         TERME 31 
          XM(ELTSEG3(IELEM)+NSEG11*(ORISEG3(IELEM)-1))
     *  = XM(ELTSEG3(IELEM)+NSEG11*(ORISEG3(IELEM)-1)) + XMT(IELEM,05)
C         TERME 21 
          XM(ELTSEG1(IELEM)+NSEG11*(2-ORISEG1(IELEM)))
     *  = XM(ELTSEG1(IELEM)+NSEG11*(2-ORISEG1(IELEM))) + XMT(IELEM,03)
C         TERME 32 
          XM(ELTSEG2(IELEM)+NSEG11*(2-ORISEG2(IELEM)))
     *  = XM(ELTSEG2(IELEM)+NSEG11*(2-ORISEG2(IELEM))) + XMT(IELEM,06)
C         TERME 13 
          XM(ELTSEG3(IELEM)+NSEG11*(2-ORISEG3(IELEM)))
     *  = XM(ELTSEG3(IELEM)+NSEG11*(2-ORISEG3(IELEM))) + XMT(IELEM,02)    
      ENDDO
C
C  ASSEMBLAGE PARTIE QUASIBULLE
C  ENTRE XM(2*NSEG11+1) ET XM(NSEG11+NSEG12)
C  VOIR DANS STOSEG COMMENT SONT CONSTRUITS ELTSEG4,5 ET 6
C
      DO IELEM = 1,NELEM
C       TERME 41
        XM(ELTSEG4(IELEM)+NSEG11) = XMT(IELEM,07)
C       TERME 42
        XM(ELTSEG5(IELEM)+NSEG11) = XMT(IELEM,08)
C       TERME 43 
        XM(ELTSEG6(IELEM)+NSEG11) = XMT(IELEM,09)
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
