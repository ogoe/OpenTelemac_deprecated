C                       *********************
                        SUBROUTINE AS3_1111_Q
C                       *********************
C
     *(XM,NSEG1,XMT,NELMAX,NELEM,ELTSEG1,ELTSEG2,ELTSEG3,
     *                           ORISEG1,ORISEG2,ORISEG3)
C
C***********************************************************************
C BIEF VERSION 5.1         30/06/99   J-M HERVOUET (LNH) 01 30 87 80 18
C***********************************************************************
C
C FONCTION : ASSEMBLING EXTRA-DIAGONAL TERMS OF MATRICES
C            IN THE CASE OF EDGE-BASED STORAGE
C
C            CASE OF LINEAR-LINEAR ELEMENT WITH NON SYMMETRICAL MATRIX
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !  XM            !<-- ! TERMES EXTRA-DIAGONAUX ASSEMBLES XA12,23,31
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
      INTEGER         , INTENT(IN)    :: NELMAX,NELEM,NSEG1
      INTEGER         , INTENT(IN)    :: ELTSEG1(NELMAX)
      INTEGER         , INTENT(IN)    :: ELTSEG2(NELMAX)
      INTEGER         , INTENT(IN)    :: ELTSEG3(NELMAX)
      INTEGER         , INTENT(IN)    :: ORISEG1(NELMAX)
      INTEGER         , INTENT(IN)    :: ORISEG2(NELMAX)
      INTEGER         , INTENT(IN)    :: ORISEG3(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: XMT(NELMAX,*)
      DOUBLE PRECISION, INTENT(INOUT) :: XM(NSEG1,2)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER ISEG,IELEM
C
C-----------------------------------------------------------------------
C
C  INITIALISATION
C
      DO ISEG = 1 , NSEG1
        XM(ISEG,1) = 0.D0
        XM(ISEG,2) = 0.D0
      ENDDO
C
C  ASSEMBLAGE
C
      DO IELEM = 1,NELEM
C         TERME 12
          XM(ELTSEG1(IELEM),  ORISEG1(IELEM)) 
     *  = XM(ELTSEG1(IELEM),  ORISEG1(IELEM)) + XMT(IELEM,1)
C         TERME 23
          XM(ELTSEG2(IELEM),  ORISEG2(IELEM))
     *  = XM(ELTSEG2(IELEM),  ORISEG2(IELEM)) + XMT(IELEM,3)
C         TERME 31 
          XM(ELTSEG3(IELEM),  ORISEG3(IELEM))
     *  = XM(ELTSEG3(IELEM),  ORISEG3(IELEM)) + XMT(IELEM,5)
C         TERME 21 
          XM(ELTSEG1(IELEM),3-ORISEG1(IELEM))
     *  = XM(ELTSEG1(IELEM),3-ORISEG1(IELEM)) + XMT(IELEM,4)
C         TERME 32 
          XM(ELTSEG2(IELEM),3-ORISEG2(IELEM))
     *  = XM(ELTSEG2(IELEM),3-ORISEG2(IELEM)) + XMT(IELEM,6)
C         TERME 13 
          XM(ELTSEG3(IELEM),3-ORISEG3(IELEM))
     *  = XM(ELTSEG3(IELEM),3-ORISEG3(IELEM)) + XMT(IELEM,2)    
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
