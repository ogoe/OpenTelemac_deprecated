C                       *******************
                        SUBROUTINE AS3_1311
C                       *******************
C
     *(XM,NSEG11,NSEG13,XMT,NELMAX,NELEM,ELTSEG,ORISEG)
C
C***********************************************************************
C BIEF VERSION 6.0       02/02/2010   J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C FONCTION : ASSEMBLING EXTRA-DIAGONAL TERMS OF MATRICES (XMT)
C            IN THE CASE OF EDGE-BASED STORAGE
C
C            CASE OF LINEAR - QUADRATIC ELEMENT
C
C            LOCAL NUMBERING OF SEGMENTS IN A TRIANGLE (SEE COMP_SEG)
C
C            01 --> 1 - 2
C            02 --> 2 - 3
C            03 --> 3 - 1
C            04 --> 1 - 4
C            05 --> 2 - 5
C            06 --> 3 - 6
C            07 --> 2 - 4
C            08 --> 3 - 5
C            09 --> 1 - 6
C            10 --> 1 - 5
C            11 --> 2 - 6
C            12 --> 3 - 4
C            13 --> 4 - 5
C            14 --> 5 - 6
C            15 --> 6 - 4
C
C            TERMS IN XMT (STOTAGE GIVEN BY ARRAY CAQ(6,3,2) IN MATRIY:
C
C            01  -->  1-2
C            02  -->  1-3
C            03  -->  2-1
C            04  -->  2-3
C            05  -->  3-1
C            06  -->  3-2
C            07  -->  4-1
C            08  -->  4-2
C            09  -->  4-3
C            10  -->  5-1
C            11  -->  5-2
C            12  -->  5-3
C            13  -->  6-1
C            14  -->  6-2
C            15  -->  6-3
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !  XM            !<-- ! TERMES EXTRA-DIAGONAUX ASSEMBLES XA12,23,31
C !  NSEG11        ! -->! NUMBER OF LINEAR SEGMENTS
C !  NSEG13        ! -->! NUMBER OF QUADRATIC SEGMENTS -
C !                !    ! THE NUMBER OF PURELY QUADRATIC SEGMENTS
C !                !    ! (THEY ARE NOT CONSIDERED IN RECTANGULAR
C !                !    !  MATRICES)
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
      INTEGER         , INTENT(IN)    :: NELMAX,NELEM,NSEG11,NSEG13
      INTEGER         , INTENT(IN)    :: ELTSEG(NELMAX,15)
      INTEGER         , INTENT(IN)    :: ORISEG(NELMAX,15)
      DOUBLE PRECISION, INTENT(IN)    :: XMT(NELMAX,*)
      DOUBLE PRECISION, INTENT(INOUT) :: XM(NSEG11+NSEG13-3*NELEM)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER ISEG,IELEM
C
C-----------------------------------------------------------------------
C
C     INITIALISATION
C
C     -6*NELEM : SEGMENTS 10,11,12 NEED NO ASSEMBLY
C                SEGMENTS 13,14,15 ARE NOT CONSIDERED
      DO ISEG = 1 , NSEG11+NSEG13-6*NELEM
        XM(ISEG) = 0.D0
      ENDDO
C
C     ASSEMBLING, LINEAR PART
C
      DO IELEM = 1,NELEM
C
C        SEGMENT 1 (TERMS 1-2 ET 2-1)
C
         XM(ELTSEG(IELEM,1)+NSEG11*(ORISEG(IELEM,1)-1)) 
     *  =XM(ELTSEG(IELEM,1)+NSEG11*(ORISEG(IELEM,1)-1))+XMT(IELEM,01)
         XM(ELTSEG(IELEM,1)+NSEG11*(2-ORISEG(IELEM,1))) 
     *  =XM(ELTSEG(IELEM,1)+NSEG11*(2-ORISEG(IELEM,1)))+XMT(IELEM,03)
C
C        SEGMENT 2 (TERMS 2-3 ET 3-2)
C 
         XM(ELTSEG(IELEM,2)+NSEG11*(ORISEG(IELEM,2)-1)) 
     *  =XM(ELTSEG(IELEM,2)+NSEG11*(ORISEG(IELEM,2)-1))+XMT(IELEM,04)
         XM(ELTSEG(IELEM,2)+NSEG11*(2-ORISEG(IELEM,2))) 
     *  =XM(ELTSEG(IELEM,2)+NSEG11*(2-ORISEG(IELEM,2)))+XMT(IELEM,06)
C
C        SEGMENT 3 (TERMS 3-1 ET 1-3)
C 
         XM(ELTSEG(IELEM,3)+NSEG11*(ORISEG(IELEM,3)-1)) 
     *  =XM(ELTSEG(IELEM,3)+NSEG11*(ORISEG(IELEM,3)-1))+XMT(IELEM,05)
         XM(ELTSEG(IELEM,3)+NSEG11*(2-ORISEG(IELEM,3))) 
     *  =XM(ELTSEG(IELEM,3)+NSEG11*(2-ORISEG(IELEM,3)))+XMT(IELEM,02)
C
      ENDDO
C
C     ASSEMBLING, SEGMENTS BETWEEN LINEAR AND QUADRATIC POINTS
C     (I.E. THE REST BUT NOT 13, 14 AND 15)
C
C     ASSEMBLING THE QUADRATIC PART
C     BETWEEN XM(2*NSEG11+1) AND XM(NSEG11+NSEG13-3*NELEM)
C     SEE IN COMP_SEG HOW ELTSEG4,5,6,7,8,9,10,11,12 ARE BUILT,
C     THEIR NUMBERING STARTS AT NSEG11+1, HENCE HERE THE STORAGE IN
C     XM STARTS AT 2*NSEG11+1
C
      DO IELEM = 1,NELEM
C       TERM OF SEGMENT 1-4
        XM(ELTSEG(IELEM,04)+NSEG11) = 
     *  XM(ELTSEG(IELEM,04)+NSEG11) + XMT(IELEM,07)
C       TERM OF SEGMENT 2-5
        XM(ELTSEG(IELEM,05)+NSEG11) = 
     *  XM(ELTSEG(IELEM,05)+NSEG11) + XMT(IELEM,11)
C       TERM OF SEGMENT 3-6 
        XM(ELTSEG(IELEM,06)+NSEG11) = 
     *  XM(ELTSEG(IELEM,06)+NSEG11) + XMT(IELEM,15)
C       TERM OF SEGMENT 2-4
        XM(ELTSEG(IELEM,07)+NSEG11) = 
     *  XM(ELTSEG(IELEM,07)+NSEG11) + XMT(IELEM,08)
C       TERM OF SEGMENT 3-5
        XM(ELTSEG(IELEM,08)+NSEG11) = 
     *  XM(ELTSEG(IELEM,08)+NSEG11) + XMT(IELEM,12)
C       TERM OF SEGMENT 1-6 
        XM(ELTSEG(IELEM,09)+NSEG11) = 
     *  XM(ELTSEG(IELEM,09)+NSEG11) + XMT(IELEM,13)
      ENDDO
C
C     THESE 3 SEGMENTS ARE NOT SHARED, NO ASSEMBLY
C
      DO IELEM = 1,NELEM
C       TERM OF SEGMENT 1-5
        XM(ELTSEG(IELEM,10)+NSEG11) = XMT(IELEM,10)
C       TERM OF SEGMENT 2-6
        XM(ELTSEG(IELEM,11)+NSEG11) = XMT(IELEM,14)
C       TERM OF SEGMENT 3-4 
        XM(ELTSEG(IELEM,12)+NSEG11) = XMT(IELEM,09)
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
