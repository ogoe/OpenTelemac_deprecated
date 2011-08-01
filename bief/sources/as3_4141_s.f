C                       *********************
                        SUBROUTINE AS3_4141_S
C                       *********************
C
     *(XM,NSEG1,XMT,DIM1XMT,DIM2XMT,STOXMT,NELMAX,NELEM,ELTSEG)
C
C***********************************************************************
C BIEF VERSION 6.0         11/08/09   J-M HERVOUET (LNHE) 01 30 87 80 18
C
C 11/08/09 JMH : CROSSED AND VERTICAL SEGMENTS SWAPPED (SEE STOSEG41)
C 14/10/09 JMH : DIM1XMT,DIM2XMT,STOXMT ADDED, + CASE STOXMT=2
C
C***********************************************************************
C
C FONCTION : ASSEMBLING EXTRA-DIAGONAL TERMS OF MATRICES
C            IN THE CASE OF EDGE-BASED STORAGE
C
C            CASE OF LINEAR-LINEAR PRISM
C            AND SYMMETRIC MATRIX
C
C            LOCAL NUMBERING OF SEGMENTS CHOSEN HERE IN A PRISM
C
C            01 : POINT 1 TO 2
C            02 : POINT 2 TO 3
C            03 : POINT 3 TO 1
C            04 : POINT 4 TO 5
C            05 : POINT 5 TO 6
C            06 : POINT 6 TO 4
C            07 : POINT 1 TO 4
C            08 : POINT 2 TO 5
C            09 : POINT 3 TO 6
C            10 : POINT 1 TO 5
C            11 : POINT 2 TO 4
C            12 : POINT 2 TO 6
C            13 : POINT 3 TO 5
C            14 : POINT 3 TO 4
C            15 : POINT 1 TO 6
C
C            LOCAL NUMBERING OF ELEMENT BY ELEMENT EXTRA-DIAGONAL TERMS
C
C            01 : POINTS 1-2
C            02 : POINTS 1-3
C            03 : POINTS 1-4
C            04 : POINTS 1-5
C            05 : POINTS 1-6
C            06 : POINTS 2-3
C            07 : POINTS 2-4
C            08 : POINTS 2-5
C            09 : POINTS 2-6
C            10 : POINTS 3-4
C            11 : POINTS 3-5
C            12 : POINTS 3-6
C            13 : POINTS 4-5
C            14 : POINTS 4-6
C            15 : POINTS 5-6
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !  XM            !<-- ! TERMES EXTRA-DIAGONAUX ASSEMBLES XA12,23,31
C !  NSEG1         ! -->! NOMBRE DE SEGMENTS
C !  XMT           ! -->! TERMES EXTRA-DIAGONAUX
C !  NELMAX        ! -->! PREMIERE DIMENSION DE IKLE ET W.
C !                !    ! (CAS D'UN MAILLAGE ADAPTATIF)
C !  NELEM         ! -->! NOMBRE D'ELEMENTS DANS LE MAILLAGE.
C !  ELTSEG        ! -->! LISTE DES ELEMENTS DE CHAQUE SEGMENT
C !________________!____!_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELANT : ASSVEC
C
C***********************************************************************
C
      USE BIEF !, EX_AS3_4141_S => AS3_4141_S
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(IN)    :: NELMAX,NELEM,NSEG1
      INTEGER         , INTENT(IN)    :: DIM1XMT,DIM2XMT,STOXMT
      INTEGER         , INTENT(IN)    :: ELTSEG(NELMAX,15)
      DOUBLE PRECISION, INTENT(IN)    :: XMT(DIM1XMT,DIM2XMT)
      DOUBLE PRECISION, INTENT(INOUT) :: XM(NSEG1)
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
        XM(ISEG) = 0.D0
      ENDDO
C
C-----------------------------------------------------------------------
C
C  ASSEMBLAGE
C
C-----------------------------------------------------------------------
C
      IF(STOXMT.EQ.1) THEN
C
      DO IELEM = 1,NELEM
        XM(ELTSEG(IELEM,01)) = XM(ELTSEG(IELEM,01)) + XMT(IELEM,01)
        XM(ELTSEG(IELEM,02)) = XM(ELTSEG(IELEM,02)) + XMT(IELEM,06)
        XM(ELTSEG(IELEM,03)) = XM(ELTSEG(IELEM,03)) + XMT(IELEM,02)
        XM(ELTSEG(IELEM,04)) = XM(ELTSEG(IELEM,04)) + XMT(IELEM,13)
        XM(ELTSEG(IELEM,05)) = XM(ELTSEG(IELEM,05)) + XMT(IELEM,15)
        XM(ELTSEG(IELEM,06)) = XM(ELTSEG(IELEM,06)) + XMT(IELEM,14)
        XM(ELTSEG(IELEM,07)) = XM(ELTSEG(IELEM,07)) + XMT(IELEM,03)
        XM(ELTSEG(IELEM,08)) = XM(ELTSEG(IELEM,08)) + XMT(IELEM,08)
        XM(ELTSEG(IELEM,09)) = XM(ELTSEG(IELEM,09)) + XMT(IELEM,12)
        XM(ELTSEG(IELEM,10)) = XM(ELTSEG(IELEM,10)) + XMT(IELEM,04)
        XM(ELTSEG(IELEM,11)) = XM(ELTSEG(IELEM,11)) + XMT(IELEM,07)
        XM(ELTSEG(IELEM,12)) = XM(ELTSEG(IELEM,12)) + XMT(IELEM,09)
        XM(ELTSEG(IELEM,13)) = XM(ELTSEG(IELEM,13)) + XMT(IELEM,11)
        XM(ELTSEG(IELEM,14)) = XM(ELTSEG(IELEM,14)) + XMT(IELEM,10)
        XM(ELTSEG(IELEM,15)) = XM(ELTSEG(IELEM,15)) + XMT(IELEM,05)
      ENDDO
C
      ELSEIF(STOXMT.EQ.2) THEN
C
      DO IELEM = 1,NELEM
        XM(ELTSEG(IELEM,01)) = XM(ELTSEG(IELEM,01)) + XMT(01,IELEM)
        XM(ELTSEG(IELEM,02)) = XM(ELTSEG(IELEM,02)) + XMT(06,IELEM)
        XM(ELTSEG(IELEM,03)) = XM(ELTSEG(IELEM,03)) + XMT(02,IELEM)
        XM(ELTSEG(IELEM,04)) = XM(ELTSEG(IELEM,04)) + XMT(13,IELEM)
        XM(ELTSEG(IELEM,05)) = XM(ELTSEG(IELEM,05)) + XMT(15,IELEM)
        XM(ELTSEG(IELEM,06)) = XM(ELTSEG(IELEM,06)) + XMT(14,IELEM)
        XM(ELTSEG(IELEM,07)) = XM(ELTSEG(IELEM,07)) + XMT(03,IELEM)
        XM(ELTSEG(IELEM,08)) = XM(ELTSEG(IELEM,08)) + XMT(08,IELEM)
        XM(ELTSEG(IELEM,09)) = XM(ELTSEG(IELEM,09)) + XMT(12,IELEM)
        XM(ELTSEG(IELEM,10)) = XM(ELTSEG(IELEM,10)) + XMT(04,IELEM)
        XM(ELTSEG(IELEM,11)) = XM(ELTSEG(IELEM,11)) + XMT(07,IELEM)
        XM(ELTSEG(IELEM,12)) = XM(ELTSEG(IELEM,12)) + XMT(09,IELEM)
        XM(ELTSEG(IELEM,13)) = XM(ELTSEG(IELEM,13)) + XMT(11,IELEM)
        XM(ELTSEG(IELEM,14)) = XM(ELTSEG(IELEM,14)) + XMT(10,IELEM)
        XM(ELTSEG(IELEM,15)) = XM(ELTSEG(IELEM,15)) + XMT(05,IELEM)
      ENDDO
C
      ELSE
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'AS3_4141_S : STOCKAGE DE XMT INCONNU : ',STOXMT
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'AS3_4141_S: UNKNOWN STORAGE OF XMT : ',STOXMT
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END