!                    *******************
                     SUBROUTINE STOSEG41
!                    *******************
!
     &(IFABOR,NELEM,NELMAX,IELM,IKLE,NBOR,NPTFR,
     & GLOSEG,MAXSEG,ELTSEG,ORISEG,NSEG,KP1BOR,NELBOR,NULONE,NELMAX2,
     & NELEM2,NPTFR2,NPOIN2,NPLAN,KNOLG,NSEG2D)
!
!***********************************************************************
! BIEF   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    BUILDS THE DATA STRUCTURE FOR EDGE-BASED STORAGE
!+                IN PRISMS.
!code
!+       LOCAL NUMBERING OF SEGMENTS CHOSEN HERE IN A PRISM
!+
!+       HORIZONTAL
!+
!+       01 : POINT 1 TO 2 (OR THE OPPOSITE DEPENDING OF ORISEG)
!+       02 : POINT 2 TO 3 (OR THE OPPOSITE DEPENDING OF ORISEG)
!+       03 : POINT 3 TO 1 (OR THE OPPOSITE DEPENDING OF ORISEG)
!+       04 : POINT 4 TO 5 (OR THE OPPOSITE DEPENDING OF ORISEG)
!+       05 : POINT 5 TO 6 (OR THE OPPOSITE DEPENDING OF ORISEG)
!+       06 : POINT 6 TO 4 (OR THE OPPOSITE DEPENDING OF ORISEG)
!+
!+       VERTICAL
!+
!+       07 : POINT 1 TO 4
!+       08 : POINT 2 TO 5
!+       09 : POINT 3 TO 6
!+
!+       CROSSED
!+
!+       10 : POINT 1 TO 5
!+       11 : POINT 2 TO 4
!+       12 : POINT 2 TO 6
!+       13 : POINT 3 TO 5
!+       14 : POINT 3 TO 4
!+       15 : POINT 1 TO 6
!
!history  JMH
!+        11/08/09
!+
!+   CROSSED AND VERTICAL SEGMENTS INVERTED IN THE NUMBERING.
!
!history  JMH
!+        16/10/09
!+
!+   NUMBERING OF CROSSED SEGMENTS CHANGED, TO FACILITATE
!
!history  J-M HERVOUET (LNHE)
!+        19/10/2009
!+        V6P0
!+
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| ELTSEG         |<--| SEGMENTS OF EVERY TRIANGLE.
!| GLOSEG         |<--| GLOBAL NUMBERS OF POINTS OF SEGMENTS.
!| IELM           |-->| 11: TRIANGLES.
!|                |   | 21: QUADRILATERALS.
!| IFABOR         |-->| ELEMENTS BEHIND THE EDGES OF A TRIANGLE
!|                |   | IF NEGATIVE OR ZERO, THE EDGE IS A LIQUID
!|                |   | BOUNDARY
!| IKLE           |-->| CONNECTIVITY TABLE.
!| KNOLG          |-->| GLOBAL NUMBER OF A LOCAL POINT IN PARALLEL
!| KP1BOR         |-->| GIVES THE NEXT BOUNDARY POINT IN A CONTOUR
!| MAXSEG         |<--| MAXIMUM NUMBER OF SEGMENTS
!| NBOR           |-->| GLOBAL NUMBERS OF BOUNDARY POINTS.
!| NELBOR         |-->| NUMBER OF ELEMENT CONTAINING SEGMENT K OF
!|                |   | THE BOUNDARY.
!| NELEM          |-->| NUMBER OF ELEMENTS IN THE MESH
!| NELEM2         |-->| NUMBER OF ELEMENTS IN 2D
!| NELMAX         |-->| MAXIMUM NUMBER OF ELEMENTS IN 3D
!| NELMAX2        |-->| MAXIMUM NUMBER OF ELEMENTS IN 2D
!| NPLAN          |-->| NUMBER OF PLANES IN THE 3D MESH OF PRISMS
!| NPOIN2         |-->| NUMBER OF POINTS IN 2D
!| NPTFR          |-->| NUMBER OF BOUNDARY POINTS
!| NPTFR2         |-->| NUMBER OF BOUNDARY POINTS IN 2D
!| NSEG           |<--| NUMBER OF SEGMENTS OF THE MESH.
!| NULONE         |-->| LOCAL NUMBER OF BOUNDARY POINTS IN A BOUNDARY
!|                |   | ELEMENT.
!| ORISEG         |<--| ORIENTATION OF SEGMENTS OF EVERY TRIANGLE.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF, EX_STOSEG41 => STOSEG41
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)    :: NELMAX,NELMAX2,NPTFR,NSEG,MAXSEG,IELM
      INTEGER, INTENT(IN)    :: NELEM,NELEM2,NPTFR2,NPOIN2,NPLAN,NSEG2D
      INTEGER, INTENT(IN)    :: NBOR(NPTFR2),KP1BOR(NPTFR2)
      INTEGER, INTENT(IN)    :: IFABOR(NELMAX2,*),IKLE(NELMAX,6)
      INTEGER, INTENT(IN)    :: NELBOR(NPTFR2),NULONE(NPTFR2)
      INTEGER, INTENT(INOUT) :: GLOSEG(MAXSEG,2)
      INTEGER, INTENT(INOUT) :: ELTSEG(NELMAX,15),ORISEG(NELMAX,15)
      INTEGER, INTENT(IN)    :: KNOLG(*)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER I1,I2,I3,IELEM,I,ONE,TWO
      INTEGER IPLAN,ISEG2D,ISEG3D,IELEM3D,NSEGH,NSEGV
!
!-----------------------------------------------------------------------
!
      IF(IELM.NE.41) THEN
        IF (LNG.EQ.1) WRITE(LU,500) IELM
        IF (LNG.EQ.2) WRITE(LU,501) IELM
500     FORMAT(1X,'STOSEG41 (BIEF) : ELEMENT NON PREVU : ',1I6)
501     FORMAT(1X,'STOSEG41 (BIEF) : UNEXPECTED ELEMENT: ',1I6)
        CALL PLANTE(1)
        STOP
      ENDIF
!
!-----------------------------------------------------------------------
!
!     BUILDS 2D SEGMENTS (THE FIRST IN THE NUMBERING)
!
!     NOTE:
!     IN PARALLEL MODE, THE NUMBER OF SEGMENTS IN THE CONTOUR LINE IS NOT NPTFR2
!     NSEG2D=(3*NELEM2+NPTFR2)/2
      NSEGH=NSEG2D*NPLAN
      NSEGV=(NPLAN-1)*NPOIN2
!
      CALL STOSEG(IFABOR,NELEM2,NELMAX,NELMAX2,11,IKLE,NBOR,NPTFR2,
     &            GLOSEG,MAXSEG,ELTSEG,ORISEG,NSEG2D,
     &            KP1BOR,NELBOR,NULONE,KNOLG)
!
!-----------------------------------------------------------------------
!
!     COMPLETES HORIZONTAL SEGMENTS (1,2,3,4,5,6)
!
      DO IPLAN=2,NPLAN
      DO ISEG2D=1,NSEG2D
        ISEG3D=ISEG2D+(IPLAN-1)*NSEG2D
        GLOSEG(ISEG3D,1)=GLOSEG(ISEG2D,1)+NPOIN2*(IPLAN-1)
        GLOSEG(ISEG3D,2)=GLOSEG(ISEG2D,2)+NPOIN2*(IPLAN-1)
      ENDDO
      ENDDO
!
!     VERTICAL SEGMENTS (7,8,9)
!
      DO IPLAN=1,NPLAN-1
      DO I=1,NPOIN2
        ISEG3D=NSEGH+NPOIN2*(IPLAN-1)+I
        GLOSEG(ISEG3D,1)=NPOIN2*(IPLAN-1)+I
        GLOSEG(ISEG3D,2)=NPOIN2*(IPLAN  )+I
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
!     ARRAY ELTSEG GIVES GLOBAL NUMBERS OF SEGMENTS IN A PRISM
!     ARRAY ORISEG GIVES ORIENTATION OF SEGMENT
!
!     COMPLETES SEGMENTS 1,2,3
!
      IF(NPLAN.GT.2) THEN
      DO IPLAN=2,NPLAN-1
        DO IELEM=1,NELEM2
          IELEM3D=IELEM+(IPLAN-1)*NELEM2
          ELTSEG(IELEM3D,1)=ELTSEG(IELEM,1)+NSEG2D*(IPLAN-1)
          ELTSEG(IELEM3D,2)=ELTSEG(IELEM,2)+NSEG2D*(IPLAN-1)
          ELTSEG(IELEM3D,3)=ELTSEG(IELEM,3)+NSEG2D*(IPLAN-1)
          ORISEG(IELEM3D,1)=ORISEG(IELEM,1)
          ORISEG(IELEM3D,2)=ORISEG(IELEM,2)
          ORISEG(IELEM3D,3)=ORISEG(IELEM,3)
        ENDDO
      ENDDO
      ENDIF
!
!     SEGMENTS 4,5,6 (=SEGMENTS 1,2,3 + NSEG2D)
!
      DO IPLAN=1,NPLAN-1
        DO IELEM=1,NELEM2
          IELEM3D=IELEM+(IPLAN-1)*NELEM2
          ELTSEG(IELEM3D,4)=ELTSEG(IELEM3D,1)+NSEG2D
          ELTSEG(IELEM3D,5)=ELTSEG(IELEM3D,2)+NSEG2D
          ELTSEG(IELEM3D,6)=ELTSEG(IELEM3D,3)+NSEG2D
          ORISEG(IELEM3D,4)=ORISEG(IELEM,1)
          ORISEG(IELEM3D,5)=ORISEG(IELEM,2)
          ORISEG(IELEM3D,6)=ORISEG(IELEM,3)
        ENDDO
      ENDDO
!
!     SEGMENTS 7,8,9 (VERTICAL SEGMENTS)  SEE SECTION 3)
!
      DO IPLAN=1,NPLAN-1
        DO IELEM=1,NELEM2
          IELEM3D=IELEM+(IPLAN-1)*NELEM2
          I1=IKLE(IELEM,1)
          I2=IKLE(IELEM,2)
          I3=IKLE(IELEM,3)
          ELTSEG(IELEM3D,7)=NSEGH+NPOIN2*(IPLAN-1)+I1
          ELTSEG(IELEM3D,8)=NSEGH+NPOIN2*(IPLAN-1)+I2
          ELTSEG(IELEM3D,9)=NSEGH+NPOIN2*(IPLAN-1)+I3
          ORISEG(IELEM3D,7)=1
          ORISEG(IELEM3D,8)=1
          ORISEG(IELEM3D,9)=1
        ENDDO
      ENDDO
!
!     SEGMENTS 10 TO 15 (CROSSED SEGMENTS)
!
!     THE PROBLEM IS TO FIND A GLOBAL NUMBERING OF CROSSED SEGMENTS
!     WHICH IS INDEPENDENT OF THE PRISM CONSIDERED, SO THAT WHEN
!     ASSEMBLING THE ISEG3D IS THE CORRECT ONE
!     HERE WE TAKE AS CROSSED SEGMENT 1 THE ONE THAT STARTS
!     FROM THE ORIGIN OF THE ORIENTED HORIZONTAL SEGMENT
!
!     FOR EVERY VERTICAL RECTANGLE IN A PRISM
!     WHEN THE HORIZONTAL SEGMENT HAS ORISEG = 1, THE CROSSED
!     SEGMENT NUMBER 1 IS THUS THE 10,12 OR 14
!     WHEN THE HORIZONTAL SEGMENT HAS ORISEG = 2, THE CROSSED
!     SEGMENT NUMBER 1 IS THUS THE 11,13 OR 15
!
!     IN THE GLOBAL NUMBERING OF SEGMENTS
!     THE ORISEG=1 SEGMENTS ARE PUT FIRST (NSEG2D PER LAYER)
!     THE ORISEG=2 SEGMENTS ARE PUT AFTER (NSEG2D PER LAYER)
!
      DO IPLAN=1,NPLAN-1
        DO IELEM=1,NELEM2
          IELEM3D=IELEM+(IPLAN-1)*NELEM2
          DO I=1,3
            ISEG2D=ELTSEG(IELEM,I)
            ISEG3D=NSEGH+NSEGV+NSEG2D*2*(IPLAN-1)+ISEG2D
            ONE=ISEG3D
            TWO=ISEG3D+NSEG2D
            IF(ORISEG(IELEM,I).EQ.1) THEN
              ELTSEG(IELEM3D,8+2*I)  =ONE
              ELTSEG(IELEM3D,8+2*I+1)=TWO
            ELSE
              ELTSEG(IELEM3D,8+2*I+1)=ONE
              ELTSEG(IELEM3D,8+2*I)  =TWO
            ENDIF
!           ONE: FOLLOWS ISEG2D
            GLOSEG(ONE,1)=GLOSEG(ISEG2D,1)+(IPLAN-1)*NPOIN2
            GLOSEG(ONE,2)=GLOSEG(ISEG2D,2)+(IPLAN  )*NPOIN2
!           BACKWARDS WITH RESPECT TO ISEG2D
            GLOSEG(TWO,1)=GLOSEG(ISEG2D,2)+(IPLAN-1)*NPOIN2
            GLOSEG(TWO,2)=GLOSEG(ISEG2D,1)+(IPLAN  )*NPOIN2
          ENDDO
          ORISEG(IELEM3D,10)=1
          ORISEG(IELEM3D,11)=1
          ORISEG(IELEM3D,12)=1
          ORISEG(IELEM3D,13)=1
          ORISEG(IELEM3D,14)=1
          ORISEG(IELEM3D,15)=1
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
      RETURN
      END
