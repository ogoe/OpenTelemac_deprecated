!                    ***************************
                     SUBROUTINE FLUSEC_TELEMAC2D
!                    ***************************
!
     &(U,V,H,IKLE,XEL,YEL,NELMAX,NELEM,X,Y,DT,NCP,CTRLSC,INFO,TPS,
     & KNOGL,MSKSEC,BM1,BM2,T1,HPROP,MESH,S,CV1,IFABOR,COMFLU,CUMFLO)
!
!***********************************************************************
! TELEMAC2D   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    COMPUTES FLUXES THROUGH CONTROL SECTIONS
!+                AND SUMS THEM UP TO REACH OSCILLATING VOLUMES.
!+
!+            GRIDS OF DIMENSION 2 AND WATER DEPTH CONSIDERED.
!
!history  J-M HERVOUET (LNHE)
!+        14/01/2008
!+        V5P8
!+
!
!history  JAJ, JACEK.JANKOWSKI@BAW.DE
!+        16/02/2010
!+        V6P0
!+   MODIFIED
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
!| BM1            |<->| WORK MATRIX
!| BM2            |<->| WORK MATRIX
!| COMFLU         |-->| KEYWORD: COMPATIBLE COMPUTATION OF FLUXES
!| CTRLSC         |-->| DATA ON CONTROL SECTIONS
!| CUMFLO         |-->| KEYWORD: PRINTING CUMULATED FLOWRATES
!| CV1            |<->| WORK ARRAY IN A BIEF_OBJ STRUCTURE
!| DT             |-->| TIME STEP
!| H              |-->| WATER DEPTH
!| HPROP          |-->| PROPAGATION DEPTH
!| IFABOR         |-->| ELEMENTS BEHIND THE EDGES OF A TRIANGLE
!|                |   | IF NEGATIVE OR ZERO, THE EDGE IS A LIQUID
!|                |   | BOUNDARY
!| IKLE           |-->| CONNECTIVITY TABLE
!| INFO           |-->| IF YES, PRINT
!| KNOGL          |-->| FROM GLOBAL TO LOCAL NUMBERING IN PARALLEL
!| MESH           |-->| MESH STRUCTURE
!| MSKSEC         |<->| BLOCK OF WORK ARRAYS, PER SECTION
!| NCP            |-->| TWO TIMES THE NUMBER OF CONTROL SECTIONS
!| NELEM          |-->| NUMBER OF ELEMENTS
!| NELMAX         |-->| MAXIMUM NUMBER OF ELEMENTS
!| S              |-->| VOID STRUCTURE
!| T1             |<->| WOR kARRAY IN A BIEF_OBJ STRUCTURE
!| TPS            |-->| TIME IN SECONDS
!| U,V            |-->| VELOCITY FIELD COMPONENTS
!| X              |-->| ABSCISSAE OF POINTS IN THE MESH
!| XEL            |-->| ABSCISSAE OF POINTS IN THE MESH, PER ELEMENT
!| Y              |-->| ORDINATES OF POINTS IN THE MESH
!| YEL            |-->| ORDINATES OF POINTS IN THE MESH, PER ELEMENT
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_TELEMAC2D, ONLY: CHAIN
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN) :: NELMAX,NELEM
      INTEGER, INTENT(INOUT) :: NCP ! CAN BE CHANGED
      INTEGER, INTENT(INOUT) :: CTRLSC(NCP) ! CAN BE CHANGED
      INTEGER, INTENT(IN) :: IKLE(NELMAX,*),KNOGL(*)
      INTEGER, INTENT(IN) :: IFABOR(NELMAX,*)
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*),TPS,DT
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,*),YEL(NELMAX,*)
      LOGICAL, INTENT(IN) :: INFO,COMFLU,CUMFLO
      TYPE(BIEF_OBJ), INTENT(IN)    :: HPROP,S,U,V,H
      TYPE(BIEF_OBJ), INTENT(INOUT) :: BM1,BM2,T1,MSKSEC,CV1
      TYPE(BIEF_MESH) :: MESH
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER NSEMAX,ERR
      PARAMETER(NSEMAX=50)
!
      INTEGER IELEM,IEL,I1,I2,I3,J1,J2,J3,ELBEST,IGBEST,ILBEST,IELMH
      INTEGER ILPREC,ISEG,ISEC,NSEC,PT,DEP,ARR,IELMU
!
      DOUBLE PRECISION DIST,DIST1,DIST2,DIST3
      DOUBLE PRECISION H1,H2,X1,Y1,X2,Y2,UN1,UN2,NX,NY,SUR6
!
      LOGICAL OK
!
      DOUBLE PRECISION, ALLOCATABLE :: FLX(:),VOLNEG(:),VOLPOS(:)
      INTEGER, ALLOCATABLE :: NSEG(:),LISTE(:,:,:)
!
      LOGICAL DEJA
      DATA DEJA/.FALSE./
      LOGICAL :: OLD_METHOD=.FALSE.
!
      SAVE LISTE,DEJA,NSEG,VOLNEG,VOLPOS,FLX,OLD_METHOD
!
!-----------------------------------------------------------------------
!
      SUR6=1.D0/6.D0
      NSEC = NCP/2
!
!  SEARCHES PATHS JOINING POINT PAIRS:
!
      IF(.NOT.DEJA) THEN
!
!     DYNAMIC ALLOCATION OF FLX, VOLNEG, VOLPOS
!
      ALLOCATE(FLX(NCP)             ,STAT=ERR)
      ALLOCATE(VOLNEG(NCP)          ,STAT=ERR)
      ALLOCATE(VOLPOS(NCP)          ,STAT=ERR)
      ALLOCATE(NSEG(NCP)            ,STAT=ERR)
      ALLOCATE(LISTE(NCP,NSEMAX,2)  ,STAT=ERR)
!
      IF(ERR.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,100) ERR
        IF(LNG.EQ.2) WRITE(LU,200) ERR
100     FORMAT(1X,'FLUSEC : ERREUR A L''ALLOCATION DE MEMOIRE : ',/,1X,
     &            'CODE D''ERREUR : ',1I6)
200     FORMAT(1X,'FLUSEC: ERROR DURING ALLOCATION OF MEMORY: ',/,1X,
     &            'ERROR CODE: ',1I6)
      ENDIF
!
      IF (.NOT.ALLOCATED(CHAIN)) OLD_METHOD=.TRUE.
      DO 60 ISEC =1,NSEC
!!JAJ #### IN THE SERIAL CASE, OR "CLASSICAL" IN PARALLEL,
!     FOLLOW THE ALGORITHM OF FINDING SEGMENT CHAINS
!
! NOTE: IF YOU CHANGE THE ALGORITHM, CHANGE IT IN PARTEL AS WELL
!
        IF (NCSIZE.LE.1 .OR. OLD_METHOD) THEN
!
        DEP = CTRLSC(1+2*(ISEC-1))
        ARR = CTRLSC(2+2*(ISEC-1))
        VOLNEG(ISEC)=0.D0
        VOLPOS(ISEC)=0.D0
        IF(NCSIZE.GT.1) THEN
          DEP=KNOGL(DEP)
          ARR=KNOGL(ARR)
          IF(DEP.EQ.0.AND.ARR.EQ.0) THEN
            NSEG(ISEC)=0
            GO TO 60
          ENDIF
          IF((DEP.EQ.0.AND.ARR.NE.0).OR.(DEP.NE.0.AND.ARR.EQ.0)) THEN
            NSEG(ISEC)=-1
            GO TO 60
          ENDIF
        ENDIF
        PT = DEP
        ISEG = 0
        DIST=(X(DEP)-X(ARR))**2+(Y(DEP)-Y(ARR))**2
10      CONTINUE
!
        DO 20 IELEM =1,NELEM
!
          I1 = IKLE(IELEM,1)
          I2 = IKLE(IELEM,2)
          I3 = IKLE(IELEM,3)
!         IF THE ELEMENT CONTAINS THE SELECTED NODE:
          IF(PT.EQ.I1.OR.PT.EQ.I2.OR.PT.EQ.I3) THEN
            DIST1 = (X(I1)-X(ARR))**2 + (Y(I1)-Y(ARR))**2
            DIST2 = (X(I2)-X(ARR))**2 + (Y(I2)-Y(ARR))**2
            DIST3 = (X(I3)-X(ARR))**2 + (Y(I3)-Y(ARR))**2
            IF(DIST1.LT.DIST) THEN
              DIST = DIST1
              ELBEST = IELEM
              IGBEST = I1
              ILBEST = 1
              IF(I1.EQ.PT) ILPREC = 1
              IF(I2.EQ.PT) ILPREC = 2
              IF(I3.EQ.PT) ILPREC = 3
            ENDIF
            IF(DIST2.LT.DIST) THEN
              DIST = DIST2
              ELBEST = IELEM
              IGBEST = I2
              ILBEST = 2
              IF(I1.EQ.PT) ILPREC = 1
              IF(I2.EQ.PT) ILPREC = 2
              IF(I3.EQ.PT) ILPREC = 3
            ENDIF
            IF(DIST3.LT.DIST) THEN
              DIST = DIST3
              ELBEST = IELEM
              IGBEST = I3
              ILBEST = 3
              IF(I1.EQ.PT) ILPREC = 1
              IF(I2.EQ.PT) ILPREC = 2
              IF(I3.EQ.PT) ILPREC = 3
            ENDIF
          ENDIF
!
20      CONTINUE
!
        IF(IGBEST.EQ.PT) THEN
          IF(LNG.EQ.1) WRITE(LU,32)
          IF(LNG.EQ.2) WRITE(LU,33)
32        FORMAT(1X,'FLUSEC_TELEMAC2D : BLOCAGE DE L''ALGORITHME')
33        FORMAT(1X,'FLUSEC_TELEMAC2D : ALGORITHM FAILED')
          CALL PLANTE(1)
          STOP
        ELSE
          PT = IGBEST
          ISEG = ISEG + 1
          IF(ISEG.GT.NSEMAX) THEN
            IF(LNG.EQ.1) THEN
              WRITE(LU,*) 'FLUSEC_TELEMAC2D : TROP DE SEGMENTS DANS UNE'
              WRITE(LU,*) '                   SECTION. AUGMENTER NSEMAX'
            ENDIF
            IF(LNG.EQ.2) THEN
              WRITE(LU,*) 'FLUSEC_TELEMAC2D: TOO MANY SEGMENTS IN A   '
              WRITE(LU,*) '                  SECTION. INCREASE  NSEMAX'
            ENDIF
            STOP
          ENDIF
          LISTE(ISEC,ISEG,1) = IKLE(ELBEST,ILPREC)
          LISTE(ISEC,ISEG,2) = IKLE(ELBEST,ILBEST)
          IF(IGBEST.NE.ARR) GO TO 10
        ENDIF
!
        NSEG(ISEC) = ISEG
!
!JAJ #### THIS PART TO BE DONE IN THE PARALLEL CASE; FILL LISTE
! WITH READY SEGMENT CHAINS PROVIDED BY PARTEL: SEE READ_SECTIONS
! NOTE: FUTURE OPTIMISATION - USE CHAIN STRUCTURE IN THE WHOLE ROUTINE
!
        ELSE
!
          FLX(ISEC)=0.0D0
          VOLNEG(ISEC)=0.0D0
          VOLPOS(ISEC)=0.0D0
          NSEG(ISEC) = CHAIN(ISEC)%NSEG
          LISTE(ISEC,:,:)=0
          DO ISEG=1,NSEG(ISEC)
            LISTE(ISEC,ISEG,1) = CHAIN(ISEC)%LISTE(ISEG,1)
            LISTE(ISEC,ISEG,2) = CHAIN(ISEC)%LISTE(ISEG,2)
          END DO
        ENDIF
!
!       IF COMPATIBLE COMPUTATION OF FLUXES=YES
        IF(COMFLU.AND.NSEG(ISEC).GE.1) THEN
!
!       LOOKING AT ALL ELEMENTS TOUCHING THE SECTION WITH 2 POINTS
!       MARKING THEM WITH 1 ON ONE SIDE AND -1 ON THE OTHER SIDE
!
        DO IEL=1,NELEM
          MSKSEC%ADR(ISEC)%P%R(IEL)=0.D0
          J1=IKLE(IEL,1)
          J2=IKLE(IEL,2)
          J3=IKLE(IEL,3)
          DO ISEG=1,NSEG(ISEC)
            I1 = LISTE(ISEC,ISEG,1)
            I2 = LISTE(ISEC,ISEG,2)
!           LEFT SIDE
            IF    ( (J1.EQ.I1.AND.J2.EQ.I2) .OR.
     &              (J2.EQ.I1.AND.J3.EQ.I2) .OR.
     &              (J3.EQ.I1.AND.J1.EQ.I2)      ) THEN
              MSKSEC%ADR(ISEC)%P%R(IEL)=1.D0
!           RIGHT SIDE
            ELSEIF( (J1.EQ.I2.AND.J2.EQ.I1) .OR.
     &              (J2.EQ.I2.AND.J3.EQ.I1) .OR.
     &              (J3.EQ.I2.AND.J1.EQ.I1)      ) THEN
              MSKSEC%ADR(ISEC)%P%R(IEL)=-1.D0
            ENDIF
          ENDDO
        ENDDO
!
!       OTHER TRIANGLES WITH ONLY 1 POINT TOUCHING THE SECTION
!       LOOKING AT NEIGHBOURS TO FIND THEIR SIDE
!
999     CONTINUE
        OK=.TRUE.
        DO IEL=1,NELEM
          J1=IKLE(IEL,1)
          J2=IKLE(IEL,2)
          J3=IKLE(IEL,3)
          DO ISEG=1,NSEG(ISEC)
            I1 = LISTE(ISEC,ISEG,1)
            I2 = LISTE(ISEC,ISEG,2)
            IF((J1.EQ.I1.OR.J2.EQ.I1.OR.J3.EQ.I1.OR.
     &          J1.EQ.I2.OR.J2.EQ.I2.OR.J3.EQ.I2).AND.
     &          ABS(MSKSEC%ADR(ISEC)%P%R(IEL)).LT.0.5D0) THEN
!             LOOKING AT NEIGHBOURS
              IF(IFABOR(IEL,1).GT.0) THEN
                IELEM=IFABOR(IEL,1)
                IF(ABS(MSKSEC%ADR(ISEC)%P%R(IELEM)).GT.0.5D0) THEN
                  MSKSEC%ADR(ISEC)%P%R(IEL)=MSKSEC%ADR(ISEC)%P%R(IELEM)
                ENDIF
              ENDIF
              IF(IFABOR(IEL,2).GT.0) THEN
                IELEM=IFABOR(IEL,2)
                IF(ABS(MSKSEC%ADR(ISEC)%P%R(IELEM)).GT.0.5D0) THEN
                  MSKSEC%ADR(ISEC)%P%R(IEL)=MSKSEC%ADR(ISEC)%P%R(IELEM)
                ENDIF
              ENDIF
              IF(IFABOR(IEL,3).GT.0) THEN
                IELEM=IFABOR(IEL,3)
                IF(ABS(MSKSEC%ADR(ISEC)%P%R(IELEM)).GT.0.5D0) THEN
                  MSKSEC%ADR(ISEC)%P%R(IEL)=MSKSEC%ADR(ISEC)%P%R(IELEM)
                ENDIF
              ENDIF
!             CORRECTION JMH 11/01/2008
!             IF(MSKSEC%ADR(ISEC)%P%R(IEL).LT.0.5D0) OK=.FALSE.
              IF(ABS(MSKSEC%ADR(ISEC)%P%R(IEL)).LT.0.5D0) OK=.FALSE.
            ENDIF
          ENDDO
        ENDDO
        IF(.NOT.OK) GO TO 999
!
!       IF(COMFLU.AND.NSEG(ISEC).GE.1) THEN
        ENDIF
!
60    CONTINUE
!
!     IF(.NOT.DEJA) THEN
      ENDIF
!
!-----------------------------------------------------------------------
!
      DEJA = .TRUE.
!
!-----------------------------------------------------------------------
!
      IELMH=HPROP%ELM
      IELMU=U%ELM
!
      DO 70 ISEC = 1 , NSEC
      FLX(ISEC) = 0.D0
      IF(NSEG(ISEC).GE.1) THEN
!
      IF(COMFLU) THEN
!
!     COMPUTING THE FLUX AS IN THE CONTINUITY EQUATION
!     (HOWEVER IMPLICITATION SHOULD PERHAPS BE ALSO USED)
!
      CALL MATRIX(BM1,'M=N     ','MATFGR         X',IELMH,IELMU,
     &            1.D0,HPROP,S,S,S,S,S,MESH,.TRUE.,MSKSEC%ADR(ISEC)%P)
      CALL MATRIX(BM2,'M=N     ','MATFGR         Y',IELMH,IELMU,
     &            1.D0,HPROP,S,S,S,S,S,MESH,.TRUE.,MSKSEC%ADR(ISEC)%P)
!
      CALL MATVEC( 'X=AY    ',CV1,BM1,U,0.D0,MESH)
      CALL MATVEC( 'X=X+AY  ',CV1,BM2,V,0.D0,MESH)
!
!     SUMMING CV1 FOR ALL POINTS OF THE SECTION, THIS IS THE FLUX !
!     (OBTAINED BY CONTINUITY EQUATION AND AN INTEGRATION BY PARTS)
!
      DO ISEG = 1 , NSEG(ISEC)
        I1   = LISTE(ISEC,ISEG,1)
        FLX(ISEC) = FLX(ISEC) + CV1%R(I1)
      ENDDO
!     LAST SEGMENT, ADDING THE LAST POINT
      I2   = LISTE(ISEC,NSEG(ISEC),2)
      FLX(ISEC) = FLX(ISEC) + CV1%R(I2)
!
!     WHEN BOTH UPWIND AND DOWNSTREAM ELEMENTS ARE TAKEN INTO ACCOUNT
!     WITH DIFFERENT SIGNS, WE GET TWICE THE FLUX
!
      FLX(ISEC)=FLX(ISEC)*0.5D0
!
      ELSE
!
!       COMPUTING THE FLUX DIRECTLY, REGARDLESS OF THE WEAK FORM
!       OF THE IMPERMEABILITY CONDITION
!
        DO ISEG = 1 , NSEG(ISEC)
          I1 = LISTE(ISEC,ISEG,1)
          I2 = LISTE(ISEC,ISEG,2)
          X1 = X(I1)
          X2 = X(I2)
          Y1 = Y(I1)
          Y2 = Y(I2)
          H1 = H%R(I1)
          H2 = H%R(I2)
          NX = Y1-Y2
          NY = X2-X1
          UN1= U%R(I1)*NX + V%R(I1)*NY
          UN2= U%R(I2)*NX + V%R(I2)*NY
          FLX(ISEC) = FLX(ISEC) + ((H1+H2)*(UN1+UN2)+H2*UN2+H1*UN1)*SUR6
        ENDDO
!
      ENDIF
!
      IF(FLX(ISEC).GT.0.D0) THEN
        VOLPOS(ISEC) = VOLPOS(ISEC) + FLX(ISEC)*DT
      ELSE
        VOLNEG(ISEC) = VOLNEG(ISEC) + FLX(ISEC)*DT
      ENDIF
!
!     IF(NSEG(ISEC).GT.1)...
      ENDIF
!
70    CONTINUE
!
!-----------------------------------------------------------------------
!
!     PRINTING THE RESULTS / !JAJ HERE ALLREDUCES FOR VALUES
!
      CALL FLUXPR_TELEMAC2D
     &  (NSEC,CTRLSC,FLX,VOLNEG,VOLPOS,INFO,TPS,
     &   NSEG,NCSIZE,CUMFLO)
!
!-----------------------------------------------------------------------
!
      RETURN
      END SUBROUTINE FLUSEC_TELEMAC2D
