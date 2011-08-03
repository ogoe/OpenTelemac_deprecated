!                    *****************
                     SUBROUTINE RESCJG
!                    *****************
!
     &(X, A,B , MESH,D,AD,AG,G,R, CFG,INFOGR,AUX)
!
!***********************************************************************
! BIEF   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    SOLVES THE LINEAR SYSTEM A X = B
!+                USING THE CONJUGATE RESIDUAL METHOD.
!code
!+-----------------------------------------------------------------------
!+                        PRECONDITIONING
!+-----------------------------------------------------------------------
!+    PRECON VALUE     I                  MEANING
!+-----------------------------------------------------------------------
!+                     I
!+        0 OR 1       I  NO PRECONDITIONING
!+                     I
!+        2            I  DIAGONAL PRECONDITIONING USING THE MATRIX
!+                     I  DIAGONAL
!+        3            I  BLOCK-DIAGONAL PRECONDITIONING
!+                     I
!+        5            I  DIAGONAL PRECONDITIONING USING THE ABSOLUTE
!+                     I  VALUE OF THE MATRIX DIAGONAL
!+                     I
!+        7            I  CROUT EBE PRECONDITIONING
!+                     I
!+       11            I  GAUSS-SEIDEL EBE PRECONDITIONING
!+                     I
!+-----------------------------------------------------------------------
!
!history  J-M HERVOUET (LNHE)
!+        27/02/06
!+        V5P6
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
!| A              |-->| MATRIX OF THE SYSTEM
!| AD             |<->| WORK ARRAY: MATRICE A MULTIPLIED BY D.
!| AG             |<->| WORK ARRAY: MATRICE A MULTIPLIED BY G.
!| AUX            |-->| MATRIX FOR PRECONDITIONING.
!| B              |-->| RIGHT-HAND SIDE OF THE SYSTEM
!| CFG            |-->| STRUCTURE OF SOLVER CONFIGURATION
!| D              |<->| WORK ARRAY: DIRECTION OF DESCENT.
!| G              |<->| DESCENT GRADIENT.
!| INFOGR         |-->| IF YES, PRINT A LOG.
!| MESH           |-->| MESH STRUCTURE.
!| R              |<->| RESIDUAL (MAY BE IN THE SAME MEMORY SPACE AS
!|                |   | GRADIENT DEPENDING ON CONDITIONING)
!| X              |<->| INITIAL VALUE, THEN SOLUTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF, EX_RESCJG => RESCJG
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      LOGICAL, INTENT(IN) :: INFOGR
!
!     STRUCTURES
!
      TYPE(BIEF_OBJ), INTENT(INOUT) :: D,AD,G,AG,R,X,B
      TYPE(SLVCFG)  , INTENT(INOUT) :: CFG
!
!     MESH STRUCTURE
!
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
!
!     MATRIX STRUCTURE
!
      TYPE(BIEF_OBJ), INTENT(IN)    :: A
      TYPE(BIEF_OBJ), INTENT(INOUT) :: AUX
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER M
!
      DOUBLE PRECISION XL,RMRM,RMDM,RMGM,TESTL,GAD
      DOUBLE PRECISION AGAD,BETA,ADAD,RO,DAD,RM1GM1,GMGM,GM1GM1,STO1
      DOUBLE PRECISION STO2,TGMTGM,C
!
      LOGICAL RELAT,PREC,CROUT,GSEB,PREBE,PRE3D
!
!-----------------------------------------------------------------------
!
      INTRINSIC SQRT
!
!-----------------------------------------------------------------------
!
!   INITIALISES
!     STO1 AND STO2 AVOID A WARNING MESSAGE FROM CRAY COMPILER
      STO1  =0.D0
      STO2  =0.D0
      TGMTGM=0.D0
      CROUT =.FALSE.
      IF(7*(CFG%PRECON/7).EQ.CFG%PRECON) CROUT=.TRUE.
      GSEB=.FALSE.
      IF(11*(CFG%PRECON/11).EQ.CFG%PRECON) GSEB=.TRUE.
      PREBE=.FALSE.
      IF(13*(CFG%PRECON/13).EQ.CFG%PRECON) PREBE=.TRUE.
      PRE3D=.FALSE.
      IF(17*(CFG%PRECON/17).EQ.CFG%PRECON) PRE3D=.TRUE.
      PREC=.FALSE.
      IF(CROUT.OR.GSEB.OR.PREBE.OR.PRE3D) PREC=.TRUE.
!
!-----------------------------------------------------------------------
!   INITIALISES
!-----------------------------------------------------------------------
!
      M   = 0
!
!  NORMALISES THE SECOND MEMBER TO COMPUTE THE RELATIVE PRECISION:
!
      XL = P_DOTS(B,B,MESH)
      IF (XL.LT.1.D0) THEN
         XL = 1.D0
         RELAT = .FALSE.
      ELSE
         RELAT = .TRUE.
      ENDIF
!
! COMPUTES THE INITIAL RESIDUAL AND POSSIBLY EXITS:
!
      CALL MATRBL( 'X=AY    ',R,A,X,  C,MESH)
!
      CALL OS( 'X=X-Y   ' , R , B , B , C )
      RMRM   = P_DOTS(R,R,MESH)
      RMGM   = RMRM
      GMGM   = RMRM
!
      IF (RMRM.LT.CFG%EPS**2*XL) GO TO 900
!
!-----------------------------------------------------------------------
! PRECONDITIONING :
!-----------------------------------------------------------------------
!
      IF(PREC) THEN
!
! COMPUTES C G0 = R
!
        IF(CROUT.OR.GSEB.OR.PREBE) THEN
          CALL DOWNUP(G, AUX , R , 'D' , MESH)
          IF(NCSIZE.GT.1) CALL PARMOY(G,MESH)
        ELSEIF(PRE3D) THEN
          CALL CPSTVC(R%ADR(1)%P,G%ADR(1)%P)
          CALL TRID3D(AUX%X%R,G%ADR(1)%P%R,R%ADR(1)%P%R,
     &                MESH%NPOIN,BIEF_NBPTS(11,MESH))
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
! COMPUTES THE DIRECTION OF INITIAL DESCENT:
!-----------------------------------------------------------------------
!
      CALL OS( 'X=Y     ' , D , G , G , C )
!
!-----------------------------------------------------------------------
! COMPUTES THE INITIAL PRODUCT A D :
!-----------------------------------------------------------------------
!
      CALL MATRBL( 'X=AY    ',AD,A,D,  C,MESH)
!
!-----------------------------------------------------------------------
!
      IF(PREC) THEN
!
!   COMPUTES  C DPRIM = AD  (DPRIM PUT IN B)
!
        IF(CROUT.OR.GSEB.OR.PREBE) THEN
          CALL DOWNUP(B, AUX , AD , 'D' , MESH)
          IF(NCSIZE.GT.1) CALL PARMOY(B,MESH)
        ELSEIF(PRE3D) THEN
          CALL CPSTVC(R%ADR(1)%P,G%ADR(1)%P)
          CALL TRID3D(AUX%X%R,B%ADR(1)%P%R,AD%ADR(1)%P%R,
     &                MESH%NPOIN,BIEF_NBPTS(11,MESH))
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
! COMPUTES INITIAL RO :
!-----------------------------------------------------------------------
!
      DAD = P_DOTS(D,AD,MESH)
      IF(PREC) THEN
       ADAD = P_DOTS(AD,B,MESH)
      ELSE
       ADAD = P_DOTS(AD,AD,MESH)
      ENDIF
      RO = DAD/ADAD
!
!-----------------------------------------------------------------------
!
! COMPUTES X1 = X0 - RO  * D
!
      CALL OS( 'X=X+CY  ' , X , D , D , -RO )
!
!-----------------------------------------------------------------------
!  ITERATIONS LOOP:
!-----------------------------------------------------------------------
!
2     M  = M  + 1
!
!-----------------------------------------------------------------------
! COMPUTES THE RESIDUAL : R(M) = R(M-1) - RO(M-1) A D(M-1)
!-----------------------------------------------------------------------
!
      CALL OS( 'X=X+CY  ' , R , AD , AD , -RO )
!
!  SOME VALUES WILL CHANGE IN CASE OF PRECONDITIONING
!
      GM1GM1 = GMGM
      RM1GM1 = RMGM
      RMRM   = P_DOTS(R,R,MESH)
      RMDM   = RMRM
      RMGM   = RMRM
      GMGM   = RMRM
!
! CHECKS END:
!
      IF (RMRM.LE.XL*CFG%EPS**2) GO TO 900
!
!-----------------------------------------------------------------------
! PRECONDITIONING : SOLVES C G = R
!-----------------------------------------------------------------------
!
      IF(PREC) THEN
!       UPDATES G BY RECURRENCE (IN B: DPRIM)
        CALL OS( 'X=X+CY  ' , G , B , B , -RO )
      ENDIF
!
!-----------------------------------------------------------------------
! COMPUTES AG :
!-----------------------------------------------------------------------
!
      CALL MATRBL( 'X=AY    ',AG,A,G,  C,MESH)
!
!-----------------------------------------------------------------------
! COMPUTES D BY RECURRENCE:
!-----------------------------------------------------------------------
!
      IF(PREC) THEN
        AGAD = P_DOTS(AG,B,MESH)
      ELSE
        AGAD = P_DOTS(AG,AD,MESH)
      ENDIF
      BETA = - AGAD / ADAD
!
      CALL OS( 'X=CX    ' , D , D , D , BETA )
      CALL OS( 'X=X+Y   ' , D , G , G , C    )
!
!-----------------------------------------------------------------------
! COMPUTES A D :
!-----------------------------------------------------------------------
!
      CALL OS( 'X=CX    ' , AD , AD , AD , BETA )
      CALL OS( 'X=X+Y   ' , AD , AG , AG  , C   )
!
      IF(PREC) THEN
!
!   COMPUTES  C DPRIM = AD  (DPRIM PUT IN B)
!
        IF(CROUT.OR.GSEB.OR.PREBE) THEN
          CALL DOWNUP(B , AUX , AD , 'D' , MESH)
          IF(NCSIZE.GT.1) CALL PARMOY(B,MESH)
        ELSEIF(PRE3D) THEN
          CALL TRID3D(AUX%X%R,B%ADR(1)%P%R,AD%ADR(1)%P%R,
     &                MESH%NPOIN,BIEF_NBPTS(11,MESH))
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
! COMPUTES RO
!-----------------------------------------------------------------------
!
      GAD = P_DOTS(G,AD,MESH)
      IF(PREC) THEN
        ADAD = P_DOTS(AD,B,MESH)
      ELSE
        ADAD = P_DOTS(AD,AD,MESH)
      ENDIF
      RO = GAD/ADAD
!
! COMPUTES X(M) = X(M-1) - RO * D
!
      CALL OS( 'X=X+CY  ' , X , D , D , -RO )
!
      IF(M.LT.CFG%NITMAX) GO TO 2
!
!-----------------------------------------------------------------------
!
      IF(INFOGR) THEN
        TESTL = SQRT( RMRM / XL )
        IF (RELAT) THEN
           IF (LNG.EQ.1) WRITE(LU,103) M,TESTL
           IF (LNG.EQ.2) WRITE(LU,104) M,TESTL
        ELSE
           IF (LNG.EQ.1) WRITE(LU,203) M,TESTL
           IF (LNG.EQ.2) WRITE(LU,204) M,TESTL
        ENDIF
      ENDIF
      GO TO 1000
!
!-----------------------------------------------------------------------
!
900   CONTINUE
!
      IF(INFOGR) THEN
        TESTL = SQRT( RMRM / XL )
        IF (RELAT) THEN
           IF (LNG.EQ.1) WRITE(LU,101) M,TESTL
           IF (LNG.EQ.2) WRITE(LU,102) M,TESTL
        ELSE
           IF (LNG.EQ.1) WRITE(LU,201) M,TESTL
           IF (LNG.EQ.2) WRITE(LU,202) M,TESTL
        ENDIF
      ENDIF
!
1000  RETURN
!
!-----------------------------------------------------------------------
!
!   FORMATS
!
101   FORMAT(1X,'RESCJG (BIEF) : ',
     &                     1I8,' ITERATIONS, PRECISION RELATIVE:',G16.7)
102   FORMAT(1X,'RESCJG (BIEF) : ',
     &                     1I8,' ITERATIONS, RELATIVE PRECISION:',G16.7)
201   FORMAT(1X,'RESCJG (BIEF) : ',
     &                     1I8,' ITERATIONS, PRECISION ABSOLUE :',G16.7)
202   FORMAT(1X,'RESCJG (BIEF) : ',
     &                     1I8,' ITERATIONS, ABSOLUTE PRECISION:',G16.7)
103   FORMAT(1X,'RESCJG (BIEF) : MAX D'' ITERATIONS ATTEINT:',
     &                     1I8,' PRECISION RELATIVE:',G16.7)
104   FORMAT(1X,'RESCJG (BIEF) : EXCEEDING MAXIMUM ITERATIONS:',
     &                     1I8,' RELATIVE PRECISION:',G16.7)
203   FORMAT(1X,'RESCJG (BIEF) : MAX D'' ITERATIONS ATTEINT:',
     &                     1I8,' PRECISION ABSOLUE :',G16.7)
204   FORMAT(1X,'RESCJG (BIEF) : EXCEEDING MAXIMUM ITERATIONS:',
     &                     1I8,' ABSOLUTE PRECISON:',G16.7)
!
!-----------------------------------------------------------------------
!
      END
