C                       ********************
                        SUBROUTINE PREVERSEG
C                       ********************
C
     *(XAUX,AD,AX,TYPDIA,TYPEXT,NPOIN,MESH,NSEG3D)
C
C***********************************************************************
C BIEF VERSION 6.0        11/08/09    J-M HERVOUET (LNHE) 01 30 87 80 18
C
C 11/08/09 JMH : CROSSED AND VERTICAL SEGMENTS SWAPPED (SEE STOSEG41)
C                                       
C***********************************************************************
C
C FONCTION : BUILDING, BY LUMPING A MATRIX DEFINED ON PRISMS,
C            TRIDIAGONAL SYSTEMS FOR EVERY VERTICAL.
C
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |                | -->| 
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : AUCUN
C
C**********************************************************************
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF, EX_PREVERSEG => PREVERSEG
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NPOIN,NSEG3D
C
      DOUBLE PRECISION, INTENT(IN)    :: AD(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: XAUX(NPOIN,*),AX(NSEG3D,2)
C
      CHARACTER(LEN=1), INTENT(IN) :: TYPDIA,TYPEXT
C
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I2,I3,NPLAN,IAN,ICOM,NPOIN2,SEGUP,SEGDOWN,NSEG2D
      INTEGER IPLAN,NSEGH
C
C-----------------------------------------------------------------------
C
C     HERE WE CONSIDER THAT NPOIN < NELMAX TO USE XAUX AS XAUX(NPOIN,3)
C
C     XAUX(I,1) IS COEFFICIENT OF POINT BELOW I IN EQUATION OF POINT I
C     XAUX(I,2) IS THE DIAGONAL
C     XAUX(I,3) IS COEFFICIENT OF POINT ABOVE I IN EQUATION OF POINT I
C
C-----------------------------------------------------------------------
C     INITIALIZING THE DIAGONAL 
C-----------------------------------------------------------------------
C
      IF(TYPDIA(1:1).EQ.'0') THEN
        CALL OV('X=C     ',XAUX(1,2),AD,AD,0.D0,NPOIN)
      ELSEIF(TYPDIA(1:1).EQ.'I') THEN
        CALL OV('X=C     ',XAUX(1,2),AD,AD,1.D0,NPOIN)
      ELSEIF(TYPDIA(1:1).EQ.'Q') THEN
        CALL OV('X=Y     ',XAUX(1,2),AD,AD,0.D0,NPOIN)
      ELSE
        WRITE(LU,*) TYPDIA
        IF(LNG.EQ.1) WRITE(LU,*) 'DIAGONALE INCONNUE DANS PREVERSEG'
        IF(LNG.EQ.2) WRITE(LU,*) 'UNKNOWN TYPE OF DIAGONAL IN PREVERSEG'
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C     LUMPING THE OFF-DIAGONAL TERMS CORRESPONDING TO VERTICAL SEGMENTS
C-----------------------------------------------------------------------
C 
      NPOIN2 = NBPTS(11)
      NPLAN  = NPOIN/NPOIN2
      NSEG2D = NBSEG(11)
      NSEGH  = NSEG2D*NPLAN
C     
      IF(TYPEXT.EQ.'Q') THEN
C       PLANE ON THE BOTTOM
        DO I2=1,NPOIN2
          SEGUP=NSEGH+I2
          XAUX(I2,1)=0.D0
          XAUX(I2,3)=AX(SEGUP,1)
        ENDDO
C       PLANE AT THE FREE SURFACE
        DO I2=1,NPOIN2
          I3=I2+(NPLAN-1)*NPOIN2
          SEGDOWN=NSEGH+NPOIN2*(NPLAN-2)+I2
          XAUX(I3,1)=AX(SEGDOWN,2)
          XAUX(I3,3)=0.D0
        ENDDO
C       OTHER PLANES
        IF(NPLAN.GT.2) THEN
        DO IPLAN=2,NPLAN-1     
          DO I2=1,NPOIN2
            I3=I2+(IPLAN-1)*NPOIN2
            SEGDOWN=NSEGH+NPOIN2*(IPLAN-2)+I2
            SEGUP  =SEGDOWN+NPOIN2
            XAUX(I3,1)=AX(SEGDOWN,2)
            XAUX(I3,3)=AX(SEGUP,1)
          ENDDO
        ENDDO
        ENDIF
      ELSEIF(TYPEXT.EQ.'S') THEN
C       PLANE ON THE BOTTOM
        DO I2=1,NPOIN2
          SEGUP=NSEGH+I2
          XAUX(I2,1)=0.D0
          XAUX(I2,3)=AX(SEGUP,1)
        ENDDO
C       PLANE AT THE FREE SURFACE
        DO I2=1,NPOIN2
          I3=I2+(NPLAN-1)*NPOIN2
          SEGDOWN=NSEGH+NPOIN2*(NPLAN-2)+I2
          XAUX(I3,1)=AX(SEGDOWN,1)
          XAUX(I3,3)=0.D0
        ENDDO
C       OTHER PLANES
        IF(NPLAN.GT.2) THEN
        DO IPLAN=2,NPLAN-1     
          DO I2=1,NPOIN2
            I3=I2+(IPLAN-1)*NPOIN2
            SEGDOWN=NSEGH+NPOIN2*(IPLAN-2)+I2
            SEGUP  =SEGDOWN+NPOIN2
            XAUX(I3,1)=AX(SEGDOWN,1)
            XAUX(I3,3)=AX(SEGUP,1)
          ENDDO
        ENDDO
        ENDIF
      ELSEIF(TYPEXT.EQ.'0') THEN
C       NOTHING TO DO (BUT WHAT IS THE USE OF AN ITERATIVE SOLVER ?)
      ELSE
        WRITE(LU,*) TYPEXT
        IF(LNG.EQ.1) WRITE(LU,*) 'TYPE DE TERMES EXTRA-DIAGONAUX'
        IF(LNG.EQ.1) WRITE(LU,*) 'INCONNUS DANS PREVERSEG'
        IF(LNG.EQ.2) WRITE(LU,*) 'UNKNOWN TYPE OF OFF-DIAGONAL TERMS'
        IF(LNG.EQ.2) WRITE(LU,*) 'IN PREVERSEG'
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C     PARALLELISM
C
      IF(NCSIZE.GT.1) THEN
        IAN    = 3
        ICOM   = 2
        CALL PARCOM2(XAUX(1,1),XAUX(1,2),XAUX(1,3),
     *               NPOIN2,NPLAN,ICOM,IAN,MESH)
      ENDIF    
C
C-----------------------------------------------------------------------
C
      RETURN
      END
