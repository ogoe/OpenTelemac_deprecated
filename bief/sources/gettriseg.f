C                       ********************
                        SUBROUTINE GETTRISEG
C                       ********************
C
     *(XAUX,AD,AX,TETA,NPOIN,MESH,NSEG3D,NSEG2D,NPLAN,NPOIN2)
C
C***********************************************************************
C BIEF VERSION 6.0        11/08/09    J-M HERVOUET (LNHE) 01 30 87 80 18
C 
C 11/08/09 JMH : CROSSED AND VERTICAL SEGMENTS SWAPPED (SEE STOSEG41)
C                                      
C***********************************************************************
C
C FONCTION: GETTING THE TRIDIAGONAL PART OF A DIFFUSION MATRIX ON PRISMS
C           REMOVING IT FROM THE INITIAL MATRIX
C
C           HERE SEGMENT STORAGE FOR MDIFF !!!!!!!!!!!!!!!!!!!!!!!!
C
C           IF MTRI IS THIS TRIDIAGONAL PART, MAUX THE RESULT AND MDIF
C           THE DIFFUSION MATRIX, WE DO HERE:
C
C           MAUX = TETA * MTRI
C           MDIF CHANGED INTO (1-TETA) * MDIF 
C            
C
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     NELEM      | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |     NELMAX     | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |     AX         | -->|  TERMES EXTRA-DIAGONAUX PAR SEGMENTS
C |                |    |  (ICI DIMENSION 1 CAR SYMETRIQUE)
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
      USE BIEF, EX_GETTRISEG => GETTRISEG
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NPOIN,NSEG3D,NSEG2D,NPLAN,NPOIN2
C
      DOUBLE PRECISION, INTENT(IN)    :: TETA
      DOUBLE PRECISION, INTENT(INOUT) :: XAUX(NPOIN,*),AX(NSEG3D)
      DOUBLE PRECISION, INTENT(INOUT) :: AD(NPOIN)
C
      TYPE(BIEF_MESH) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I2,I3,IPLAN,IAN,ICOM,SEGUP,SEGDOWN,NSEGH,NSEGV
C
C-----------------------------------------------------------------------
C
      NSEGH=NSEG2D*NPLAN
      NSEGV=NPOIN2*(NPLAN-1)
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
C     INITIALIZING THE DIAGONAL TERMS
C-----------------------------------------------------------------------
C
      CALL OV('X=CY    ',XAUX(1,2),AD,AD,TETA,NPOIN)       
      CALL OV('X=CX    ',AD,AD,AD,1.D0-TETA,NPOIN)
C
C-----------------------------------------------------------------------
C     TRIDIAGONAL TERMS    
C-----------------------------------------------------------------------
C
C     PLANE ON THE BOTTOM
C
      DO I2=1,NPOIN2
        SEGUP=NSEGH+I2
        XAUX(I2,1)=0.D0
        XAUX(I2,3)=TETA*AX(SEGUP)
      ENDDO
C
C     PLANE AT THE FREE SURFACE
C
      DO I2=1,NPOIN2
        I3=I2+(NPLAN-1)*NPOIN2
        SEGDOWN=NSEGH+NPOIN2*(NPLAN-2)+I2
        XAUX(I3,1)=TETA*AX(SEGDOWN)
        XAUX(I3,3)=0.D0
      ENDDO
C
C     OTHER PLANES
C
      IF(NPLAN.GT.2) THEN
      DO IPLAN=2,NPLAN-1     
        DO I2=1,NPOIN2
          I3=I2+(IPLAN-1)*NPOIN2
          SEGDOWN=NSEGH+NPOIN2*(IPLAN-2)+I2
          SEGUP  =SEGDOWN+NPOIN2
          XAUX(I3,1)=TETA*AX(SEGDOWN)
          XAUX(I3,3)=TETA*AX(SEGUP)
        ENDDO
      ENDDO
      ENDIF
C
      CALL OV('X=CX    ',AX(NSEGH+1:NSEGH+NSEGV),AX,AX,1.D0-TETA,NSEGV)
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
