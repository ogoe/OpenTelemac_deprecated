C                       ********************
                        SUBROUTINE GETTRIEBE
C                       ********************
C
     *(XAUX,AD,AX,TETA,IKLE,NPOIN,NELEM,NELMAX,MESH)
C
C***********************************************************************
C BIEF VERSION 5.9        13/08/08    J-M HERVOUET (LNHE) 01 30 87 80 18
C                                       
C***********************************************************************
C
C FONCTION : GETTING THE TRIDIAGONAL PART OF A DIFFUSION MATRIX ON PRISMS
C            REMOVING IT FROM THE INITIAL MATRIX
C
C            IF MTRI IS THIS TRIDIAGONAL PART, MAUX THE RESULT AND MDIF
C            THE DIFFUSION MATRIX, WE DO HERE:
C
C            MAUX = TETA * MTRI
C            MDIF CHANGED INTO (1-TETA) * MDIF 
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
      USE BIEF, EX_GETTRIEBE => GETTRIEBE
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX,NPOIN
      INTEGER, INTENT(IN) :: IKLE(NELMAX,6)
C
      DOUBLE PRECISION, INTENT(IN)    :: TETA
      DOUBLE PRECISION, INTENT(INOUT) :: XAUX(NPOIN,*),AX(NELMAX,*)
      DOUBLE PRECISION, INTENT(INOUT) :: AD(NPOIN)
C
      TYPE(BIEF_MESH) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I1,I2,I3,I4,I5,I6,IELEM,NPLAN,IAN,ICOM,NPOIN2
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
C     INITIALIZING THE DIAGONAL AND OFF-DIAGONAL TERMS
C-----------------------------------------------------------------------
C
      CALL OV('X=C     ',XAUX(1,1),AD,AD,0.D0,NPOIN)
      CALL OV('X=CY    ',XAUX(1,2),AD,AD,TETA,NPOIN)
      CALL OV('X=C     ',XAUX(1,3),AD,AD,0.D0,NPOIN)        
C
      CALL OV('X=CX    ',AD,AD,AD,1.D0-TETA,NPOIN)
C
C-----------------------------------------------------------------------
C     ADDING TRIDIAGONAL TERMS    
C-----------------------------------------------------------------------
C     
      DO IELEM=1,NELEM
C
        I1=IKLE(IELEM,1)
        I2=IKLE(IELEM,2)
        I3=IKLE(IELEM,3)
        I4=IKLE(IELEM,4)
        I5=IKLE(IELEM,5)
        I6=IKLE(IELEM,6)
        XAUX(I1,3)=XAUX(I1,3)+TETA*AX(IELEM,03) ! TERME 1-4
        XAUX(I2,3)=XAUX(I2,3)+TETA*AX(IELEM,08) ! TERME 2-5
        XAUX(I3,3)=XAUX(I3,3)+TETA*AX(IELEM,12) ! TERME 3-6
        XAUX(I4,1)=XAUX(I4,1)+TETA*AX(IELEM,03) ! TERME 4-1
        XAUX(I5,1)=XAUX(I5,1)+TETA*AX(IELEM,08) ! TERME 5-2
        XAUX(I6,1)=XAUX(I6,1)+TETA*AX(IELEM,12) ! TERME 6-3
C
        AX(IELEM,03)=AX(IELEM,03)*(1.D0-TETA)
        AX(IELEM,08)=AX(IELEM,08)*(1.D0-TETA)
        AX(IELEM,12)=AX(IELEM,12)*(1.D0-TETA)
C
      ENDDO
C
C-----------------------------------------------------------------------
C
C     PARALLELISM
C
      IF(NCSIZE.GT.1) THEN
        IAN    = 3
        ICOM   = 2
        NPOIN2 = NBPTS(11)
        NPLAN=NPOIN/NPOIN2
        CALL PARCOM2(XAUX(1,1),XAUX(1,2),XAUX(1,3),
     *               NPOIN2,NPLAN,ICOM,IAN,MESH)
      ENDIF    
C
C-----------------------------------------------------------------------
C
      RETURN
      END
