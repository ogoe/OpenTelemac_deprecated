C                       ********************
                        SUBROUTINE PREVEREBE
C                       ********************
C
     *(XAUX,AD,AX,TYPDIA,TYPEXT,IKLE,NPOIN,NELEM,NELMAX,MESH)
C
C***********************************************************************
C BIEF VERSION 5.9        02/06/08    J-M HERVOUET (LNHE) 01 30 87 80 18
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
      USE BIEF, EX_PREVEREBE => PREVEREBE
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
      DOUBLE PRECISION, INTENT(IN) :: AD(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: XAUX(NPOIN,*),AX(NELMAX,*)
C
      CHARACTER(LEN=1), INTENT(IN) :: TYPDIA,TYPEXT
C
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
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
      CALL OV('X=C     ',XAUX(1,3),AD,AD,0.D0,NPOIN)
C
      IF(TYPDIA(1:1).EQ.'0') THEN
        CALL OV('X=C     ',XAUX(1,2),AD,AD,0.D0,NPOIN)
      ELSEIF(TYPDIA(1:1).EQ.'I') THEN
        CALL OV('X=C     ',XAUX(1,2),AD,AD,1.D0,NPOIN)
      ELSEIF(TYPDIA(1:1).EQ.'Q') THEN
        CALL OV('X=Y     ',XAUX(1,2),AD,AD,0.D0,NPOIN)
      ELSE
       WRITE(LU,*) TYPDIA
       IF(LNG.EQ.1) WRITE(LU,*) 'DIAGONALE INCONNUE DANS PREVEREBE'
       IF(LNG.EQ.2) WRITE(LU,*) 'UNKNOWN TYPE OF DIAGONAL IN PREVEREBE'
       CALL PLANTE(1)
       STOP
      ENDIF
C
C-----------------------------------------------------------------------
C     LUMPING THE OFF-DIAGONAL TERMS
C-----------------------------------------------------------------------
C      
      IF(TYPEXT.EQ.'Q') THEN
        DO IELEM=1,NELEM
          I1=IKLE(IELEM,1)
          I2=IKLE(IELEM,2)
          I3=IKLE(IELEM,3)
          I4=IKLE(IELEM,4)
          I5=IKLE(IELEM,5)
          I6=IKLE(IELEM,6)
          XAUX(I1,3)=XAUX(I1,3)+AX(IELEM,03) ! TERME 1-4
          XAUX(I2,3)=XAUX(I2,3)+AX(IELEM,08) ! TERME 2-5
          XAUX(I3,3)=XAUX(I3,3)+AX(IELEM,12) ! TERME 3-6
          XAUX(I4,1)=XAUX(I4,1)+AX(IELEM,18) ! TERME 4-1
          XAUX(I5,1)=XAUX(I5,1)+AX(IELEM,23) ! TERME 5-2
          XAUX(I6,1)=XAUX(I6,1)+AX(IELEM,27) ! TERME 6-3
        ENDDO
      ELSEIF(TYPEXT.EQ.'S') THEN
        DO IELEM=1,NELEM
          I1=IKLE(IELEM,1)
          I2=IKLE(IELEM,2)
          I3=IKLE(IELEM,3)
          I4=IKLE(IELEM,4)
          I5=IKLE(IELEM,5)
          I6=IKLE(IELEM,6)
          XAUX(I1,3)=XAUX(I1,3)+AX(IELEM,03) ! TERME 1-4
          XAUX(I2,3)=XAUX(I2,3)+AX(IELEM,08) ! TERME 2-5
          XAUX(I3,3)=XAUX(I3,3)+AX(IELEM,12) ! TERME 3-6
          XAUX(I4,1)=XAUX(I4,1)+AX(IELEM,03) ! TERME 4-1
          XAUX(I5,1)=XAUX(I5,1)+AX(IELEM,08) ! TERME 5-2
          XAUX(I6,1)=XAUX(I6,1)+AX(IELEM,12) ! TERME 6-3
        ENDDO
      ELSEIF(TYPEXT.EQ.'0') THEN
C       NOTHING TO DO (BUT WHAT IS THE USE OF AN ITERATIVE SOLVER ?)
      ELSE
        WRITE(LU,*) TYPEXT
        IF(LNG.EQ.1) WRITE(LU,*) 'TYPE DE TERMES EXTRA-DIAGONAUX'
        IF(LNG.EQ.1) WRITE(LU,*) 'INCONNUS DANS PREVEREBE'
        IF(LNG.EQ.2) WRITE(LU,*) 'UNKNOWN TYPE OF OFF-DIAGONAL TERMS'
        IF(LNG.EQ.2) WRITE(LU,*) 'IN PREVEREBE'
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
