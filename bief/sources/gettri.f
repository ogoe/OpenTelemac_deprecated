C                       *****************
                        SUBROUTINE GETTRI
C                       *****************
C
     *(M,MDIFF,TETA,MESH3D,NPLAN,NPOIN2,NSEG2D)
C
C***********************************************************************
C BIEF VERSION 5.6         16/06/05    J-M HERVOUET (LNH) 01 30 87 80 18
C                                       
C***********************************************************************
C
C FONCTION : GETTING THE TRIDIAGONAL PART OF A DIFFUSION MATRIX ON PRISMS
C            REMOVING IT FROM THE INITIAL MATRIX
C
C            IF MTRI IS THIS TRIDIAGONAL PART, M THE RESULT AND MDIF
C            THE DIFFUSION MATRIX, WE DO HERE:
C
C            M = TETA * MTRI
C            MDIF CHANGED INTO (1-TETA) * MDIF             
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
      USE BIEF, EX_GETTRI => GETTRI 
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C-----------------------------------------------------------------------
C
      INTEGER, INTENT(IN) :: NPLAN,NPOIN2,NSEG2D
C
      DOUBLE PRECISION, INTENT(IN)    :: TETA
      DOUBLE PRECISION, INTENT(INOUT) :: M(NPOIN2,*)
C
      TYPE(BIEF_OBJ),  INTENT(INOUT)  :: MDIFF
      TYPE(BIEF_MESH), INTENT(IN)     :: MESH3D
C
C-----------------------------------------------------------------------
C
      IF(MDIFF%STO.EQ.1) THEN
C
        CALL GETTRIEBE(M,MDIFF%D%R,MDIFF%X%R,TETA,
     *                 MESH3D%IKLE%I,
     *                 MESH3D%NPOIN,MESH3D%NELEM,MESH3D%NELMAX,MESH3D)
C
      ELSEIF(MDIFF%STO.EQ.3) THEN
C
        CALL GETTRISEG(M,MDIFF%D%R,MDIFF%X%R,TETA,
     *                 MESH3D%NPOIN,MESH3D,
     *                 MESH3D%NSEG,NSEG2D,NPLAN,NPOIN2)
C
      ELSE
C
        WRITE(LU,*) 'UNKNOWN STORAGE FOR MDIFF IN GETTRI'
        CALL PLANTE(1)
        STOP
C
      ENDIF      
C
C-----------------------------------------------------------------------
C
      RETURN
      END
