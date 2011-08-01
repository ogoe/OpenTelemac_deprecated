C                       ****************
                        SUBROUTINE OVD_2
C                       ****************
C
     * ( OP , X , DIMX , Y , DIMY , Z , DIMZ , C , DIM1 , NPOIN ,
     *   IOPT , INFINI, ZERO )
C
C***********************************************************************
C BIEF VERSION 5.2           29/11/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C  FONCTION : BETWEEN OS AND OVD WHEN 2-DIMENSION VECTORS ARE INVOLVED
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |      OP        | -->| CHAINE DE CARACTERES INDIQUANT L'OPERATION
C |                |   >| A EFFECTUER.
C |      X , DIMX  |<-- | STRUCTURE RESULTAT ET DIMENSION A TRAITER
C |      Y , DIMY  | -->| STRUCTURE OPERANDE ...
C |      Z , DIMZ  | -->| STRUCTURE OPERANDE ...
C |      C         | -->| CONSTANTE DONNEE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER,          INTENT(IN)    :: DIMX,DIMY,DIMZ,DIM1,NPOIN,IOPT
      DOUBLE PRECISION, INTENT(IN)    :: C,INFINI,ZERO
      CHARACTER(LEN=8), INTENT(IN)    :: OP
      DOUBLE PRECISION, INTENT(INOUT) :: X(DIM1,*)
      DOUBLE PRECISION, INTENT(IN)    :: Y(DIM1,*),Z(DIM1,*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CALL OVD( OP , X(1,DIMX) , Y(1,DIMY) , Z(1,DIMZ) , C , NPOIN ,
     *          IOPT , INFINI , ZERO )
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
