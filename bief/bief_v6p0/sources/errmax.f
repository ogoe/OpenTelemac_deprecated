C                       *****************                                  
                        SUBROUTINE ERRMAX                                        
C                       *****************                                       
C     
     *(X1,X2,ERR,IERR)
C     
C***********************************************************************
C BIEF VERSION 5.2                               27/04/93    E. BARROS
C                     UPGRADE TO 5.1    02/10/00    A. LEOPARDI (UNINA)                                                                     
C***********************************************************************
C     
C     FUNCTION: COMPUTE MAX DIFFERENCES BETWEEN COMPUTED 2 ARRAYS      
C     
C-----------------------------------------------------------------------
C     ARGUMENTS                                         
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    X1,X2       | -->| ARRAYS TO COMPARE 
C |    ERR         |<-- | MAXIMUM ABSOLUTE DIFFERENCE
C |    IERR        |<-- | POINT WHERE THE DIFFERENCE OCCURS
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C     
C     APPELE PAR :            HOMERE_PIT
C     
C     SOUS-PROGRAMME APPELE : RIEN
C     
C***********************************************************************
C     
      USE BIEF_DEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER          , INTENT(OUT) :: IERR
      DOUBLE PRECISION , INTENT(OUT) :: ERR
      TYPE (BIEF_OBJ)  , INTENT(IN)  :: X1,X2    
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C
      INTRINSIC ABS
C     
C---------------------------------------------------------------------
C
      IERR=1
      ERR=-1.D0    
      DO I=1,X1%DIM1                              
C   
        IF(ABS(X1%R(I)-X2%R(I)).GT.ERR) THEN      
           ERR=ABS(X1%R(I)-X2%R(I))
           IERR=I
        ENDIF
C     
      ENDDO     
C     
C---------------------------------------------------------------------
C 
      RETURN
      END
