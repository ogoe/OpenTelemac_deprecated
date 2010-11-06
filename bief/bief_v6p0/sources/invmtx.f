C                           *****************      
                            SUBROUTINE INVMTX
C                           *****************
C
     *(AM,BM,NP) 
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.7    28/07/2006                        Chun WANG
C
C***********************************************************************
C
C      FONCTION: This subroutine inverts a matrix of np by np. BM is
C                the inversion of AM. 
C      NOTE::    This subroutine calls LUDCMP and LUBKSB, which were
C                copied from "Numeric recipes"-- a well-known book.     
C
C      This subroutine is called by SPECTRE
C
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER LNG,LU             
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NP
      DOUBLE PRECISION, INTENT(INOUT) :: AM(NP,NP),BM(NP,NP)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C       
      INTEGER INDX(500),N,I,J
C
C-----------------------------------------------------------------------
C
      IF(NP.GT.500) THEN
        WRITE(LU,*) 'NP MUST BE LESS THAN 500 IN INVMTX'
        CALL PLANTE(1)
        STOP
      ENDIF
C       
      N = NP
      DO I=1,N !Set up identity matrix. 
        DO J=1,N 
          BM(I,J)=0.D0 
        ENDDO 
        BM(I,I)=1.D0 
      ENDDO
C
C     Decompose the matrix just once.
C             
      CALL LUDCMP(AM,N,NP,INDX)
C
C     Find inverse by columns.
C 
      DO J=1,N  
        CALL LUBKSB(AM,N,NP,INDX,BM(1,J)) 
      ENDDO 
C
C-----------------------------------------------------------------------
C
      RETURN
      END
