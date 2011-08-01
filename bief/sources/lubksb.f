C                            *****************
                             SUBROUTINE LUBKSB 
C                            *****************  
C
     *(A,N,NP,INDX,B)
C
C************************************************************************
C  TELEMAC 2D VERSION 5.7    28/07/2006                         Chun WANG
C************************************************************************
C
C      FONCTION:  Solves the set of n linear equations A \Delta X = B. 
C                 Here a is input, not as the matrix A but rather as 
C                 its LU decomposition, determined by the routine ludcmp. 
C                 indx is input as the permutation vector returned by 
C                 ludcmp. b(1:n) is input as the right­hand side vector B, 
C                 and returns with the solution vector X. a, n, np, and 
C                 indx are not modified by this routine C and can be left 
C                 in place for successive calls with different right­hand 
C                 sides b. This routine takes into account the possibility 
C                 that b will begin with many zero elements, so it is 
C                 efficient for use in matrix inversion. 
C
C      This subroutine is called by INVMTX
C
C-----------------------------------------------------------------------     
C
      IMPLICIT NONE
      INTEGER LNG,LU             
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: N,NP
      INTEGER, INTENT(IN) :: INDX(N)
      DOUBLE PRECISION, INTENT(INOUT) :: B(N)
      DOUBLE PRECISION, INTENT(IN)    :: A(NP,NP)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C    
      INTEGER I,II,J,LL 
      DOUBLE PRECISION XSOM 
C
C-----------------------------------------------------------------------     
C
      II=0 !When ii is set to a positive value, it will become the in­ 
           !dex of the first nonvanishing element of b. We now do 
           !the forward substitution, equation (2.3.6). The only new 
           !wrinkle is to unscramble the permutation as we go. 
C
      DO I=1,N 
        LL=INDX(I) 
        XSOM=B(LL) 
        B(LL)=B(I) 
        IF(II.NE.0) THEN 
          DO J=II,I-1 
           XSOM=XSOM-A(I,J)*B(J) 
          ENDDO  
        ELSEIF(XSOM.NE.0.) THEN 
          II=I !A nonzero element was encountered, so from now on we will 
               !have to do the sums in the loop above. 
        ENDIF 
        B(I)=XSOM 
      ENDDO 
      DO I=N,1,-1 !Now we do the backsubstitution, equation (2.3.7). 
        XSOM=B(I) 
        DO J=I+1,N 
          XSOM=XSOM-A(I,J)*B(J) 
        ENDDO 
        B(I)=XSOM/A(I,I) !Store a component of the solution vector X. 
      ENDDO
C
C-----------------------------------------------------------------------     
C  
      RETURN 
      END
