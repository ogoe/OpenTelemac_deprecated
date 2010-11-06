C                          *****************     
                           SUBROUTINE LUDCMP 
C                          *****************
C
     *(A,N,NP,INDX)
C
C************************************************************************
C  TELEMAC 2D VERSION 5.7    28/07/2006                         Chun WANG
C
C************************************************************************
C
C      FONCTION:  Given a matrix a(1:n,1:n), with physical dimension np 
C                 by np, this routine replaces it by the LU decomposition 
C                 of a rowwise permutation of itself. a and n are input. 
C                 a is output, arranged as in equation (2.3.14) above; 
C                 indx(1:n) is an output vector that records the row 
C                 permutation effected by the partial pivoting; d is 
C                 output as \Sigma1 depending on whether the number of 
C                 row interchanges was even or odd, respectively. This 
C                 routine is used in combination with lubksb to solve 
C                 linear equations or invert a matrix. 
C
C      This subroutine is called by INVMTX
C
C------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER LNG,LU             
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: N,NP
      INTEGER, INTENT(INOUT)          :: INDX(N)
      DOUBLE PRECISION, INTENT(INOUT) :: A(NP,NP)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C  
      DOUBLE PRECISION D 
      INTEGER I,IMAX,J,K 
      DOUBLE PRECISION AAMAX,DUM,XSOM,VV(500) 
C
C------------------------------------------------------------------------
C
      D=1.D0 !No row interchanges yet. 
C
C     Loop over rows to get the implicit scaling information 
C
      DO I=1,N 
        AAMAX=0.D0 
        DO J=1,N 
          IF(ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J)) 
        ENDDO  
        IF(AAMAX.LT.1.D-20) THEN
          WRITE(LU,*) 'Singular matrix in ludcmp'
          CALL PLANTE(1)
          STOP 
        ENDIF                                                     
        VV(I)=1.D0/AAMAX !Save the scaling. 
      ENDDO
C
C     This is the loop over columns of Crout's method
C  
      DO J=1,N  
         DO I=1,J-1 !This is equation (2.3.12) except for i = j. 
            XSOM=A(I,J) 
            DO K=1,I-1 
               XSOM=XSOM-A(I,K)*A(K,J) 
            ENDDO  
            A(I,J)=XSOM 
         ENDDO 
         AAMAX=0.D0 !Initialize for the search for largest pivot element. 
         DO I=J,N !This is i = j of equation (2.3.12) and 
                  ! i = j +1 : : : N of equation (2.3.13). 
            XSOM=A(I,J) 
            DO K=1,J-1 
               XSOM=XSOM-A(I,K)*A(K,J) 
            ENDDO  
            A(I,J)=XSOM 
            DUM=VV(I)*ABS(XSOM) !Figure of merit for the pivot. 
            IF (DUM.GE.AAMAX) THEN !Is it better than the best so far? 
               IMAX=I 
               AAMAX=DUM 
            ENDIF 
         ENDDO 
         IF (J.NE.IMAX) THEN !Do we need to interchange rows? 
            DO K=1,N  
              DUM=A(IMAX,K) 
              A(IMAX,K)=A(J,K) 
              A(J,K)=DUM 
            ENDDO 
            D=-D !...and change the parity of d. 
            VV(IMAX)=VV(J) !Also interchange the scale factor. 
         ENDIF 
         INDX(J)=IMAX 
         A(J,J)=MAX(A(J,J),1.D-20) 
C
C  If the pivot element is zero the matrix is singular (at least to the 
C  precision of the algorithm).
C
         IF(J.NE.N) THEN !Now, finally, divide by the pivot element. 
           DUM=1.D0/A(J,J) 
           DO I=J+1,N 
             A(I,J)=A(I,J)*DUM 
           ENDDO  
         ENDIF
C 
      ENDDO 
C
C------------------------------------------------------------------------
C 
      RETURN 
      END
