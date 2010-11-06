C                    *****************************
                     SUBROUTINE BUILD_GLOBAL_BOUND
C                    *****************************
C
     *(KNOLG,NPOIN,NPOIN_TOT,NPTFR,NPTFR_TOT,
     * X,Y,K,C,CG,LIHBOR,XT,YT,KT,CTT,CGT,LIHBORT,NBOR,NBOR_TOT)
C
C***********************************************************************
C
C  PARALLEL VERSION 6.0   16/01/2010    C. DENIS (SINETICS)
C
C***********************************************************************
C
C      FONCTION: RECONSTRUIT LA FRONTIERE GLOBAL DU MAILLAGE 
C      ATTENTION A UTILISER UNIQUEMENT SI ON EST EN TRAITEMENT PARALLELE
C
C
C-----------------------------------------------------------------------
C
C APPELE PAR : 
C
C***********************************************************************
C      
      IMPLICIT NONE
      INCLUDE 'mpif.h'
C
      INTEGER, INTENT(IN) :: NPOIN_TOT
      INTEGER, INTENT(IN) :: NPOIN
      INTEGER, INTENT(IN) :: NPTFR
      INTEGER, INTENT(IN) :: NPTFR_TOT
C      
      INTEGER, INTENT(IN), DIMENSION(NPOIN) :: KNOLG
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NPOIN)  :: X, Y, K, C,CG
      DOUBLE PRECISION,INTENT(OUT),DIMENSION(NPOIN_TOT)::XT,YT,
     *                                                   KT,CTT,CGT 
      INTEGER, INTENT(IN) :: LIHBOR(NPTFR)
      INTEGER, INTENT(OUT) :: LIHBORT(NPTFR_TOT)
      INTEGER, INTENT(IN) :: NBOR(NPTFR), NBOR_TOT(NPTFR_TOT)
      INTEGER :: I,J
      INTEGER :: IER
      INTEGER, ALLOCATABLE :: TEMP1(:)
      INTEGER, ALLOCATABLE :: TEMP2(:)
      DOUBLE PRECISION, ALLOCATABLE :: TEMP3(:)
      DOUBLE PRECISION :: TMP
      ALLOCATE(TEMP3(NPOIN_TOT))
      YT(:)=0.0
      XT(:)=0.0
      CTT(:)=0.0
      CGT(:)=0.0
      KT(:)=0.0
C     XT MERGING      
      TEMP3(:)=-HUGE(TMP)
      DO I=1,NPOIN
         TEMP3(KNOLG(I))=X(I)
      END DO
C
      CALL MPI_ALLREDUCE(TEMP3,XT,NPOIN_TOT,MPI_REAL8,MPI_MAX,
     *                   MPI_COMM_WORLD,ier)
      WHERE(XT .EQ. -HUGE(TMP))
         XT=0.0
      END WHERE
C     YT MERGING
      TEMP3(:)=-HUGE(TMP)
      DO I=1,NPOIN
         TEMP3(KNOLG(I))=Y(I)
      END DO
      CALL MPI_ALLREDUCE(TEMP3,YT, NPOIN_TOT,MPI_REAL8,MPI_MAX,
     *                   MPI_COMM_WORLD,ier)     
      WHERE(YT .EQ. -HUGE(TMP))
         YT=0.0
      END WHERE
C     CT MERGING
      TEMP3(:)=-HUGE(TMP)
      DO I=1,NPOIN
         TEMP3(KNOLG(I))=C(I)
      END DO
      CALL MPI_ALLREDUCE(TEMP3,CTT, NPOIN_TOT,MPI_REAL8,MPI_MAX,
     *                   MPI_COMM_WORLD,ier)
      WHERE(CTT .EQ. -HUGE(TMP))
         CTT=0.0
      END WHERE
C     CGT MERGING
      TEMP3(:)=-HUGE(TMP)
      DO I=1,NPOIN
         TEMP3(KNOLG(I))=CG(I)
      END DO
      CALL MPI_ALLREDUCE(TEMP3,CGT, NPOIN_TOT,MPI_REAL8,MPI_MAX,
     *                   MPI_COMM_WORLD,ier)
      WHERE(CGT .EQ.-HUGE(TMP))
         CGT=0.0
      END WHERE
C     KT MERGING
      TEMP3(:)=-HUGE(TMP)
      DO I=1,NPOIN
         TEMP3(KNOLG(I))=K(I)
      END DO
      CALL MPI_ALLREDUCE(TEMP3,KT, NPOIN_TOT,MPI_REAL8,MPI_MAX,
     *                   MPI_COMM_WORLD,ier)
      WHERE(KT .EQ. -HUGE(TMP))
         KT=0.0
      END WHERE
      DEALLOCATE(TEMP3)
      ALLOCATE(TEMP2(NPTFR_TOT))
      TEMP2(:)=-HUGE(I)
      DO I=1,NPTFR_TOT
         DO J=1,NPTFR
            IF (NBOR_TOT(I) .EQ. KNOLG(NBOR(J))) THEN
               TEMP2(I)=LIHBOR(J)
            END IF
         END DO
      END DO
      CALL MPI_ALLREDUCE(TEMP2,LIHBORT, NPTFR_TOT,MPI_INTEGER,MPI_MAX,
     *                   MPI_COMM_WORLD,ier)
      
      WHERE(LIHBORT .EQ. -HUGE(I))
         LIHBORT=0.0
      END WHERE
      DEALLOCATE(TEMP2)
C
C-----------------------------------------------------------------------
C
      RETURN      
      END
