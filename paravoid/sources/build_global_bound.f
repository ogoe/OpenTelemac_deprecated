C                  *****************************      
                   SUBROUTINE BUILD_GLOBAL_BOUND
C                  *****************************
C
     *(KNOLG,NPOIN,NPOIN_TOT,NPTFR,NPTFR_TOT,
     * X,Y,K,C,CG,LIHBOR,XT,YT,KT,CTT,CGT,LIHBORT,NBOR,NBOR_TOT)
C
C***********************************************************************
C
C  PARAVOID VERSION 6.0   16/01/2010    C. DENIS (SINETICS)
C
C***********************************************************************
C
C      FONCTION: RECONSTRUIT LA FRONTIERE GLOBAL DU MAILLAGE 
C      ATTENTION A UTILISER UNIQUEMENT SI ON EST EN TRAITEMENT PARALLELE
C
C      COQUILLE VIDE EN SEQUENTIEL 
C-----------------------------------------------------------------------
C
C APPELE PAR : 
C
C***********************************************************************
C     
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER, INTENT(IN) :: NPOIN_TOT
      INTEGER, INTENT(IN) :: NPOIN
      INTEGER, INTENT(IN) :: NPTFR
      INTEGER, INTENT(IN) :: NPTFR_TOT
C      
      INTEGER, INTENT(IN), DIMENSION(NPOIN) :: KNOLG
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NPOIN)  :: X, Y, K, C,CG
      DOUBLE PRECISION , INTENT(IN), DIMENSION(NPOIN_TOT) :: XT,YT,KT,
     *                                                       CTT,CGT 
      INTEGER, INTENT(IN) :: LIHBOR(NPTFR)
      INTEGER, INTENT(IN) :: LIHBORT(NPTFR_TOT)
      INTEGER, INTENT(IN) :: NBOR(NPTFR), NBOR_TOT(NPTFR_TOT)
C
C-----------------------------------------------------------------------
C    
      IF(LNG.EQ.1) WRITE(LU,*) 'BUILD_GLOBAL_BOUND VERSION VIDE'
      IF(LNG.EQ.2) WRITE(LU,*) 'BUILD_GLOBAL_BOUND IN ITS VOID VERSION'
C
C-----------------------------------------------------------------------
C
      RETURN   
      END
