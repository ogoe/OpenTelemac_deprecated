C                       ****************************
                        SUBROUTINE BEDLOAD_ENGEL_OLD
C                       ****************************
C
     *(TETAP,CF,NPOIN,GRAV,DM,DENS,TETA,QSC)
C
C***********************************************************************
C SISYPHE VERSION 5.8  11/07/2007 J-M HERVOUET (SUPPRESSION DES OS)    
C SISYPHE VERSION 5.4  --/11/2003   C.VILLARET                         
C SISYPHE VERSION 5.1  11/09/1995  E. PELTIER                          
C SISYPHE VERSION 5.1  11/09/1995  C. LENORMANT                        
C SISYPHE VERSION 5.1  11/09/1995  J.-M. HERVOUET 
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT
C***********************************************************************
C
C  FONCTION  : Bed-load transport formula of Engelund-Hansen
C              Formulation differente de bedload_engel
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |  UCMOY         | -->| NORM OF VELOCITY 
C |  QSC           |<-- | 
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : AUCUN
C
C***********************************************************************
C
      USE INTERFACE_SISYPHE,EX_BEDLOAD_ENGEL_OLD => BEDLOAD_ENGEL_OLD
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ),   INTENT(IN)    :: TETAP,CF
      INTEGER,          INTENT(IN)    :: NPOIN
      DOUBLE PRECISION, INTENT(IN)    :: GRAV, DM, DENS
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: TETA! work array T1
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: QSC
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER          :: I
      DOUBLE PRECISION :: CENGEL
!
      INTRINSIC SQRT
!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!
!     CONTRAINTE DE PEAU ADIMENSIONNELLE: TETAP
!
!     CONTRAINTE TOTALE ADIMENSIONNELLE
!
      DO I = 1, NPOIN
        IF(TETAP%R(I) <= 0.06D0) THEN
          TETA%R(I) = 0.D0
        ELSEIF(TETAP%R(I) <  0.384D0) THEN
          TETA%R(I) = SQRT( 2.5D0 * (TETAP%R(I) - 0.06D0))
        ELSEIF(TETAP%R(I) <  1.080D0) THEN
          TETA%R(I) = 1.066D0 * TETAP%R(I)**0.176D0
        ELSE 
          TETA%R(I) = TETAP%R(I)
        ENDIF
      ENDDO
!
!     TRANSPORT PAR CHARRIAGE
! 
      CENGEL = 0.1D0*SQRT(DENS*GRAV*DM**3)
      DO I=1,NPOIN
        QSC%R(I)=CENGEL*SQRT(TETA%R(I)**5)/MAX(CF%R(I),1.D-6)
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
