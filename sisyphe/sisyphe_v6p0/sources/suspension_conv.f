!                      **************************
                       SUBROUTINE SUSPENSION_CONV
!                      **************************
!
     *(TOB, XMVE, KSR, NPOIN, ZREF, U2D, V2D, HN, HMIN,
     * UCONV, VCONV, KARMAN, ZERO, XWC,T1,ALPHA)
!
!***********************************************************************
! SISYPHE VERSION 6.0  01/08/06        C. VILLARET (LNHE) 01 30 87 83 28
!
! 02/05/2008 : ARGUMENT ALPHA AJOUTE, POUR CONSERVER LE COEFFICIENT
!              DE CORRECTION
!  
!***********************************************************************
!
!      FONCTION:
!      =========
!
!      Correction des vitesses U2D, V2D 
!      called by suspension_computation and init_transport
!
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________.
! !  NOM           !MODE!                  ROLE                        !
! !________________!____!______________________________________________!
! !  TOB           ! -->! SHEAR STRESS
! !  XMVE          ! -->! WATER DENSITY
! !  CF            ! -->! FRICTION COEFFICIENT
! !  NPOIN         ! -->! NUMBER OF POINTS
! !________________!____!______________________________________________!
! MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
!
!***********************************************************************
!
      USE BIEF
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU         
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C 
      TYPE (BIEF_OBJ),  INTENT(IN)    :: HN,U2D,V2D,ZREF,KSR
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: UCONV,VCONV,T1,ALPHA
      TYPE (BIEF_OBJ),  INTENT(IN)    :: TOB
      INTEGER,          INTENT(IN)    :: NPOIN
      DOUBLE PRECISION, INTENT(IN)    :: ZERO,XWC,HMIN
      DOUBLE PRECISION, INTENT(IN)    :: KARMAN,XMVE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C  
      DOUBLE PRECISION :: R1,I1,I2,A,B,AUX,LAUX,LL,USTAR,ROUSE
!
      INTEGER :: I
!
      CALL OS('X=N(Y,Z)', X=T1, Y=U2D, Z=V2D)
!
      LL=LOG(30.D0)
!
      DO I = 1, NPOIN
!
        IF(TOB%R(I).GT.ZERO) THEN                           
!
          USTAR = SQRT(TOB%R(I)/XMVE)             
!
!         B --> KS/H
!
!         AUX = 1.D0 + KARMAN*SQRT(2.D0/MAX(CF%R(I),ZERO))
!         B = 30.D0*EXP(-AUX)
!
          B = KSR%R(I) /MAX(HN%R(I),1.1D0*KSR%R(I))            
          A = ZREF%R(I)/MAX(HN%R(I),1.1D0*ZREF%R(I))
!        
! take max value of A = ZREF/H and B=KSR/H
          A=MAX(A,B)           
!
! version simplifiée 
          ROUSE=MIN(XWC/MAX(USTAR,ZERO),1.D0)/KARMAN
          R1=  1.D0-ROUSE
          LAUX=LOG(A)
! 
          IF(ABS(R1).LT.1.D-8) THEN
            I1= -LAUX
            I2= -LAUX**2/2.D0
          ELSE 
            AUX=A**R1
            I1=(1.D0-AUX)/R1
            I2=-(I1+LAUX*AUX)/R1                     
          ENDIF
!
!         AUX=LOG(A/30.D0)
          AUX=LAUX - LL
          ALPHA%R(I)=-(I2-AUX*I1)/(I1*(AUX+1.D0))
!
        ELSE     
!
          ALPHA%R(I)=1.D0
!
        ENDIF
!           
!       CHECK 0<ALPHA<1
!  
        ALPHA%R(I)=MIN(ALPHA%R(I),1.D0)
        ALPHA%R(I)=MAX(ALPHA%R(I),0.D0) 
!      
        UCONV%R(I) = ALPHA%R(I)*U2D%R(I)
        VCONV%R(I) = ALPHA%R(I)*V2D%R(I)
!  
      ENDDO
C 
C-----------------------------------------------------------------------
C         
      RETURN
      END
