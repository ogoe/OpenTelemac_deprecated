C                       *****************
                        SUBROUTINE DISPER
C                       *****************
C
     *( VISC , U , V , H , CF , ELDER , PROPNU )
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.6      26/05/06     C MOULIN (LNH) 30 87 83 81
C                                             + MODIFS JMH
C***********************************************************************
C
C     FONCTION  : CALCUL DES COEFFICIENTS DE DISPERSION TENSORIELS
C                 EN FONCTION DES COEFFICIENTS LONGITUDINAL ET
C                 TRANSVERSAL.
C
C----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !   VISC         !<-- ! COEFF DU TENSEUR  DE DISPERSION (DIM. NPOIN)
C !   U,V          ! -->! COMPOSANTES DE LA VITESSE
C !   H            ! -->! HAUTEUR D'EAU
C !   CF           ! -->! COEFFICIENT DE FROTTEMENT
C !   ELDER        ! -->! COEFFICIENTS ADIMENSIONNELS DE DISPERSION
C !   PROPNU       ! -->! VISCOSITE LAMINAIRE
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C SOUS-PROGRAMME APPELANT : TELMAC
C SOUS-PROGRAMMES APPELES :
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(IN)  :: ELDER(2),PROPNU
      DOUBLE PRECISION, INTENT(IN)  :: H(*),CF(*),U(*),V(*)
      TYPE(BIEF_OBJ), INTENT(INOUT) :: VISC      
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,NPOIN,NPX
C
      DOUBLE PRECISION KL,KT,COST,SINT,NORMV,USTAR          
C
      INTRINSIC SQRT,MAX
C
C-----------------------------------------------------------------------
C CALCUL DES COEFFICIENTS DE DISPERSION
C-----------------------------------------------------------------------
C
      NPOIN = VISC%DIM1
      NPX   = VISC%MAXDIM1
C
      DO 20 I=1,NPOIN
C
         NORMV = MAX(SQRT(U(I)**2+V(I)**2),1.D-6)
         COST = U(I)/NORMV
         SINT = V(I)/NORMV
         USTAR = SQRT( 0.5D0 * CF(I) * ( U(I)**2 + V(I)**2 ) )
         KL = ELDER(1) * USTAR * H(I)
         KT = ELDER(2) * USTAR * H(I)
         VISC%R(I      ) = PROPNU + ( KL - KT ) * COST**2    + KT
         VISC%R(I+NPX  ) = PROPNU + ( KT - KL ) * COST**2    + KL
         VISC%R(I+2*NPX) = PROPNU + ( KL - KT ) * COST*SINT
C
20    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
