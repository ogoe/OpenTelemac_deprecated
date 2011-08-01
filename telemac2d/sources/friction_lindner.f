C                       ***************************
                        SUBROUTINE FRICTION_LINDNER
C                       ***************************
C
     & (VA,HA,CF,VK,G,DP,SP,CP)
C
C***********************************************************************
C  TELEMAC-2D VERSION 5.5                 J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C 20/04/04 : subroutine written by F. Huvelin
C (Written from the C++ program RISMO2D of the BAW)
C
C
          ! ---------------------------------------------- !
          !        Compute friction coefficient for        !
          !    non-submerged vegetation from parameters    !
          ! ---------------------------------------------- !
!
! The algorithm was developed by LINDNER (1982) and PASCHE (1986).
! I have made some changes to the computation of depth ratio (1990)
! and to the wake length equation (1992): the slope of energy line
! is no longer requested (this essential for use in 2D methods).
! Michael Schroeder in November 1992, BAW
C
C
C               TTTTT EEEEE L     EEEEE M   M   AA  CCCCC
C                 T   E     L     E     MM MM  A  A C
C                 T   EEE   L     EEE   M M M  AAAA C
C                 T   E     L     E     M   M  A  A C
C                 T   EEEEE LLLLL EEEEE M   M  A  A CCCCC
C
C
C----------------------------------------------------------------------C
C                             ARGUMENTS                                C
C .________________.____.______________________________________________C
C |      NOM       |MODE|                   ROLE                       C
C |________________|____|______________________________________________C
C | VA             | => | Velocity                                     C
C | HA             | => | Flow depth                                   C
C | CF             | => | Friction coefficient for bottom roughness    C
C | VK             | => | kinemtic viscosity                           C
C | G              | => | gravity acceleration                         C
C | DP             | => | Diameter of roughness element                C
C | SP             | => | Spacing of roughhness element                C
C | CP             | <= | Friction coeff for non-submerged vegetation  C
C |________________|____|______________________________________________C
C                    <=  input value                                   C
C                    =>  output value                                  C 
C ---------------------------------------------------------------------C
!
!=======================================================================!
!=======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS                !
!=======================================================================!
!=======================================================================!
!
      IMPLICIT NONE      
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(IN)  :: VA,HA,CF,VK,G,DP,SP
      DOUBLE PRECISION, INTENT(OUT) :: CP
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER,          PARAMETER :: kMaxIter = 200
      DOUBLE PRECISION, PARAMETER :: kPrecision = 1.0D-3
      INTEGER                     :: cWRav
      INTEGER                     :: cWRmax
      INTEGER                     :: cWRcount
      INTEGER                     :: aNLav
      INTEGER                     :: aNLmax
      INTEGER                     :: aNLcount
      INTEGER                     :: itErr
!
      INTEGER :: i, J
      INTEGER :: icWR               ! iteration counter: cWR
      INTEGER :: iaNL               ! iteration counter: aNL
      INTEGER :: realRoots
      LOGICAL :: lcWR
!  
      DOUBLE PRECISION :: cW, cWR, RcWR, dcWR, aNL, aNB, RaNL
      DOUBLE PRECISION :: cWR1, cWR2, dcWR1, dcWR2
      DOUBLE PRECISION :: lambda, Fr
      DOUBLE PRECISION :: x(3), vRatio, hRatio
      DOUBLE PRECISION :: alfa, aCof, bCof, cCof, dCof
!
      DOUBLE PRECISION :: tmp1, tmp2
!
!=======================================================================!
!=======================================================================!
!                               PROGRAMME                               !
!=======================================================================!
!=======================================================================!
!
      lcWR     = .TRUE.
      cWRav    = 0
      cWRmax   = 0
      cWRcount = 0
      aNLav    = 0
      aNLmax   = 0
      aNLcount = 0
      itErr    = 0

      IF ((dp < 1.0e-3)     .OR.(sp < 1.0e-2).OR.
     &    (ABS(va) < 1.0e-3).OR.(ha < 1.0e-3)     ) THEN

         CP = 0.D0

      ELSE

         ! Initialization
         ! --------------
         cWR   = 1.0      ! drag coefficient
         aNL   = sp/2.0   ! wake length of a cylinder
         cWR1  = 1.0 
         cWR2  = 1.0
         dcWR1 = 0.0 
         dcWR2 = 0.0

         ! Start of iteration for cWR
         ! --------------------------
         DO icWR = 1, kMaxIter

            ! superposed friction coefficient
            ! -------------------------------
            lambda = 8.0D0*cf  +  4.0D0*cWR*ha*dp/sp/sp

            ! drag coefficient cW for one cylinder
            ! ------------------------------------ 
            CALL dragCoeff(va, dp, vk, cW)

            ! wake length of a cylinder (iterative computation)
            ! -------------------------------------------------
            DO J=1, kMaxIter

               tmp1 = 1.0D0  +  aNL*lambda/4.0D0/ha
               tmp2 = 30.0D0/ABS(tmp1)**(1.5)
               RaNL = cW*dp*ABS(tmp2)**(1.429)

               ! test for convergence
               ! --------------------
               IF (ABS((RaNL-aNL)/RaNL) < kPrecision) THEN
                  aNL  = RaNL
                  iaNL = -1*J
                  EXIT
               ENDIF

               aNL = 0.5 * (RaNL + aNL)
            ENDDO

            ! statistics of cWR iteration
            ! ---------------------------
            IF ( iaNL > 0 ) THEN
               aNL = sp/2.0D0
            ELSE
               iaNL = ABS(iaNL)
               aNLcount = aNLcount + 1
               aNLav = iaNL + aNLav
               IF (iaNL > aNLmax) aNLmax = iaNL
            ENDIF

            ! wake width
            ! ----------
            aNB = 0.24 * ABS(aNL)**(0.59) * ABS(cW*dp)**(0.41)

            ! ratio of velocity in front of and behind cylinder
            ! -------------------------------------------------
            vRatio = 1.151 * ABS(aNL/sp)**(-0.483)
     &             +   0.5 * ABS(aNB/sp)**(1.1)

            ! ratio of flow depth
            ! -------------------
            Fr = va / SQRT( g * ha ) ! Froude number
    
            alfa = dp / sp
            aCof =  Fr * Fr * (1.0D0 - alfa * cWR/2.0D0)
            bCof = -Fr * Fr - (1.0D0 - alfa) / 2.0D0
            cCof =  0.0D0
            dCof = (1.0D0 - alfa) / 2.0D0
            hRatio = 1.0D0

            IF (ABS(aCof) < 1.0e-10) THEN
               hRatio = SQRT( -dCof / bCof)

            ELSE
               CALL cubeEquation(aCof, bCof, cCof, dCof, realRoots, x)

               DO i = 1, RealRoots
                  IF (x(i) > 0.0  .AND.  x(i) < 1.0)  THEN
                     hRatio = x(i)
                     EXIT
                  ENDIF
               ENDDO
            ENDIF

            ! revise drag coefficient cWR
            ! ---------------------------
            RcWR = 1.3124D0*cW*vRatio + 2.0D0*(1.0D0-hRatio)/Fr/Fr

            ! test for convergence
            ! --------------------
            IF ( ABS((RcWR-cWR)/RcWR) < kPrecision ) THEN
               lcWR = .FALSE.
!               icWR = -1/icWR
               EXIT
            ENDIF

            ! use PEGASUS algorithm for cWR iteration
            ! ---------------------------------------
            dcWR = RcWR - cWR
    
            IF ((icWR >= 3) .AND. (dcWR1*dcWR2 < 0.0D0)) THEN

               IF (dcWR2*dcWR < 0.0D0) THEN
                  dcWR1 = dcWR2/(dcWR2+dcWR)*dcWR1

               ELSE
                  cWR1  = cWR2
                  dcWR1 = dcWR2
               ENDIF
               cWR2  = cWR
               dcWR2 = dcWR
               cWR   = cWR2 - dcWR2*(cWR2-cWR1)/(dcWR2-dcWR1)

            ELSE
               cWR1 = cWR2
               dcWR1 = dcWR2
               cWR2 = cWR
               dcWR2 = dcWR

               IF ((icWR >= 2) .AND. (dcWR1*dcWR2 < 0.0 )) THEN
                  cWR = cWR2 - dcWR2*(cWR2-cWR1)/(dcWR2-dcWR1) 
               ELSE
                  cWR = RcWR
               ENDIF
            ENDIF

         ENDDO !icWR = 1, kMaxIter


         IF (lcWR) THEN
            itErr = itErr + 1
            CP = -1.D0

         ELSE
            ! statistics of cWR iteration
            ! ---------------------------
            icWR = -1/icWR ! as the program RISMO2d from the BAW
            icWR = -icWR
            cWRcount = cWRcount + 1
            cWRav = icWR + cWRav

            IF (icWR > cWRmax) cWRmax = icWR

            CP = lambda/8.0D0 - cf
         ENDIF

      ENDIF
!
!=======================================================================!
!=======================================================================!
!
      RETURN
      END SUBROUTINE FRICTION_LINDNER
