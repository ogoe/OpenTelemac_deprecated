! Subroutine to calculate the new tau from secondary currents
!************************************************
       SUBROUTINE BEDLOAD_SECCURRENT(IELMU)
!************************************************
!
!
       USE DECLARATIONS_SISYPHE
       USE BIEF
       IMPLICIT NONE
!
       INTEGER LNG,LU
       COMMON/INFO/LNG,LU
!
       INTEGER I, IELMU
       DOUBLE PRECISION C, ALPHA
!
! remember: QU = U_tel*H_tel, QV=V_tel*H_tel
!
!
!
! calculation of Pi
!       PI = ACOS(-1.D0)
!
!RK change for secondary currents
! calculating the gradient of the free surface in x-direction
       CALL VECTOR(T5,'=','GRADF          X',IELMU,
     &      1.D0,Z,S,S,S,S,S,MESH,MSK,MASKEL)
! for parallel computing
       IF (NCSIZE.GT.1) CALL PARCOM (T5, 2, MESH)
! calculating the gradient of the free surface in y-direction
       CALL VECTOR(T6,'=','GRADF          Y',IELMU,
     &      1.D0,Z,S,S,S,S,S,MESH,MSK,MASKEL)
! for parallel computing
       IF (NCSIZE.GT.1) CALL PARCOM (T6, 2, MESH)
! calculating the mass-matrix
      CALL VECTOR(T4,'=','MASBAS          ',IELMU,
     &      1.D0,S,S,S,S,S,S,MESH,MSK,MASKEL)
! for parallel computing
       IF (NCSIZE.GT.1) CALL PARCOM (T4, 2, MESH)
! for the weak formulation in FEM, there must be divided by the mass-matrix
       CALL OS ('X=Y/Z   ', T5,T5,T4,C,2,0.D0,1.D-12)
       CALL OS ('X=Y/Z   ', T6,T6,T4,C,2,0.D0,1.D-12)
!
!
! calculation of the x- and y-parts of the secondary current according to Engelung
! tau_x_sec = C*QV, tau_y_sec = C*QU
!
! at the moment alpha must be set here  (0,75 for very rough bottoms, 1 for smooth ones)
! attention: the variable alpha is more than the alpha from the theory
      ALPHA = 1.0D0
      ALPHA = 7.D0 / ALPHA * XMVE *GRAV
!     WRITE(LU,*)'ALPHA',1.D0/ALPHA*7.D0*GRAV*XMVE
!
!
      CALL OS( 'X=YZ    ' , T1 , T6      , QU   , C   ) ! dzsdy*QU
      CALL OS( 'X=Y/Z   ' , T1 , T1      , HN   , C   ) ! dzsdy*QU/HN
      CALL OS( 'X=YZ    ' , T2 , T5      , QV   , C   ) ! dzsdx*QV
      CALL OS( 'X=Y/Z   ' , T2 , T2      , HN   , C   ) ! dzsdx*QV/HN
      CALL OS( 'X=-Y    ' , T2 , T2      , T3   , C   )
      CALL OS( 'X=X+Y   ' , T1 , T2      , T3   , C   ) ! QU*dzsdy - QV*dzsdx
!
      CALL OS( 'X=YZ    ' , T2 , QU      , QU   , C   ) ! QU**2
      CALL OS( 'X=Y/Z   ' , T2 , T2      , HN   , C   ) ! QU**2/HN
      CALL OS( 'X=Y/Z   ' , T2 , T2      , HN   , C   ) ! QU**2/HN**2
      CALL OS( 'X=YZ    ' , T3 , QV      , QV   , C   ) ! QV**2
      CALL OS( 'X=Y/Z   ' , T3 , T3      , HN   , C   ) ! QV**2/HN
      CALL OS( 'X=Y/Z   ' , T3 , T3      , HN   , C   ) ! QV**2/HN**2
      CALL OS( 'X=X+Y   ' , T2 , T3      , T3   , C   ) ! QU**2+QV**2
!
      CALL OS('X=Y/Z   ' , T1 , T1 , T2, C ,2 , 0.D0,1.D-12) !(QU*dzsdy - QV*dzsdx)/(QU**2+QV**2)
!
      CALL OS( 'X=CX    ' , T1 , T2      , T3   , ALPHA   ) ! T1 * 7/Alpha*XMVE*GRAV
      CALL OS( 'X=XY    ' , T1 , HN      , T3   , C   ) ! T1*HN
!
! only for Strickler roughness
! T4: chestr as kstr
!      CALL OS( 'X=C     ' , T4 , T2      , T3   ,71.2D0   ) ! set of kstr

!      CALL OS( 'X=XC    ' , T1 , HN      , T3   , GRAV   ) ! T1*HN*GRAV
!      CALL OS( 'X=YZ    ' , T2 , T4  , T4  , C   ) ! Chestr**2
!      CALL OS( 'X=Y/Z   ' , T1 , T1 , T2, C ,2 , 0.D0,1.D-12) ! T1 / chestr**2
!      C = 1.D0/3.D0
!      CALL OS( 'X=Y**C   ' , T2 , HN   , HN   , C   ) ! HN**1/3
!      CALL OS( 'X=Y/Z   ' , T1 , T1 , T2, C ,2 , 0.D0,1.D-12) ! T1 / HN**1/3
!
!
! for all roughness laws 
      CALL OS( 'X=CXY   ' , T1 , CF      , T3   , 0.5D0   )!T1*CF/2
!
! tau_x_sek = -c*qv : T5
! tau_y_sek = c*qu : T6
      CALL OS('X=YZ    ' , T5 , T1    , QV,  C ) ! c*qv
      CALL OS('X=Y/Z   ' , T5 , T5    , HN,  C ) ! c*qv/HN
      CALL OS('X=YZ    ' , T6 , T1    , QU,  C ) ! c*qu
      CALL OS('X=Y/Z   ' , T6 , T6    , HN,  C ) ! c*qu/HN
      CALL OS('X=-Y    ' , T6 , T6    , QV,   C ) ! -c*qu
! sqrt(tau_x_sek**2+tau_y_sek**2) : T3
      CALL OS('X=YZ    ' , T2 , T5    , T5,  C ) ! T2 = (c*qv)**2
      CALL OS('X=YZ    ' , T3 , T6    , T6,  C ) ! T3 = (c*qu)**2
      CALL OS('X=X+Y   ' , T2 , T3    , T3,  C ) ! T2 = (c*qv)**2+(c*qu)**2
      CALL OS('X=SQR(Y)' , T3 , T2    , T3,  C ) ! T3 = sqrt((c*qu)**2+(c*qv)**2
!      print*,'taux',T5%r(1061),T6%r(1061)
!
! tau_x_ges = tob*effpnt*calfa + tau_x_sek : T1
! tau_y_ges = tob*effpnt*salfa + tau_y_sek : T2
      CALL OS( 'X=YZ    ' , T1 , TOB      , COEFPN   , C   ) ! tob*effpnt
      CALL OS( 'X=YZ    ' , T2 , T1      ,  SALFA   , C   ) ! tob*effpnt*salfa
      CALL OS( 'X=YZ    ' , T1 , T1      , CALFA   , C   ) ! tob*effpnt*calfa
      CALL OS('X=X+Y   ' , T1 , T5    , T3,  C ) ! tau_x_ges = tob*calfa+tau_x_sek
      CALL OS('X=X+Y   ' , T2 , T6    , T3,  C ) ! tau_y_ges = tob*salfa+tau_y_sek
!tau_ges=sqrt(tau_x_ges**2+tau_y_ges**2)
      CALL OS( 'X=YZ    ' , T3 , T1      , T1   , C   ) ! tau_x_ges**2
      CALL OS( 'X=YZ    ' , T4 , T2      , T2   , C   ) ! tau_y_ges**2
      CALL OS('X=X+Y   ' , T4 , T3    , T3,  C ) !tau_x_ges**2+tau_y_ges**2
      CALL OS('X=SQR(Y)' , T4 , T4    , T3,  C ) ! sqrt(tau_x_ges**2+tau_y_ges**2)
!
!
! new angle
! calfa_new = cos(tau_x_ges/tau_ges)
! salfa_new = sin(tau_y_ges/tau_ges)
      CALL OS('X=Y/Z   ' , T1 , T1 , T4, C ,2 , 0.D0,1.D-12) !tau_x_ges/tau_ges
      CALL OS('X=Y/Z   ' , T2 , T2 , T4, C ,2 , 0.D0,1.D-12) !tau_y_ges/tau_ges
!
! taken from effpnt ueber
! to be sure, that tau_x_ges/tau_ges are between (-1,1)
       DO i=1,NPOIN
         if(T1%R(i).lt.-1.D0.or.T1%R(i).gt.1.D0.or.
     &      T2%R(i).lt.-1.D0.or.T2%R(i).gt.1.D0) THEN
            print*,'not acceptable border crossing',i
         ENDIF
         T1%R(i) = MIN(T1%R(I),1.D0)
         T1%R(i) = MAX(T1%R(I),-1.D0)  
         T2%R(i) = MIN(T2%R(I),1.D0)
         T2%R(i) = MAX(T2%R(i),-1.D0)
       ENDDO
!
      CALL OS( 'X=Y     ' ,X=CALFA ,Y=T1 ) ! (tau_x_ges/tau_ges)
      CALL OS( 'X=Y     ' ,X=SALFA ,Y=T2 ) ! (tau_y_ges/tau_ges)
!
! coefpn_new = tau_ges / TOB
      CALL OS('X=Y/Z   ' , COEFPN , T4 , TOB, C ,2 , 0.D0,1.D-12) !coefpn=tau_ges/tob
!
!from effpnt
C   TRAITEMENT DES POINTS FRONTIERES A DEBIT IMPOSE :
C  PAS DE MODIFICATION DU QS QUAND IL EST DETERMINIE PAR UTILISATEUR
      DO 10 I = 1 , NPTFR
        IF (LIQBOR%I(I).EQ.5) THEN
          COEFPN%R(MESH%NBOR%I(I)) = 1.D0
        ENDIF
10    CONTINUE
!
      RETURN
      END


