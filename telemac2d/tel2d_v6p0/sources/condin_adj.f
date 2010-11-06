C                       *********************
                        SUBROUTINE CONDIN_ADJ
C                       *********************
C
     *(ALIRE,NRES,TROUVE)
C
C***********************************************************************
C TELEMAC 2D VERSION 6.0    24/04/2009  J-M HERVOUET TEL: 01 30 87 80 18
C
C***********************************************************************
C
C     FONCTION  : INITIALISATION DES GRANDEURS PHYSIQUES POUR LE DEBUT
C                 DU CALCUL EN MODE ADJOINT
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |                | -- |  
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC2D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: ALIRE(*),NRES
      INTEGER, INTENT(INOUT) :: TROUVE(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,ITER
      DOUBLE PRECISION HIST(1),AT1
C
C-----------------------------------------------------------------------
C
C     CONDIN FOR ADJOINT PROBLEM: T=T(N+1) : P=0,Q=0,R=0
C
      CALL OS( 'X=C     ',PP,PP,PP,0.D0)
      CALL OS( 'X=C     ',QQ,QQ,QQ,0.D0)
      CALL OS( 'X=C     ',RR,RR,RR,0.D0)
C
C     JUST IN CASE CV1,.. IS WRITTEN IN THE RESULT FILE
C
      CALL OS( 'X=C     ',CV1,CV1,CV1,0.D0)
      CALL OS( 'X=C     ',CV2,CV2,CV2,0.D0)
      CALL OS( 'X=C     ',CV3,CV3,CV3,0.D0)
C
C     READING THE LAST TIME IN THE TELEMAC RESULTS FILE (NRES)
C     INITIALISE U,V AND H
C
      REWIND NRES
C
      CALL BIEF_SUITE(VARSOR,VARCL,ITER,NRES,'SERAFIN ',
     *           HIST,0,NPOIN,AT,TEXTE,VARCLA,
     *           NVARCL,TROUVE,ALIRE,LISTIN,.TRUE.,MAXVAR)
C
C     GIVING MEASUREMENTS HD,UD AND VD AT THE LAST TIME STEP
C     (ITER AND AT GIVEN BY THE PREVIOUS CALL TO SUITE)
C
      CALL MESURES(ITER,AT)      
C     INITIALISATION DE HH, UU, VV
C
      CALL OS( 'X=Y     ' , HH   , H , H , 0.D0 )
      CALL OS( 'X=Y     ' , UU   , U , U , 0.D0 )
      CALL OS( 'X=Y     ' , VV   , V , V , 0.D0 )
      CALL OS( 'X=C     ' , HIT1 , HIT1 , HIT1 , 0.D0 )
      CALL OS( 'X=C     ' , UIT1 , UIT1 , UIT1 , 0.D0 )
      CALL OS( 'X=C     ' , VIT1 , VIT1 , VIT1 , 0.D0 )
C         
C     READING OF TELEMAC2D RESULTS (RESULTS FILE - UNIT NRES)
C     THIS IS TO HAVE UN OF THE LAST TIME STEP INTO U.
C
C     ATTENTION : SUPPOSE QUE NVARRES A ETE CALCULE AVANT
C
      DO I=1,2*(NVARRES+1)
        BACKSPACE NRES
      ENDDO 
      CALL LITENR(VARSOR,VARCL,NRES,'STD',HIST,0,NPOIN,AT1,TEXTE,
     *           TEXRES,NVARRES,VARCLA,0,TROUVE,ALIRE,W,.FALSE.,MAXVAR)
C
C-----------------------------------------------------------------------
C
      AT = AT + DT
C
C-----------------------------------------------------------------------
C
      RETURN
      END
