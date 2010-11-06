C                       *****************
                        SUBROUTINE INCIDE
C                       *****************
C
     *(COTOND,H,C0,PATMOS,ATMOS,ZF,MESH,LT,AT,GRAV,ROEAU,PRIVE)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2    17/08/94    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C     FONCTION  :  CALCUL DE L'ONDE INCIDENTE IMPOSEE AUX BORDS
C
C     EN CHAQUE POINT D'ONDE INCIDENTE, IL FAUT DONNER :
C
C
C     ( 1 - NINC . NBOR ) A COS ( PHI - OMEGA T)
C
C     AVEC :
C
C     NINC  : DIRECTION DE L'ONDE INCIDENTE.
C     NBOR  : NORMALE A LA PAROI (XSGBOR ET YSGBOR).
C     A     : AMPLITUDE DE L'ONDE.
C     OMEGA : PULSATION DE L'ONDE.
C     PHI   : PHASE DE L'ONDE.
C     T     : TEMPS.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  COTOND        |<-- | ONDE RESULTAT.
C |  H             | -->| HAUTEUR D'EAU.
C |  C0            | -->| CELERITE DE REFERENCE
C |  PATMOS        | -->| PRESSION ATMOSPHERIQUE
C |  ATMOS         | -->| LOGIQUE INDIQUANT SI PATMOS EST REMPLI.
C |  ZF            | -->| FOND
C |  MESH          | -->| STRUCTURE DU MAILLAGE
C |  LT,AT         | -->| NUMERO DE L'ITERATION,TEMPS
C |  GRAV          | -->| PESANTEUR
C |  ROEAU         | -->| MASSE VOLUMIQUE DE L'EAU
C |  PRIVE         | -->| TABLEAU DE TRAVAIL DEFINI DANS PRINCI
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PROPAG
C
C SOUS-PROGRAMME APPELE : OV
C
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
      INTEGER, INTENT(IN)            :: LT
      LOGICAL, INTENT(IN)            :: ATMOS 
      DOUBLE PRECISION, INTENT(IN)   :: AT,GRAV,ROEAU
      TYPE(BIEF_OBJ), INTENT(IN)     :: PATMOS,H,C0,ZF,PRIVE
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: COTOND
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION PI                 
C
C-----------------------------------------------------------------------
C  
      CALL OS( 'X=C     ' , X=COTOND , C=0.D0 )
C
C     PI = 3.141592653589D0
C     T=200.D0
C     W=2.*PI/T
C     A=0.25
C
C      DO 10 K=261,271
C       COTOND%R(K) = 2.*A*SIN(W*AT)
C10    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
