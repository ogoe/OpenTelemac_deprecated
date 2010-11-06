C                       *****************
                        SUBROUTINE ANACOS
C                       *****************
C
     *( UC    , VC    , X     , Y     , NPOIN2 ) 
C
C***********************************************************************
C  TOMAWAC VERSION 5.2    07/06/01       
C***********************************************************************
C
C     FONCTION  : PERMET LA SPECIFICATION D'UN COURANT ANALYTIQUE 
C                 (! STATIONNAIRE !)
C
C     FUNCTION  : SPECIFICATION OF AN ANALYTICAL CURRENT 
C                 (! STATIONNARY !)
C                 
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    UC,VC       !<-- ! COMPOSANTES DU CHAMP DE COURANT              !
C !    X,Y         ! -->! COORDONNEES DES POINTS DU MAILLAGE 2D        !
C !    NPOIN2      ! -->! NOMBRE DE POINTS 2D                          !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : CONDIW
C
C  SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NPOIN2
      DOUBLE PRECISION X (NPOIN2) , Y (NPOIN2)
      DOUBLE PRECISION UC(NPOIN2) , VC(NPOIN2)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  IP
      DOUBLE PRECISION UCONST, VCONST
C
C
      UCONST=0.D0
      VCONST=0.D0
C
      DO 100 IP=1,NPOIN2
        UC(IP)=UCONST
        VC(IP)=VCONST
  100 CONTINUE
C
      RETURN
      END
