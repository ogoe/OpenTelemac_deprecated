C                       *****************
                        SUBROUTINE ANAVEN
C                       *****************
C
     *( UV    , VV    , X     , Y     , NPOIN2, AT    , DDC   , VX_CTE,
     *  VY_CTE) 
C
C***********************************************************************
C  TOMAWAC VERSION 1.0    07/06/95       M. BENOIT (LNH) 30 87 72 66
C***********************************************************************
C
C     FONCTION  : PERMET LA SPECIFICATION D'UN VENT ANALYTIQUE
C                 (EVENTUELLEMENT VARIABLE EN TEMPS) 
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    UV,VV       !<-- ! COMPOSANTES DU CHAMP DE VENT INITIAL         !
C !    X,Y         ! -->! COORDONNEES DES POINTS DU MAILLAGE 2D        !
C !    NPOIN2      ! -->! NOMBRE DE POINTS 2D                          !
C !    AT          ! -->! TEMPS DU CALCUL                              !
C !    DDC         ! -->! DATE DE DEBUT DU CALCUL                      !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : CONDIW,SEMIMP
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
      DOUBLE PRECISION AT    , DDC   , VX_CTE, VY_CTE
      DOUBLE PRECISION X (NPOIN2)    , Y (NPOIN2)
      DOUBLE PRECISION UV(NPOIN2)    , VV(NPOIN2)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  IP
C
C
      DO IP=1,NPOIN2
        UV(IP)=VX_CTE
        VV(IP)=VY_CTE
      ENDDO
C
      RETURN
      END
