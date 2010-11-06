C                       *****************
                        SUBROUTINE ANAMAR
C                       *****************
C
     *( UC  , VC  , ZM  , ZM1 , ZM2 , DZHDT , X  , Y  , NPOIN2 ,
     *  AT  , DDC , LT  )
C
C***********************************************************************
C  TOMAWAC VERSION 5.0
C***********************************************************************
C
C     FONCTION  : PERMET LA SPECIFICATION D'UNE MAREE ANALYTIQUE :
C                 HAUTEUR, COURANT VARIABLES EN TEMPS 
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    UC,VC       !<-- ! COMPOSANTES DU CHAMP DE COURANT DE LA MAREE  !
C !    ZM1         !<-- ! HAUTEUR DE LA MAREE PAR RAPPORT A ZREPOS A T !
C !    ZM2         !<-- ! HAUTEUR DE LA MAREE PAR RAPPORT A ZREPOS A T2!
C !    DZHDT       !<-- ! VARIATION TEMPORELLE DE LA HAUTEUR DE MAREE  !
C !    X,Y         ! -->! COORDONNEES DES POINTS DU MAILLAGE 2D        !
C !    NPOIN2      ! -->! NOMBRE DE POINTS 2D                          !
C !    AT          ! -->! TEMPS DU CALCUL                              !
C !    DDC         ! -->! DATE DE DEBUT DU CALCUL                      !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : CONDIW,WAC
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
      INTEGER  LT
      DOUBLE PRECISION AT    , DDC
      DOUBLE PRECISION X (NPOIN2), Y (NPOIN2), ZM1(NPOIN2), ZM2(NPOIN2)
      DOUBLE PRECISION UC(NPOIN2), VC(NPOIN2), DZHDT(NPOIN2),ZM(NPOIN2)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER          IP, I, J
      DOUBLE PRECISION UCONST, VCONST
C
C-----------------------------------------------------------------------
C     EXEMPLE 1
C-----------------------------------------------------------------------
C
      UCONST=0.D0
      VCONST=0.D0
C
      DO 100 IP=1,NPOIN2
        UC(IP)   = UCONST
        VC(IP)   = VCONST
        ZM(IP)   = 0.D0
        DZHDT(IP)= 0.D0
  100 CONTINUE
C
C
      RETURN
      END
