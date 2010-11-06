C                       *****************
                        SUBROUTINE MARUTI
C                       *****************
C
     *(X,Y,NPOIN,NMAR, BINMAR,NBOR,NPTFR,AT,DDC,TV1,TV2,
     * NP,XRELV,YRELV,ZR,Z1,Z2,NPMAX)
C
C***********************************************************************
C  TOMAWAC VERSION 1.0    30/08/95        F.MARCOS     (LNH) 30 87 72 66
C***********************************************************************
C
C   FONCTION : CE SOUS-PROGRAMME PERMET LA LECTURE DES MAREES
C               A UN FORMAT DEFINI PAR L'UTILISATEUR
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    X,Y         ! -->!  COORDONNEES DU MAILLAGE                     !
C !    NPOIN       ! -->!  NOMBRE DE POINTS DU MAILLAGE                !
C !    NMAR        ! -->!  NUMERO D'UNITE LOGIQUE DU FICHIER DES MAREES!
C !    BINMAR      ! -->!  BINAIRE DU FICHIER DES MAREES               !
C !    NBOR        ! -->!  NUMERO GLOBAUX DES POINTS DE BORD           !
C !    NPTFR       ! -->!  NOMBRE DE POINTS DE BORD                    !
C !    AT          ! -->!  TEMPS                                       !
C !    DDC         ! -->!  DATE DE DEBUT DU CALCUL                     !
C !    TV1         ! -->!  DATE CORRESPONDANT A LA MAREE Z1            !
C !    TV2         ! -->!  DATE CORRESPONDANT A LA MAREE Z2            !
C !    NP          ! -->!  NOMBRE DE POINTS RELEVES                    !
C !    XRELV       !<-->!  TABLEAU DES ABSCISSES DES POINTS RELEVES    !
C !    YRELV       !<-->!  TABLEAU DES ORDONNEES DES POINTS RELEVES    !
C !    ZR          !<-->!  TABLEAU DES MAREES RELEVEES                 !
C !    Z1          !<-->!  TABLEAU DES MAREES RELEVEES AU TEMPS 1      !
C !    Z2          !<-->!  TABLEAU DES MAREES RELEVEES AU TEMPS 2      !
C !    NPMAX       ! -->!  NOMBRE DE POINTS RELEVES MAXIMUM            !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : LECHAM,NOUMAR
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      USE DECLARATIONS_TOMAWAC ,ONLY : MESH
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
      INTEGER NMAR,NPOIN,NPMAX,NP,NPTFR,NBOR(NPTFR,2)
C
      DOUBLE PRECISION X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION XRELV(NPMAX),YRELV(NPMAX), ZR(NPMAX)
      DOUBLE PRECISION Z1(NPMAX),Z2(NPMAX)
      DOUBLE PRECISION AT,DDC,TV1,TV2
C
      CHARACTER*3 BINMAR
C-----------------------------------------------------------------------
C
C PROGRAMMER ICI LA LECTURE DE VOTRE FICHIER DE MAREE
C
C     POUR LE PREMIER PASSAGE VOUS DEVEZ IDENTIFIER LES TEMPS TV1 ET TV2
C     QUI ENCADRENT LE PREMIER PAS DE TEMPS. ENSUITE, A L'AIDE DES 
C     TABLEAUX XRELV,YRELV,UR,VR, LUS OU DEDUIS DU FICHIER DE MAREE ,
C     VOUS DEVEZ (EVENTUELLEMENT) INTERPOLER LES MAREES LUES SUR LE FICHIER
C     DANS LES TABLEAUX U1,V1 U2,V2. 
C     PROCEDURE D'INTERPOLATION FASP :
C     CALL FASP(X,Y,Z1,NPOIN,XRELV,YRELV,ZR,NP,NBOR,MESH%KP1BOR%I,NPTFR,0.D0)
C     LE CODE INTERPOLERA AUTOMATIQUEMENT LA MAREE ENTRE LES 2 PAS DE 
C     TEMPS
C
C     LES AUTRES PASSAGES SONT EFFECTUES LORSQU'ON A BESOIN D'UN NOUVEL 
C     ENREGISTREMENT (AT>TV2).IL FAUT UNIQUEMENT ICI CALCULER TV2,Z2
C-----------------------------------------------------------------------
C
      WRITE(LU,*) '*********************************************'
      WRITE(LU,*) '  VOUS FAITES APPEL A LA PROCEDURE MARUTI    '
      WRITE(LU,*) '    (FORMAT DU FICHIER DES MAREES = 3)     '
      WRITE(LU,*) '     MAIS VOUS NE L''AVEZ PAS MODIFIEE       '
      WRITE(LU,*) '*********************************************'
      CALL PLANTE(0)
      
      RETURN
      END
