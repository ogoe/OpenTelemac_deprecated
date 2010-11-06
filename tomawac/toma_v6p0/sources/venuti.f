C                       *****************
                        SUBROUTINE VENUTI
C                       *****************
C
     *(X,Y,NPOIN,NVEN, BINVEN,NBOR,NPTFR,AT,DDC,TV1,TV2,
     * NP,XRELV,YRELV,UR,VR,U1,V1,U2,V2,NPMAX)
C
C***********************************************************************
C  TOMAWAC VERSION 1.0    30/08/95        F.MARCOS     (LNH) 30 87 72 66
C***********************************************************************
C
C   FONCTION : CE SOUS-PROGRAMME PERMET LA LECTURE DES VENTS
C               A UN FORMAT DEFINI PAR L'UTILISATEUR
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    X,Y         ! -->!  COORDONNEES DU MAILLAGE                     !
C !    NPOIN       ! -->!  NOMBRE DE POINTS DU MAILLAGE                !
C !    NVEN        ! -->!  NUMERO D'UNITE LOGIQUE DU FICHIER DES VENTS !
C !    BINVEN      ! -->!  BINAIRE DU FICHIER DES VENTS                !
C !    NBOR        ! -->!  NUMERO GLOBAUX DES POINTS DE BORD           !
C !    NPTFR       ! -->!  NOMBRE DE POINTS DE BORD                    !
C !    AT          ! -->!  TEMPS                                       !
C !    DDC         ! -->!  DATE DE DEBUT DU CALCUL                     !
C !    TV1         ! -->!  DATE CORRESPONDANT AU CHAMP DE VENT U1,V1   !
C !    TV2         ! -->!  DATE CORRESPONDANT AU CHAMP DE VENT U2,V2   !
C !    NP          ! -->!  NOMBRE DE POINTS RELEVES                    !
C !    XRELV       !<-->!  TABLEAU DES ABSCISSES DES POINTS RELEVES    !
C !    YRELV       !<-->!  TABLEAU DES ORDONNEES DES POINTS RELEVES    !
C !    UR,VR       !<-->!  TABLEAU DES VENTS RELEVES                   !
C !    U1,V1       !<-->!  TABLEAU DES VENTS RELEVES AU TEMPS 1        !
C !    U2,V2       !<-->!  TABLEAU DES VENTS RELEVES AU TEMPS 2        !
C !    NPMAX       ! -->!  NOMBRE DE POINTS RELEVES MAXIMUM            !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : INIVEN,NOUVEN
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
      INTEGER NVEN,NPOIN,NPMAX,NP,NPTFR,NBOR(NPTFR,2)
C
      DOUBLE PRECISION X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION XRELV(NPMAX),YRELV(NPMAX), UR(NPMAX),VR(NPMAX)
      DOUBLE PRECISION U1(NPMAX),V1(NPMAX),U2(NPMAX),V2(NPMAX)
      DOUBLE PRECISION AT,DDC,TV1,TV2
C
      CHARACTER*3 BINVEN
C-----------------------------------------------------------------------
C
C PROGRAMMER ICI LA LECTURE DE VOTRE FICHIER DE VENT
C
C     POUR LE PREMIER PASSAGE VOUS DEVEZ IDENTIFIER LES TEMPS TV1 ET TV2
C     QUI ENCADRENT LE PREMIER PAS DE TEMPS. ENSUITE, A L'AIDE DES 
C     TABLEAUX XRELV,YRELV,UR,VR, LUS OU DEDUIS DU FICHIER DE VENT ,
C     VOUS DEVEZ (EVENTUELLEMENT) INTERPOLER LES VENTS LUS SUR LE FICHIER
C     DANS LES TABLEAUX U1,V1 U2,V2. 
C     PROCEDURE D'INTERPOLATION FASP :
C     CALL FASP(X,Y,U1,NPOIN,XRELV,YRELV,UR,NP,NBOR,MESH%KP1BOR%I,NPTFR,0.D0)
C     CALL FASP(X,Y,V1,NPOIN,XRELV,YRELV,VR,NP,NBOR,MESH%KP1BOR%I,NPTFR,0.D0)
C     LE CODE INTERPOLERA AUTOMATIQUEMENT LE VENT ENTRE LES 2 PAS DE 
C     TEMPS
C
C     LES AUTRES PASSAGES SONT EFFECTUES LORSQU'ON A BESOIN D'UN NOUVEL 
C     ENREGISTREMENT (AT>TV2).IL FAUT UNIQUEMENT ICI CALCULER TV2,U2,V2
C-----------------------------------------------------------------------
C
      WRITE(LU,*) '*********************************************'
      WRITE(LU,*) '  VOUS FAITES APPEL A LA PROCEDURE VENUTI    '
      WRITE(LU,*) '    (FORMAT DU FICHIER DES VENTS = 4)        '
      WRITE(LU,*) '     MAIS VOUS NE L''AVEZ PAS MODIFIEE       '
      WRITE(LU,*) '*********************************************'
      CALL PLANTE(0)
      
      RETURN
      END
