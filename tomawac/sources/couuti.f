C                       *****************
                        SUBROUTINE COUUTI
C                       *****************
C
     *(X,Y,NPOIN,NCOU, BINCOU,NBOR,NPTFR,AT,DDC,TC1,TC2,
     * NP,XRELC,YRELC,UR,VR,UC1,VC1,UC2,VC2,NPMAX)
C
C***********************************************************************
C  TOMAWAC VERSION 1.0    30/08/95        F.MARCOS     (LNH) 30 87 72 66
C***********************************************************************
C
C   FONCTION : CE SOUS-PROGRAMME PERMET LA LECTURE DES COURANTS
C               A UN FORMAT DEFINI PAR L'UTILISATEUR
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    X,Y         ! -->!  COORDONNEES DU MAILLAGE                     !
C !    NPOIN       ! -->!  NOMBRE DE POINTS DU MAILLAGE                !
C !    NCOU        ! -->!  No D'UNITE LOGIQUE DU FICHIER DES COURANTS  !
C !    BINCOU      ! -->!  BINAIRE DU FICHIER DES COURANTS             !
C !    NBOR        ! -->!  NUMERO GLOBAUX DES POINTS DE BORD           !
C !    NPTFR       ! -->!  NOMBRE DE POINTS DE BORD                    !
C !    AT          ! -->!  TEMPS                                       !
C !    DDC         ! -->!  DATE DE DEBUT DU CALCUL                     !
C !    TC1         ! -->!  DATE CORRESPONDANT AU COURANT (UC1,VC1)     !
C !    TC2         ! -->!  DATE CORRESPONDANT AU COURANT (UC2,VC2)     !
C !    NP          ! -->!  NOMBRE DE POINTS RELEVES                    !
C !    XRELC       !<-->!  TABLEAU DES ABSCISSES DES POINTS RELEVES    !
C !    YRELC       !<-->!  TABLEAU DES ORDONNEES DES POINTS RELEVES    !
C !    UR,VR       !<-->!  TABLEAU DES COURANTS RELEVES                   !
C !    UC1,VC1     !<-->!  TABLEAU DES COURANTS RELEVES AU TEMPS 1        !
C !    UC2,VC2     !<-->!  TABLEAU DES COURANTS RELEVES AU TEMPS 2        !
C !    NPMAX       ! -->!  NOMBRE DE POINTS RELEVES MAXIMUM            !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : 
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
      INTEGER NCOU,NPOIN,NPMAX,NP,NPTFR,NBOR(NPTFR,2)
C
      DOUBLE PRECISION X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION XRELC(NPMAX),YRELC(NPMAX), UR(NPMAX),VR(NPMAX)
      DOUBLE PRECISION UC1(NPMAX),VC1(NPMAX),UC2(NPMAX),VC2(NPMAX)
      DOUBLE PRECISION AT,DDC,TC1,TC2
C
      CHARACTER*3 BINCOU
C-----------------------------------------------------------------------
C
C PROGRAMMER ICI LA LECTURE DE VOTRE FICHIER DE COURANT
C
C     POUR LE PREMIER PASSAGE VOUS DEVEZ IDENTIFIER LES TEMPS TC1 ET TC2
C     QUI ENCADRENT LE PREMIER PAS DE TEMPS. ENSUITE, A L'AIDE DES 
C     TABLEAUX XRELC,YRELC,UR,VR, LUS OU DEDUIS DU FICHIER DE COURANT ,
C     VOUS DEVEZ (EVENTUELLEMENT) INTERPOLER LES COURANTS LUS SUR LE FICHIER
C     DANS LES TABLEAUX UC1,VC1 UC2,VC2. 
C     PROCEDURE D'INTERPOLATION FASP :
C     CALL FASP(X,Y,UC1,NPOIN,XRELC,YRELC,UR,NP,NBOR,MESH%KP1BOR%I,NPTFR,0.D0)
C     CALL FASP(X,Y,VC1,NPOIN,XRELC,YRELC,VR,NP,NBOR,MESH%KP1BOR%I,NPTFR,0.D0)
C     LE CODE INTERPOLERA AUTOMATIQUEMENT LE COURANT ENTRE LES 2 PAS DE 
C     TEMPS
C
C     LES AUTRES PASSAGES SONT EFFECTUES LORSQU'ON A BESOIN D'UN NOUVEL 
C     ENREGISTREMENT (AT>TC2).IL FAUT UNIQUEMENT ICI CALCULER TC2,UC2,VC2
C
C
C PROGRAM HERE THE READING OF THE CURRENTS FILE
C
C     DURING THE FIRST PASSAGE YOU MUST IDENTIFY THE TIMES TC1 AND TC2
C     WHICH SURROUND THE FIRST TIME STEP. NEXT, WITH THE  
C     TABLES XRELC,YRELC,UR,VR, READ OR WITH THE CURRENT FILE ,
C     YOU MUST (EVENTUALLY) INTERPOLATE THE CURRENTS READ ON THE FILE
C     INTO THE TABLES UC1,VC1 UC2,VC2. 
C     INTERPOLATION SUBROUTINE FASP :
C     CALL FASP(X,Y,UC1,NPOIN,XRELC,YRELC,UR,NP,NBOR,MESH%KP1BOR%I,NPTFR,0.D0)
C     CALL FASP(X,Y,VC1,NPOIN,XRELC,YRELC,VR,NP,NBOR,MESH%KP1BOR%I,NPTFR,0.D0)
C     THE CODE WILL INTERPOLATE AUTOMATICALLY THE CURRENT BETWEEN THE 
C     2 TIME STEPS
C
C     THE OTHER PASSAGES ARE MADE WHEN A NEW RECORDING (AT>TC2) IS 
C     REQUIED.IT IS JUST NECESSARY TO CALCULATE HERE TC2,UC2,VC2
C-----------------------------------------------------------------------
C
      IF(LNG.EQ.1) THEN
         WRITE(LU,*) '*********************************************'
         WRITE(LU,*) '  VOUS FAITES APPEL A LA PROCEDURE COUUTI    '
         WRITE(LU,*) '    (FORMAT DU FICHIER DES COURANTS = 3)     '
         WRITE(LU,*) '     MAIS VOUS NE L''AVEZ PAS MODIFIEE       '
         WRITE(LU,*) '*********************************************'
      ELSE
         WRITE(LU,*) '*********************************************'
         WRITE(LU,*) '       YOU CALL THE SUBROUTINE COUUTI        '
         WRITE(LU,*) '        (CURRENTS FILE FORMAT = 3)           '
         WRITE(LU,*) '       BUT YOU DID NOT MODIFIED IT           '
         WRITE(LU,*) '*********************************************'
      ENDIF
      CALL PLANTE(0)
      
      RETURN
      END
