C                       *****************
                        SUBROUTINE RADIAT
C                       *****************
     *( FX    , FY    , SXX   , SXY   , SYY   , XK1   , FS    , CG1   ,
     *  DEPTH1, CGSUC1, DSXXDX, DSXYDX, DSXYDY, DSYYDY, NPOIN2_DIM    )
C
C**********************************************************************
C  TOMAWAC - V1.0    M. BENOIT               (EDF/DER/LNH)  -  13/12/95
C**********************************************************************
C
C  FONCTION : CALCUL DES CONTRAINTES DE RADIATION ET DES FORCES 
C  ********   MOTRICES POUR LA GENERATION DES COURANTS DE HOULE
C             (VOIR REMARQUES POUR MODE DE CALCUL ET UNITE DES FORCES)
C
C  ARGUMENTS :
C  ***********
C  +-------------+----+--------------------------------------------+
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C  +-------------+----+--------------------------------------------+
C  ! FX(-)       !<-- ! COMPOSANTE EN X DE LA FORCE MOTRICE        !
C  ! FY(-)       !<-- ! COMPOSANTE EN Y DE LA FORCE MOTRICE        !
C  ! SXX(-)      !<-- ! CONTRAINTE DE RADIATION SXX                !
C  ! SXY(-)      !<-- ! CONTRAINTE DE RADIATION SXY                !
C  ! SYY(-)      !<-- ! CONTRAINTE DE RADIATION SYY                !
C  ! XK(-,-)     ! -->! TABLEAU DES NOMBRES D'ONDE                 !
C  ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL                       !
C  ! FREQ(-)     ! -->! TABLEAU DES FREQUENCES DE DISCRETISATION   !
C  ! CG(-,-)     ! -->! TABLEAU DES VITESSES DE GROUPE             !
C  ! DFREQ(-)    ! -->! TABLEAU DES PAS DE FREQUENCE               !
C  ! SINTET(-)   ! -->! TABLEAU DES   SINUS DES DIRECTIONS         !
C  ! COSTET(-)   ! -->! TABLEAU DES COSINUS DES DIRECTIONS         !
C  ! DEPTH(-)    ! -->! TABLEAU DES PROFONDEURS                    !
C  ! DEUPI       ! -->! 2.PI                                       !
C  ! GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR               !
C  ! ROEAU       ! -->! MASSE VOLUMIQUE DE L'EAU                   !
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !
C  ! NELEM2      ! -->! NOMBRE D'ELEMENTS DU MAILLAGE SPATIAL      !
C  ! IELM2       ! -->! TYPE D'ELEMENT                             !
C  ! MESH,XMESH  ! -->! STRUCTURE DE MAILLAGE                      !
C  ! X(-)        ! -->! ABSCISSES DES POINTS DU MAILLAGE SPATIAL   !
C  ! Y(-)        ! -->! ORDONNEES DES POINTS DU MAILLAGE SPATIAL   !
C  ! W1(-,-)     !<-->! TABLEAU DE TRAVAIL (DIMENSION NELEM2*4)    !
C  ! CGSURC(-,-) !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2*NF)   !
C  ! DSXXDX(-)   !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C  ! DSXYDX(-)   !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C  ! DSXYDY(-)   !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C  ! DSYYDY(-)   !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C  ! T1(-)       !<-->! TABLEAU DE TRAVAIL STRUCTURE (NPOIN2)      !
C  ! T2(-)       !<-->! TABLEAU DE TRAVAIL STRUCTURE (NPOIN2)      !
C  ! T3(-)       !<-->! TABLEAU DE TRAVAIL STRUCTURE (NPOIN2)      !
C  ! T4(-)       !<-->! TABLEAU DE TRAVAIL STRUCTURE (NPOIN2)      !
C  +-------------+----+--------------------------------------------+
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C  +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  PRE2D
C  ********    - PROGRAMME(S) APPELE(S) :  GRADF
C
C  REMARQUES :
C  ***********
C  - LE CALCUL EST EFFECTUE SELON LA FOMULATION "THEORIQUE", AVEC
C    CALCUL DE TERMES DU TENSEUR DES CONTRAINTES DE RADIATION, PUIS
C    DE LEUR GRADIENTS SPATIAUX.
C
C  - LE RESULATS DE CE CALCUL EST DONNE SOUS LA FORME :
C       FI = - 1/D D( SIJ )/D( XJ )    UNITE : M/S**2
C
C**********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TOMAWAC
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      DOUBLE PRECISION DEUPI
C
      INTEGER NPOIN2_DIM
      DOUBLE PRECISION FS(NPOIN2_DIM,NPLAN,NF),CG1(NPOIN2_DIM,NF)
      DOUBLE PRECISION DEPTH1(NPOIN2_DIM), CGSUC1(NPOIN2_DIM,NF)
      DOUBLE PRECISION XK1(NPOIN2_DIM,NF)
      DOUBLE PRECISION DSXXDX(NPOIN2_DIM),DSXYDX(NPOIN2_DIM),
     *                 FX(NPOIN2_DIM)
      DOUBLE PRECISION DSXYDY(NPOIN2_DIM),DSYYDY(NPOIN2_DIM),
     *                 FY(NPOIN2_DIM)
      DOUBLE PRECISION SXX(NPOIN2_DIM),SXY(NPOIN2_DIM),SYY(NPOIN2_DIM)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  JP    , JF    , IP
      DOUBLE PRECISION COEF  , COCO  , SISI  , SICO  , OMEGA , ROGER
      DOUBLE PRECISION DTETAR, C
C
C
      DEUPI=6.283185307D0
      ROGER=ROEAU*GRAVIT
      DTETAR=DEUPI/DBLE(NPLAN)
C
      DO IP=1,NPOIN2
        SXX(IP) = 0.D0
        SXY(IP) = 0.D0
        SYY(IP) = 0.D0
      ENDDO
C
C.....CALCUL DU TABLEAU (DE TRAVAIL) DES N = CG/C.
C     """"""""""""""""""""""""""""""""""""""""""""
      DO JF = 1,NF
        OMEGA=DEUPI*FREQ(JF)
        DO IP=1,NPOIN2
          CGSUC1(IP,JF)=CG1(IP,JF)*XK1(IP,JF)/OMEGA
        ENDDO
      ENDDO
C
C.....CALCUL DES CONTRAINTES DE RADIATION INTEGREES SUR LE SPECTRE.
C.....(SOMMATIONS SUR LA PARTIE DISCRETISEE DU SPECTRE)
C     """""""""""""""""""""""""""""""""""""""""""""""""
      DO JP=1,NPLAN
        COCO=COSTET(JP)*COSTET(JP)
        SICO=SINTET(JP)*COSTET(JP)
        SISI=SINTET(JP)*SINTET(JP)
        DO JF=1,NF
          COEF=GRAVIT*DFREQ(JF)*DTETAR
          DO IP=1,NPOIN2
            SXX(IP)=SXX(IP)
     &             +(CGSUC1(IP,JF)*(1.D0+SISI)-0.5D0)*FS(IP,JP,JF)*COEF
            SXY(IP)=SXY(IP)
     &             +(CGSUC1(IP,JF)*SICO)*FS(IP,JP,JF)*COEF
            SYY(IP)=SYY(IP)
     &             +(CGSUC1(IP,JF)*(1.D0+COCO)-0.5D0)*FS(IP,JP,JF)*COEF
          ENDDO
        ENDDO
      ENDDO
C
C.....CALCUL DES GRADIENTS SPATIAUX DES CONTRAINTES DE RADIATION
C     """"""""""""""""""""""""""""""""""""""""""""""""""""""""""
C.....PAS DE MASQUAGE D ELEMENTS.
      DO IP=1,NELEM2
        W1(IP)=1.D0
      ENDDO
C
C.....DERIVEES EN X
      CALL OV('X=Y     ',T4,SXX,T3,C,NPOIN2)
      CALL VECTOR
     * (ST1,'=','GRADF          X',IELM2,1.D0,ST4,
     *  ST3,ST3,ST3,ST3,ST3,MESH,.FALSE.,SW1)
C
      CALL OV('X=Y     ',T4,SXY,T3,C,NPOIN2)
      CALL VECTOR
     * (ST2,'=','GRADF          X',IELM2,1.D0,ST4,
     *  ST3,ST3,ST3,ST3,ST3,MESH,.FALSE.,SW1)
C
      CALL VECTOR
     * (ST3,'=','GRADF          X',IELM2,1.D0,MESH%X,
     *  ST3,ST3,ST3,ST3,ST3,MESH,.FALSE.,SW1)
C
!BD_INCKA modif //
       IF(NCSIZE.GT.1) THEN
          CALL PARCOM(ST1,2,MESH)
          CALL PARCOM(ST2,2,MESH) 
          CALL PARCOM(ST3,2,MESH)
       ENDIF
!BD_INCKA fin modif //
      CALL OV('X=Y/Z   ',DSXXDX,T1,T3,C,NPOIN2)
      CALL OV('X=Y/Z   ',DSXYDX,T2,T3,C,NPOIN2)
C
C.....DERIVEES EN Y.
      CALL OV('X=Y     ',T4,SYY,T3,C,NPOIN2)
      CALL VECTOR
     * (ST1,'=','GRADF          Y',IELM2,1.D0,ST4,
     *  ST3,ST3,ST3,ST3,ST3,MESH,.FALSE.,SW1)
C
      CALL OV('X=Y     ',T4,SXY,T3,C,NPOIN2)
      CALL VECTOR
     * (ST2,'=','GRADF          Y',IELM2,1.D0,ST4,
     *  ST3,ST3,ST3,ST3,ST3,MESH,.FALSE.,SW1)
C
      CALL VECTOR
     * (ST3,'=','GRADF          Y',IELM2,1.D0,MESH%Y,
     *  ST3,ST3,ST3,ST3,ST3,MESH,.FALSE.,SW1)
C
!BD_INCKA modif //
       IF(NCSIZE.GT.1) THEN
          CALL PARCOM(ST1,2,MESH)
          CALL PARCOM(ST2,2,MESH) 
          CALL PARCOM(ST3,2,MESH)
       ENDIF
!BD_INCKA fin modif //
      CALL OV('X=Y/Z   ',DSYYDY,T1,T3,C,NPOIN2)
      CALL OV('X=Y/Z   ',DSXYDY,T2,T3,C,NPOIN2)
C
C
C.....CALCUL DES FORCES MOTRICES POUR LES COURANTS DE HOULE
C     """""""""""""""""""""""""""""""""""""""""""""""""""""
      DO IP=1,NPOIN2
        FX(IP)= - (DSXXDX(IP)+DSXYDY(IP))/DEPTH1(IP)
        FY(IP)= - (DSXYDX(IP)+DSYYDY(IP))/DEPTH1(IP)
      ENDDO
C
      RETURN
      END
