C                       *****************
                        SUBROUTINE USTAR1
C                       *****************
C
     *( USTAR , Z0    , TAUWAV, UV    , VV    , CDRAG , ALPHA , XKAPPA,
     *  ZVENT , GRAVIT, NPOIN2)
C
C**********************************************************************
C  TOMAWAC - V1.0    M. BENOIT               (EDF/DER/LNH)  -  25/04/95
C**********************************************************************
C
C  FONCTION : CALCUL DE LA VITESSE DE FROTTEMENT ET DE LA LONGUEUR DE
C  ********** RUGOSITE EN TOUS LES POINTS DU MAILLAGE SPATIAL.
C             UTILISATION DE LA THEORIE DE JANSSEN (1989, 1991).
C
C  ARGUMENTS :
C  ***********
C +-------------+----+--------------------------------------------+
C ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C +-------------+----+--------------------------------------------+
C ! USTAR(-)    !<-- ! TABLEAU DES VITESSES DE FROTTEMENT         !
C ! Z0(-)       !<-- ! TABLEAU DES LONGUEURS DE RUGOSITE          !
C ! TAUWAV(-)   ! -->! TABLEAU DES CONTRAINTES DUES A LA HOULE    !
C ! UV(-)       ! -->! TABLEAU DES COMPOSANTES OUEST-EST DU VENT  !
C ! VV(-)       ! -->! TABLEAU DES COMPOSANTES SUD-NORD  DU VENT  !
C ! CDRAG       ! -->! COEFFICIENT DE TRAINEE                     !
C ! ALPHA       ! -->! CONSTANTE DE LA LOI DE CHARNOCK            !
C ! XKAPPA      ! -->! CONSTANTE DE VON KARMAN                    !
C ! ZVENT       ! -->! COTE A LAQUELLE EST MESURE LE VENT (M)     !
C ! GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR               !
C ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !
C +-------------+----+--------------------------------------------+
C ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  SEMIMP
C  ********    - PROGRAMME(S) APPELE(S) :  TAUTOT
C
C  REMARQUES :
C  ***********
C  - LA DETERMINATION DE TAUT A PARTIR DE UVENT ET TAUW SE FAIT
C    PAR APPEL A TAUTOT.
C
C  REFRENCES : - JANSSEN P.A.E.M (1989) : WIND-INDUCED STRESS AND THE
C              DRAG OF AIR FLOW OVER SEA WAVES. JPO, VOl 19, PP 745-754
C              - JANSSEN P.A.E.M (1991) : QUASI-LINEAR THEORY OF WIND-
C              WAVE GENERATION APPLIED TO WAVE FORECASTING. JPO, VOL 21
C              PP 1631-1642.
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NPOIN2
      DOUBLE PRECISION USTAR(NPOIN2) , Z0(NPOIN2) , TAUWAV(NPOIN2)
      DOUBLE PRECISION    UV(NPOIN2) , VV(NPOIN2)
      DOUBLE PRECISION CDRAG , ALPHA , XKAPPA , ZVENT, GRAVIT
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  ITRMIN, ITRMAX, ITR   , IP
      DOUBLE PRECISION TAUT  , UVENT , TAUW  , USMIN , SEUIL , X
      DOUBLE PRECISION USTEMP
C
C
      USMIN =1.D-6
      SEUIL =1.D-7
      ITRMIN=1
      ITRMAX=15
C
C.....BOUCLE PRINCIPALE SUR LES POINTS DU MAILLAGE SPATIAL.
C     """""""""""""""""""""""""""""""""""""""""""""""""""""
      DO IP=1,NPOIN2
C
C.......CALCUL DE LA CONTRAINTE TOTALE.
C       """""""""""""""""""""""""""""""
        UVENT=SQRT(UV(IP)**2+VV(IP)**2)
        TAUW =TAUWAV(IP)
        CALL TAUTOT
     *( TAUT  , UVENT , TAUW  , CDRAG , ALPHA , XKAPPA, ZVENT , SEUIL ,
     *  GRAVIT, ITR   , ITRMIN, ITRMAX)
C
C.......CALCUL DE LA VITESSE DE FROTTEMENT.
C       """""""""""""""""""""""""""""""""""
        USTAR(IP)=SQRT(TAUT)
C
C.......CALCUL DE LA LONGUEUR DE RUGOSITE.
C       """"""""""""""""""""""""""""""""""
        USTEMP=MAX(USTAR(IP),USMIN)
        X     =MIN(TAUWAV(IP)/USTEMP**2,0.999D0)
        Z0(IP)=ALPHA*USTEMP**2/(GRAVIT*SQRT(1.D0-X))
C
      ENDDO
C
      RETURN
      END
