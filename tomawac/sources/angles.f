C                       *****************
                        SUBROUTINE ANGLES
C                       *****************
C
     *( XLAMD , DTPLUS, DTMOIN)
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 1.2 ----                         C
C                                                                     C
C  ANGLES :  CALCUL DES ANGLES DES VECTEURS RESONNANTS POUR LE CAS DE C
C  ********  DE LA CONFIGURATION D'INTERACTION STANDARD DE LA METHODE C
C            "DISCRETE INTERACTION APPROXIMATION (DIA)" PROPOSEE PAR  C
C            HASSELMANN ET HASSELMANN (1985).                         C
C            PROCEDURE SPECIFIQUE AU CAS OU LES DIRECTIONS SONT       C
C            REPARTIES REGULIEREMENT SUR [0;2.PI].                    C
C                                                                     C
C   - CREE POUR VERSION 1.2  LE 26/06/96 PAR M. BENOIT                C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  ARGUMENTS D'APPEL :                                                C
C  *******************                                                C
C  +-------------+----+--------------------------------------------+  C
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! XLAMD       ! -->! COEFFICIENT LAMBDA DE LA CONFIGUARTION STD !  C
C  ! DTPLUS      !<-- ! ANGLE ASSOCIE A LA FREQUENCE (1+XLAMD).FREQ!  C
C  ! DTMOIN      !<-- ! ANGLE ASSOCIE A LA FREQUENCE (1-XLAMD).FREQ!  C
C  +-------------+----+--------------------------------------------+  C
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !  C
C  +---------------------------------------------------------------+  C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  PRENL1                     C
C  ********    - PROGRAMME(S) APPELE(S) :    -                        C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  REMARQUES :                                                        C
C  ***********                                                        C
C                                                                     C
C  - LES ANGLES DTPLUS ET DTMOIN SONT DONNES EN DEGRES ET SONT PAR    C
C    CONVENTION TOUS DEUX POSITIFS.                                   C
C                                                                     C
C  - REFERENCE PRINCIPALE :                                           C
C         * HASSELMANN S., HASSELMANN K. ET AL.(1985) : COMPUTATIONS  C
C               AND PARAMETERIZATIONS OF THE NONLINEAR ENERGY         C
C               TRANSFER IN GRAVITY-WAVE SPECTRUM. PART2 : PARAME-    C
C               TERIZATIONS OF THE NONLINEAR ENERGY TRANSFER FOR      C
C               APPLICATION IN WAVE MODELS. JPO, VOL 15, PP 1378-1391 C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
C.....VARIABLES TRANSMISES.
C     """""""""""""""""""""
      DOUBLE PRECISION XLAMD , DTPLUS, DTMOIN
C
C.....VARIABLES LOCALES.
C     """"""""""""""""""
      DOUBLE PRECISION CNVD  , AUX
C
C
C.....1. VERIFICATION QUE LAMBDA EST COMPRIS ENTRE 0 ET 0.5.
C     """"""""""""""""""""""""""""""""""""""""""""""""""""""
      IF ((XLAMD.LT.0.D0).OR.(XLAMD.GT.0.5D0)) THEN
        IF(LNG.EQ.1) THEN
           WRITE(LU,1001) XLAMD
        ELSE
           WRITE(LU,1002) XLAMD
        ENDIF
        STOP
      ENDIF
C
C.....2. CALCUL DE TETA_PLUS (DTPLUS) ET TETA_MOINS (DTMOIN).
C     """""""""""""""""""""""""""""""""""""""""""""""""""""""
      CNVD=180.D0/3.141592654D0
      AUX=2.D0*XLAMD*(1.D0+XLAMD*XLAMD)
      DTPLUS=ACOS( (1.D0+AUX)/(1.D0+XLAMD)**2 )*CNVD
      DTMOIN=ACOS( (1.D0-AUX)/(1.D0-XLAMD)**2 )*CNVD
C
C
 1001 FORMAT('/!/-----------------------------------------------/!/'/
     *       '/!/  ARRET DU PROGRAMME DANS SUBROUTINE ANGLES    /!/'/
     *       '/!/-----------------------------------------------/!/'/
     *       '/!/  LA VALEUR CHOISIE POUR LE PARAMETRE LAMBDA   /!/'/
     *       '/!/  INTERVENANT DANS LA METHODE DISCRETE INTER-  /!/'/
     *       '/!/  ACTION APPROXIMATION DOIT ETRE INCLUSE DANS  /!/'/
     *       '/!/  L INTERVALLE [ 0. ; 0.5 ].                   /!/'/
     *       '/!/  OR LA VALEUR UTILISEE EST : ', G10.4,'       /!/'/
     *       '/!/-----------------------------------------------/!/')
C
 1002 FORMAT('/!/-----------------------------------------------/!/'/
     *       '/!/        PROGRAM STOP IN SUBROUTINE ANGLES      /!/'/
     *       '/!/-----------------------------------------------/!/'/
     *       '/!/       THE VALUE OF THE LAMBDA PARAMETER       /!/'/
     *       '/!/             USED IN THE DIA METHOD            /!/'/
     *       '/!/  MUST BE INCLUDED IN THE INTERVAL [ 0. ; 0.5 ]/!/'/
     *       '/!/  THE VALUE HERE IMPOSED IS : ', G10.4,'       /!/'/
     *       '/!/-----------------------------------------------/!/')
C
      RETURN
      END
