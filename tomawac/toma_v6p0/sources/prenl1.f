C                       *****************
                        SUBROUTINE PRENL1
C                       *****************
C
     *( IANGNL, COEFNL, NPLAN , NF    , RAISF , XLAMD )
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 1.2 ----                         C
C                                                                     C
C  PRENL1 :  PREPARATION DU CALCUL DU TERME SOURCE D'INTERACTIONS     C
C  ********  NON-LINEAIRES ENTRE QUADRUPLETS DE FREQUENCES A L'AIDE   C
C            DE LA METHODE "DISCRETE INTERACTION APPROXIMATION (DIA)" C
C            PROPOSEE PAR HASSELMANN ET HASSELMANN (1985).            C
C            PROCEDURE SPECIFIQUE AU CAS OU LES FREQUENCES SONT EN    C
C            PROGRESSION GEOMETRIQUE ET LES DIRECTIONS REGULIEREMENT  C
C            ESPACEES SUR [0;2.PI].                                   C
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
C  ! IANGNL(-,-) !<-- ! TABLEAU DES INDICES ANGULAIRES POUR DIA    !  C
C  ! COEFNL(-)   !<-- ! VECTEUR DES COEFFICIENTS DE CALCUL POUR DIA!  C
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !  C
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !  C
C  ! RAISF       ! -->! RAISON FREQUENTIELLE DE DISCRETISATION     !  C
C  ! XLAMD       ! -->! COEFFICIENT LAMBDA DE LA CONFIGUARTION STD !  C
C  +-------------+----+--------------------------------------------+  C
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !  C
C  +---------------------------------------------------------------+  C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  WAC                        C
C  ********    - PROGRAMME(S) APPELE(S) :  ANGLES, INTANG             C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  REMARQUES :                                                        C
C  ***********                                                        C
C                                                                     C
C  - CETTE PROCEDURE EST DESTINEE A ETRE UTILISEE EN CONJONCTION AVEC C
C    LE SOUS-PROGRAMME QNLIN1 DONT ELLE OPTIMISE LE FONCTIONNEMENT.   C
C                                                                     C
C  - REFERENCES PRINCIPALES SUR LA METHODE DIA :                      C
C         * HASSELMANN S., HASSELMANN K. (1985) : COMPUTATIONS        C
C               AND PARAMETERIZATIONS OF THE NONLINEAR ENERGY         C
C               TRANSFER IN GRAVITY-WAVE SPECTRUM. PART1 : A NEW      C
C               METHOD FOR EFFICIENT COMPUTATION OF THE EXACT NON-    C
C               LINEAR TRANSFER INTEGRAL. JPO, VOL 15, PP 1369-1377.  C
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
C.....VARIABLES TRANSMISES.
C     """""""""""""""""""""
      INTEGER  NPLAN , NF
      INTEGER  IANGNL(NPLAN,8)
      DOUBLE PRECISION RAISF , XLAMD
      DOUBLE PRECISION COEFNL(16)
C
C.....VARIABLES LOCALES.
C     """"""""""""""""""
      INTEGER  JP
      DOUBLE PRECISION DELTA1, DELTA2, DTMOIN, DTPLUS, DTETAD
      DOUBLE PRECISION APLUS , AMOIN , BPLUS , BMOIN , FPLUS , FMOIN
C
C
C=====C---------------------------------------------------C
C  1  C CALCULS LIES A L'INTERPOLATION ANGULAIRE          C
C=====C---------------------------------------------------C
C
C.....1.1 DETERMINATION DES DIRECTIONS RESONNANTES.
C         (AVEC LA CONVENTION  0 < DTPLUS < DTMOIN)
C     """""""""""""""""""""""""""""""""""""""""""""
      CALL  ANGLES( XLAMD , DTPLUS, DTMOIN)
C
C.....1.2 INDICES ANGULAIRES POUR LA CONFIGURATION 'STANDARD'.
C         (CORRESPONDANT AU COUPLE (-DTPLUS,DTMOIN))
C     """"""""""""""""""""""""""""""""""""""""""""""""""""""""
      DELTA1=-DTPLUS
      DELTA2= DTMOIN
      DO 110 JP=1,NPLAN
        CALL INTANG( IANGNL(JP,2) , IANGNL(JP,1) , JP , NPLAN , DELTA1)
        CALL INTANG( IANGNL(JP,3) , IANGNL(JP,4) , JP , NPLAN , DELTA2)
  110 CONTINUE
C
C.....1.3 INDICES ANGULAIRES POUR LA CONFIGURATION 'IMAGE'.
C         (CORRESPONDANT AU COUPLE (DTPLUS,-DTMOIN))
C     """""""""""""""""""""""""""""""""""""""""""""""""""""
      DELTA1= DTPLUS
      DELTA2=-DTMOIN
      DO 120 JP=1,NPLAN
        CALL INTANG( IANGNL(JP,5) , IANGNL(JP,6) , JP , NPLAN , DELTA1)
        CALL INTANG( IANGNL(JP,8) , IANGNL(JP,7) , JP , NPLAN , DELTA2)
  120 CONTINUE
C
C.....1.4 COEFFICIENTS D'INTERPOLATION ANGULAIRE.
C     """""""""""""""""""""""""""""""""""""""""""
      DTETAD=360.D0/DBLE(NPLAN)
      APLUS=DTPLUS/DTETAD-DBLE(INT(DTPLUS/DTETAD))
      AMOIN=DTMOIN/DTETAD-DBLE(INT(DTMOIN/DTETAD))
C
C
C=====C---------------------------------------------------C
C  2  C CALCULS LIES A L'INTERPOLATION FREQUENTIELLE      C
C=====C---------------------------------------------------C
      FPLUS=DLOG(1.D0+XLAMD)/DLOG(RAISF)
      FMOIN=DLOG(1.D0-XLAMD)/DLOG(RAISF)
      BPLUS=(RAISF**(FPLUS-IDINT(FPLUS)     )-1.D0)/(RAISF-1.D0)
      BMOIN=(RAISF**(FMOIN-IDINT(FMOIN)+1.D0)-1.D0)/(RAISF-1.D0)
C
C
C=====C---------------------------------------------------C
C  3  C AFFECTATION DES COEFFICIENTS POUR QNLIN1          C
C=====C---------------------------------------------------C
      COEFNL( 1)=(1.D0-APLUS) * (1.D0-BPLUS)
      COEFNL( 2)=      APLUS  * (1.D0-BPLUS)
      COEFNL( 3)=(1.D0-APLUS) *       BPLUS
      COEFNL( 4)=      APLUS  *       BPLUS
      COEFNL( 5)=(1.D0-AMOIN) * (1.D0-BMOIN)
      COEFNL( 6)=      AMOIN  * (1.D0-BMOIN)
      COEFNL( 7)=(1.D0-AMOIN) *       BMOIN
      COEFNL( 8)=      AMOIN  *       BMOIN
      COEFNL( 9)=FPLUS
      COEFNL(10)=FMOIN
      COEFNL(11)=1.D0/(1.D0+XLAMD)**4
      COEFNL(12)=1.D0/(1.D0-XLAMD)**4
      COEFNL(13)=DBLE(1)
      COEFNL(14)=DBLE(NF+IDINT(1.D0-FMOIN))
C
      RETURN
      END
