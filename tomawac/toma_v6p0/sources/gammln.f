C                       ***************
                        FUNCTION GAMMLN
C                       ***************
C
     *( XX    , DEUPI )
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 1.2 ----                         C
C                                                                     C
C  DELFRA : CALCUL DU LOGARITHME NEPERIEN DE LA FONCTION GAMMA        C
C  ********   (FONCTION EULERIENNE DE DEUXIEME ESPECE)                C
C                                                                     C
C   - CREE POUR VERSION 1.0  LE 15/11/95 PAR M. BENOIT                C
C   - MOD. POUR VERSION 1.2  LE 07/11/96 PAR M. BENOIT                C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  ARGUMENTS D'APPEL :                                                C
C  *******************                                                C
C  +-------------+----+--------------------------------------------+  C
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! GAMMLN      !<-- ! VALEUR DE LA FONCTION GAMMA                !  C
C  ! XX          ! -->! VALEUR EN LAQUELLE LOG(GAMMA) EST CALCULE  !  C
C  ! DEUPI       ! -->! 2.PI                                       !  C
C  +-------------+----+--------------------------------------------+  C
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !  C
C  +---------------------------------------------------------------+  C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  DELFRA                     C
C  ********    - PROGRAMME(S) APPELE(S) :                             C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  REMARQUES :                                                        C
C  ***********                                                        C
C  - SI XX EST UN ENTIER NOTE N, GAMMA(N) = (N-1)!                    C
C  - REFERENCE : PRESS ET AL. (1989) : NUMERICAL RECIPES. THE ART OF  C
C                SCIENTIFIC COMPUTING. (CF. PP 156-157).              C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      DOUBLE PRECISION GAMMLN, XX    , DEUPI
C
C.....VARIABLES LOCALES 
C     """""""""""""""""
      INTEGER  J
      DOUBLE PRECISION STP   , X     , XC    , TMP   , SER   , AUX
      DOUBLE PRECISION COF(6)
C
C
      COF(1)= 76.180091730D0
      COF(2)=-86.505320330D0
      COF(3)= 24.014098220D0
      COF(4)= -1.231739516D0
      COF(5)=  0.001208580D0
      COF(6)= -0.000005364D0
      STP   =  2.506628275D0
C
      IF (XX.LT.1.D0) THEN
        XC=2.D0-XX
      ELSE
        XC=XX
      ENDIF
      X=XC-1.D0
      TMP=X+5.5D0
      TMP=(X+0.5D0)*DLOG(TMP)-TMP
      SER=1.D0
      DO 11 J=1,6
        X=X+1.D0
        SER=SER+COF(J)/X
   11 CONTINUE
      GAMMLN=TMP+DLOG(STP*SER)
      IF (XX.LT.1D0) THEN
        AUX=0.5D0*DEUPI*(1.D0-XX)
        GAMMLN=DLOG(AUX/SIN(AUX))-GAMMLN
      ENDIF
C
      RETURN
      END
