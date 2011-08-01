C                       *****************
                        SUBROUTINE CQUEUE
C                       *****************
C
     *( NF    , RAISF , TAILF , JFRE  , JBIS  , COEF1 )
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 1.2 ----                         C
C                                                                     C
C  CQUEUE :  AJUSTEMENT DES INDICES FREQUENTIELS ET CALCUL DES COEFFI-C
C  ********  CIENTS DE QUEUE SUR LA PARTIE NON-DISCRETISEE DU SPECTRE C
C            POUR LA METHODE "DISCRETE INTERACTION APPROXIMATION"     C
C            PROPOSEE PAR HASSELMANN ET HASSELMANN (1985).            C
C            PROCEDURE SPECIFIQUE AU CAS OU LES FREQUENCES SONT EN    C
C            PROGRESSION GEOMETRIQUE.                                 C
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
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !  C
C  ! RAISF       ! -->! RAISON FREQUENTIELLE DE DISCRETISATION     !  C
C  ! TAILF       ! -->! FACTEUR DE QUEUE                           !  C
C  ! JFRE        ! -->! INDICE FREQUENTIEL DE LA COMPOSANTE FREQ.  !  C
C  ! JBIS        !<-- ! INDICE AJUSTE DANS L'INTERVALLE [1;NF]     !  C
C  ! COEF1       !<-- ! COEF. MULTIPLICATIF  F(JFRE)=COEF1*F(JBIS) !  C
C  +-------------+----+--------------------------------------------+  C
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !  C
C  +---------------------------------------------------------------+  C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  QNLIN1                     C
C  ********    - PROGRAMME(S) APPELE(S) :    -                        C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  REMARQUES :                                                        C
C  ***********                                                        C
C                                                                     C
C  - LE SPECTRE EST SUPPOSE NUL EN TOUTES LES FREQUENCES INFERIEURES  C
C    A LA PREMIERE FREQUENCES DE DISCRETISATION.                      C
C                                                                     C
C  - AU-DELA DE LA DERNIERE FREQUENCE DE DISCRETISATION LE SPECTRE    C
C    EST SUPPOSE EN DECROISSANCE AVEC UNE FORME EN FREQ**(-TAILF).    C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES.
C     """""""""""""""""""""
      INTEGER  NF    , JFRE  , JBIS
      DOUBLE PRECISION RAISF , TAILF , COEF1
C
C
      IF (JFRE.GT.NF) THEN
        JBIS = NF
        COEF1= 1.D0/RAISF**(DBLE(JFRE-NF)*TAILF)
      ELSEIF (JFRE.LT.1) THEN
        JBIS = 1
        COEF1= 0.D0
      ELSE
        JBIS = JFRE
        COEF1= 1.D0
      ENDIF
C
      RETURN
      END
