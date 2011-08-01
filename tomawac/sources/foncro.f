C                       ***************
                        FUNCTION FONCRO
C                       ***************
C
     *( X     , B     , N     , A     , XM    )
C
C**********************************************************************
C  TOMAWAC - V1.1    F. BECQ                 (EDF/DER/LNH)  -  26/03/96
C**********************************************************************
C
C  FONCTION : CALCUL DE LA VALEUR DE LA FONCTION A INTEGRER POUR LE
C  ********** DEFERLEMENT PAR LA METHODE DE ROELVINK (1993).
C
C  ARGUMENTS :
C  ***********
C +-------------+----+--------------------------------------------+
C ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C +-------------+----+--------------------------------------------+
C ! FONCRO      !<-- ! VALEUR DE LA FONCTION                      !
C ! X           ! -->! VALEUR A LAQUELLE LA FONCTION EST EVALUEE  !
C ! B           ! -->! PARAMETRE B DE LA FONCTION A EVALUER       !
C ! N           ! -->! EXPOSANT N  DE LA FONCTION A EVALUER       !
C ! A           ! -->! PARAMETRE A DE LA FONCTION A EVALUER       !
C ! XM          ! -->! PARAMETRE M DE LA FONCTION A EVALUER       !
C +-------------+----+--------------------------------------------+
C ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  QGAUSS, BORNES
C  ********    - PROGRAMME(S) APPELE(S) :  - 
C
C  REMARQUES :
C  ***********
C
C  REFERENCES :
C  ************
C
C**********************************************************************
C
      IMPLICIT NONE
C
C     VARIABLES TRANSMISES.
C     """""""""""""""""""""
      INTEGER  N
      DOUBLE PRECISION X      , B     , A     , XM    , FONCRO
C
C     VARIABLES LOCALES.
C     """"""""""""""""""
      DOUBLE PRECISION AUX
C
C
      AUX   = A*X**XM
      FONCRO= XM*AUX*DEXP(-AUX)*(1.D0-DEXP(-(B*X)**N))
C
      RETURN
      END
