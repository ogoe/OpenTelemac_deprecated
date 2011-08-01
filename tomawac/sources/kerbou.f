C                       ***************
                        FUNCTION KERBOU
C                       ***************
C
     *( XK1   , XK2   , FREQ1  , FREQ2  , DEPTH , TETA1 , TETA2 )
C
C**********************************************************************
C  TOMAWAC - V1.1                           (EDF/DER/LNH)  -  11/06/98
C**********************************************************************
C
C  FONCTION : CALCUL DU COEFFICIENT DE COUPLAGE.
C  **********
C
C  ARGUMENTS :
C  ***********
C  +-------------+----+--------------------------------------------+
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C  +-------------+----+--------------------------------------------+
C  ! XK1         ! -->! NOMBRE D'ONDE DE 1                         !
C  ! XK2         ! -->! NOMBRE D'ONDE DE 2                         !
C  ! FREQ1       ! -->! FREQUENCE DE L'ONDE 1                      !
C  ! FREQ2       ! -->! FREQUENCE DE L'ONDE 2                      !
C  ! TETA1       ! -->! PLAN DE L'ONDE 1                           !
C  ! TETA2       ! -->! PLAN DE L'ONDE 2                           !
C  ! DEPTH       ! -->! PROFONDEUR                                 !
C  +-------------+----+--------------------------------------------+
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C  +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  QTRIA2
C  ********    - PROGRAMME(S) APPELE(S) :
C
C  REMARQUES :
C  ***********
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      DOUBLE PRECISION  XK1, XK2, FREQ1, FREQ2 , TETA1 , TETA2
      DOUBLE PRECISION  DEPTH
      DOUBLE PRECISION  KERBOU
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      DOUBLE PRECISION  VAR1 , VAR2 , VAR3 , VAR4 , DANG
      DOUBLE PRECISION  PI , DEUPI , GRAVIT
      PARAMETER(PI=3.141592654D0, DEUPI=2.D0*PI )
      PARAMETER(GRAVIT=9.81D0)
C
C
      VAR1  = XK1*XK1
      VAR2  = XK2*XK2
      VAR3  = XK1*XK2
      VAR4  = DEUPI*DEUPI*FREQ1*FREQ2
      DANG  = DCOS(TETA1-TETA2)
C
      KERBOU = GRAVIT*0.5D0*(VAR1+VAR2+2.D0*VAR3*DANG) +
     *         (VAR4/VAR3)*((VAR1+VAR2)*DANG+VAR3*(1+DANG*DANG))/
     *         DEPTH
C
      RETURN
      END
