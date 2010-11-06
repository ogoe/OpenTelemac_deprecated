C                       ***************
                        FUNCTION QGAUSS
C                       ***************
C
     *( B     , N     , A     , XM    )
C
C**********************************************************************
C  TOMAWAC - V1.1    F. BECQ                 (EDF/DER/LNH)  -  26/03/96
C**********************************************************************
C
C  FONCTION : CALCULE L'INTEGRALE DE 0 A INFINI DE LA FONCTION DONNEE 
C  ********** PAR "FONCRO", A L'AIDE DES QUADRATURES DE GAUSS.
C
C  ARGUMENTS :
C  ***********
C +-------------+----+--------------------------------------------+
C ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C +-------------+----+--------------------------------------------+
C ! QGAUSS      !<-- ! VALEUR DE L'INTEGRALE                      !
C ! B           ! -->! PARAMETRE B DE LA FONCTION A INTEGRER      !
C ! N           ! -->! EXPOSANT N  DE LA FONCTION A INTEGRER      !
C ! A           ! -->! PARAMETRE A DE LA FONCTION A INTEGRER      !
C ! XM          ! -->! PARAMETRE M DE LA FONCTION A INTEGRER      !
C +-------------+----+--------------------------------------------+
C ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  QBREK3
C  ********    - PROGRAMME(S) APPELE(S) :  BORNES, FONCRO
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
      DOUBLE PRECISION QGAUSS, B     , A     , XM
C
C     VARIABLES LOCALES.
C     """"""""""""""""""
      INTEGER  J     , I     , NFOIS
      DOUBLE PRECISION XB    , XR    , DX    , DA    , SS    , W(5)
      DOUBLE PRECISION A1    , A2    , A3    , Y2    , X(5)
C
C     FONCTIONS EXTERNES.
C     """""""""""""""""""
      DOUBLE PRECISION  FONCRO
      EXTERNAL          FONCRO
C
      DATA X/.1488743389D0,.4333953941D0,.6794095682D0,
     *       .8650633666D0,.9739065285D0/
      DATA W/.2955242247D0,.2692667193D0,.2190863625D0,
     *       .1494513491D0,.0666713443D0/
C
C
      NFOIS = 1
C
      CALL BORNES
     *( B     , N     , A     , XM    , A2    , A3    )
      QGAUSS = 0.D0
      DA = (A3-A2)/DBLE(NFOIS)
C
      DO I=1,NFOIS
         A1 = A2
         A2 = A2+DA
         XB = 0.5D0*(A1+A2)
         XR = 0.5D0*(A2-A1)
         SS = 0.D0
         DO J=1,5
            DX = XR*X(J)
            SS = SS + W(J)*(FONCRO(XB+DX,B,N,A,XM)
     *                     +FONCRO(XB-DX,B,N,A,XM))
         ENDDO
         Y2 = XR*SS
         QGAUSS = QGAUSS + Y2
      ENDDO
C
      RETURN
      END
