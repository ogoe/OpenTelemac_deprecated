C                       *****************
                        SUBROUTINE BORNES
C                       *****************
C
     *( B     , N     , A     , XM    , X0    , X1    )
C
C**********************************************************************
C  TOMAWAC - V1.1    F. BECQ                 (EDF/DER/LNH)  -  26/03/96
C**********************************************************************
C
C  FONCTION : CALCUL DES BORNES D'INTEGRATION POUR L'INTEGRATION DE LA
C  ********** FONCTION "FONCRO", A L'AIDE DES QUADRATURES DE GAUSS.
C
C  ARGUMENTS :
C  ***********
C +-------------+----+--------------------------------------------+
C ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C +-------------+----+--------------------------------------------+
C ! B           ! -->! PARAMETRE B DE LA FONCTION A INTEGRER      !
C ! N           ! -->! EXPOSANT N  DE LA FONCTION A INTEGRER      !
C ! A           ! -->! PARAMETRE A DE LA FONCTION A INTEGRER      !
C ! XM          ! -->! PARAMETRE M DE LA FONCTION A INTEGRER      !
C ! X0          !<-- ! BORNE INFERIEURE DE L'INTERVALLE           !
C ! X1          !<-- ! BORNE SUPERIEURE DE L'INTERVALLE           !
C +-------------+----+--------------------------------------------+
C ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  QGAUSS
C  ********    - PROGRAMME(S) APPELE(S) :  FONCRO
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
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
C     VARIABLES TRANSMISES.
C     """""""""""""""""""""
      INTEGER  N
      DOUBLE PRECISION B     , A     , XM    , X0    , X1
C
C     VARIABLES LOCALES.
C     """"""""""""""""""
      INTEGER  I0    , I1    , II    , JJ    , IMAX  , INP
      DOUBLE PRECISION X(11) , Y(11) , EPS   , EPS1  , DX
C
C.....FONCTIONS EXTERNES
C     """"""""""""""""""
      DOUBLE PRECISION  FONCRO
      EXTERNAL          FONCRO
C
C
      I1  = 11
      I0  = 1
      X(I0)= 0.D0
      X(I1)= 20.D0
      Y(1) = 0.D0
      EPS1 = 0.01D0
      EPS  = 0.0001D0
      INP  = 0
C
      DO 10 II=1,20
         DX = (X(I1)-X(I0))/10.D0
         X(1) = X(I0)
         IMAX = 0
         I0   = 1
         I1   = 11
         DO JJ=2,11
            X(JJ)=X(JJ-1)+DX
            Y(JJ)=FONCRO(X(JJ),B,N,A,XM)
            IF(Y(JJ).EQ.0.D0.AND.JJ.EQ.2.AND.INP.EQ.0D0) THEN
               X(I1) = X(I1)/10.D0
               INP   = 1
               GOTO 10
            END IF
            IF(Y(JJ).LT.Y(JJ-1)) THEN
               IF(IMAX.EQ.0) THEN
                  IMAX = JJ-1
                  EPS  = EPS1*Y(IMAX)
               END IF
               IF (Y(JJ).LT.EPS) THEN
                  I1 = JJ
                  GOTO 30
               END IF
            ELSEIF(IMAX.EQ.0.AND.Y(JJ).LT.EPS.AND.JJ.NE.2) THEN
               I0 = JJ
            END IF
         END DO
   30    CONTINUE
         IF((I1-I0).LE.2) THEN
            GOTO 10
         ELSE
            GOTO 20
         END IF
   10 CONTINUE
C
   20 CONTINUE
C
      X0 = X(I0)
      X1 = X(I1)
C
      RETURN
      END
