C                       *****************
                        SUBROUTINE PERALE
C                       *****************
C
     *(PALE,GAMMA,PERPIC,NPALE,TRA01,NPOIN,PRIVE,NPRIV,PMIN,PMAX)
C
C***********************************************************************
C
C  ARTEMIS VERSION 5.1   02/06/99    D. AELBRECHT (LNH) 01 30 87 74 12
C
C  LINKED TO BIEF VERS. 5.0          J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C      FONCTION:    DISCRETISE UN SPECTRE D'ENERGIE EN NPALE BANDES
C                   D'EGALE ENERGIE. LE RESULTAT EST LA DONNEE DES
C                   PERIODES DE CHACUNE DES BANDES.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   PALE         |<-- |  PERIODES DE DISCRETISATION DU SPECTRE       |
C |   GAMMA        | -->|  COEFFICIENT DANS LA FORMULE DU SPECTRE      |
C |   PERPIC       | -->|  PERIODE DU PIC DU SPECTRE                   |
C |   NPALE        | -->|  NOMBRE DE PERIODES DE DISCRETISATION        |
C |   TRA01        |<-->|  TABLEAU DE TRAVAIL                          |
C |   NPOIN        | -->|  NOMBRE DE POINTS DU MAILLAGE                |
C |   PRIVE        |<-->|  TABLEAU PRIVE DE L'UTILISATEUR              |
C |   NPRIV        |<-->|  NOMBRE DE TABLEAUX PRIVES                   |
C |   PMIN         | -->|  FREQUENCE MINIMUM DU SPECTRE                |
C |   PMAX         | -->|  FREQUENCE MAXIMUM DU SPECTRE                |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : ARTEMI
C
C***********************************************************************
C
      USE INTERFACE_ARTEMIS, EX_PERALE=> PERALE                 
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER NPALE,NPOIN,NPRIV,NPAS,I,K
C
      DOUBLE PRECISION PALE(NPALE),TRA01(NPOIN)
      DOUBLE PRECISION PERPIC,GAMMA ,SUMB  ,SUMICI   ,DF    ,VAR
      DOUBLE PRECISION PMIN  ,PMAX  ,FMIN  ,FMAX  
      DOUBLE PRECISION B     ,B1    ,B2
C
      TYPE(BIEF_OBJ) :: PRIVE
C
      DOUBLE PRECISION FP,GAM,DELTA
      COMMON /COEFHE/ FP,GAM,DELTA
C
C      DOUBLE PRECISION SPE
C      EXTERNAL SPE
C
      INTRINSIC LOG,FLOAT
C
C-----------------------------------------------------------------------
C
C FREQUENCE DE PIC
      FP   = 1.D0 / PERPIC
      FMIN = 1.D0 / PMAX
      FMAX = 1.D0 / PMIN
      IF (FMAX.GE.99.D0) THEN
         FMAX = 2.5D0 * FP
      ENDIF
C
C ON PASSE GAMMA PAR COMMON DANS LA FONCTION SPE (ON NE PEUT PAS
C L'APPELER GAMMA CAR C'EST ICI UN ARGUMENT DU SOUS-PROGRAMME)
      GAM  = GAMMA
C
C-----------------------------------------------------------------------
C
      IF (GAMMA.GT.0.99D0 .AND. GAMMA.LT.1.01D0) THEN
C
C
C        SPECTRE DE PIERSON-MOSKOWITZ
C        ----------------------------
C
         B1 = EXP(-1.25D0 * (FP/FMAX)**4)
         B2 = EXP(-1.25D0 * (FP/FMIN)**4)
         B  = B1 - B2
         DO 10 I=1,NPALE
            PALE(NPALE-I+1) = PERPIC *
     *    (-0.8D0*LOG( B2 + B*FLOAT(2*I-1)/FLOAT(2*NPALE) ))**(0.25D0) 
10       CONTINUE
C
      ELSE
C
C
C        SPECTRE DE JONSWAP
C        ------------------
C
C        LES FREQUENCES LIMITANT LE SPECTRE A GAUCHE ET A DROITE
C        SONT DONNEES PAR MOTS CLES DANS LES ARGUMENTS
C
      IF (FMAX.LE.FP) THEN
         FMAX = 2.5D0 * FP
         WRITE(LU,110) FMAX
 110     FORMAT(/,1X,'(PERALE) : FMAX < FP ??? =>',1X, 
     *          'CORRECTION : FMAX =',1X,F5.3,' Hz',/)
      ENDIF
C
C        NOMBRE D'INTERVALLES D'INTEGRATION POUR LA METHODE DES TRAPEZES
C
         NPAS = 2000*NPALE
C
C        LONGUEUR D'UN INTERVALLE D'INTEGRATION
C
         DF = (FMAX-FMIN)/FLOAT(NPAS)
C
C        COEFFICIENT POUR LA FONCTION DU SPECTRE (CALCULE ICI POUR NE
C        PAS LE RECALCULER A CHAQUE APPEL DE LA FONCTION SPE)
C
         DELTA = 0.0624D0 * FP**4 /
     *           ( 0.230D0 + 0.0336D0*GAM - 0.185D0 / (1.9D0+GAM) )
C
C        SURFACE DU SPECTRE PAR METHODE DES TRAPEZES
C
         SUMB = (SPE(FMIN) + SPE(FMAX))/2.D0
         DO 20 I = 2,NPAS-1
            SUMB = SUMB + SPE(FMIN+FLOAT(I)*DF)
20       CONTINUE
C
C        ON DIVISE EN 2*NPALE BANDES D'EGALES ENERGIES
C
         SUMB = SUMB/FLOAT(2*NPALE)
C
C        ON RECHERCHE LES FREQUENCES TOUS LES (2*I-1)*SUMB (I=1,NPALE)
C
         SUMICI = SPE(FMIN)/2.D0
         I   = 1
         DO 30 K=1,NPAS
            VAR = SPE(FMIN+DF*FLOAT(K))
            SUMICI = SUMICI + VAR/2.D0
            IF (SUMICI.GT.SUMB*FLOAT(2*I-1)) THEN
               PALE(NPALE-I+1) = 1.D0 / ( FMIN + DF*(FLOAT(K)-0.5D0) )
               I = I + 1
               IF (I.GT.NPALE) RETURN
            ENDIF
            SUMICI = SUMICI + VAR/2.D0
30       CONTINUE
C
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
