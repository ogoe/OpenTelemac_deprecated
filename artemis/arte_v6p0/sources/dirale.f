C                       *****************
                        SUBROUTINE DIRALE
C                       *****************
C
     *(DALE,EXPOS,TETAH,TETMIN,TETMAX,NDALE,TRA01,NPOIN,PRIVE,NPRIV)
C
C***********************************************************************
C
C  ARTEMIS VERSION 5.1   02/06/99    D. AELBRECHT (LNH) 01 30 87 74 12
C
C  LINKED TO BIEF VERS. 5.0          J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C      FONCTION:    DISCRETISE UN SPECTRE D'ENERGIE DIRECTIONNEL EN
C                   EN NDALE BANDES D'EGALE ENERGIE. LE RESULTAT
C                   EST LA DONNEE DES DIRECTIONS CHACUNE DES BANDES.
C
C      ON UTILISE LA FORMULE DONNEE PAR GODA DANS 'RANDOM SEAS AND
C      DESIGN OF MARITIME STRUCTURES' - UNIVERSITY OF TOKYO PRESS
C
C      G = ( COS( (TETA-TETAH))/2 ) )**EXPOS
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   DALE         |<-- |  DIRECTIONS DE DISCRETISATION DU SPECTRE     |
C |   EXPOS        | -->|  COEFFICIENT DANS LA FORMULE DU SPECTRE      |
C |   TETAH        | -->|  DIRECTION PRINCIPALE DE PROPAGATION         |
C |   TETMIN       | -->|  VALEUR MINIMUM DE L'ANGLE DE PROPAGATION    |
C |   TETMAX       | -->|  VALEUR MAXIMUM DE L'ANGLE DE PROPAGATION    |
C |   NDALE        | -->|  NOMBRE DE BANDES DE DISCRETISATION          |
C |   TRA01        |<-->|  TABLEAU DE TRAVAIL                          |
C |   NPOIN        | -->|  NOMBRE DE POINTS DU MAILLAGE                |
C |   PRIVE        |<-->|  TABLEAU PRIVE DE L'UTILISATEUR              |
C |   NPRIV        |<-->|  NOMBRE DE TABLEAUX PRIVES                   |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : ARTEMI
C
C***********************************************************************
C
      USE BIEF
      USE INTERFACE_ARTEMIS, EX_DIRALE => DIRALE 
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER NDALE,NPOIN,NPAS,I,K,NPRIV
C
      DOUBLE PRECISION DALE(NDALE),TRA01(NPOIN)
      DOUBLE PRECISION EXPOS,TETMIN,TETMAX,TETAH,DTETA,SUMB,VAR,SUMICI
C
      TYPE(BIEF_OBJ) :: PRIVE
C
      DOUBLE PRECISION EXPO
      COMMON /COEFHD/ EXPO
C
c      DOUBLE PRECISION SPD
c      EXTERNAL SPD
C
C-----------------------------------------------------------------------
C
C ON PASSE EXPOS PAR COMMON DANS LA FONCTION SPD (ON NE PEUT PAS
C L'APPELER EXPOS CAR C'EST ICI UN ARGUMENT DU SOUS-PROGRAMME)
      EXPO = EXPOS
C
C-----------------------------------------------------------------------
C
C     NOMBRE D'INTERVALLES D'INTEGRATION POUR LA METHODE DES TRAPEZES
      NPAS = 2000*NDALE
C
C     LONGUEUR D'UN INTERVALLE D'INTEGRATION
      DTETA = (TETMAX-TETMIN)/FLOAT(NPAS)
C
C     SURFACE DU SPECTRE
      SUMB = (SPD(TETMIN-TETAH) + SPD(TETMAX-TETAH))/2.D0
      DO 20 I = 2,NPAS-1
         SUMB = SUMB + SPD(TETMIN-TETAH+FLOAT(I)*DTETA)
20    CONTINUE
C
C     ON DIVISE EN 2*NDALE BANDES D'EGALES ENERGIES
      SUMB = SUMB/FLOAT(2*NDALE)
C
C     ON RECHERCHE LES ANGLES TOUS LES (2*I-1)*SUMB (I=1,NDALE)
      SUMICI = SPD(TETMIN-TETAH)/2.D0
      I   = 1
      DO 30 K=1,NPAS
         VAR = SPD(TETMIN-TETAH+DTETA*FLOAT(K))
         SUMICI = SUMICI + VAR/2.D0
         IF (SUMICI.GE.SUMB*FLOAT(2*I-1)) THEN
            DALE(I) =  TETMIN + DTETA*(FLOAT(K)-0.5D0)
            I = I + 1
            IF (I.GT.NDALE) RETURN
         ENDIF
         SUMICI = SUMICI + VAR/2.D0
30    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
