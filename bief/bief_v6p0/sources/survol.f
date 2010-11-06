C                       *****************
                        SUBROUTINE SURVOL
C                       *****************
C
     *(SURFAC, XEL,YEL,ZEL,NELEM,NELMAX,IELM)
C
C***********************************************************************
C BIEF VERSION 5.1           05/02/91    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : CALCUL DE LA SURFACE (VOLUME) DES ELEMENTS D'UN MAILLAGE
C
C-----------------------------------------------------------------------
C  SIGNIFICATION DE IELM :
C
C  TYPE D'ELEMENT      NOMBRE DE POINTS          PROGRAMME ICI
C
C  11 : TRIANGLE P1            3                       OUI
C  21 : QUADRILATERE Q1        4                       OUI
C  41 : PRISMES TELEMAC-3D     6                       NON
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      SURFAC    |<-- |  SURFACE OU VOLUME DES ELEMENTS.
C |      XEL,YEL,..| -->|  COORDONNEES DES POINTS PAR ELEMENT.
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |      IELM      | -->|  TYPE D'ELEMENT
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : SURV11 , SURV21 , (SURV41) , PLANTE
C
C**********************************************************************
C
      USE BIEF, EX_SURVOL => SURVOL
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER, INTENT(IN) :: IELM,NELEM,NELMAX
C
      DOUBLE PRECISION, INTENT(INOUT) :: SURFAC(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,*)
      DOUBLE PRECISION, INTENT(IN) :: YEL(NELMAX,*)
      DOUBLE PRECISION, INTENT(IN) :: ZEL(NELMAX,*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      IF(IELM.EQ.11) THEN
C
        CALL SURV11(SURFAC, XEL,YEL,NELEM,NELMAX)
C
C     ELSEIF(IELM.EQ.21) THEN
C
C       CALL SURV21(SURFAC, XEL,YEL,NELEM,NELMAX)
C
C     ELSEIF(IELM.EQ.41) THEN
C
C       CALL SURV41(SURFAC, XEL,YEL,ZEL,NELEM,NELMAX)
C
C  IELM NON PREVU : ERREUR
C
      ELSE
         IF (LNG.EQ.1) WRITE(LU,100) IELM
         IF (LNG.EQ.2) WRITE(LU,101) IELM
100      FORMAT(1X,'SURVOL (BIEF) : IELM = ',1I6,' ELEMENT NON PREVU')
101      FORMAT(1X,
     *   'SURVOL (BIEF) : IELM = ',1I6,' ELEMENT NOT AVAILABLE')
         CALL PLANTE(1)
         STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
