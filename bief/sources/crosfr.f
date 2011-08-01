C                       *****************
                        SUBROUTINE CROSFR
C                       *****************
C
     *(X,Y,XR,YR,XMAIL,YMAIL,NPMAX,NBOR,KP1BOR,NPTFR,DM,OK)
C
C***********************************************************************
C BIEF VERSION 5.9        20/03/08    J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C  FONCTION : ON SOUHAITE INTERPOLER LES FONDS POUR UN POINT DE
C             COORDONNEES X ET Y. UN POINT QUI VA SERVIR A CETTE
C             INTERPOLATION A POUR COORDONNEES XR ET YR. ON VERIFIE
C             ICI QUE CE POINT N'EST PAS EN DEHORS DU DOMAINE. POUR
C             CE FAIRE, ON VERIFIE QUE LE SEGMENT JOIGNANT (X,Y) ET
C             (XR,YR) NE COUPE PAS LA FRONTIERE DU DOMAINE.
C
C  NOTE JMH : NE MARCHE PAS EN PARALLELE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |    X,Y         | -->|  COORDONNEES DU POINT OU L'ON VEUT INTERPOLER
C |    XR,YR       | -->|  COORDONNEES DU POINT SITUE DANS LE CADRAN R
C |    XMAIL,YMAIL | -->|  COORDONNEES DES POINTS DU MAILLAGE
C |    NPMAX       | -->|  NOMBRE MAX DE POINTS DU MAILLAGE
C |    NBOR        | -->|  NUMEROTATION DES ELEMENTS DE BORD
C |    NPTRF       | -->|  NOMBRE DE POINTS FRONTIERE
C |    DM          | -->|  DISTANCE MINIMALE A LA FRONTIERE
C |    OK          | <--|  .TRUE. SI ON TROUVE UN POINT DANS LE CADRAN R
C |                |    |  .FALSE: SINON
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : FASP
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(IN) :: X,Y,XR,YR,DM
      INTEGER, INTENT(IN)          :: NPTFR,NPMAX
      DOUBLE PRECISION, INTENT(IN) :: XMAIL(NPMAX),YMAIL(NPMAX)
      INTEGER, INTENT(IN)          :: NBOR(NPTFR),KP1BOR(NPTFR)
      LOGICAL, INTENT(INOUT)       :: OK
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER KA
C
      DOUBLE PRECISION DM2,XA,YA,XB,YB,DET,ALFA,BETA,EPS,DISTA2,DISTB2
C
C-----------------------------------------------------------------------
C
C     ON REFUSE DE PRENDRE DES POINTS TROP PRES DU BORD
C     DM     : DISTANCE MINIMUM
      DM2 = DM**2
C
      DO 10 KA=1,NPTFR
C
C INTERSECTION D'UN SEGMENT DE LA FRONTIERE ET DU SEGMENT
C FORME PAR LES POINTS (X,Y) ET (XR,YR)
C
        XA = XMAIL(NBOR(KA))
        YA = YMAIL(NBOR(KA))
        XB = XMAIL(NBOR(KP1BOR(KA)))
        YB = YMAIL(NBOR(KP1BOR(KA)))
C
        DET = (XR-X)*(YA-YB) - (YR-Y)*(XA-XB)
C
        IF(ABS(DET).LT.1.D-6) GO TO 10
C
        ALFA = ( (XA-X)*(YA-YB) - (YA-Y)*(XA-XB) ) / DET
        BETA = ( (XR-X)*(YA-Y ) - (YR-Y)*(XA-X ) ) / DET
C
        EPS=0.05D0
        IF(ALFA.GE.EPS.AND.ALFA.LE.1.D0-EPS.AND.
     *     BETA.GE.EPS.AND.BETA.LE.1.D0-EPS) THEN
          OK = .FALSE.
          GO TO 1000
        ENDIF
C
C ON ELIMINE AUSSI LES POINTS TROP PROCHES DE LA FRONTIERE
C
        DISTA2 = (XR-XA)**2 + (YR-YA)**2
        DISTB2 = (XR-XB)**2 + (YR-YB)**2
        IF(DISTA2.LT.DM2.OR.DISTB2.LT.DM2) THEN
          OK = .FALSE.
          GO TO 1000
        ENDIF
C
10    CONTINUE
C
C-----------------------------------------------------------------------
C
1000  RETURN
      END
