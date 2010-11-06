C                       *****************
                        SUBROUTINE PROXIM
C                       *****************
C
     *(IP,XP,YP,X,Y,NP,NPOIN,IKLE,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 6.0      03/07/2009    J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C  FONCTION : RECHERCHE DES POINTS DU MAILLAGE LES PLUS PROCHES
C             POUR UN ENSEMBLE DE POINTS DONNES.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |    IP          |<-- |  ADRESSES DES POINTS TROUVES.
C |    XP,YP       | -->|  COORDONNEES DES POINTS DONNES.
C |    X,Y         | -->|  COORDONNEES DES POINTS DU MAILLAGE
C |    NP          | -->|  NOMBRE DE POINTS DONNES.
C |    NPOIN       | -->|  NOMBRE DE POINTS DU MAILLAGE.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : TELMAC
C
C***********************************************************************
C
      USE BIEF, EX_PROXIM => PROXIM
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NP,NPOIN,NELEM,NELMAX
      INTEGER, INTENT(INOUT) :: IP(NP)
      INTEGER, INTENT(IN)    :: IKLE(NELMAX,3)
C      
      DOUBLE PRECISION, INTENT(IN) :: XP(NP),YP(NP),X(NPOIN),Y(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,K,IELEM
      DOUBLE PRECISION X1,Y1,X2,Y2,X3,Y3,A31,A12,A23,DIST2,D2,ALERT
      DOUBLE PRECISION XX,YY
C
      INTRINSIC SQRT
C
      DOUBLE PRECISION P_DSUM,P_DMAX
      EXTERNAL         P_DSUM,P_DMAX
C
C-----------------------------------------------------------------------
C
      DO 10 K=1,NP
        IP(K)=0
        DIST2=1.D10
        ALERT=0.D0
        XX=-1.D10
        YY=-1.D10
C
C       BOUCLE SUR LES TRIANGLES :
C
        DO 20 IELEM=1,NELEM
          X1=X(IKLE(IELEM,1))
          X2=X(IKLE(IELEM,2))
          X3=X(IKLE(IELEM,3))
          Y1=Y(IKLE(IELEM,1))
          Y2=Y(IKLE(IELEM,2))
          Y3=Y(IKLE(IELEM,3))
          A31=XP(K)*Y3-YP(K)*X3+X3*Y1-X1*Y3+X1*YP(K)-XP(K)*Y1
          A12=XP(K)*Y1-YP(K)*X1+X1*Y2-X2*Y1+X2*YP(K)-XP(K)*Y2
          A23=XP(K)*Y2-YP(K)*X2+X2*Y3-X3*Y2+X3*YP(K)-XP(K)*Y3
          IF(A31.GT.-1.D-6.AND.A12.GT.-1.D-6.AND.A23.GT.-1.D-6) THEN
C           ON PREND LE POINT LE PLUS PROCHE
            DO I=1,3
              D2=(XP(K)-X(IKLE(IELEM,I)))**2+(YP(K)-Y(IKLE(IELEM,I)))**2
              IF(D2.LT.DIST2) THEN
                IP(K)=IKLE(IELEM,I)
                DIST2=D2
              ENDIF
            ENDDO
          ENDIF
20      CONTINUE
        IF(IP(K).EQ.0) THEN
          IF(LNG.EQ.1) WRITE(LU,*) 'POINT SOURCE ',K,' HORS DOMAINE'
          IF(LNG.EQ.2) WRITE(LU,*) 'SOURCE POINT ',K,' OUTSIDE DOMAIN'
          IF(NCSIZE.LE.1) THEN
            WRITE(LU,*) ' '
            IF(LNG.EQ.1) WRITE(LU,*) 'INTERDIT EN MODE SCALAIRE'
            IF(LNG.EQ.2) WRITE(LU,*) 'NOT ALLOWED IN SCALAR MODE'
            CALL PLANTE(1)
            STOP
          ENDIF
        ELSE
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'POINT SOURCE ',K,' ASSIMILE AU POINT ',IP(K)
            WRITE(LU,*) 'SITUE A ',SQRT(DIST2),' METRES'
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'SOURCE POINT ',K,' PUT ON POINT ',IP(K)
            WRITE(LU,*) 'LOCATED AT ',SQRT(DIST2),' METRES'
          ENDIF
          IF(SQRT(DIST2).GT.1.D-6.AND.NCSIZE.GT.1) THEN
            XX=X(IP(K))
            YY=Y(IP(K))
            ALERT=1.D0
          ENDIF
        ENDIF
        IF(NCSIZE.GT.1) THEN
          XX=P_DMAX(XX)
          YY=P_DMAX(YY)
          ALERT=P_DSUM(ALERT)
        ENDIF
        IF(ALERT.GT.0.5D0) THEN
          WRITE(LU,*) ' '
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'EN MODE PARALLELE LES SOURCES PONCTUELLES'
            WRITE(LU,*) 'DOIVENT COINCIDER EXACTEMENT AVEC DES POINTS'
            WRITE(LU,*) 'DU MAILLAGE, POUR LA SOURCE ',K,' CHOISIR :'
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'IN PARALLEL SOURCES MUST COINCIDE WITH'
            WRITE(LU,*) 'NODES IN THE MESH, FOR SOURCE',K,' CHOOSE:'
          ENDIF
          WRITE(LU,*) 'X=',XX,' Y=',YY
          CALL PLANTE(1)
          STOP
        ENDIF
10    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
