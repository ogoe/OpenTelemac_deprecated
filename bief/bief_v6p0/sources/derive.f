C                       *****************
                        SUBROUTINE DERIVE
C                       *****************
C
     *( U , V , DT , X , Y , IKLE , IFABOR , LT , IELM , NDP , NPOIN ,
     *  NELEM , NELMAX , SURDET , XFLOT , YFLOT ,
     *  SHPFLO , DEBFLO , FINFLO , ELTFLO , NFLOT , NITFLO,FLOPRD,T8)
C
C***********************************************************************
C BIEF VERSION 5.1           18/08/94       J-M JANIN (LNH) 30 87 72 84
C
C***********************************************************************
C
C      FONCTION:
C
C   - CALCULE, AU MOMENT DU LARGAGE D'UN FLOTTEUR, LES COORDONNEES BARY-
C     CENTRIQUES DE CELUI-CI DANS LE MAILLAGE.
C
C   - CALCULE, AUX PAS DE TEMPS SUIVANTS, LES POSITIONS SUCCESSIVES DE
C     CE FLOTTEUR TRANSPORTE SANS FROTTEMENT PAR LE COURANT.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    U,V         | -->| COMPOSANTE DE LA VITESSE   
C |    DT          | -->| PAS DE TEMPS.   
C |    X,Y         | -->| COORDONNEES DES POINTS DU MAILLAGE. 
C |    IKLE        | -->| TRANSITION ENTRE LES NUMEROTATIONS LOCALE  
C |                |    | ET GLOBALE. 
C |    IFABOR      | -->| NUMEROS DES ELEMENTS AYANT UNE FACE COMMUNE 
C |                |    | AVEC L'ELEMENT .  SI IFABOR<0 OU NUL  
C |                |    | ON A UNE FACE LIQUIDE,SOLIDE,OU PERIODIQUE  
C |    LT          | -->| NUMERO DU PAS DE TEMPS  
C |    IELM        | -->| TYPE D'ELEMENT.   
C |    NDP         | -->| NOMBRE DE POINTS PAR ELEMENT   
C |    NPOIN       | -->| NOMBRE DE POINTS DU MAILLAGE.   
C |    NELEM       | -->| NOMBRE D'ELEMENTS.  
C |    NELMAX      | -->| NOMBRE MAXIMAL D'ELEMENTS DANS LE MAILLAGE 2D 
C |    SURDET      | -->| VARIABLE UTILISEE PAR LA TRANSFORMEE ISOPARAM.
C |  XFLOT,YFLOT   |<-->| POSITIONS SUCCESSIVES DES FLOTTEURS.  
C |    SHPFLO      |<-->| COORDONNEES BARYCENTRIQUES INSTANTANNEES DES 
C |                |    | FLOTTEURS DANS LEURS ELEMENTS RESPECTIFS.  
C |    DEBFLO      | -->| NUMEROS DES PAS DE TEMPS DE LARGAGE DE  
C |                |    | CHAQUE FLOTTEUR.  
C |    FINFLO      |<-->| NUMEROS DES PAS DE TEMPS DE FIN DE CALCUL DE 
C |                |    | DERIVE POUR CHAQUE FLOTTEUR.  
C |                |    | FORCE ICI SI UN FLOTTEUR SORT PAR UNE FR. LIQ.
C |    ELTFLO      |<-->| NUMEROS DES ELEMENTS DANS LESQUELS SE TROUVE 
C |                |    | A CET INSTANT CHACUN DES FLOTTEURS.        
C |    NFLOT       | -->| NOMBRE DE FLOTTEURS.                        
C |    NITFLO      | -->| NOMBRE MAXIMAL D'ENREGISTREMENTS DES        
C |                |    | POSITIONS SUCCESSIVES DES FLOTTEURS.       
C |    FLOPRD      | -->| NOMBRE DE PAS DE TEMPS ENTRE 2 ENREGITREMENTS
C |                |    | DES POSITIONS SUCCESSIVES DES FLOTTEURS.     
C |    T8          | -- | TABLEAU DE TRAVAIL  
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMME APPELE : DERIQ , DERIT
C
C***********************************************************************
C
      USE BIEF, EX_DERIVE => DERIVE
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(IN)    :: NPOIN,LT,IELM,NDP,NELEM
      INTEGER         , INTENT(IN)    :: NITFLO,FLOPRD,NELMAX,NFLOT
      DOUBLE PRECISION, INTENT(IN)    :: U(NPOIN),V(NPOIN),DT
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN),Y(NPOIN)
      INTEGER         , INTENT(IN)    :: IKLE(NELMAX,NDP)
      INTEGER         , INTENT(IN)    :: IFABOR(NELMAX,NDP)
      DOUBLE PRECISION, INTENT(IN)    :: SURDET(NELEM)
      DOUBLE PRECISION, INTENT(INOUT) :: XFLOT(NITFLO,NFLOT)
      DOUBLE PRECISION, INTENT(INOUT) :: YFLOT(NITFLO,NFLOT)
      INTEGER         , INTENT(INOUT) :: DEBFLO(NFLOT),FINFLO(NFLOT)
      INTEGER         , INTENT(INOUT) :: ELTFLO(NFLOT)
      DOUBLE PRECISION, INTENT(INOUT) :: SHPFLO(NDP,NFLOT)
      DOUBLE PRECISION, INTENT(INOUT) :: T8(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER LTP,IFLOT,N1,N2,N3,NRK,IELEM,LTT,NSP(1)
C
      DOUBLE PRECISION DET1,DET2,DET3,DX(1),DY(1)
C
C-----------------------------------------------------------------------
C
      LTT=(LT-1)/FLOPRD + 1
      LTP=(LT-2)/FLOPRD + 1
C
      DO 10 IFLOT=1,NFLOT
C
        IF(LT.EQ.DEBFLO(IFLOT)) THEN
C
C-----------------------------------------------------------------------
C
C   - CALCULE, AU MOMENT DU LARGAGE D'UN FLOTTEUR, LES COORDONNEES BARY-
C     CENTRIQUES DE CELUI-CI DANS LE MAILLAGE
C
C-----------------------------------------------------------------------
C
          XFLOT(LTT,IFLOT) = XFLOT(1,IFLOT)
          YFLOT(LTT,IFLOT) = YFLOT(1,IFLOT)
C
          IF(IELM.EQ.11) THEN
C
C MAILLAGE DE TRIANGLES P1
C ========================
C
            DO 20 IELEM=1,NELEM
              N1=IKLE(IELEM,1)
              N2=IKLE(IELEM,2)
              N3=IKLE(IELEM,3)
C
C DET1 = (N2N3,N2FLOT)  DET2 = (N3N1,N3FLOT)  DET3 = (N1N2,N1FLOT)
C ----------------------------------------------------------------
C
              DET1=(X(N3)-X(N2))*(YFLOT(LTT,IFLOT)-Y(N2))
     *            -(Y(N3)-Y(N2))*(XFLOT(LTT,IFLOT)-X(N2))
              DET2=(X(N1)-X(N3))*(YFLOT(LTT,IFLOT)-Y(N3))
     *            -(Y(N1)-Y(N3))*(XFLOT(LTT,IFLOT)-X(N3))
              DET3=(X(N2)-X(N1))*(YFLOT(LTT,IFLOT)-Y(N1))
     *            -(Y(N2)-Y(N1))*(XFLOT(LTT,IFLOT)-X(N1))
              IF(DET1.GE.0.D0.AND.DET2.GE.0.D0.AND.DET3.GE.0.D0) GOTO 30
C
20          CONTINUE
C
            IF(LNG.EQ.1) WRITE(LU,33) IFLOT
            IF(LNG.EQ.2) WRITE(LU,34) IFLOT
33          FORMAT(1X,'ERREUR D''INTERPOLATION DANS DERIVE :',/,
     *             1X,'LARGAGE DU FLOTTEUR',I6,/,
     *             1X,'EN DEHORS DU DOMAINE DE CALCUL')
34          FORMAT(1X,'INTERPOLATION ERROR IN DERIVE :',/,
     *             1X,'DROP POINT OF FLOAT',I6,/,
     *             1X,'OUT OF THE DOMAIN')
            STOP
C
C ELEMENT CONTENANT LE POINT DE LARGAGE TROUVE, CALCUL DES SHPFLO
C ---------------------------------------------------------------
C
30          CONTINUE
            SHPFLO(1,IFLOT) = DET1*SURDET(IELEM)
            SHPFLO(2,IFLOT) = DET2*SURDET(IELEM)
            SHPFLO(3,IFLOT) = DET3*SURDET(IELEM)
            ELTFLO (IFLOT)  = IELEM
C
          ELSE
            IF(LNG.EQ.1) WRITE(LU,123) IELM
            IF(LNG.EQ.2) WRITE(LU,124) IELM
123         FORMAT(1X,'DERIVE : TYPE D''ELEMENT NON PREVU : ',1I6)
124         FORMAT(1X,'DERIVE : UNEXPECTED TYPE OF ELEMENT: ',1I6)
            CALL PLANTE(1)
            STOP
          ENDIF
C
        ELSEIF(LT.GT.DEBFLO(IFLOT).AND.LT.LE.FINFLO(IFLOT)) THEN
C
C-----------------------------------------------------------------------
C
C   - CALCUL, AUX PAS DE TEMPS SUIVANTS, LES POSITIONS SUCCESSIVES DE
C     CE FLOTTEUR TRANSPORTE SANS FROTTEMENT PAR LE COURANT
C
C-----------------------------------------------------------------------
C
C NOMBRE DE SOUS-PAS DE RUNGE-KUTTA PAR ELEMENT TRAVERSE
C ======================================================
C
          NRK = 3
C
          XFLOT(LTT,IFLOT) = XFLOT(LTP,IFLOT)
          YFLOT(LTT,IFLOT) = YFLOT(LTP,IFLOT)
C
          IF(IELM.EQ.11) THEN
C
C  TRIANGLES P1
C  ============
C
            CALL CHAR11( U , V , DT , NRK , X , Y , IKLE , IFABOR ,
     *                   XFLOT(LTT,IFLOT) , YFLOT(LTT,IFLOT) , DX , DY ,
     *                   SHPFLO(1,IFLOT) , ELTFLO(IFLOT) , NSP ,
     *                   1 , NPOIN , NELEM , NELMAX , SURDET ,  1,T8 )
C
          ELSE
C
            IF(LNG.EQ.1) WRITE(LU,123) IELM
            IF(LNG.EQ.2) WRITE(LU,124) IELM
            STOP
C
          ENDIF
C
C  CAS DES FLOTTEURS PERDUS
C  ========================
C
          IF(ELTFLO(IFLOT).LE.0) FINFLO(IFLOT) = LT
C
        ENDIF
C
10    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
