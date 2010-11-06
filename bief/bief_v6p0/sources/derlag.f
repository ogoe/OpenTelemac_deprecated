C                       *****************
                        SUBROUTINE DERLAG
C                       *****************
C
     *( U , V , DT , X , Y , IKLE , IFABOR , LT , IELM , NDP , NPOIN ,
     *  NELEM , NELMAX , SURDET , XLAG , YLAG , DX , DY ,
     *  NSP , SHPLAG , DEBLAG , FINLAG , ELTLAG , NLAG , RESUX , RESUY ,
     *  NBOR , NELBOR , NULONE , NPTFR , MSK,MASKEL,MASKPT,T8)
C
C***********************************************************************
C BIEF VERSION 5.9           02/09/08       J-M JANIN (LNH) 30 87 72 84
C
C JMH 02/09/2008 : APPEl DE GTSH11 A LA PLACE DE GTSHP11
C
C***********************************************************************
C
C      FONCTION:
C
C   - FIXE, AU DEBUT DU CALCUL DE CHAQUE DERIVE, LES COORDONNEES BARY-
C     CENTRIQUES DANS LE MAILLAGE
C
C   - CALCUL, AUX PAS DE TEMPS SUIVANTS, LES POSITIONS SUCCESSIVES DE
C     LA DERIVE LAGRANGIENNE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    U,V         | -->| COMPOSANTE DE LA VITESSE                     |
C |    DT          | -->| PAS DE TEMPS.                                |
C |    X,Y         | -->| COORDONNEES DES POINTS DU MAILLAGE.          |
C |    IKLE        | -->| TRANSITION ENTRE LES NUMEROTATIONS LOCALE    |
C |                |    | ET GLOBALE.                                  |
C |    IFABOR      | -->| NUMEROS DES ELEMENTS AYANT UNE FACE COMMUNE  |
C |                |    | AVEC L'ELEMENT .  SI IFABOR<0 OU NUL         |
C |                |    | ON A UNE FACE LIQUIDE,SOLIDE,OU PERIODIQUE   |
C |    LT          | -->| NUMERO DU PAS DE TEMPS                       |
C |    IELM        | -->| TYPE DE MAILLAGE.                            |
C |    NDP         | -->| NOMBRE DE POINTS PAR ELEMENT                 |
C |    NPOIN       | -->| NOMBRE DE POINTS DU MAILLAGE.                |
C |    NELEM       | -->| NOMBRE D'ELEMENTS.                           |
C |    NELMAX      | -->| NOMBRE MAXIMAL D'ELEMENTS DANS LE MAILLAGE 2D|
C |    SURDET      | -->| VARIABLE UTILISEE PAR LA TRANSFORMEE ISOPARAM.
C |   XLAG,YLAG    |<-->| POSITIONS INSTANTANNEES DES DERIVES.         |
C |    DX,DY       | -- | STOCKAGE DES SOUS-PAS .
C |    NSP         | -- | NOMBRE DE SOUS PAS DE RUNGE KUTTA.
C |    SHPLAG      |<-->| COORDONNEES BARYCENTRIQUES INSTANTANNEES DES |
C |                |    | DERIVES DANS LEURS ELEMENTS RESPECTIFS.      |
C |    DEBLAG      | -->| NUMEROS DES PAS DE TEMPS DE DEBUT DE CALCUL  |
C |                |    | DES DERIVES.                                 |
C |    FINLAG      | -->| NUMEROS DES PAS DE TEMPS DE FIN DE CALCUL DES|
C |                |    | DERIVES.                                     |
C |    ELTLAG      |<-->| NUMEROS DES ELEMENTS DANS LESQUELS SE TROUVE |
C |                |    | A CET INSTANT CHACUNE DES DERVIES.           |
C |    NLAG        | -->| NOMBRE DE DERIVES.                           |
C |  RESUX,RESUY   |<-- | RESULTAT POUR ECRITURE SUR FICHIER DE LA     |
C |                |    | DERNIERE DERIVE ACHEVEE.                     |
C |    NBOR        | -->| NUMEROS GLOBAUX DES POINTS DE BORD.
C |    NELBOR      | -->| NUMERO DE L'ELEMENT ADJACENT AU K IEME       |
C |                |    | SEGMENT DE BORD.                             |
C |    NULONE      | -->| NUMERO LOCAL D'UN POINT DE BORD DANS         |
C |                |    | L'ELEMENT ADJACENT DONNE PAR NELBOR.         |
C |    NPTFR       | -->| NOMBRE DE POINTS FRONTIERES.                 |
C |    MSK         | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |    MASKEL      | -->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |    MASKPT      | -->|  TABLEAU DE MASQUAGE DES POINTS.
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      USE BIEF   !, EX_DERLAG => DERLAG
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(IN)    :: NPOIN,LT,IELM,NDP,NELEM,NLAG
      INTEGER         , INTENT(IN)    :: NPTFR,NELMAX
      DOUBLE PRECISION, INTENT(IN)    :: U(NPOIN),V(NPOIN),DT
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN),Y(NPOIN)
      INTEGER         , INTENT(IN)    :: IKLE(NELMAX,NDP)
      INTEGER         , INTENT(IN)    :: IFABOR(NELMAX,NDP)
      DOUBLE PRECISION, INTENT(IN)    :: SURDET(NELEM)
      DOUBLE PRECISION, INTENT(INOUT) :: XLAG(NPOIN,NLAG)
      DOUBLE PRECISION, INTENT(INOUT) :: YLAG(NPOIN,NLAG)
      INTEGER         , INTENT(INOUT) :: DEBLAG(NLAG),FINLAG(NLAG)
      INTEGER         , INTENT(INOUT) :: ELTLAG(NPOIN,NLAG)
      DOUBLE PRECISION, INTENT(INOUT) :: T8(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: DX(NPOIN),DY(NPOIN)
      INTEGER         , INTENT(INOUT) :: NSP(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: RESUX(NPOIN),RESUY(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: SHPLAG(NDP,NPOIN,NLAG)
      INTEGER         , INTENT(IN)    :: NBOR(NPTFR),NELBOR(NPTFR)
      INTEGER         , INTENT(IN)    :: NULONE(NPTFR)
      LOGICAL         , INTENT(IN)    :: MSK
      DOUBLE PRECISION, INTENT(IN)    :: MASKEL(NELMAX),MASKPT(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER ILAG,JLAG,LTT,NRK,IPOIN
C
      DOUBLE PRECISION Z(1),C
C
C-----------------------------------------------------------------------
C
      DO 10 ILAG=1,NLAG
C
        IF(LT.EQ.DEBLAG(ILAG)) THEN
C
C-----------------------------------------------------------------------
C
C   - FIXE, AU DEBUT DU CALCUL DE CHAQUE DERIVE, LES COORDONNEES BARY-
C     CENTRIQUES DANS LE MAILLAGE
C
C-----------------------------------------------------------------------
C
          CALL OV( 'X=CY    ' , XLAG(1,ILAG) , U , Z , -1.D0 , NPOIN )
          CALL OV( 'X=CY    ' , YLAG(1,ILAG) , V , Z , -1.D0 , NPOIN )
C
          IF(IELM.EQ.11) THEN
C
C  TRIANGLES P1
C  ============
C
C      REMPLISSAGE DES SHP ET DES ELT OPTIMISE
C
            CALL GTSH11(XLAG(1,ILAG),YLAG(1,ILAG),X,Y,SHPLAG(1,1,ILAG),
     *                  ELTLAG(1,ILAG),IKLE,NSP,NSP,NPOIN,
     *                  NELEM,NELMAX,1,MSK,MASKEL)
C
          ELSE
           IF(LNG.EQ.1) THEN 
             WRITE(LU,*) IELM,' : ELEMENT NON PREVU DANS DERLAG'
           ENDIF
           IF(LNG.EQ.2) THEN 
             WRITE(LU,*) IELM,': ELEMENT NOT IMPLEMENTED IN DERLAG'
           ENDIF
           STOP 
          ENDIF
C
          CALL OV( 'X=Y     ' , XLAG(1,ILAG) , X , Z , C , NPOIN )
          CALL OV( 'X=Y     ' , YLAG(1,ILAG) , Y , Z , C , NPOIN )
C
        ELSEIF(LT.GT.DEBLAG(ILAG).AND.LT.LE.FINLAG(ILAG)) THEN
C
C-----------------------------------------------------------------------
C
C   - CALCULE, AUX PAS DE TEMPS SUIVANTS, LES POSITIONS SUCCESSIVES DE
C     LA DERIVE LAGRANGIENNE
C
C-----------------------------------------------------------------------
C
C NOMBRE DE SOUS-PAS DE RUNGE-KUTTA PAR ELEMENT TRAVERSE
C ======================================================
C
          NRK     =  3
C
C  TRIANGLES P1
C  ============
C
          CALL CHAR11( U , V , DT , NRK , X , Y , IKLE , IFABOR ,
     *                 XLAG(1,ILAG) , YLAG(1,ILAG) , DX , DY ,
     *                 SHPLAG(1,1,ILAG) , ELTLAG(1,ILAG) , NSP ,
     *                 NPOIN , NPOIN , NELEM , NELMAX , SURDET , 1 ,T8)
C
        ENDIF
C
C-----------------------------------------------------------------------
C
C   - ANNULATION DES DERIVES SORTANT DU DOMAINE
C
C-----------------------------------------------------------------------
C
        IF(LT.EQ.FINLAG(ILAG)) THEN
          DO IPOIN=1,NPOIN
            IF(ELTLAG(IPOIN,ILAG).LE.0) THEN
              XLAG(IPOIN,ILAG) = X(IPOIN)
              YLAG(IPOIN,ILAG) = Y(IPOIN)
            ENDIF
          ENDDO
        ENDIF
C
10    CONTINUE
C
C-----------------------------------------------------------------------
C
C   - STOCKAGE POUR SORTIE DES RESULTATS DE LA DERNIERE DERIVE CALCULEE
C
C-----------------------------------------------------------------------
C
      CALL OV( 'X=C     ' , RESUX , Y , Z , 0.D0 , NPOIN )
      CALL OV( 'X=C     ' , RESUY , Y , Z , 0.D0 , NPOIN )
      LTT=0
      JLAG=1
      DO ILAG=1,NLAG
        IF(FINLAG(ILAG).GT.LTT.AND.FINLAG(ILAG).LE.LT) THEN
          LTT=FINLAG(ILAG)
          JLAG=ILAG
        ENDIF
      ENDDO
      IF(LTT.NE.0) THEN
        CALL OV( 'X=Y-Z   ' , RESUX , XLAG(1,JLAG) , X , C , NPOIN )
        CALL OV( 'X=Y-Z   ' , RESUY , YLAG(1,JLAG) , Y , C , NPOIN )
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
