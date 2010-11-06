C                       ****************
                        SUBROUTINE ASTRO
C                       ****************
C
     *(YEAR,MONTH,DAY,HOUR,MIN,SEC,AT,ARL,ARS,DL,DS,AL,AS)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.2    01/03/94      E. DAVID    (LHF) 76 33 42 08
C                                         F LEPEINTRE (LNH) 30 87 78 54
C                                         J-M JANIN   (LNH) 30 87 72 84
C***********************************************************************
C
C      FONCTION:
C      =========
C
C  CALCULE LES TERMES ASTRONOMIQUES NECESSAIRES POUR LE CALCUL
C  DES TERMES DE FORCAGE DE LA MAREE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C | YEAR,MONTH,DAY | -->| DATE DU CALCUL DES TERMES ASTROS             |
C |  HOUR,MIN,SEC  | -->| TEMPS DU CALCUL EN TEMPS UNIVERSEL (TU) ||   |
C |  AT            | -->| TEMPS
C |  L,S           |    | DESIGNE RESPECTIVEMENT LA LUNE ET LE SOLEIL  |
C |  AR L,S        |<-- | RAPPORT RAYON TERRE / DISTANCE A L'ASTRE     |
C |  D  L,S        |<-- | DECLINAISON DE L'ASTRE                       |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE PAR : MAREE
C SOUS-PROGRAMMES APPELES   : JULTIM, DMO
C
C***********************************************************************
C
C
C DESCRIPTION DES PARAMETRES ASTRONOMIQUES :
C
C  NOTE : "PERIODE" EST LA PERIODE DU PARAMETRE EXPRIME EN JOUR/ANNEE/
C          SIECLE JULIENS
C
C .____.__________________________________________.________.___________.
C |CODE|                   ROLE                   |NOTATION| PERIODE   |
C |____|__________________________________________|________|___________|
C |    |             FONCTION DU TEMPS            |        |           |
C | H  |   LONGITUDE SOLAIRE MOYENNE              |   H,L  | 365,25  J |
C | S  |   LONGITUDE LUNAIRE MOYENNE              |    S   | 27,32   J |
C | P  |   LONGITUDE DU PERIGEE LUNAIRE MOYENNE   |    P   | 8,85    A |
C | O  |   LONGITUDE DU NOEUD LUNAIRE ASCENDANT   | N,OMEGA| 18,61   A |
C | OM |   INCLINAISON DE L'EQUATEUR SUR          |  OMEGA | 27665,7 S |
C |    |     L'ECLIPTIQUE                         |        |           |
C | L  |   LONGITUDE VRAIE DE LA LUNE             |    L   |           |
C | CR |   PARALLAXE VRAIE DE LA LUNE             |   C/R  |           |
C | DL |   DECLINAISON DE LA LUNE                 |  DELTA |           |
C | AL |   ASCENSION DROITE DE LA LUNE            |  ALPHA |           |
C | NU |   ASCENSION DROITE DE G (EQUATEUR-ORBITE)|    NU  |           |
C | ET |   EXCENTRICITE DE L'ORBITE TERRESTRE     |    E   |           |
C | MA |   ANOMALIE MOYENNE DU SOLEIL             |    M   | 365,26  J |
C | EA |   ANOMALIE EXCENTRIQUE DU SOLEIL         |    E   |           |
C | TS |   DISTANCE TERRE-SOLEIL (UA)             |    R   |           |
C | TL |   DISTANCE TERRE-LUNE   (KM)             |    R   |           |
C | ARL|   RAPPORT RT / TL                        |   A/R  |           |
C | ARS|   RAPPORT RT / TS                        |   A/R  |           |
C | VS |   ANOMALIE VRAIE DU SOLEIL               |    V   |           |
C | LS |   LONGITUDE VRAIE DU SOLEIL              |  TETA  |           |
C | AS |   ASCENSION DROITE DU SOLEIL             |  ALPHA |           |
C | DS |   DECLINAISON DU SOLEIL                  |  DELTA |           |
C |    |                                          |        |           |
C |    |               CONSTANTES                 |        |           |
C | I0 |   INCLINAISON DE L'ORBITE LUNAIRE        |   I    |           |
C |    |     SUR L'ECLIPTIQUE                     |        |           |
C | E  |   EXENTRICITE DE L'ORBITE LUNAIRE        |   E    |           |
C | M  |   RAPPORT DU MOUVEMENT MOYEN DU SOLEIL   |   M    |           |
C |    |     AU MOUVEMENT MOYEN DE LA LUNE        |        |           |
C | UA |   UNITE ASTRONOMIQUE EN KM               |   UA   |           |
C | RT |   RAYON MOYEN DE LA TERRE                |   A    |           |
C | C  |   DEMI GRAND AXE DE L'ORBITE LUNAIRE (KM)|  C,C   |           |
C | AC |   RAPPORT RT / C                         |  A/C   |           |
C |____|__________________________________________|________|___________|
C
C AUTRE CONSTANTES :
C
C   PI  :  VALEUR DE PI
C   API :  TRANSFORME DEGRES EN RADIANS
C
C VARIABLES INTERMEDIAIRES
C
C   KSI,I,X,EA1 : VARIABLES INTERMEDIAIRES UTILISES POUR LE CALCUL DE
C                 L'ASCENSION.
C                 DROITE DE LA LUNE (A) ET DE LA DECLINAISON LUNAIRE (D)
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: YEAR,MONTH,DAY,HOUR,MIN,SEC
      DOUBLE PRECISION, INTENT(IN)    :: AT
      DOUBLE PRECISION, INTENT(INOUT) :: ARL,ARS,DL,DS,AL,AS
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      DOUBLE PRECISION T,H,S,P,O,OM,L,CR,ET,MA,EA,TS,VS,LS,JULTIM
      DOUBLE PRECISION I0,E,M,UA,RT,C,AC,PI,API,EA1,KSI,I,X,NU
C
      INTRINSIC ACOS,ASIN,ATAN,COS,SIN,SQRT,TAN,ABS,MOD
C
      DOUBLE PRECISION DMO,ATANC
      EXTERNAL         DMO,ATANC
C
C-----------------------------------------------------------------------
C
      PI  = ACOS (-1.D0)
      API = PI / 180.D0
      I0  = 5.145576994D0 * API
      E   = 0.05490D0
      M   = 0.074804D0
      UA  = 149503899.D0
      RT  = 6378.D0
      C   = 384403.D0
      AC  = RT / C
C
      T   = JULTIM(YEAR,MONTH,DAY,HOUR,MIN,SEC,AT)
C
      H   = DMO ( 279.69668D0 + 36000.76892D0       * T
     *                        + 0.0003025D0         * T * T )
      S   = DMO ( 270.434164D0 + 481267.8831D0      * T
     *                         - 0.001133D0         * T * T
     *                         + 0.0000019D0        * T * T * T )
      P   = DMO ( 334.328019444D0 + 4069.03220556D0 * T
     *                             - 0.01034D0      * T * T
     *                             + 0.0000125D0    * T * T * T )
      O   = DMO ( 259.183275D0 - 1934.1420D0        * T
     *                      + 0.002078D0            * T * T
     *                      + 0.0000022D0           * T * T * T )
      OM  = DMO ( 23.452294D0 - 0.0130125D0         * T
     *                        - 0.00000164D0        * T * T
     *                        + 0.000000503D0       * T * T * T )
      L   = S + 2*E*SIN(S-P) +  5.D0/4.D0 *E*E*SIN(2*(S-P)) +
     *                         15.D0/4.D0 *M*E*SIN(S-2*H+P) +
     *                         11.D0/8.D0 *M*M*SIN(2*(S-H))
      CR  = 1.D0 + E*COS(S-P) +            E*E*COS(2*(S-P)) +
     *                         15.D0/8.D0 *M*E*COS(S-2*H+P) +
     *                                     M*M*COS(2*(S-H))
      KSI=MOD(O-ATAN(SIN(O)/(SIN(I0)/TAN(OM)+COS(I0)*COS(O))),PI)
C
C KSI VARIE DE -12 A +12 DEGRES EN 18.7 ANS
C
      IF (KSI.GT.+API*13.D0) KSI=KSI-PI
      IF (KSI.LT.-API*13.D0) KSI=KSI+PI
C
C CALCUL DE I
C
      I   = ACOS( COS(OM)*COS(I0) - SIN(OM)*SIN(I0)*COS(O) )
C
C CALCUL DE X : TAN(X) = TAN(L-KSI) * COS(I)
C     AVEC  X ENTRE 0 ET 2 PI
C
      X   = ATANC( TAN(L-KSI) * COS(I) , L )
C
      ET  = 0.01675104D0 - 0.0000418D0 * T - 0.000000126D0 * T * T
      MA  = DMO ( 358.47583D0 + 35999.04975D0 * T
     *                        - 0.00015D0     * T * T
     *                        + 0.0000033D0   * T * T * T )
      EA1 = MA
10    CONTINUE
      EA  = MA + ET * SIN (EA1)
      IF (ABS(EA-EA1).GT.1.D-12) THEN
        EA1=EA
        GOTO 10
      ENDIF
      TS  = 1.0000002D0 * ( 1.D0-ET*COS(EA) ) * UA
      VS  = 2.D0 * ATAN( SQRT((1.D0+ET)/(1.D0-ET)) * TAN(EA/2.D0) )
      LS  = H + VS - MA
C
C PARAMETRES DE SORTIE
C
      ARL = AC * CR
      ARS = RT / TS
      DL  = ASIN ( SIN(L-KSI) * SIN(I) )
      DS  = ASIN ( SIN(OM)*SIN(LS) )
      NU  = ATANC(SIN(O)/(SIN(OM)/TAN(I0)+COS(OM)*COS(O)),O)
      AL  = X + NU
      AS  = ATAN ( (COS(OM)*SIN(LS) / COS(LS)) )
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
