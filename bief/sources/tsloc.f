C                       *******************************
                        DOUBLE PRECISION FUNCTION TSLOC
C                       *******************************
C
     * (YEAR,MONTH,DAY,HOUR,MINUTE,SEC,AT)
C
C***********************************************************************
C BIEF VERSION 5.9          02/06/08    E. DAVID (LHF) 76 33 42 08
C
C***********************************************************************
C
C     FONCTIONS :  CALCUL DU TEMPS SIDERAL LOCAL EN RADIAN POUR
C                  LA DATE DONNEE COMPTEE EN TEMPS UNIVERSEL
C
C     ATTENTION |
C
C     DANS CETTE FONCTION, LE TEMPS DU CALCUL DOIT ETRE SEPARE EN :
C
C     - TEMPS JUSQU'A 0H DE LA MEME JOURNEE,
C     - TEMPS RESTANT JUSQU'A LA DATE PRECISE EN SECONDES.
C
C     CECI NECESSITE UNE RECOMBINAISON ENTRE LA DATE DE DEPART
C     ET LE TEMPS DE CALCUL AT :
C
C     CALCUL DU TEMPS EN SIECLE JULIEN A 0H DE LA MEME DATE
C
C     AT1 : NOMBRE DE JOURS EN SECONDES A AJOUTER A LA DATE
C           DE REFERENCE (YEAR, MONTH, DAY) POUR ARRIVER A
C           LA DATE DU CALCUL A L'HEURE (0H, 0MIN, 0SEC).
C     ATR : NOMBRE DE SECONDES ENTRE LA DATE REPRESENTEE PAR AT1
C           ET LA DATE EXACTE DU CALCUL. EN FAIT, ATR CORRESPOND
C           A L'HEURE TU EXPRIME EN SECONDES.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C º      NOM       ºMODEº                   ROLE
C º________________º____º______________________________________________
C º   YEAR         º -->º ANNEE.
C º   MONTH        º -->º MOIS.
C º   DAY          º -->º JOUR.
C º   HOUR         º -->º HEURE EN TEMPS UNIVERSEL.
C º   MIN          º -->º MINUTE EN TEMPS UNIVERSEL.
C º   SEC          º -->º SECONDE EN TEMPS UNIVERSEL.
C º   AT           º -->º TEMPS
C º________________º____º______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
      USE BIEF, EX_TSLOC => TSLOC
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER, INTENT(IN)          :: MONTH,DAY,HOUR,MINUTE,SEC
      INTEGER, INTENT(INOUT)       :: YEAR
      DOUBLE PRECISION, INTENT(IN) :: AT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      DOUBLE PRECISION T,TETA,TETA0,UT
      DOUBLE PRECISION AT1,ATR
C
      INTRINSIC ACOS,INT
C
C-----------------------------------------------------------------------
C
      ATR = AT + ( HOUR * 60.D0 + MINUTE ) * 60.D0 + SEC
      AT1 = INT ( ATR / ( 24.D0 * 3600.D0 ) ) * ( 24.D0 * 3600.D0 )
      ATR = ATR - AT1
      T = JULTIM(YEAR,MONTH,DAY,0,0,0,AT1)
C
C CALCUL DU TEMPS SIDERAL A GREENWICH A 0 H (EN HEURES)
C
      TETA0 = 6.6460656D0 + 2400.051262D0 * T + 0.00002581D0 * T**2
C
C CALCUL DU TEMPS SIDERAL A GREENWICH A L'HEURE TU
C
      UT = ATR / 3600.D0
      TETA = TETA0 + UT*1.002737908D0
C
C CALCUL DU TEMPS SIDERAL LOCAL EN RADIANS SANS TENIR COMPTE DE LA
C LONGITUDE.
C
      TSLOC = TETA * ACOS(-1.D0) / 12.D0
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
