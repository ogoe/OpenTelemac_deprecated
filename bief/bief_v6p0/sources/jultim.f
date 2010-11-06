C             ********************************
              DOUBLE PRECISION FUNCTION JULTIM
C             ********************************
C
     *(YEAR,MONTH,DAY,HOUR,MIN,SEC,AT)
C
C***********************************************************************
C BIEF VERSION 5.1          12/07/95    E. DAVID (LHF) 76 33 42 08
C
C***********************************************************************
C
C     FONCTIONS :  CALCUL DU TEMPS ECOULE DEPUIS LE 31/12/1899
C                  EXPRIME EN SIECLES JULIEN
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
C PROGRAMMES APPELES :
C
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER,          INTENT(IN)    :: MONTH,DAY,HOUR,MIN,SEC
      INTEGER,          INTENT(INOUT) :: YEAR
      DOUBLE PRECISION, INTENT(IN)    :: AT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C      
      INTEGER GREG,Y,M
      DOUBLE PRECISION J
C
      INTRINSIC INT
C
      PARAMETER (GREG=15+31*(10+12*1582))
C
C-----------------------------------------------------------------------
C
      IF(YEAR.EQ.0) THEN
        IF (LNG.EQ.1) WRITE (LU,100)
        IF(LNG.EQ.2)  WRITE (LU,101)
        STOP
      ENDIF
      IF(YEAR.LT.0) YEAR=YEAR+1
C
      IF (MONTH.GT.2) THEN
       Y=YEAR
       M=MONTH+1
      ELSE
       Y=YEAR-1
       M=MONTH+13
      ENDIF
C
      J=INT(365.25D0*Y)+INT(30.6001D0*M)+DAY+1720995.D0
      IF(DAY+31*(MONTH+12*YEAR).GE.GREG) THEN
        J=J+2-INT(0.01D0*Y)+INT(0.25D0*INT(0.01D0*Y))
      ENDIF
      J=J-2415020.5D0
      JULTIM=(J+(HOUR+(MIN+(SEC+AT)/60.D0)/60.D0)/24.D0)/36525.D0
C
C---------------------------------------------------------------
C
100   FORMAT (//,10X,'**********************************',
     *         /,10X,'       FONCTION JULTIM',
     *         /,10X,' LA VALEUR DE L''ANNEE EST NULLE',
     *         /,10X,' CALCUL IMPOSSIBLE ...'
     *         /,10X,'**********************************')
101   FORMAT (//,10X,'**********************************',
     *         /,10X,'       JULTIM FUNCTION',
     *         /,10X,' THE VALUE FOR THE YEAR IS ZERO',
     *         /,10X,' COMPUTATION NOT POSSIBLE ...'
     *         /,10X,'**********************************')
C
C---------------------------------------------------------------
C
      RETURN
      END
 
