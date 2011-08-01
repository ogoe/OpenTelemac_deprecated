C                       ***************************
                        SUBROUTINE NOMVAR_TELEMAC2D
C                       ***************************
C
     *(TEXTE,TEXTPR,MNEMO,NPERIAF,NTRAC,NAMETRAC)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.8   31/08/07   J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C FONCTION  :  FIXE LES NOMS DES VARIABLES DU CODE POUR LES FICHIERS
C              DE RESULTAT ET DE GEOMETRIE (TEXTE) ET POUR LE FICHIER
C              DE RESULTATS DU CALCUL PRECEDENT (TEXTPR)
C
C              EN GENERAL TEXTE ET TEXTPR SONT EGAUX SAUF SI ON FAIT
C              UNE SUITE A PARTIR D'UN AUTRE LOGICIEL.
C
C-----------------------------------------------------------------------
C
C FUNCTION  :  GIVE THE NAMES OF VARIABLES AS UNDERSTOOD IN
C              SELAFIN FORMAT
C
C              FOR RESULTS FILE AND GEOMETRY FILE: IN TEXTE
C
C              FOR PREVIOUS COMPUTATION FILE: IN TEXTPR
C
C              TEXTE AND TEXTPR ARE GENERALLY THE SAME, EXCEPT TO
C              RECOVER FILES IN ANOTHER LANGUAGE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________|
C |   TEXTE        |<-- | SEE ABOVE
C |   TEXTPR       |<-- | SEE ABOVE
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PREDON
C
C SOUS-PROGAMME APPELE : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER(LEN=32), INTENT(INOUT) :: TEXTE(*),TEXTPR(*)
      CHARACTER(LEN=8),  INTENT(INOUT) :: MNEMO(*)
      INTEGER, INTENT(IN)              :: NPERIAF,NTRAC
      CHARACTER(LEN=32), INTENT(IN)    :: NAMETRAC(32)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER(LEN=2) I_IN_2_LETTERS(32)
      DATA I_IN_2_LETTERS /'1 ','2 ','3 ','4 ','5 ','6 ','7 ','8 ','9 ',
     *                     '10','11','12','13','14','15','16','17','18',
     *                     '19','20','21','22','23','24','25','26','27',
     *                     '28','29','30','31','32'/
      INTEGER I
C
C-----------------------------------------------------------------------
C
C  ENGLISH
C
      IF(LNG.EQ.2) THEN
C
      TEXTE (1 ) = 'VELOCITY U      M/S             '
      TEXTE (2 ) = 'VELOCITY V      M/S             '
      TEXTE (3 ) = 'CELERITY        M/S             '
      TEXTE (4 ) = 'WATER DEPTH     M               '
      TEXTE (5 ) = 'FREE SURFACE    M               '
      TEXTE (6 ) = 'BOTTOM          M               '
      TEXTE (7 ) = 'FROUDE NUMBER                   '
      TEXTE (8 ) = 'SCALAR FLOWRATE M2/S            '
      TEXTE (9 ) = 'EX TRACER                       '
      TEXTE (10) = 'TURBULENT ENERG.JOULE/KG        '
      TEXTE (11) = 'DISSIPATION     WATT/KG         '
      TEXTE (12) = 'VISCOSITY       M2/S            '
      TEXTE (13) = 'FLOWRATE ALONG XM2/S            '
      TEXTE (14) = 'FLOWRATE ALONG YM2/S            '
      TEXTE (15) = 'SCALAR VELOCITY M/S             '
      TEXTE (16) = 'WIND ALONG X    M/S             '
      TEXTE (17) = 'WIND ALONG Y    M/S             '
      TEXTE (18) = 'AIR PRESSURE    PASCAL          '
      TEXTE (19) = 'BOTTOM FRICTION                 '
      TEXTE (20) = 'DRIFT ALONG X   M               '
      TEXTE (21) = 'DRIFT ALONG Y   M               '
      TEXTE (22) = 'COURANT NUMBER                  '
      TEXTE (23) = 'VARIABLE 23     UNIT   ??       '
      TEXTE (24) = 'VARIABLE 24     UNIT   ??       '
      TEXTE (25) = 'VARIABLE 25     UNIT   ??       '
      TEXTE (26) = 'VARIABLE 26     UNIT   ??       '
      TEXTE (27) = 'HIGH WATER MARK M               '
      TEXTE (28) = 'HIGH WATER TIME S               '
      TEXTE (29) = 'HIGHEST VELOCITYM/S             '
      TEXTE (30) = 'TIME OF HIGH VELS               '
      TEXTE (31) = 'FRICTION VEL.   M/S             '
C
C TEXTPR IS USED FOR READING PREVIOUS COMPUTATION FILES.
C IN GENERAL TEXTPR=TEXTE BUT YOU CAN FOLLOW UP A COMPUTATION
C FROM ANOTHER CODE WITH DIFFERENT NAMES THAT YOU HAVE TO
C WRITE HERE.
C
      TEXTPR (1 ) = 'VELOCITY U      M/S             '
      TEXTPR (2 ) = 'VELOCITY V      M/S             '
      TEXTPR (3 ) = 'CELERITY        M/S             '
      TEXTPR (4 ) = 'WATER DEPTH     M               '
      TEXTPR (5 ) = 'FREE SURFACE    M               '
      TEXTPR (6 ) = 'BOTTOM          M               '
      TEXTPR (7 ) = 'FROUDE NUMBER                   '
      TEXTPR (8 ) = 'SCALAR FLOWRATE M2/S            '
      TEXTPR (9 ) = 'EX TRACER                       '
      TEXTPR (10) = 'TURBULENT ENERG.JOULE/KG        '
      TEXTPR (11) = 'DISSIPATION     WATT/KG         '
      TEXTPR (12) = 'VISCOSITY       M2/S            '
      TEXTPR (13) = 'FLOWRATE ALONG XM2/S            '
      TEXTPR (14) = 'FLOWRATE ALONG YM2/S            '
      TEXTPR (15) = 'SCALAR VELOCITY M/S             '
      TEXTPR (16) = 'WIND ALONG X    M/S             '
      TEXTPR (17) = 'WIND ALONG Y    M/S             '
      TEXTPR (18) = 'AIR PRESSURE    PASCAL          '
      TEXTPR (19) = 'BOTTOM FRICTION                 '
      TEXTPR (20) = 'DRIFT ALONG X   M               '
      TEXTPR (21) = 'DRIFT ALONG Y   M               '
      TEXTPR (22) = 'COURANT NUMBER                  '
      TEXTPR (23) = 'VARIABLE 23     UNIT   ??       '
      TEXTPR (24) = 'VARIABLE 24     UNIT   ??       '
      TEXTPR (25) = 'VARIABLE 25     UNIT   ??       '
      TEXTPR (26) = 'VARIABLE 26     UNIT   ??       '
      TEXTPR (27) = 'HIGH WATER MARK M               '
      TEXTPR (28) = 'HIGH WATER TIME S               '
      TEXTPR (29) = 'HIGHEST VELOCITYM/S             '
      TEXTPR (30) = 'TIME OF HIGH VELS               '
      TEXTPR (31) = 'FRICTION VEL.   M/S             '
C
C-----------------------------------------------------------------------
C
C  FRANCAIS OU AUTRE
C
      ELSE
C
      TEXTE (1 ) = 'VITESSE U       M/S             '
      TEXTE (2 ) = 'VITESSE V       M/S             '
      TEXTE (3 ) = 'CELERITE        M/S             '
      TEXTE (4 ) = 'HAUTEUR D''EAU   M               '
      TEXTE (5 ) = 'SURFACE LIBRE   M               '
      TEXTE (6 ) = 'FOND            M               '
      TEXTE (7 ) = 'FROUDE                          '
      TEXTE (8 ) = 'DEBIT SCALAIRE  M2/S            '
      TEXTE (9 ) = 'EX TRACEUR                      '
      TEXTE (10) = 'ENERGIE TURBUL. JOULE/KG        '
      TEXTE (11) = 'DISSIPATION     WATT/KG         '
      TEXTE (12) = 'VISCOSITE TURB. M2/S            '
      TEXTE (13) = 'DEBIT SUIVANT X M2/S            '
      TEXTE (14) = 'DEBIT SUIVANT Y M2/S            '
      TEXTE (15) = 'VITESSE SCALAIREM/S             '
      TEXTE (16) = 'VENT X          M/S             '
      TEXTE (17) = 'VENT Y          M/S             '
      TEXTE (18) = 'PRESSION ATMOS. PASCAL          '
      TEXTE (19) = 'FROTTEMENT                      '
      TEXTE (20) = 'DERIVE EN X     M               '
      TEXTE (21) = 'DERIVE EN Y     M               '
      TEXTE (22) = 'NBRE DE COURANT                 '
      TEXTE (23) = 'VARIABLE 23     UNITES ??       '
      TEXTE (24) = 'VARIABLE 24     UNITES ??       '
      TEXTE (25) = 'VARIABLE 25     UNITES ??       '
      TEXTE (26) = 'VARIABLE 26     UNITES ??       '
      TEXTE (27) = 'COTE MAXIMUM    M               '
      TEXTE (28) = 'TEMPS COTE MAXI S               '
      TEXTE (29) = 'VITESSE MAXIMUM M/S             '
      TEXTE (30) = 'T VITESSE MAXI  S               '
      TEXTE (31) = 'VITESSE DE FROT.M/S             '
C
C TEXTPR SERT A LA LECTURE DES FICHIERS DE CALCULS PRECEDENTS
C A PRIORI TEXTPR=TEXTE MAIS ON PEUT ESSAYER DE FAIRE UNE SUITE
C DE CALCUL A PARTIR D'UN AUTRE CODE.
C
      TEXTPR (1 ) = 'VITESSE U       M/S             '
      TEXTPR (2 ) = 'VITESSE V       M/S             '
      TEXTPR (3 ) = 'CELERITE        M/S             '
      TEXTPR (4 ) = 'HAUTEUR D''EAU   M               '
      TEXTPR (5 ) = 'SURFACE LIBRE   M               '
      TEXTPR (6 ) = 'FOND            M               '
      TEXTPR (7 ) = 'FROUDE                          '
      TEXTPR (8 ) = 'DEBIT SCALAIRE  M2/S            '
      TEXTPR (9 ) = 'EX TRACEUR                      '
      TEXTPR (10) = 'ENERGIE TURBUL. JOULE/KG        '
      TEXTPR (11) = 'DISSIPATION     WATT/KG         '
      TEXTPR (12) = 'VISCOSITE TURB. M2/S            '
      TEXTPR (13) = 'DEBIT SUIVANT X M2/S            '
      TEXTPR (14) = 'DEBIT SUIVANT Y M2/S            '
      TEXTPR (15) = 'VITESSE SCALAIREM/S             '
      TEXTPR (16) = 'VENT X          M/S             '
      TEXTPR (17) = 'VENT Y          M/S             '
      TEXTPR (18) = 'PRESSION ATMOS. PASCAL          '
      TEXTPR (19) = 'FROTTEMENT                      '
      TEXTPR (20) = 'DERIVE EN X     M               '
      TEXTPR (21) = 'DERIVE EN Y     M               '
      TEXTPR (22) = 'NBRE DE COURANT                 '
      TEXTPR (23) = 'VARIABLE 23     UNITES ??       '
      TEXTPR (24) = 'VARIABLE 24     UNITES ??       '
      TEXTPR (25) = 'VARIABLE 25     UNITES ??       '
      TEXTPR (26) = 'VARIABLE 26     UNITES ??       '
      TEXTPR (27) = 'COTE MAXIMUM    M               '
      TEXTPR (28) = 'TEMPS COTE MAXI S               '
      TEXTPR (29) = 'VITESSE MAXIMUM M/S             '
      TEXTPR (30) = 'T VITESSE MAXI  S               '
      TEXTPR (31) = 'VITESSE DE FROT.M/S             '
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C   ALIAS DES NOMS DE VARIABLES POUR LE FICHIER DES PARAMETRES
C
C     UVCHSBFQTKEDIJMXYPWAGLNORZ
C     VITESSE U
      MNEMO(1)   = 'U       '
C     VITESSE V
      MNEMO(2)   = 'V       '
C     CELERITE
      MNEMO(3)   = 'C       '
C     HAUTEUR D'EAU
      MNEMO(4)   = 'H       '
C     SURFACE LIBRE
      MNEMO(5)   = 'S       '
C     FOND
      MNEMO(6)   = 'B       '
C     FROUDE
      MNEMO(7)   = 'F       '
C     DEBIT SCALAIRE
      MNEMO(8)   = 'Q       '
C     EX TRACEUR
      MNEMO(9)   = '?       '
C     ENERGIE TURBUL.
      MNEMO(10)   = 'K       '
C     DISSIPATION
      MNEMO(11)   = 'E       '
C     VISCOSITE TURB.
      MNEMO(12)   = 'D       '
C     DEBIT SUIVANT X
      MNEMO(13)   = 'I       '
C     DEBIT SUIVANT Y
      MNEMO(14)   = 'J       '
C     VITESSE SCALAIRE
      MNEMO(15)   = 'M       '
C     VENT X
      MNEMO(16)   = 'X       '
C     VENT Y
      MNEMO(17)   = 'Y       '
C     PRESSION ATMOS.
      MNEMO(18)   = 'P       '
C     FROTTEMENT
      MNEMO(19)   = 'W       '
C     DERIVE EN X
      MNEMO(20)   = 'A       '
C     DERIVE EN Y
      MNEMO(21)   = 'G       '
C     NBRE DE COURANT
      MNEMO(22)   = 'L       '
C     VARIABLE 23
      MNEMO(23)   = 'N       '
C     VARIABLE 24
      MNEMO(24)   = 'O       '
C     VARIABLE 25
      MNEMO(25)   = 'R       '
C     VARIABLE 26
      MNEMO(26)   = 'Z       '
C     VARIABLE 27
      MNEMO(27)   = 'MAXZ    '
C     VARIABLE 28
      MNEMO(28)   = 'TMXZ    '
C     VARIABLE 29
      MNEMO(29)   = 'MAXV    '
C     VARIABLE 30
      MNEMO(30)   = 'TMXV    '
C     VARIABLE 31
      MNEMO(31)   = 'US      '
C
C-----------------------------------------------------------------------
C
C     ANALYSES DE FOURIERS
C
      IF(NPERIAF.GT.0) THEN
        DO I=1,NPERIAF
          IF(LNG.EQ.1) THEN
            TEXTE(32+NTRAC+2*(I-1)) =  'AMPLI PERIODE '
     *                         //I_IN_2_LETTERS(I)
     *                         //'M               '
            TEXTE(33+NTRAC+2*(I-1)) =  'PHASE PERIODE '
     *                         //I_IN_2_LETTERS(I)
     *                         //'DEGRES          '
            TEXTPR(32+NTRAC+2*(I-1)) =  'AMPLI PERIODE '
     *                         //I_IN_2_LETTERS(I)
     *                         //'M               '
            TEXTPR(33+NTRAC+2*(I-1)) =  'PHASE PERIODE '
     *                         //I_IN_2_LETTERS(I)
     *                         //'DEGRES          '
          ELSE
            TEXTE(32+NTRAC+2*(I-1)) =  'AMPLI PERIOD  '
     *                         //I_IN_2_LETTERS(I)
     *                         //'M               '
            TEXTE(33+NTRAC+2*(I-1)) =  'PHASE PERIOD  '
     *                         //I_IN_2_LETTERS(I)
     *                         //'DEGRES          '
            TEXTPR(32+NTRAC+2*(I-1)) =  'AMPLI PERIOD  '
     *                         //I_IN_2_LETTERS(I)
     *                         //'M               '
            TEXTPR(33+NTRAC+2*(I-1)) =  'PHASE PERIOD  '
     *                         //I_IN_2_LETTERS(I)
     *                         //'DEGRES          '
          ENDIF
          MNEMO(32+NTRAC+2*(I-1)) = 'AMPL'//I_IN_2_LETTERS(I)//'  '
          MNEMO(33+NTRAC+2*(I-1)) = 'PHAS'//I_IN_2_LETTERS(I)//'  '
        ENDDO 
      ENDIF
C
C-----------------------------------------------------------------------
C
C     TRACEURS
C
      IF(NTRAC.GT.0) THEN
        DO I=1,NTRAC
          TEXTE(31+I)  = NAMETRAC(I)
          TEXTPR(31+I) = NAMETRAC(I)
          MNEMO(31+I)  = 'T'//I_IN_2_LETTERS(I)//'   '
        ENDDO 
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
