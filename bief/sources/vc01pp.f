C                       *****************
                        SUBROUTINE VC01PP
C                       *****************
C
     *( XMUL,SF,F,Z,SURFAC,
     *  IKLE1,IKLE2,IKLE3,IKLE4,IKLE5,IKLE6,NELEM,NELMAX,
     *  W1,W2,W3,W4,W5,W6)
C
C***********************************************************************
C BIEF VERSION 5.1           09/12/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C                    /
C    VEC(I) = XMUL  /    PSI(I) * F  D(OMEGA)
C                  /OMEGA
C
C
C    PSI(I) EST UNE BASE DE TYPE PRISME P1
C
C    F EST UN VECTEUR DE TYPE IELMF
C
C
C    ATTENTION : LE RESULTAT EST DANS W  SOUS FORME NON ASSEMBLEE.
C
C                MAILLAGE REEL
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      XMUL      | -->|  COEFFICIENT MULTIPLICATEUR.
C |      SF,SG,SH  | -->|  STRUCTURES DES FONCTIONS F,G ET H
C |      SU,SV,SW  | -->|  STRUCTURES DES FONCTIONS U,V ET W
C |      F,G,H     | -->|  FONCTIONS INTERVENANT DANS LA FORMULE.
C |      U,V,W     | -->|  COMPOSANTES D'UN VECTEUR
C |                |    |  INTERVENANT DANS LA FORMULE.
C |      XEL,YEL,..| -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |      SURFAC    | -->|  SURFACE DES ELEMENTS.
C |      IKLE1,... | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE.
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE.
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE.
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |      W1,2,3    |<-- |  VECTEUR RESULTAT SOUS FORME NON ASSEMBLEE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C  PROGRAMMES APPELES : ASSVEC
C
C-----------------------------------------------------------------------
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C***********************************************************************
C
      USE BIEF, EX_VC01PP => VC01PP
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX),IKLE3(NELMAX)
      INTEGER, INTENT(IN) :: IKLE4(NELMAX),IKLE5(NELMAX),IKLE6(NELMAX)
C
      DOUBLE PRECISION, INTENT(IN) :: Z(*)
      DOUBLE PRECISION, INTENT(IN) :: SURFAC(NELMAX)
      DOUBLE PRECISION,INTENT(INOUT)::W1(NELMAX),W2(NELMAX),W3(NELMAX)
      DOUBLE PRECISION,INTENT(INOUT)::W4(NELMAX),W5(NELMAX),W6(NELMAX)
      DOUBLE PRECISION, INTENT(IN) :: XMUL
C
C     STRUCTURE DE F ET DONNEES REELLES
C
      TYPE(BIEF_OBJ),   INTENT(IN) :: SF
      DOUBLE PRECISION, INTENT(IN) :: F(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER IELEM,IELMF
      DOUBLE PRECISION SUR360,COEF,H1,H2,H3,SHT,SH1,SH2,SH3
      DOUBLE PRECISION F1,F2,F3,F4,F5,F6,SFI,SFS,SF1,SF2,SF3,SF4,SF5,SF6
      DOUBLE PRECISION HF1,HF2,HF3,HF4,HF5,HF6,SHFI,SHFS
      DOUBLE PRECISION SHF1,SHF2,SHF3,SHF4,SHF5,SHF6
C
C***********************************************************************
C
      IELMF=SF%ELM
C
C-----------------------------------------------------------------------
C
C   F LINEAIRE
C
      IF(IELMF.EQ.41) THEN
C
         SUR360 = XMUL / 360.D0
C
         DO 3 IELEM = 1 , NELEM
C
            COEF = SUR360 * SURFAC(IELEM)
C
            H1  = COEF * (Z(IKLE4(IELEM)) - Z(IKLE1(IELEM)))
            H2  = COEF * (Z(IKLE5(IELEM)) - Z(IKLE2(IELEM)))
            H3  = COEF * (Z(IKLE6(IELEM)) - Z(IKLE3(IELEM)))
            SHT = H1 + H2 + H3
            SH1 = H1 + SHT
            SH2 = H2 + SHT
            SH3 = H3 + SHT
C
            F1  = F(IKLE1(IELEM))
            F2  = F(IKLE2(IELEM))
            F3  = F(IKLE3(IELEM))
            F4  = F(IKLE4(IELEM))
            F5  = F(IKLE5(IELEM))
            F6  = F(IKLE6(IELEM))
            SFI = F1 + F2 + F3
            SFS = F4 + F5 + F6
            SF1 = F1 + SFI
            SF2 = F2 + SFI
            SF3 = F3 + SFI
            SF4 = F4 + SFS
            SF5 = F5 + SFS
            SF6 = F6 + SFS
C
            HF1  = H1 * F1
            HF2  = H2 * F2
            HF3  = H3 * F3
            HF4  = H1 * F4
            HF5  = H2 * F5
            HF6  = H3 * F6
            SHFI = HF1 + HF2 + HF3
            SHFS = HF4 + HF5 + HF6
            SHF1 = HF1 + SHFI
            SHF2 = HF2 + SHFI
            SHF3 = HF3 + SHFI
            SHF4 = HF4 + SHFS
            SHF5 = HF5 + SHFS
            SHF6 = HF6 + SHFS
C
            W1(IELEM) = SH1 * (SF1+SF1+SF4) + SHF1 + SHF1 + SHF4
            W2(IELEM) = SH2 * (SF2+SF2+SF5) + SHF2 + SHF2 + SHF5
            W3(IELEM) = SH3 * (SF3+SF3+SF6) + SHF3 + SHF3 + SHF6
            W4(IELEM) = SH1 * (SF1+SF4+SF4) + SHF1 + SHF4 + SHF4
            W5(IELEM) = SH2 * (SF2+SF5+SF5) + SHF2 + SHF5 + SHF5
            W6(IELEM) = SH3 * (SF3+SF6+SF6) + SHF3 + SHF6 + SHF6
C
3        CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSE
C
       IF (LNG.EQ.1) WRITE(LU,101) IELMF,SF%NAME
       IF (LNG.EQ.2) WRITE(LU,102) IELMF,SF%NAME
101    FORMAT(1X,'VC01PP (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,' CAS NON PREVU',/,
     *        1X,'NOM REEL DE F : ',A6)
102    FORMAT(1X,'VC01PP (BIEF) :',/,
     *        1X,'DISCRETISATION OF F : ',1I6,' NOT IMPLEMENTED',/,
     *        1X,'REAL NAME OF F: ',A6)
       CALL PLANTE(1)
       STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
