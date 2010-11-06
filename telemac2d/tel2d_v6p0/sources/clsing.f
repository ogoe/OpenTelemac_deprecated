C                       *****************
                        SUBROUTINE CLSING
C                       *****************
C
     *(NWEIRS,NPSING,NPSMAX,NUMDIG,X,Y,ZF,CHESTR,NKFROT,KARMAN,
     * ZDIG,PHIDIG,NBOR,H,T,NTRAC,IOPTAN,UNORM,
     * UBOR,VBOR,TBOR,LIHBOR,LIUBOR,LIVBOR,LITBOR,GRAV)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.6  19/04/96     V. GUINOT   (LHF)
C              MODIFIE LE  23/11/05  J.-M. HERVOUET (LNH) 01 30 87 80 18
C***********************************************************************
C
C      FONCTION:   GESTION DU CALCUL DES DEBITS ET REMPLISSAGE DES
C      =========   CONDITIONS AUX LIMITES
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   NWEIRS       | -->| NOMBRE DE SINGULARITES LINEIQUES.
C |   NPSING(N)    | -->| NOMBRE DE POINTS DE CHAQUE COTE DE LA
C |                |    | SINGULARITE N.
C |   NPSMAX       | -->| NOMBRE MAXIMUM DE POINTS POUR UN COTE D'UNE
C |                |    | SINGULARITE.
C |   NUMDIG(K,N,I)| -->| NUMERO DES POINTS DES DIGUES
C |                |    | DANS LA NUMEROTATION DES POINTS DE BORD
C |                |    | DES CONDITIONS AUX LIMITES) DU I-EME
C |                |    | POINT SUR LE COTE K DE L'OUVRAGE N
C |   X,Y          | -->| COORDONNEES DES NOUEDS.
C |   ZF           | -->| COTE DU FOND.
C |   CHESTR       | -->| COEFFICIENTS DE FROTTEMENT SUR LE FOND.
C |   KFROT        | -->| LOI DE FROTTEMENT SUR LE FOND.
C |   KARMAN       | -->| CONSTANTE DE KARMAN.
C |   ZDIG(N,I)    | -->| COTE DE LA DIGUE AU I-EME POINT DE
C |                |    | LA N-IEME SINGULARITE
C |   PHIDIG(N,I)  | -->| COEFFICIENT DE DEBIT AU I-EME POINT DE
C |                |    | LA N-IEME SINGULARITE
C |   NBOR         | -->| NUMEROTATION GLOBALE DES POINTS DE BORD.
C |   H            | -->| HAUTEUR AU PAS DE TEMPS COURANT.
C |   T            | -->| TRACEUR AU PAS DE TEMPS COURANT.
C |   TRAC         | -->| SI OUI, IL Y A UN TRACEUR.
C |   IOPTAN       | -->| OPTION DE CALCUL DES VITESSES TANGENTIELLES.
C |   UNORM (J)    | -->| VITESSE NORMALE A LA SOUS-ITERATION COURANTE
C |   UBOR(J)      |<-- | VALEUR DE LA CONDITION EN VIESSE U AU
C |                |    | J-EME POINT LIMITE
C |   VBOR(J)      |<-- | VALEUR DE LA CONDITION EN VITESSE V A
C |                |    | J-EME POINT LIMITE
C |   TBOR(J)      |<-- | VALEUR DE LA CONDITION EN TRACEUR AU
C |                |    | J-EME POINT LIMITE
C |   LIHBOR(J)    |<-- | TYPE DE LA CONDITION EN HAUTEUR AU
C |                |    | J-EME POINT LIMITE
C |   LIUBOR(J)    |<-- | TYPE DE LA CONDITION EN VITESSE U AU
C |                |    | J-EME POINT LIMITE
C |   LIVBOR(J)    |<-- | TYPE DE LA CONDITION EN VITESSE V AU
C |                |    | J-EME POINT LIMITE
C |   LITBOR(J)    |<-- | TYPE DE LA CONDITION EN TRACEUR AU
C |                |    | J-EME POINT LIMITE
C |   GRAV         | -->| GRAVITE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : CLSING
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C***********************************************************************
C
      USE BIEF
      USE INTERFACE_TELEMAC2D, EX_CLSING => CLSING
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NWEIRS,NPSMAX,IOPTAN
      INTEGER, INTENT(IN) :: NKFROT(*),NBOR(*)
      INTEGER, INTENT(IN) :: NPSING(NWEIRS),NUMDIG(2,NWEIRS,NPSMAX)
      INTEGER, INTENT(INOUT) :: LIUBOR(*),LIVBOR(*),LIHBOR(*) 
      INTEGER, INTENT(IN) :: NTRAC
      DOUBLE PRECISION, INTENT(IN) :: PHIDIG(NWEIRS,NPSMAX)
      DOUBLE PRECISION, INTENT(IN) :: ZDIG(NWEIRS,NPSMAX),H(*)
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*),ZF(*),CHESTR(*)
      DOUBLE PRECISION, INTENT(IN) :: KARMAN,GRAV
      DOUBLE PRECISION, INTENT(INOUT) :: UBOR(*),VBOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: UNORM(*)
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: TBOR,LITBOR
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: T
C  
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,N,IA,IB,NA,NB                      
C
      DOUBLE PRECISION HMIN,PHI,QAB,YAA,YBB,YDEN,YS 
C
C-----------------------------------------------------------------------
C
      HMIN=1.D-3
C
C     CALCUL DES DEBITS UNITAIRES
C
      DO 10 N=1,NWEIRS
      DO 20 I=1,NPSING(N)
        IA=NUMDIG(1,N,I)
        IB=NUMDIG(2,N,I)
        NA=NBOR(IA)
        NB=NBOR(IB)
        YAA=H(NA)+ZF(NA)
        YBB=H(NB)+ZF(NB)
        YS=ZDIG(N,I)
        PHI=PHIDIG(N,I)
C
        IF(YAA.GT.YBB) THEN
C         CAS OU L'AMONT EST EN A
          YDEN=YS/3.D0+2.D0*YAA/3.D0
          IF(YBB.LT.YDEN) THEN
            CALL LOIDEN(YAA,YS,PHI,QAB,GRAV)
          ELSE
            CALL LOINOY(YAA,YBB,YS,PHI,QAB,GRAV)
          ENDIF
        ELSE
C         CAS OU L'AMONT EST EN B
          YDEN=YS/3.D0+2.D0*YBB/3.D0
          IF(YAA.LT.YDEN) THEN
            CALL LOIDEN(YBB,YS,PHI,QAB,GRAV)
          ELSE
            CALL LOINOY(YBB,YAA,YS,PHI,QAB,GRAV)
          ENDIF
          QAB=-QAB
        ENDIF
C
C CALCUL DU DEBIT NORMAL
C
        IF(H(NA).LE.HMIN) THEN
          UNORM(IA)=0.D0
        ELSE
          UNORM(IA)=-QAB/H(NA)
        ENDIF
C
        IF(H(NB).LE.HMIN) THEN
          UNORM(IB)=0.D0
        ELSE
          UNORM(IB)=-QAB/H(NB)
        ENDIF
C
20    CONTINUE
10    CONTINUE
C
C     DETERMINATION DE LA VALEUR NUMERIQUE
C     DES CONDITIONS AUX LIMITES :
C
      CALL CLHUVT(NWEIRS,NPSING,NPSMAX,NUMDIG,ZDIG,X,Y,ZF,
     *            IOPTAN,UNORM,CHESTR,NKFROT,KARMAN,T,NTRAC,H,
     *            UBOR,VBOR,TBOR,NBOR,LIHBOR,LIUBOR,LIVBOR,LITBOR)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
