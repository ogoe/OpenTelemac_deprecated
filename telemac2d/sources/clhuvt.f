C                       *****************
                        SUBROUTINE CLHUVT
C                       *****************
C
     *(NWEIRS,NPSING,NPSMAX,NUMDIG,ZDIG,
     * X,Y,ZF,IOPTAN,UNORM,CHESTR,
     * NKFROT,KARMAN,T,NTRAC,H,UBOR,VBOR,TBOR,NBOR,
     * LIHBOR,LIUBOR,LIVBOR,LITBOR)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.6  19/04/96     V. GUINOT   (LHF)
C              MODIFIE LE  23/11/05  J.-M. HERVOUET (LNH) 01 30 87 80 18
C***********************************************************************
C
C      FONCTION:   DETERMINATION DES HAUTEURS, VITESSES... A IMPOSER AU
C      =========
C                  NIVEAU DES POINTS, A PARTIR DES HAUTEURS ET DEBITS
C
C                  MOYENS SUR LES SEGMENTS CONSTITUANT LA SINGULARITE.
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
C |   ZDIG(N,I)    | -->| COTE DE LA DIGUE AU I-EME POINT DE
C |                |    | LA N-IEME SINGULARITE
C |   X,Y          | -->| COORDONNEES DES NOUEDS.
C |   ZF           | -->| COTE DU FOND.
C |   IOPTAN       | -->| OPTION DE CALCUL DES VITESSES TANGENTIELLES.
C |   UNORM(J)     | -->| VITESSE NORMALE.
C |   CHESTR       | -->| COEFFICIENT DE FROTTEMENT.
C |   KFROT        | -->| LOI DE FROTTEMENT SUR LE FOND.
C |   KARMAN       | -->| CONSTANTE DE KARMAN.
C |   T            | -->| TRACEUR.
C |   TRAC         | -->| SI OUI, IL Y A UN TRACEUR.
C |   H            | -->| HAUTEUR.
C |   UBOR(J)      |<-- | VALEUR DE LA CONDITION EN VIESSE U AU
C |                |    | J-EME POINT LIMITE
C |   VBOR(J)      |<-- | VALEUR DE LA CONDITION EN VITESSE V A
C |                |    | J-EME POINT LIMITE
C |   TBOR(J)      |<-- | VALEUR DE LA CONDITION EN TRACEUR AU
C |                |    | J-EME POINT LIMITE
C |   NBOR         | -->| NUMEROTATION GLOBALE DES POINTS DE BORD.
C |   LIHBOR(J)    |<-- | TYPE DE LA CONDITION EN HAUTEUR AU
C |                |    | J-EME POINT LIMITE
C |   LIUBOR(J)    |<-- | TYPE DE LA CONDITION EN VITESSE U AU
C |                |    | J-EME POINT LIMITE
C |   LIVBOR(J)    |<-- | TYPE DE LA CONDITION EN VITESSE V AU
C |                |    | J-EME POINT LIMITE
C |   LITBOR(J)    |<-- | TYPE DE LA CONDITION EN TRACEUR AU
C |                |    | J-EME POINT LIMITE
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
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C  
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NWEIRS,NPSMAX,IOPTAN
      INTEGER, INTENT(IN)    :: NPSING(NWEIRS),NUMDIG(2,NWEIRS,NPSMAX)
      INTEGER, INTENT(IN)             :: NBOR(*),NKFROT(*)
      INTEGER, INTENT(INOUT)          :: LIUBOR(*),LIHBOR(*),LIVBOR(*)
      INTEGER, INTENT(IN)             :: NTRAC
      DOUBLE PRECISION, INTENT(IN)    :: ZDIG(NWEIRS,NPSMAX)
      DOUBLE PRECISION, INTENT(IN)    :: X(*),Y(*),ZF(*),CHESTR(*),H(*)
      DOUBLE PRECISION, INTENT(IN)    :: UNORM(*)
      DOUBLE PRECISION, INTENT(INOUT) :: UBOR(*),VBOR(*)
      DOUBLE PRECISION, INTENT(IN)    :: KARMAN
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: TBOR,LITBOR
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: T
C    
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,I1,I2,K,N,N0,N1,N2,ITRAC       
C
      DOUBLE PRECISION DL,NX,NY,PENTE,CZ,HH,TX,TY,UTAN,XX
C
      INTRINSIC ABS,SQRT,SIGN
C
C-----------------------------------------------------------------------
C
C     BOUCLE SUR LES OUVRAGES
C
      DO 10 N=1,NWEIRS
C
C       BOUCLE SUR LES COTES DE L'OUVRAGE
C
        DO 20 K=1,2
C
C         BOUCLE SUR LES POINTS DE CHAQUE COTE
C
          DO 30 I=1,NPSING(N)
C
          I1=NUMDIG(K,N,I)
          N1=NBOR(I1)
C
          IF(I.EQ.1) THEN
            N0=N1
            N2=NBOR(NUMDIG(K,N,I+1))
            XX=0.D0
          ELSEIF(I.LT.NPSING(N)) THEN
            N0=NBOR(NUMDIG(K,N,I-1))
            N2=NBOR(NUMDIG(K,N,I+1))
            XX=1.D0
          ELSE
            N0=NBOR(NUMDIG(K,N,I-1))
            N2=N1
            XX=0.D0
          ENDIF
C
C         CALCUL DU VECTEUR NORMAL SORTANT COTE 1, RENTRANT COTE 2
C
          TX=X(N2)-X(N0)
          TY=Y(N2)-Y(N0)
          DL=SQRT(TX*TX+TY*TY)
          TX=TX/DL
          TY=TY/DL
          NX=-TY
          NY=TX
C
C         CALCUL DE LA VITESSE TANGENTIELLE
C
          IF (IOPTAN.EQ.0) THEN
C
             UTAN=0.D0
C
          ELSEIF(IOPTAN.EQ.1) THEN
C
C            ON PREND LA HAUTEUR SUR LE SEUIL (A DEBATTRE)
C            HH = H(N1)
             HH = H(N1)+ZF(N1)-ZDIG(N,I)
C            LINE ADDED ON 23/11/2005 BY JMH (HH MAY BE NEGATIVE)
             HH=MAX(HH,0.D0)
             PENTE=(H(N0)-H(N2)+ZF(N0)-ZF(N2))/DL
C
             IF (NKFROT(N1).EQ.2) THEN
                UTAN = CHESTR(N1)*SQRT(HH*ABS(PENTE))*SIGN(1.D0,PENTE)
             ELSEIF (NKFROT(N1).EQ.3) THEN
                UTAN = CHESTR(N1)*HH**(2.D0/3.D0)*SQRT(ABS(PENTE))
     *                                           *SIGN(1.D0,PENTE)
             ELSEIF (NKFROT(N1).EQ.4) THEN
                UTAN = HH**(2.D0/3.D0)*SQRT(ABS(PENTE))
     *                                *SIGN(1.D0,PENTE)/CHESTR(N1)
             ELSEIF (NKFROT(N1).EQ.5) THEN
                HH   = MAX(HH,1.D-9)
              CZ = MAX(1.D-9,LOG(11.D0*HH/MAX(CHESTR(N1),1.D-9))/KARMAN)
                UTAN = CZ*SQRT(HH*ABS(PENTE))*SIGN(1.D0,PENTE)
             ELSE
                IF (LNG.EQ.1) THEN
                   WRITE(LU,*)'CLHUVT : OPTION INCONNUE :',NKFROT(N1)
                   WRITE(LU,*)'         POUR LA LOI DE FROTTEMENT'
                ELSEIF(LNG.EQ.2) THEN
                   WRITE(LU,*)'CLHUVT : UNKNOWN OPTION:',NKFROT(N1)
                   WRITE(LU,*)'         FOR THE FRICTION LAW'
                ENDIF
                CALL PLANTE(1)
                STOP
             ENDIF
C            POUR AVOIR DES VITESSES TANGENTIELLES NULLES DANS LES COINS
             UTAN = XX*UTAN
          ELSE
             IF (LNG.EQ.1) THEN
                WRITE(LU,*)'CLHUVT : OPTION INCONNUE :',IOPTAN
                WRITE(LU,*)'         POUR LES VITESSES TANGENTIELLES'
             ELSEIF(LNG.EQ.2) THEN
                WRITE(LU,*)'CLHUVT : UNKNOWN OPTION:',IOPTAN
                WRITE(LU,*)'         FOR THE TANGENTIAL VELOCITY'
             ENDIF
             CALL PLANTE(1)
             STOP
          ENDIF
C
C         ON CALCULE LES VITESSES U ET V
C         DANS LE REPERE (X,Y) ORDINAIRE.
C
          UBOR(I1)=UTAN*TX+UNORM(I1)*NX
          VBOR(I1)=UTAN*TY+UNORM(I1)*NY
C
30    CONTINUE
20    CONTINUE
C
C
C-----------------------------------------------------------------------
C
C  TYPES DE CONDITIONS POUR LA HAUTEUR ET LES VITESSES :
C
      DO 40 I=1,NPSING(N)
C
        I1=NUMDIG(1,N,I)
        I2=NUMDIG(2,N,I)
        LIHBOR(I1)=4
        LIUBOR(I1)=6
        LIVBOR(I1)=6
        LIHBOR(I2)=4
        LIUBOR(I2)=6
        LIVBOR(I2)=6
C
C       CORRECTION : TYPE PAROI SI LA VITESSE NORMALE EST NULLE
C
        IF(ABS(UNORM(I1)).LT.1.D-10) THEN
          LIHBOR(I1)=2
          LIUBOR(I1)=2
          LIVBOR(I1)=2
          LIHBOR(I2)=2
          LIUBOR(I2)=2
          LIVBOR(I2)=2
          IF(NTRAC.GT.0) THEN
            DO ITRAC=1,NTRAC
              LITBOR%ADR(ITRAC)%P%I(I1)=2
              LITBOR%ADR(ITRAC)%P%I(I2)=2
            ENDDO
          ENDIF
        ENDIF
C
40    CONTINUE
C
C-----------------------------------------------------------------------
C
C  TYPES DE CONDITIONS POUR LE TRACEUR ET VALEURS A IMPOSER
C
      IF(NTRAC.GT.0) THEN
C
        DO ITRAC=1,NTRAC
        DO I=1,NPSING(N)
C
          I1=NUMDIG(1,N,I)
          I2=NUMDIG(2,N,I)
C
          IF(UNORM(I1).LT.-1.D-8) THEN
C           VITESSE SORTANTE EN 1, ENTRANTE EN 2
            LITBOR%ADR(ITRAC)%P%I(I1)=4
            LITBOR%ADR(ITRAC)%P%I(I2)=5
            TBOR%ADR(ITRAC)%P%R(I2)=T%ADR(ITRAC)%P%R(NBOR(I1))
          ELSEIF(UNORM(I1).GT.1.D-8) THEN
C           VITESSE SORTANTE EN 2, ENTRANTE EN 1
            LITBOR%ADR(ITRAC)%P%I(I1)=5
            TBOR%ADR(ITRAC)%P%R(I1)=T%ADR(ITRAC)%P%R(NBOR(I2))
            LITBOR%ADR(ITRAC)%P%I(I2)=4
          ELSE
C           VITESSE NULLE
            LITBOR%ADR(ITRAC)%P%I(I1)=2
            LITBOR%ADR(ITRAC)%P%I(I2)=2
          ENDIF
C
        ENDDO
        ENDDO
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C FIN DE LA BOUCLE SUR LES SEUILS
C
10    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
