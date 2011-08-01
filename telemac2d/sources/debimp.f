C                       *****************
                        SUBROUTINE DEBIMP
C                       *****************
C
     *(Q,UBOR,VBOR,U,V,H,NUMLIQ,IFRLIQ,WORK1,WORK2,NPTFR,MASK,MESH,
     * KP1BOR,EQUA)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2    24/04/97    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C      FONCTION:  IMPOSE DES CONDITIONS AUX LIMITES DE DEBIT ,
C                 AVEC L'HYPOTHESE D'AFFINITE DES PROFILS DE
C                 VITESSE A L'ENTREE.
C
C      PRINCIPE:  1) CALCUL DU DEBIT SUR LA FRONTIERE COMPRISE ENTRE
C                    LES POINTS DE BORD NDEB ET NFIN.
C
C                 2) CORRECTION SUR LES VITESSES IMPOSEES TROUVEES ENTRE
C                    NDEB ET NFIN POUR QUE LE DEBIT DEVIENNE EGAL AU
C                    DEBIT VOULU. CETTE CORRECTION EST UNE SIMPLE REGLE
C                    DE TROIS ET NE MODIFIE DONC PAS LE PROFIL DES
C                    VITESSES U,V DONNE OU CALCULE AUPARAVANT.
C
C      LE RESULTAT EST MIS DANS LES TABLEAUX UBOR ET VBOR POUR ETRE
C      IMPOSE A U ET V AU COURS DES DIFFERENTES ETAPES DE CALCUL.
C
C      SI LE DEBIT A TRAVERS UNE FRONTIERE EST NUL, ON NE PEUT PAS
C      FAIRE DE REGLE DE TROIS, ON REMPLACE ALORS UBOR ET VBOR PAR
C      U ET V OBTENUS AU PAS DE TEMPS N.
C
C      REMARQUE : LE CAS OU LE SEGMENT (NDEB,NFIN) PRESENTE UNE DIS-
C                 CONTINUITE DANS LA NUMEROTATION (CAS OU CE SEGMENT
C                 CONTIENT LE POINT NUMERO 1) EST TRAITE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   Q            | -->|  VALEUR DU DEBIT IMPOSE                      |
C |   NDEB,NFIN    | -->|  NUMEROS DES POINTS DE BORD ENTRE LESQUELS ON|
C |                |    |  VEUT IMPOSER LE DEBIT                       |
C |   UBOR,VBOR    |<-- |  VALEURS DE U ET V DIRICHLET A L'ENTREE      |
C |   U , V , H    | -->|  VALEURS DE U,V,H AU TEMPS N
C |   WORK1,WORK2  | -->|  TABLEAUX DE TRAVAIL.                        |
C |   NPTFR        | -->|  NOMBRE DE POINTS FRONTIERE.                 |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : BORD (SI ON IMPOSE UN DEBIT)
C
C SOUS-PROGRAMME APPELE : FLUBOR
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
      INTEGER, INTENT(IN)             :: NPTFR,IFRLIQ
      INTEGER, INTENT(IN)             :: NUMLIQ(NPTFR),KP1BOR(NPTFR,2)
      CHARACTER(LEN=20), INTENT(IN)   :: EQUA
      DOUBLE PRECISION, INTENT(INOUT) :: UBOR(NPTFR),VBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: MASK(*),Q
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
      TYPE(BIEF_OBJ), INTENT(IN)      :: H,U,V
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: WORK1,WORK2
C       
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K,IELM      
C     
      DOUBLE PRECISION Q1
C
      INTRINSIC ABS
C
      DOUBLE PRECISION P_DSUM
      EXTERNAL         P_DSUM
C
C=======================================================================
C     CALCUL DU FLUX
C=======================================================================
C
C  DANS LA BOUCLE SUIVANTE ON RESTREINT LE MASQUE DES SEGMENTS DIRICHLETS
C  A CEUX DE LA FRONTIERE LIQUIDE NUMERO IFRLIQ. COMME NUMLIQ EST
C  DONNE PAR POINTS, ON RISQUE UNE ERREUR POUR LE SEGMENT SUIVANT
C  LE DERNIER POINT DE LA FRONTIERE. EN FAIT CE SEGMENT SERA SOLIDE
C  ET AURA UN MASQUE DEJA NUL.
C
      IF(EQUA(1:15).NE.'SAINT-VENANT VF') THEN
C
      DO K=1,NPTFR
        IF(NUMLIQ(K).EQ.IFRLIQ) THEN
          WORK1%R(K)=MASK(K)
        ELSE
          WORK1%R(K)=0.D0
        ENDIF
      ENDDO
C
      ELSE
C
C     LES VOLUMES FINIS COMPTENT LES SEGMENTS SOLIDES VOISINS
C     DES FRONTIERES LIQUIDES
C
      DO K=1,NPTFR
        IF(NUMLIQ(K).EQ.IFRLIQ) THEN
          WORK1%R(K)          =MASK(K)
          WORK1%R(KP1BOR(K,1))=MASK(K)
          WORK1%R(KP1BOR(K,2))=MASK(K)
        ELSE
          WORK1%R(K)=0.D0
        ENDIF
      ENDDO
C
      ENDIF
C
      IELM=11
      CALL VECTOR(WORK2,'=','FLUBDF          ',IELBOR(IELM,1),
     *            1.D0,H,H,H,U,V,V,MESH,.TRUE.,WORK1)
C     CONVENTION DE SIGNE INVERSEE ENTRE UTILISATEUR ET CODE
C     POUR L'UTILISATEUR : DEBIT POSITIF = ENTRANT
C     POUR LE CODE : U.N < 0 = ENTRANT
      Q1 = - BIEF_SUM(WORK2)
      IF(NCSIZE.GT.1) Q1 = P_DSUM(Q1)
C
      IF(ABS(Q1).LT.1.D-10) THEN
C
C DEBIT NUL : MESSAGE D'AVERTISSEMENT
C
        IF(ABS(Q).GT.1.D-10) THEN
          IF(LNG.EQ.1) WRITE(LU,30) IFRLIQ
          IF(LNG.EQ.2) WRITE(LU,31) IFRLIQ
30        FORMAT(1X,'DEBIMP : PROBLEME SUR LA FRONTIERE ',1I6,/,1X,
     *     '         DONNER UN PROFIL DE VITESSES        ',/,1X,
     *     '         DANS LE FICHIER DES CONDITIONS AUX LIMITES',/,1X,
     *     '         OU VERIFIER LES HAUTEURS D''EAU')
31        FORMAT(1X,'DEBIMP : PROBLEM ON BOUNDARY NUMBER ',1I6,/,1X,
     *     '         GIVE A VELOCITY PROFILE  ',/,1X,
     *     '         IN THE BOUNDARY CONDITIONS FILE',/,1X,
     *     '         OR CHECK THE WATER DEPTHS')
          CALL PLANTE(1)
          STOP
        ELSE
          Q1 = 1.D0
        ENDIF
C
      ENDIF
C
C=======================================================================
C   CALCUL DE UBOR ET VBOR
C=======================================================================
C
      DO 40 K=1,NPTFR
C
        IF(NUMLIQ(K).EQ.IFRLIQ) THEN
          UBOR(K) = UBOR(K) * Q / Q1
          VBOR(K) = VBOR(K) * Q / Q1
        ENDIF
C
40    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
