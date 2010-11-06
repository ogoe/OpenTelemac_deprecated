C                       *****************
                        SUBROUTINE TRACVF
C                       *****************
C
     *(F,FN,FSCEXP,H,HN,FXMAT,FXMATPAR,
     * V2DPAR,UNSV2D,DDT,FXBOR,FBOR,SMH,YASMH,T1,T2,T4,T5,T6,T7,T8,
     * MESH,LIMTRA,KDIR,KDDL,OPTSOU,IOPT2,FLBORTRA,MSK,IT,DT,TDT)
C
C***********************************************************************
C BIEF VERSION 6.0           06/02/09     C-T PHAM (LNHE) 01 30 87 85 93
C***********************************************************************
C
C  FONCTION  : CALCUL DU TRACEUR POUR SCHEMA VOLUMES FINIS
C              A COMPLETER
C
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|_______________________________________________
C |    F           |<-- | VALEURS DU TRACEUR A L'ETAPE N+1.
C |    FN          | -->| VALEURS DU TRACEUR A L'ETAPE N.
C |    H           | -->| VALEURS DE LA HAUTEUR D'EAU A L'ETAPE N+1.
C |    HN          | -->| VALEURS DE LA HAUTEUR D'EAU A L'ETAPE N.
C |    FXMAT       | -->| MATRICE DE STOCKAGE DES FLUX.
C |    MAS         | -->| VECTEUR MASS ASSEMBLE LUMPE.
C |    DT          | -->| PAS DE TEMPS.
C |    FXBOR       | -->| MATRICE DES FLUX SUR LE BORD.
C |    FBOR        | -->| VALEURS DU TRACEUR SUR LE BORD.
C |    SM          | -->| TERMES SOURCES.
C |    SMH         | -->| TERME SOURCE DE L'EQUATION DE CONTINUITE.
C |    MESH        | -->| STRUCTURE DE MAILLAGE.
C |    IOPT2       | -->| 0 : UCONV RESPECTE LA CONTINUITE
C |                |    | 1 : UCONV NE RESPECTE PAS LA CONTINUITE
C |    MSK         | -->| IF YES, MASKING OF DRY ELEMENTS
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : CVDFTR
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_TRACVF => TRACVF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)           :: KDIR,KDDL,OPTSOU,LIMTRA(*)
      INTEGER, INTENT(IN)           :: IOPT2,IT
      DOUBLE PRECISION, INTENT(IN)  :: DDT,DT,TDT
      TYPE(BIEF_OBJ), INTENT(INOUT) :: F,T1,T2,T4,T5,T6,T7,T8,FLBORTRA
      TYPE(BIEF_OBJ), INTENT(IN)    :: FN,H,HN,V2DPAR,SMH,FBOR,FSCEXP
      TYPE(BIEF_OBJ), INTENT(IN)    :: FXBOR,UNSV2D
      DOUBLE PRECISION, INTENT(IN)  :: FXMAT(*),FXMATPAR(*)
      TYPE(BIEF_MESH), INTENT(INOUT):: MESH
      LOGICAL, INTENT(IN)           :: YASMH,MSK
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C
C-----------------------------------------------------------------------
C
      IF(IOPT2.EQ.0) THEN
C
C-----------------------------------------------------------------------
C
C     CAS OU LE CHAMP CONVECTEUR RESPECTE L'EQUATION DE CONTINUITE
C     (LA HAUTEUR POURRAIT ETRE CALCULEE PAR INTERPOLATION EN TEMPS)
C
C     T4 VA PRENDRE LES VALEURS SUCCESSIVES DE F (INITIALISE DANS CVTRVF)
C
      CALL TVF(F%R,FN%R,T4%R,T5%R,FXMAT,FXMATPAR,UNSV2D%R,DDT,
     *         FXBOR%R,T7%R,T8,FBOR%R,SMH%R,YASMH,FSCEXP%R,
     *         MESH%NSEG,MESH%NPOIN,MESH%NPTFR,
     *         MESH%GLOSEG%I,MESH%GLOSEG%DIM1,
     *         MESH%NBOR%I,LIMTRA,KDIR,KDDL,
     *         OPTSOU,T5%R,IOPT2,FLBORTRA%R,DDT/DT,MESH,F)
C
C-----------------------------------------------------------------------
C
C     CAS OU LE CHAMP CONVECTEUR NE RESPECTE PAS L'EQUATION DE CONTINUITE
C
      ELSEIF(IOPT2.EQ.1) THEN
C
C     T1 VA PRENDRE LES VALEURS SUCCESSIVES DE HN CALCULE AVEC LA CONTINUITE
C     T2 VA PRENDRE LES VALEURS SUCCESSIVES DE H  CALCULE AVEC LA CONTINUITE
C     T4 VA PRENDRE LES VALEURS SUCCESSIVES DE F
C     T5 VA PRENDRE LES VALEURS SUCCESSIVES DE LA VRAIE HAUTEUR
C
C     H2 HAUTEUR PAR EQUATION DE CONTINUITE
C
      CALL HVF(T2%R,T1%R,FXMAT,UNSV2D%R,DDT,T7%R,SMH%R,
     *         YASMH,MESH%NSEG,MESH%NPOIN,MESH%NPTFR,
     *         MESH%GLOSEG%I,MESH%GLOSEG%DIM1,MESH%NBOR%I,OPTSOU,
     *         T8,MESH,MSK)
C
      CALL TVF(F%R,FN%R,T4%R,T2%R,FXMAT,FXMATPAR,UNSV2D%R,DDT,
     *         FXBOR%R,T7%R,T8,FBOR%R,SMH%R,YASMH,FSCEXP%R,
     *         MESH%NSEG,MESH%NPOIN,MESH%NPTFR,
     *         MESH%GLOSEG%I,MESH%GLOSEG%DIM1,
     *         MESH%NBOR%I,LIMTRA,KDIR,KDDL,
     *         OPTSOU,T5%R,IOPT2,FLBORTRA%R,DDT/DT,MESH,F)
C
C-----------------------------------------------------------------------
C
      ELSE
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'TRACVF : OPTION INCONNUE'
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'TRACVF: UNKNOWN OPTION'
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
