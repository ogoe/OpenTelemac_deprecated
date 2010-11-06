C                       **************
                        SUBROUTINE TVF
C                       **************
C
     *(F,FN,FC,H,FXMAT,FXMATPAR,
     * UNSV2D,DT,FXBOR,FXBORPAR,T7,FBOR,SMH,YASMH,FSCEXP,
     * NSEG,NPOIN,NPTFR,GLOSEG,SIZGLO,NBOR,LIMTRA,KDIR,KDDL,OPTSOU,HLIN,
     * IOPT2,FLBORTRA,SURNIT,MESH,SF)
C
C***********************************************************************
C BIEF VERSION 5.9           27/02/09     C-T PHAM (LNHE) 01 30 87 85 93
C
C
C  27/02/2009 JMH : DISTINCTION ENTRE FXBOR ET FXBORTRA
C
C***********************************************************************
C
C  FONCTION  : CALCUL DU TRACEUR POUR SCHEMA VOLUMES FINIS
C              A COMPLETER
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|_______________________________________________
C |    F           |<-- | VALEURS DU TRACEUR A L'ETAPE N+1 
C |                |    | DE LA SOUS-ITERATION
C |    FN          | -->| VALEURS DU TRACEUR A L'ETAPE N 
C |    FC          | -->| VALEURS DU TRACEUR A L'ETAPE N
C |                |    | DE LA SOUS-ITERATION
C |    H           | -->| VALEURS DE LA HAUTEUR D'EAU A L'ETAPE N+1.
C |                |    | EN SUPPOSANT LA CONTINUITE RESOLUE
C |    FXMAT       | -->| MATRICE DE STOCKAGE DES FLUX.
C |    FXMATPAR    | -->| IDEM, ASSEMBLE EN PARALLELE.
C |    MAS         | -->| VECTEUR MASS ASSEMBLE LUMPE.
C |    DT          | -->| PAS DE TEMPS.
C |    FXBOR       | -->| FLUX SUR LE BORD (DEFINI SUR LE BORD)
C |                |    | NON ASSEMBLE
C |    FXBORPAR    | -->| FLUX SUR LE BORD (DEFINI SUR TOUT LE DOMAINE
C |                |    | ET ASSEMBLE EN PARALLELE)
C |    FBOR        | -->| VALEURS DU TRACEUR SUR LE BORD.
C |    SM          | -->| TERMES SOURCES.
C |    SMH         | -->| TERME SOURCE DE L'EQUATION DE CONTINUITE.
C |    YASMH       | -->| IF YES, SOURCE TERMS IN SMH
C |    FSCEXP      | -->| FSCE-(1-TETAT)*FN, SEE DIFSOU
C |                |    | SO HERE FSCE-FN, THIS IS NOT VERY CONVENIENT
C |                |    | AS WE NEED HERE FSCE-FC (LOOK UNDER IF(YASMH))
C |    NSEG        | -->| NOMBRE DE SEGMENTS DANS LE MAILLAGE.
C |    NELEM       | -->| NOMBRE D'ELEMENTS DANS LE MAILLAGE.
C |    NPOIN       | -->| NOMBRE DE NOEUDS DANS LE MAILLAGE.
C |    NPTFR       | -->| NOMBRE DE NOEUDS SUR LA FRONTIERE.
C |    GLOSEG      | -->| GLOBAL NUMBER OF THE 2 POINTS OF A SEGMENT.
C |    NBOR        | -->| TABLEAU D'INDICES DE NOEUDS SUR LE BORD.
C |    HLIN        | -->| VALEURS DE LA HAUTEUR D'EAU A L'ETAPE N+1.
C |                |    | AVEC INTERPOLATION LINEAIRE EN TEMPS
C |                |    | ENTRE HN ET H
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : CVDFTR
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_TVF => TVF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NSEG,NPOIN,NPTFR,KDIR,KDDL
      INTEGER, INTENT(IN)             :: SIZGLO,OPTSOU,IOPT2
      INTEGER, INTENT(IN)             :: GLOSEG(SIZGLO,2)
      INTEGER, INTENT(IN)             :: NBOR(NPTFR),LIMTRA(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: DT,SURNIT
      DOUBLE PRECISION, INTENT(INOUT) :: FLBORTRA(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: F(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: FXBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: FC(NPOIN),H(NPOIN),HLIN(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: SMH(NPOIN),UNSV2D(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: FSCEXP(NPOIN),FN(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: FBOR(NPTFR),FXBORPAR(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: FXMAT(NSEG),FXMATPAR(NSEG)
      LOGICAL, INTENT(IN)             :: YASMH
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: T7,SF
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,N
C
C-----------------------------------------------------------------------
C
      IF(IOPT2.EQ.0) THEN
C       CHAMP CONVECTEUR CONSERVATIF
        DO I = 1,NPOIN
          F(I) = FC(I)
        ENDDO
      ELSEIF(IOPT2.EQ.1) THEN
C       CHAMP CONVECTEUR NON CONSERVATIF
        DO I = 1,NPOIN
          F(I) = FC(I)*MAX(H(I),1.D-8)/MAX(HLIN(I),1.D-8)
        ENDDO
      ELSE
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'TVF : OPTION INCONNUE'
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'TVF: UNKNOWN OPTION'
        ENDIF
        CALL PLANTE(1)
        STOP      
      ENDIF
C
      IF(NCSIZE.GT.1) THEN
C       THE CONTRIBUTION OF FLUXES IS BUILT APART FOR
C       PRELIMINARY PARALLEL ASSEMBLING BEFORE ADDING ON F
        DO I = 1,NPOIN
          T7%R(I) = 0.D0
        ENDDO
        DO I = 1,NSEG
          IF(FXMATPAR(I).LT.0.D0) THEN
            T7%R(GLOSEG(I,1)) = T7%R(GLOSEG(I,1))
     *      - DT/HLIN(GLOSEG(I,1))*UNSV2D(GLOSEG(I,1))
     *      *FXMAT(I)*(FC(GLOSEG(I,2))-FC(GLOSEG(I,1)))
          ELSEIF(FXMATPAR(I).GT.0.D0) THEN
            T7%R(GLOSEG(I,2)) = T7%R(GLOSEG(I,2))
     *      + DT/HLIN(GLOSEG(I,2))*UNSV2D(GLOSEG(I,2))
     *      *FXMAT(I)*(FC(GLOSEG(I,1))-FC(GLOSEG(I,2)))
          ENDIF
        ENDDO
        CALL PARCOM(T7,2,MESH)
        DO I = 1,NPOIN
          F(I) = F(I)+T7%R(I)
        ENDDO        
      ELSE
        DO I = 1,NSEG
          IF(FXMATPAR(I).LT.0.D0) THEN
            F(GLOSEG(I,1)) = F(GLOSEG(I,1))
     *      - DT/HLIN(GLOSEG(I,1))*UNSV2D(GLOSEG(I,1))
     *      *FXMAT(I)*(FC(GLOSEG(I,2))-FC(GLOSEG(I,1)))
          ELSEIF(FXMATPAR(I).GT.0.D0) THEN
            F(GLOSEG(I,2)) = F(GLOSEG(I,2))
     *      + DT/HLIN(GLOSEG(I,2))*UNSV2D(GLOSEG(I,2))
     *      *FXMAT(I)*(FC(GLOSEG(I,1))-FC(GLOSEG(I,2)))
          ENDIF
        ENDDO
      ENDIF
C
C     SOURCE TERMS
C
      IF(YASMH) THEN
        IF(OPTSOU.EQ.1) THEN
          DO I=1,NPOIN
            F(I)=F(I)+DT/HLIN(I)*SMH(I)*(FSCEXP(I)+FN(I)-FC(I))
          ENDDO
        ELSEIF(OPTSOU.EQ.2) THEN
          DO I=1,NPOIN
         F(I)=F(I)+DT/HLIN(I)*UNSV2D(I)*SMH(I)*(FSCEXP(I)+FN(I)-FC(I))
          ENDDO
        ENDIF
      ENDIF
C
C ON THE DIRICHLET BOUNDARIES, FLUX TERMS TAKEN INTO ACCOUNT
C ON OTHERS, FBOR IS TAKEN AS FN, SO NO CONTRIBUTION
C
      DO I=1,NPTFR
        IF(LIMTRA(I).EQ.KDIR) THEN
          N=NBOR(I)
          F(N)=F(N)-DT/HLIN(N)*UNSV2D(N)*FXBORPAR(N)*(FBOR(I)-FC(N))
        ELSEIF(LIMTRA(I).EQ.KDDL) THEN
          N=NBOR(I)
          FLBORTRA(I)=FLBORTRA(I)+FXBOR(I)*FC(N)*SURNIT
        ENDIF
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
