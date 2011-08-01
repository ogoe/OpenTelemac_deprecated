C                           **************
                            SUBROUTINE HVF
C                           **************
C
     *(H,HN,FXMAT,UNSV2D,DT,FXBOR,SMH,YASMH,NSEG,NPOIN,NPTFR,GLOSEG,
     * SIZGLO,NBOR,OPTSOU,T7,MESH,MSK)
C
C***********************************************************************
C BIEF VERSION 5.9      09/02/09     CHI-TUAN PHAM (LNHE) 01 30 87 85 93
C
C
C 09/02/2009 JMH : SEQUENCE IF(MSK) : AVOIDING NEGATIVE DEPTHS
C                  
C
C***********************************************************************
C
C  FONCTION  : RECALCULER UNE HAUTEUR INTERMEDIAIRE SI IL Y A DES
C              SOUS-ITERATIONS
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|_______________________________________________
C |    H           |<-- | VALEURS DE LA HAUTEUR D'EAU A L'ETAPE N+1.
C |    HN          | -->| VALEURS DE LA HAUTEUR D'EAU A L'ETAPE N.
C |    FXMAT       | -->| MATRICE DE STOCKAGE DES FLUX.
C |    MAS         | -->| VECTEUR MASS ASSEMBLE LUMPE.
C |    DT          | -->| PAS DE TEMPS.
C |    FXBOR       | -->| FLUX AU BORD (DEFINI SUR TOUT LE DOMAINE
C |                |    | ET ASSEMBLE EN PARALLELE)
C |    SMH         | -->| TERME SOURCE DE L'EQUATION DE CONTINUITE.
C |    YASMH       | -->| IF YES, SMH MUST BE TAKEN INTO ACCOUNT
C |    NSEG        | -->| NOMBRE DE SEGMENTS DANS LE MAILLAGE.
C |    NELEM       | -->| NOMBRE D'ELEMENTS DANS LE MAILLAGE.
C |    NPOIN       | -->| NOMBRE DE NOEUDS DANS LE MAILLAGE.
C |    NPTFR       | -->| NOMBRE DE NOEUDS SUR LA FRONTIERE.
C |    GLOSEG      | -->| GLOBAL NUMBER OF THE 2 POINTS OF A SEGMENT
C |    NBOR        | -->| TABLEAU D'INDICES DE NOEUDS SUR LE BORD.
C |    OPTSOU      | -->| OPTION FOR THE TREATMENT OF SOURCES
C |                |    | 1: NORMAL  2: DIRAC
C |                |    | SEE PROPAG IN TELEMAC-2D
C |    MSK         | -->| MSK : IF YES, MASKING OF DRY ELEMENTS
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : CVDFTR
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_HVF => HVF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NSEG,NPOIN,NPTFR,OPTSOU,SIZGLO
      INTEGER, INTENT(IN)             :: GLOSEG(SIZGLO,2)
      INTEGER, INTENT(IN)             :: NBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: DT
      DOUBLE PRECISION, INTENT(INOUT) :: H(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: HN(NPOIN),UNSV2D(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: FXBOR(NPOIN),SMH(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: FXMAT(NSEG*2)
      LOGICAL, INTENT(IN)             :: YASMH,MSK
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: T7
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,N
C
C-----------------------------------------------------------------------
C
      DO I = 1,NPOIN
        H(I) = HN(I) 
      ENDDO
C
C     TERMES SOURCES
C
      IF(YASMH) THEN
        IF(OPTSOU.EQ.1) THEN
          DO I = 1,NPOIN
            H(I) = H(I) + DT*SMH(I)
          ENDDO
        ELSEIF(OPTSOU.EQ.2) THEN
          DO I = 1,NPOIN
            H(I) = H(I) + DT*UNSV2D(I)*SMH(I)
          ENDDO
        ENDIF
      ENDIF
C
      IF(NCSIZE.GT.1) THEN
        DO I = 1,NPOIN
          T7%R(I) = 0.D0
        ENDDO
        DO I = 1,NSEG
          T7%R(GLOSEG(I,1))=T7%R(GLOSEG(I,1))
     *                     -DT*UNSV2D(GLOSEG(I,1))*FXMAT(I)
          T7%R(GLOSEG(I,2))=T7%R(GLOSEG(I,2))
     *                     +DT*UNSV2D(GLOSEG(I,2))*FXMAT(I)
        ENDDO
        CALL PARCOM(T7,2,MESH)
        DO I = 1,NPOIN
          H(I) = H(I) + T7%R(I)
        ENDDO
      ELSE
        DO I = 1,NSEG
          H(GLOSEG(I,1))=H(GLOSEG(I,1))-DT*UNSV2D(GLOSEG(I,1))*FXMAT(I)
          H(GLOSEG(I,2))=H(GLOSEG(I,2))+DT*UNSV2D(GLOSEG(I,2))*FXMAT(I)
        ENDDO
      ENDIF
C
C     ON THE BOUNDARIES : BOUNDARY FLUX TERMS
C
      DO I=1,NPTFR
        N=NBOR(I)
        H(N) = H(N) - DT*UNSV2D(N)*FXBOR(N)
      ENDDO
C
C-----------------------------------------------------------------------
C
C     WHEN NEGATIVE DEPTHS APPEAR WHILE COMPUTING H, THE PREVIOUS
C     VALUE OF H IS KEPT
C
      IF(MSK) THEN
        DO I = 1,NPOIN
          IF(H(I).LT.0.D0) H(I) = MAX(1.D-2,HN(I))
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
