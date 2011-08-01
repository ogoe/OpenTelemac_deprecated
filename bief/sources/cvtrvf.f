C                       *****************
                        SUBROUTINE CVTRVF
C                       *****************
C
     *(F,FN,FSCEXP,DIFT,CONV,H,HN,HPROP,UCONV,VCONV,DM1,ZCONV,SOLSYS,
     * VISC,VISC_S,SM,SMH,YASMH,SMI,YASMI,FBOR,MASKTR,MESH,
     * T1,T2,T3,T4,T5,T6,T7,T8,HNT,HT,AGGLOH,TE1,DT,ENTET,BILAN,
     * OPDTRA,MSK,MASKEL,S,MASSOU,OPTSOU,LIMTRA,KDIR,KDDL,NPTFR,FLBOR,
     * YAFLBOR,V2DPAR,UNSV2D,IOPT,FLBORTRA,MASKPT)
C
C***********************************************************************
C  BIEF VERSION 6.0   09/10/09   CHI-TUAN PHAM     (LNHE) 01 30 87 ?? ??
C                                           
C***********************************************************************
C
C  FONCTION : CONVECTEUR EN VOLUMES FINIS, UPWIND, EXPLICITE
C
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   F            |<-- |  VALEURS A L' ETAPE N+1.
C |   FN           | -->|  VALEURS A L' ETAPE N.
C |   FSCEXP       | -->|  PARTIE EXPLICITE DU TERME SOURCE
C |                |    |  EGALE A ZERO PARTOUT SAUF POUR LES POINTS
C |                |    |  SOURCES OU IL Y A FSCE - (1-TETAT) FN
C |                |    |  VOIR DIFSOU
C |   DIFT         | -->|  LOGIQUE INDIQUANT S'IL Y A DIFFUSION DE F
C |   CONV         | -->|  LOGIQUE INDIQUANT S'IL Y A CONVECTION DE F
C |   H , HN       | -->|  VALEURS DE LA HAUTEUR A L' ETAPE N+1 ET N
C |   HPROP        | -->|  HAUTEUR DE PROPAGATION (FAITE DANS CVDFTR).
C |   U,V,UN,VN    | -->|  VITESSES A T(N+1) ET T(N)
C |   UCONV,VCONV  | -->|  TABLEAUX DE TRAVAIL.
C |   TETAU        | -->|  IMPLICITATION SUR U
C |   VISC         | -->|  COEFFICIENTS DE VISCOSITE SUIVANT X,Y ET Z .
C |                |    |  SI P0 : VISCOSITE DONNEE PAR ELEMENT
C |                |    |  SINON : VISCOSITE DONNEE PAR POINT
C |   SM           | -->|  TERMES SOURCES .
C |   SMH          | -->|  TERME SOURCE DE L'EQUATION DE CONTINUITE
C |   YASMH        | -->|  LOGIQUE INDIQUANT DE PRENDRE EN COMPTE SMH
C |   FBOR         | -->|  CONDITIONS DE DIRICHLET SUR F.
C |   MASKTR(1,1)  | -->|  MASQUE VALANT 1. POUR LES SEGMENTS DIRICHLET
C |   MASKTR(1,2)  | -->|  MASQUE VALANT 1. POUR LES SEGMENTS DDL
C |   MASKTR(1,3)  | -->|  MASQUE VALANT 1. POUR LES SEGMENTS NEUMANN
C |                |    |  (ET ZERO SINON)
C |   MESH         | -->|  BLOC DES ENTIERS DU MAILLAGE.
C |   T1......T7   |<-->|  TABLEAUX DE TRAVAIL (T1 PAS UTILISE)
C |   HNT,HT       |<-- |  TABLEAUX DE TRAVAIL (HAUTEURS MODIFIEES POUR
C |                |    |  TENIR COMPTE DU MASS-LUMPING)
C |   AGGLOH       | -->|  MASS-LUMPING UTILISE DANS L'EQUATION DE CONTINUITE
C |   TE1          |<-->|  TABLEAU DE TRAVAIL SUR LES ELEMENTS
C |   KNEU         | -->|  CONVENTION POUR LES CONDITIONS DE NEUMANN
C |   KDDL         | -->|  CONVENTION POUR LES DEGRES DE LIBERTE
C |   DT           | -->|  PAS DE TEMPS
C |   ENTET        | -->|  LOGIQUE INDIQUANT SI ON IMPRIME DES INFOS
C |                |    |  SUR LE BILAN DE MASSE DE TRACEUR
C |   BILAN        | -->|  LOGIQUE INDIQUANT SI ON DOIT FAIRE UN BILAN
C |                |    |  DE MASSE. DANS CE CAS IL FAUT RETOURNER LA
C |                |    |  VALEUR DE L'APPORT DES TERMES SOURCES.
C |   OPDTRA       | -->|  MOT-CLE : OPTION POUR LA DIFFUSION DU TRACEUR
C |   MSK          | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |   MASKEL       | -->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |   S            | -->|  STRUCTURE BIDON
C |   MASSOU       | -->|  MASSE DE TRACEUR AJOUTEE PAR TERME SOURCE
C |                |    |  VOIR DIFSOU
C |   OPTSOU       | -->|  OPTION DE TRAITEMENT DES TERMES SOURCES.
C |                |    |  1 : NORMAL
C |                |    |  2 : DIRAC
C |   IOPT         |    |  OPTIONS DE CALCUL
C |                |    |  CHIFFRE DES DIZAINES (IOPT2):
C |                |    |  0 : UCONV RESPECTE L'EQUATION DE CONTINUITE
C |                |    |  1 : UCONV NE RESPECTE PAS LA CONTINUITE
C |                |    |  CHIFFRE DES UNITES (IOPT1):
C |                |    |  0 : CONSTANTE PAR ELEMENT NULLE
C |                |    |  1 : CONSTANTE DE CHI-TUAN PHAM
C |                |    |  2 : SCHEMA N
C |                |    |  3 : SCHEMA PSI
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : OV , MATMAS , MATDIF , MATBDF , MATVEC
C                      MATMAT , INTBDF , SOLV01
C
C***********************************************************************
C
C   ATTENTION : POUR LES ELEMENTS DE BORD OU IL N'Y A PAS DE FROTTEMENT
C               AFBOR ET BFBOR DOIVENT ETRE NULS.
C
C             : ATTENTION A LA DISCRETISATION DE VISC
C
C**********************************************************************
C
      USE BIEF, EX_CVTRVF => CVTRVF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: OPDTRA,OPTSOU,KDIR,NPTFR,SOLSYS
      INTEGER, INTENT(IN)             :: LIMTRA(NPTFR),KDDL,IOPT
      DOUBLE PRECISION, INTENT(IN)    :: DT,AGGLOH
      DOUBLE PRECISION, INTENT(INOUT) :: MASSOU
      LOGICAL, INTENT(IN)             :: BILAN,CONV,YASMH,YAFLBOR
      LOGICAL, INTENT(IN)             :: DIFT,MSK,ENTET,YASMI
      TYPE(BIEF_OBJ), INTENT(IN)      :: MASKEL,H,HN,DM1,ZCONV,MASKPT
      TYPE(BIEF_OBJ), INTENT(IN)      :: V2DPAR,UNSV2D,HPROP
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: F,SM,HNT,HT
      TYPE(BIEF_OBJ), INTENT(IN)      :: FBOR,UCONV,VCONV,FN,SMI,SMH
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: TE1,FLBORTRA
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: T1,T2,T3,T4,T5,T6,T7,T8
      TYPE(BIEF_OBJ), INTENT(IN)      :: FSCEXP,S,MASKTR,FLBOR
      TYPE(BIEF_OBJ), INTENT(IN)      :: VISC_S,VISC 
      TYPE(BIEF_MESH) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELMF,I,IOPT1,IOPT2
      LOGICAL YACSTE
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION MASSET,MASSETN,TSOU,DTMAX,DT_REMAIN,DDT,TDT
      DOUBLE PRECISION FXT2
      CHARACTER(LEN=16) FORMUL
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: SAVE_HT,SAVE_HNT
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: FXMAT,FXMATPAR
C
      DOUBLE PRECISION P_DMIN,P_DSUM
      EXTERNAL         P_DMIN,P_DSUM
C
      INTEGER NITMAX,NIT
      DATA NITMAX/200/
C
C-----------------------------------------------------------------------
C
      SAVE_HT =>HT%R
      SAVE_HNT=>HNT%R
      FXMAT=>MESH%MSEG%X%R(1:MESH%NSEG)
C     IN PARALLELISM, ASSEMBLED AND NON ASSEMBLED VERSION DIFFERENT
      IF(NCSIZE.GT.1) THEN
        FXMATPAR=>MESH%MSEG%X%R(MESH%NSEG+1:2*MESH%NSEG)
      ELSE
        FXMATPAR=>MESH%MSEG%X%R(1:MESH%NSEG)
      ENDIF
C
C-----------------------------------------------------------------------
C
C     EXTRACTING OPTIONS
C
      IOPT2=IOPT/10
      IOPT1=IOPT-10*IOPT2
C
C-----------------------------------------------------------------------
C
!     IELMF = F%ELM
!     FORCE A LINEAIRE
      IELMF=11
C
C     TAKING INTO ACCOUNT MASS-LUMPING IN THE CONTINUITY EQUATION
C
      IF(ABS(1.D0-AGGLOH).GT.1.D-8) THEN
        CALL VECTOR(HT ,'=','MASVEC          ',IELMF,
     *              1.D0-AGGLOH,H ,S,S,S,S,S,MESH,MSK,MASKEL)
        CALL VECTOR(HNT,'=','MASVEC          ',IELMF,
     *              1.D0-AGGLOH,HN,S,S,S,S,S,MESH,MSK,MASKEL)
        IF(NCSIZE.GT.1) THEN
          CALL PARCOM(HT ,2,MESH)
          CALL PARCOM(HNT,2,MESH)
        ENDIF
        CALL OS('X=YZ    ',X=HT ,Y=HT ,Z=UNSV2D)
        CALL OS('X=YZ    ',X=HNT,Y=HNT,Z=UNSV2D)
        CALL OS('X=X+CY  ',X=HT ,Y=H  ,C=AGGLOH)
        CALL OS('X=X+CY  ',X=HNT,Y=HN ,C=AGGLOH)
      ELSE
!       CALL OS('X=Y     ',X=HT ,Y=H )
!       CALL OS('X=Y     ',X=HNT,Y=HN)
        HT%R =>H%R
        HNT%R=>HN%R
      ENDIF
C
C     INITIALISATION DU FLUX DE TRACEUR AUX FRONTIERES
C
      DO I=1,MESH%NPTFR
        IF(LIMTRA(I).EQ.KDIR) THEN
!         FLBOR IS NOT ASSEMBLED IN PARALLEL
          FLBORTRA%R(I)=FLBOR%R(I)*FBOR%R(I)
        ELSE
C         FOR KDDL, WILL BE DONE IN TVF
          FLBORTRA%R(I)=0.D0
        ENDIF
      ENDDO
C
C     COMPUTATION OF THE FLUXES PHIij = FXMAT
C
      FORMUL='HUGRADP         '
      IF(SOLSYS.EQ.2) FORMUL(8:8)='2'
      CALL VECTOR(T2,'=',FORMUL,IELMF,-1.D0,
     *            HPROP,DM1,ZCONV,UCONV,VCONV,VCONV,MESH,MSK,MASKEL)
C                 T2 AS HUGRADP IS NOT USED AS AN ASSEMBLED VECTOR
C                 BUT TO GET THE NON ASSEMBLED FORM MESH%W
      NIT=0
      DT_REMAIN=DT
      TDT=0.D0
      CALL CPSTVC(HN,T7)
      CALL CPSTVC(H ,T5)
      CALL CPSTVC(H ,T4)
      CALL CPSTVC(F,T8)
C
C     T4 WILL BE F PROGRESSIVELY UPDATED
C     T5 WILL BE THE DEPTH AT THE END OF THE SUB TIME STEP
C     (INITIALISED HERE FOR CALLING CFLVF)
C
      DO I=1,HN%DIM1
        T4%R(I)=FN%R(I)
        T5%R(I)=HNT%R(I)
      ENDDO
C
C     T1 WILL BE THE DEPTH ACCORDING TO THE CONTINUITY EQUATION
C
      IF(IOPT2.EQ.1) THEN
        DO I=1,HN%DIM1
          T1%R(I)=HNT%R(I)
        ENDDO
      ENDIF
C
      IF(.NOT.YAFLBOR) THEN
C       MASK=8 FOR LIQUID BOUNDARIES
        CALL VECTOR(T3,'=','FLUBDF          ',1,1.D0,HPROP,HPROP,HPROP,
     *              UCONV,VCONV,VCONV,MESH,.TRUE.,MASKTR%ADR(8)%P)
      ENDIF
C
100   CONTINUE
      NIT=NIT+1
C
C----------------------------------------
C DIFFERENT OPTIONS TO COMPUTE THE FLUXES
C----------------------------------------
C 
      IF(NIT.EQ.1.OR.IOPT1.EQ.3) THEN
        CALL FLUX_EF_VF(FXMAT,MESH%W%R,MESH%NSEG,MESH%NELEM,
     *                  MESH%ELTSEG%I,MESH%ORISEG%I,
     *                  MESH%IKLE%I,.TRUE.,IOPT1,T4)
C       CANCELLING FLUXES TO AND FROM MASKED POINTS
        IF(MSK) THEN
          CALL FLUX_MASK(FXMAT,MESH%NSEG,
     *                   MESH%GLOSEG%I,MESH%GLOSEG%DIM1,MASKPT%R)
        ENDIF
C       ASSEMBLING THE FLUXES AT INTERFACES IN PARALLEL, THIS IS
C       FOR UPWINDING (STORED IN SECOND DIMENSION OF MESH%MSEG)
        IF(NCSIZE.GT.1) THEN
          CALL OV('X=Y     ',FXMATPAR,FXMAT,FXMAT,0.D0,MESH%NSEG)
          CALL PARCOM2_SEG(FXMATPAR,FXMATPAR,FXMATPAR,
     *                     MESH%NSEG,1,2,1,MESH,1)
        ENDIF
      ENDIF
C
C---------------------------------------------
C DETERMINING THE LARGEST ADMISSIBLE TIME STEP
C---------------------------------------------  
C
C     THIS COULD BE PUT OUTSIDE THE LOOP, BUT T7 USED LATER IN THE LOOP...
C
C     IN CFLVF, T7 WILL BE FLBOR WITH A DIMENSION NPOIN
      CALL OS('X=0     ',X=T7)
      IF(YAFLBOR) THEN
        CALL OSDB('X=Y     ',T7,FLBOR,FLBOR,0.D0,MESH)
      ELSE
        CALL OSDB('X=Y     ',T7,T3,T3,0.D0,MESH)
      ENDIF
      IF(NCSIZE.GT.1) CALL PARCOM(T7,2,MESH)
C
C     MASQUAGE EVENTUEL DE FLBOR
C
      IF(MSK) CALL OS('X=XY    ',X=T7,Y=MASKPT)
C 
C     COMPUTING THE MAXIMUM TIME STEP ENSURING MONOTONICITY
C                                                   
      CALL CFLVF(DTMAX,T5%R,HT%R,FXMAT,FXMATPAR,
C                                   FLBOR%R(NPOIN)
     *           V2DPAR%R,DT_REMAIN,T7%R   ,SMH%R,
     *           YASMH,T8,MESH%NSEG,MESH%NPOIN,MESH%NPTFR,
     *           MESH%GLOSEG%I,MESH%GLOSEG%DIM1,MESH,MSK,MASKPT)
      IF(NCSIZE.GT.1) DTMAX=P_DMIN(DTMAX)
C
      DDT=MIN(DT_REMAIN,DTMAX)
      TDT=TDT+DDT
C
C     T5 VA PRENDRE LES VALEURS SUCCESSIVES DE H
C     AUX FINS DE SOUS-PAS DE TEMPS
C
      DO I=1,HN%DIM1
        T5%R(I)=HNT%R(I)+TDT*(HT%R(I)-HNT%R(I))/DT
      ENDDO
C
C     IN TVF FACTOR HT/HLIN MAY TRIGGER DIVERGENCE ON DRY POINTS
C
      IF(MSK) THEN
        DO I=1,HN%DIM1
          IF(MASKPT%R(I).LT.0.5D0) T5%R(I)=HT%R(I)
        ENDDO
      ENDIF
C
C-----------------
C FINAL RESOLUTION
C-----------------
C
      IF(YAFLBOR) THEN
        CALL TRACVF(F,FN,FSCEXP,HT,HNT,FXMAT,FXMATPAR,V2DPAR,UNSV2D,
     *              DDT,FLBOR,FBOR,SMH,YASMH,T1,T2,T4,T5,T6,T7,T8,
     *              MESH,LIMTRA,KDIR,KDDL,OPTSOU,IOPT2,FLBORTRA,MSK,
     *              NIT,DT,TDT)
      ELSE
        CALL TRACVF(F,FN,FSCEXP,HT,HNT,FXMAT,FXMATPAR,V2DPAR,UNSV2D,
     *              DDT,T3,FBOR,SMH,YASMH,T1,T2,T4,T5,T6,T7,T8,MESH,
     *              LIMTRA,KDIR,KDDL,OPTSOU,IOPT2,FLBORTRA,MSK,
     *              NIT,DT,TDT)
      ENDIF
C
      DO I=1,HN%DIM1
        T4%R(I)=F%R(I)
      ENDDO    
      IF(IOPT2.EQ.1) THEN
        DO I=1,HN%DIM1
          T1%R(I)=T2%R(I)
        ENDDO
      ENDIF
C
      DT_REMAIN=DT_REMAIN-DDT
C
      IF(DT_REMAIN.NE.0.D0.AND.NIT.LT.NITMAX) GO TO 100
C
      IF(NIT.GE.NITMAX) THEN
        IF(LNG.EQ.1) WRITE(LU,900) NIT
        IF(LNG.EQ.2) WRITE(LU,901) NIT
900     FORMAT(1X,'CVTRVF : ',1I6,' SOUS-ITERATIONS DEMANDEES POUR LE'
     *   ,/,1X,   '         SCHEMA VF. DIMINUER LE PAS DE TEMPS')
901     FORMAT(1X,'CVTRVF: ',1I6,' SUB-ITERATIONS REQUIRED FOR THE'
     *   ,/,1X,   '         VF SCHEME. DECREASE THE TIME-STEP')
        CALL PLANTE(1)
        STOP
      ELSEIF(ENTET) THEN
        IF(LNG.EQ.1) WRITE(LU,902) NIT
        IF(LNG.EQ.2) WRITE(LU,903) NIT
902     FORMAT(1X,'CVTRVF (BIEF) : ',1I6,' SOUS-ITERATIONS')
903     FORMAT(1X,'CVTRVF (BIEF): ',1I6,' SUB-ITERATIONS')
      ENDIF
C
C-----------------------------------------------------------------------
C
C     EXPLICIT SOURCE TERM
C
      DO I = 1,MESH%NPOIN
        F%R(I) = F%R(I)+DT*SM%R(I)
      ENDDO
C
C     IMPLICIT SOURCE TERM
C
      IF(YASMI) THEN
        DO I = 1,MESH%NPOIN
          F%R(I) = F%R(I)/(1.D0-DT*SMI%R(I)/MAX(H%R(I),1.D-4))
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------
C
C     LOCAL BALANCE OF MASS (FOR CHECKING)
C
!     CALL OS('X=Y-Z   ',X=T7,Y=T5,Z=HT)
!     PRINT*,'DIFFERENCE ENTRE H RECALCULE ET H : ',DOTS(T7,T7)
C     VERIFICATION DE L'EQUATION DU TRACEUR
!     CALL CPSTVC(FBOR,T4)
C     T4 : F AT THE BOUNDARIES AS TAKEN FOR THE BOUNDARY FLUXES
!     DO I=1,NPTFR
!       IF(LIMTRA(I).EQ.KDIR) THEN
!         T4%R(I)=FBOR%R(I)
!       ELSE
!         T4%R(I)=FN%R(MESH%NBOR%I(I))
!       ENDIF
!     ENDDO
!     CALL OS('X=YZ    ',X=T6,Y=FN,Z=HNT)
!     CALL OS('X=YZ    ',X=T7,Y=F ,Z=HT )
!     MASSETN=P_DOTS(V2DPAR,T6,MESH)
!     MASSET =P_DOTS(V2DPAR,T7,MESH)
!     FXT2   =P_DOTS(FLBOR,T4,MESH)
!     PRINT*,'MASSE INIT: ',MASSETN,' MASSE FINALE: ',MASSET
!     PRINT*,'FLUX: ',FXT2
!     MASSETN = MASSETN - FXT2*DT 
!     TSOU=0.D0
!     IF(YASMH) THEN
!       IF(OPTSOU.EQ.1) THEN
!         DO I=1,MESH%NPOIN
!           MASSETN=MASSETN
!    *             +DT*V2DPAR%R(I)*SMH%R(I)*(FSCEXP%R(I)+FN%R(I))
!           TSOU=TSOU+DT*V2DPAR%R(I)*SMH%R(I)*(FSCEXP%R(I)+FN%R(I))
!         ENDDO
!       ELSEIF(OPTSOU.EQ.2) THEN
!         DO I=1,MESH%NPOIN
!           MASSETN=MASSETN
!    *             +DT*SMH%R(I)*(FSCEXP%R(I)+FN%R(I))
!           TSOU=TSOU+DT*SMH%R(I)*(FSCEXP%R(I)+FN%R(I))
!         ENDDO
!       ENDIF
!     ENDIF
!     PRINT*,'CREATION PAR SOURCE : ',TSOU     
!     PRINT*,'ERREUR DE MASSE DE TRACEUR VF : ',MASSETN-MASSET
C     VERIFICATION DE L'EQUATION DE CONTINUITE
!     DO I = 1,MESH%NPOIN
!       T5%R(I)=V2DPAR%R(I)*(HT%R(I)-HNT%R(I))
!     ENDDO
!     DO I = 1,MESH%NSEG
!       T5%R(MESH%GLOSEG%I(I)) = 
!    *  T5%R(MESH%GLOSEG%I(I)) + DT*MESH%MSEG%X%R(I)
!       T5%R(MESH%GLOSEG%I(I+MESH%NSEG)) = 
!    *  T5%R(MESH%GLOSEG%I(I+MESH%NSEG)) - DT*MESH%MSEG%X%R(I)               
!     ENDDO
!     DO I = 1,MESH%NPTFR
!       T5%R(MESH%NBOR%I(I))=T5%R(MESH%NBOR%I(I))+DT*FLBOR%R(I)                
!     ENDDO
!     IF(YASMH) THEN
!       IF(OPTSOU.EQ.1) THEN
!         DO I = 1,MESH%NPOIN
!           T5%R(I)=T5%R(I)-DT*V2DPAR%R(I)*SMH%R(I)
!         ENDDO
!       ELSEIF(OPTSOU.EQ.2) THEN
!         DO I = 1,MESH%NPOIN
!           T5%R(I)=T5%R(I)-DT*SMH%R(I)
!         ENDDO
!       ENDIF
!     ENDIF
!     MASSET=0.D0
!     MASSETN = 0.D0
!     DO I = 1,MESH%NPOIN
!       MASSET=MASSET+T5%R(I)
!       MASSETN=MAX(MASSETN,ABS(T5%R(I)))
!     ENDDO
!     PRINT*,'ERREUR DE MASSE GLOBALE : ',MASSET,' LOCALE : ',MASSETN
C
C-----------------------------------------------------------------------
C
C     RETABLISSEMENT DES POINTEURS DE HT ET HNT
C
      HT%R =>SAVE_HT
      HNT%R=>SAVE_HNT
C
C-----------------------------------------------------------------------
C
      RETURN
      END
