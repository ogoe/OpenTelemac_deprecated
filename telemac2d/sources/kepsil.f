C                       *****************
                        SUBROUTINE KEPSIL
C                       *****************
C
     *(AK,EP,AKTILD,EPTILD,AKN,EPN,VISC,CF,U,V,HN,UCONV,VCONV,
     * KBOR,EBOR,LIMKEP,IELMK,IELME,
     * SMK,SME,TM1,MAK,MAE,CM2,TE1,TE2,NPTFR,DT,
     * MESH,T1,T2,T3,TB,CMU,C1,C2,SIGMAK,SIGMAE,ESTAR,SCHMIT,
     * KMIN,KMAX,EMIN,EMAX,
     * INFOKE,KDIR,MSK,MASKEL,MASKPT,S,SLVK,SLVEP,ICONV,OPTSUP)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.5    27/11/92    J-M HERVOUET (LNH) 30 87 80 18
C                            30/05/94    L. VAN HAREN (LNH) 30 87 84 14
C***********************************************************************
C
C  FONCTION  : ETAPE DE DIFFUSION-TERMES SOURCES DU MODELE K-EPSILON
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      AK        |<-- | ENERGIE TURBULENTE AU TEMPS T(N+1)
C |      EP        |<-- | DISSIPATION TURBULENTE AU TEMPS T(N+1)
C |      AKTILD    | -->| ENERGIE TURBULENTE APRES CONVECTION
C |      EPTILD    | -->| DISSIPATION TURBULENTE APRES CONVECTION
C |      AKN       | -->| ENERGIE TURBULENTE AU TEMPS T(N)
C |      EPN       | -->| DISSIPATION TURBULENTE AU TEMPS T(N)
C |     VISC       | -->| DIFFUSION TURBULENTE
C |      CF        | -->| COEFFICIENT DE FROTTEMENT POUR K-EPSILON
C |     U , V      | -->| COMPOSANTES DE LA VITESSE
C |      HN        | -->| HAUTEUR AU TEMPS N
C |   KBOR,EBOR    | -->| K ET EPSILON IMPOSES AU BORD
C |    LIMKEP      | -->| CONDITIONS AUX LIMITES SUR K ET EPSILON
C |      SMK       | -- | SECOND MEMBRE DU SYSTEME A RESOUDRE POUR K
C |      SME       | -- | SECOND MEMBRE DU SYSTEME A RESOUDRE POUR E
C |      TM1       | -- | MATRICE DE DIFFUSION
C |      MAK       | -- | MATRICE DU SYSTEME A RESOUDRE POUR K
C |      MAE       | -- | MATRICE DU SYSTEME A RESOUDRE POUR E
C |      CM2       | -- | MATRIX
C |     NPTFR      | -->| NOMBRE DE POINTS FRONTIERES
C |       DT       | -->| PAS DE TEMPS
C |   T1,2,3,4     | -- | TABLEAUX DE TRAVAIL
C |    KARMAN      | -->| CONSTANTE DE KARMAN
C |     CMU        | -->| CONSTANTE DU MODELE K-EPSILON
C |    C1,C2       | -->| CONSTANTES DU MODELE K-EPSILON
C |    SIGMAK      | -->| CONSTANTE DU MODELE K-EPSILON
C |    SIGMAE      | -->| CONSTANTE DU MODELE K-EPSILON
C |    ESTAR       | -->| CONSTANTE DU MODELE K-EPSILON
C |    SCHMIT      | -->| CONSTANTE DU MODELE K-EPSILON
C |   KMIN,KMAX    | -->| K MINIMUM ET MAXIMUM EN CAS DE CLIPPING
C |   EMIN,EMAX    | -->| EPSILON MINIMUM ET MAXIMUM EN CAS DE CLIPPING
C |     INFOKE     | -->| LOGIQUE INDIQUANT SI LES INFORMATIONS SUR LE
C |                |    | SOLVEUR SONT A RESTITUER
C |     KNEU       | -->| CONDITION A LA LIMITE DE TYPE NEUMANN
C |     KDIR       | -->| CONDITION A LA LIMITE DE TYPE DIRICHLET
C |     KDDL       | -->| CONDITION A LA LIMITE DE TYPE DEGRE DE LIBERTE
C |     MSK        | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |     MASKEL     | -->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |     MASKPT     | -->|  MASQUES PAR POINTS.
C |     S          | -->|  STRUCTURE BIDON
C |     SLVK       | -->|  STRUCTURE WITH SOLVER OPTIONS FOR K
C |     SLVEP      | -->|  STRUCTURE WITH SOLVER OPTIONS FOR E
C |     ICONV      | -->|  TYPE OF ADVECTION ON K AND EPSILON
C |                |    |  1 : CHARACTERISTICS
C |                |    |  2 : SUPG
C |     OPTSUP     | -->|  SUPG OPTION
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMMES APPELES : CLIP , GRADF , LUMPIN , MATDIF , MATMAS ,
C                           MATMAT , MATVEC , OV , PRIDIR , SOLV01
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
      TYPE(SLVCFG), INTENT(INOUT)  :: SLVK,SLVEP
      INTEGER, INTENT(IN)          :: ICONV,NPTFR,KDIR,LIMKEP(NPTFR,2)
      INTEGER, INTENT(IN)          :: OPTSUP,IELMK,IELME
      LOGICAL, INTENT(IN)          :: INFOKE,MSK
      DOUBLE PRECISION, INTENT(IN) :: KMIN,KMAX,EMIN,EMAX,SCHMIT
      DOUBLE PRECISION, INTENT(IN) :: CMU,C1,C2,SIGMAK,SIGMAE,ESTAR
      DOUBLE PRECISION, INTENT(IN) :: DT
C     MATRIX STRUCTURES
      TYPE(BIEF_OBJ), INTENT(INOUT) :: TM1,MAK,MAE,CM2
C     VECTOR STRUCTURES
      TYPE(BIEF_OBJ), INTENT(IN)    :: UCONV,VCONV,AKN,EPN,AKTILD,EPTILD
      TYPE(BIEF_OBJ), INTENT(IN)    :: HN,VISC,U,V,MASKEL,S,MASKPT,CF
      TYPE(BIEF_OBJ), INTENT(IN)    :: KBOR,EBOR
      TYPE(BIEF_OBJ), INTENT(INOUT) :: T1,T2,T3,AK,EP,SMK,SME,TE1,TE2
C     MESH STRUCTURE
      TYPE(BIEF_MESH) :: MESH
C     BLOCK STRUCTURE
      TYPE(BIEF_OBJ) :: TB
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C      
      DOUBLE PRECISION C,SL1,CEPS,USTAR,AGGLOK,AGGLOE,TETAK
C
      INTEGER N    
C
C-----------------------------------------------------------------------
C
      INTRINSIC SQRT,MAX
C
C-----------------------------------------------------------------------
C
C   CALCULS DE LA MATRICE DE MASSE ET DE LA MATRICE DE DIFFUSION       *
C
      SL1 = 1.D0/DT
C
C     -----------------------------
C     CALCUL DES MATRICES DE MASSE
C     -----------------------------
C
      CALL MATRIX(MAK,'M=N     ','MATMAS          ',IELMK,IELMK,
     *            SL1,S,S,S,S,S,S,MESH,MSK,MASKEL)
      CALL MATRIX(MAE,'M=N     ','MATMAS          ',IELME,IELME,
     *            SL1,S,S,S,S,S,S,MESH,MSK,MASKEL)
C
C     ESSAI DE MASS-LUMPING
C
      AGGLOK = 1.D0
      AGGLOE = 1.D0
      IF(AGGLOK.GT.0.001D0) THEN
        CALL LUMP(T1,MAK,MESH,AGGLOK)
        CALL OM( 'M=CN    ' , MAK , MAK , S  , 1.D0-AGGLOK , MESH )
        CALL OM( 'M=M+D   ' , MAK , MAK , T1 , C           , MESH )
      ENDIF
      IF(AGGLOE.GT.0.001D0) THEN
        CALL LUMP(T1,MAE,MESH,AGGLOE)
        CALL OM( 'M=CN    ' , MAE , MAE , S  , 1.D0-AGGLOE , MESH )
        CALL OM( 'M=M+D   ' , MAE , MAE , T1 , C           , MESH )
      ENDIF
C
C     --------------------------------------------------
C     AGGLOMEREE DE LA MATRICE DE MASSE: MISE DANS T3
C     --------------------------------------------------
C
      CALL LUMP(T3,MAK,MESH,DT)
C
C     ---------------------
C     MATRICE DE DIFFFUSION
C     ---------------------
C
      CALL MATRIX(TM1,'M=N     ','MATDIF          ',IELMK,IELMK,
     *            1.D0,S,S,S,VISC,VISC,VISC,MESH,MSK,MASKEL)
C
C***********************************************************************
C
C     TERMES SOURCES EXPLICITES : T1 POUR K , T2 POUR EPSILON    *
C                                                                      *
C     TERME EXPLICITE SUR K :                                          *
C                                            3                         *
C                               N           U                          *
C                              K             *
C                              --   +  C  * --  +  PROD
C                              DT       K   H
C
C
C     TERME EXPLICITE SUR EPSILON:
C
C                                            4
C                                N          U              N
C                              EP            *           EP
C                              --   +  C  * --  +  C   * -- * PROD
C                              DT       E    2      E1    N
C                                           H            K
C
C
C                     2        2           2
C                  DU       DV     DU   DV           N
C      PROD = ( 2*(--) + 2*(--) + (-- + --)  ) * VISC
C                  DX       DY     DY   DX
C
C
C                           N
C                         EP
C      LE TERME  +  C1  * -- * PROD   EST MIS SOUS LA FORME :
C                          N
C                         K
C
C                               N
C                   C1 * CMU * K * PROD / VISC
C
C***********************************************************************
C
C     --------------------------------
C     PRISE EN COMPTE DE LA CONVECTION
C     --------------------------------
C
      IF(ICONV.EQ.1) THEN
C
        CALL MATVEC('X=AY    ',SMK,MAK,AKTILD,C,MESH)
        CALL MATVEC('X=AY    ',SME,MAE,EPTILD,C,MESH)
C
      ELSEIF(ICONV.EQ.2) THEN
C
        CALL MATVEC('X=AY    ',SMK,MAK,AKN,C,MESH)
        CALL MATVEC('X=AY    ',SME,MAE,EPN,C,MESH)
C       TERME DE CONVECTION CENTRE SEMI-IMPLICITE : MATRICE
        CALL MATRIX(CM2,'M=N     ','MATVGR          ',IELMK,IELMK,
     *              1.D0,S,S,S,UCONV,VCONV,VCONV,MESH,MSK,MASKEL)
C       SUPG CONTRIBUTION
        IF(OPTSUP.EQ.1) THEN
C         SUPG CLASSIQUE
          CALL KSUPG(TE1,TE2,1.D0,UCONV,VCONV,MESH)
          CALL MATRIX(CM2,'M=M+N   ','MASUPG          ',IELMK,IELMK,
     *                1.D0,TE1,TE2,S,UCONV,VCONV,VCONV,
     *                MESH,MSK,MASKEL)
        ELSEIF(OPTSUP.EQ.2) THEN
C         SUPG MODIFIE
          CALL MATRIX(CM2,'M=M+N   ','MAUGUG          ',IELMK,IELMK,
     *                0.5D0*DT,S,S,S,UCONV,VCONV,VCONV,
     *                MESH,MSK,MASKEL)
        ENDIF
C       END OF SUPG CONTRIBUTION
C       EXPLICIT RIGHT HAND SIDES
        TETAK=0.6
        CALL MATVEC( 'X=X+CAY ',SMK,CM2,AKN,TETAK-1.D0,MESH)
        CALL MATVEC( 'X=X+CAY ',SME,CM2,EPN,TETAK-1.D0,MESH)
C       ADDING SUPG MATRIX TO MAK AND MAE
        CALL OM( 'M=X(M)  ' , MAK , MAK , S , C , MESH )
        CALL OM( 'M=M+CN  ' , MAK , CM2 , S , TETAK , MESH )
        CALL OM( 'M=X(M)  ' , MAE , MAE , S , C , MESH )
        CALL OM( 'M=M+CN  ' , MAE , CM2 , S , TETAK , MESH )
C
      ELSE
C
        IF(LNG.EQ.1) WRITE(LU,100) ICONV
        IF(LNG.EQ.2) WRITE(LU,101) ICONV
100     FORMAT(1X,'KEPSIL : FORME DE LA CONVECTION INCONNUE : ',1I4)
101     FORMAT(1X,'KEPSIL: UNKNOWN TYPE OF ADVECTION:',1I4)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C     ------------------------------------------------
C     TERME DE CREATION PAR FROTTEMENT SUR LE FOND (P)
C     ------------------------------------------------
C
C     ATTENTION IL MANQUE 1.D0/(0.5D0*CF)**0.75 AU VRAI CEPS
C     (PRIS EN COMPTE APRES)
C
      CEPS = C2 * SQRT(CMU) / SQRT( ESTAR*SCHMIT )
C
      CALL CPSTVC(SMK,T1)
      CALL CPSTVC(SMK,T2)
      DO N=1,T1%DIM1
         USTAR = SQRT(0.5D0*CF%R(N)*(U%R(N)**2+V%R(N)**2))
C        T1 : TERME Pkv DE L'EQUATION DE K
         T1%R(N)=USTAR**3/MAX(SQRT(0.5D0*CF%R(N))*HN%R(N),1.D-6)
C        ON LIMITE LA CROISSANCE DUE AU FROTTEMENT
C        SUR LE FOND DE K A 50% PAR PAS DE TEMPS
         T1%R(N) = MIN (T1%R(N) , EP%R(N) + 0.5D0*AK%R(N)/DT )
C        T2 : TERME Pev DE L'EQUATION DE EPSILON
         T2%R(N) = CEPS * USTAR**4 /
     *           MAX(((0.5D0*CF%R(N))**0.75D0)*HN%R(N)**2,1.D-6)
C        ON LIMITE LA CROISSANCE DE E A 50% PAR PAS DE TEMPS
         T2%R(N) = MIN(T2%R(N),C2*EP%R(N)**2/MAX(AK%R(N),KMIN)
     *           +0.5D0*EP%R(N)/DT )
      ENDDO
C
      CALL OS( 'X=XY    ' , T1  , T3 , T3 , C )
      CALL OS( 'X=X+Y   ' , SMK , T1 , T1 , C )
      CALL OS( 'X=XY    ' , T2  , T3 , T3 , C )
      CALL OS( 'X=X+Y   ' , SME , T2 , T2 , C )
C
C     -----------------------------------
C     TERMES DE CREATION PAR CISAILLEMENT
C     -----------------------------------
C
      CALL VECTOR(T1,'=','PRODF           ',IELMK,
     *            1.D0,S,S,S,U,V,S,MESH,MSK,MASKEL)
      CALL OS( 'X=XY    ' , T1  , VISC , VISC , C )
C
C ESSAI JMH : LIMITATION DE LA TURBULENCE SUR LES PETITS FONDS ( < 2CM )
C
      DO N=1,SMK%DIM1
        IF(HN%R(N).LT.0.02D0) T1%R(N)=0.D0
      ENDDO
C
C FIN ESSAI
C
      CALL OS( 'X=X+Y   ' , SMK , T1   , T1   , C )
C
      CALL VECTOR(T2,'=','PRODF           ',IELMK,
     *            CMU*C1,S,S,S,U,V,S,MESH,MSK,MASKEL)
      CALL OS( 'X=XY    ' , T2  , AK , AK , C )
      CALL OS( 'X=X+Y   ' , SME , T2 , T2 , C )
C
C ESSAI JMH : LIMITATION DE LA TURBULENCE SUR LES PETITS FONDS ( < 2CM )
C
      DO N=1,SMK%DIM1
        SMK%R(N) = SMK%R(N) * (MIN(HN%R(N),0.02D0)/0.02D0)**2
      ENDDO
C
C FIN ESSAI
C
C***********************************************************************
C     TERMES SOURCES IMPLICITES : T1 POUR K , T2 POUR EPSILON          *
C                                                                      *
C     TERME IMPLICITE SUR K :           +      EP(N)/K(N) * K (N+1)    *
C     TERME IMPLICITE SUR EPSILON:      + C2 * EP(N)/K(N) * EP(N+1)    *
C***********************************************************************
C
      CALL OS( 'X=Y/Z   ',T1,EP,AK,C  ,IOPT=2,INFINI=0.D0,ZERO=KMIN)
      CALL OS( 'X=CY    ',T2,T1,T1,C2 )
C
C     ---------------------------------------------------
C     INTEGRATION DE CES TERMES SOURCES DANS LES MATRICES
C     ---------------------------------------------------
C
      CALL OS( 'X=XY    ' , T1 , T3 , T3 , C )
      CALL OS( 'X=XY    ' , T2 , T3 , T3 , C )
C
C     -------------------------------------------
C     AJOUT A LA DIAGONALE DE LA MATRICE DE MASSE
C     -------------------------------------------
C
      CALL OM( 'M=M+D   ' , MAK , MAK , T1 , C , MESH )
      CALL OM( 'M=M+D   ' , MAE , MAE , T2 , C , MESH )
C
C***********************************************************************
C
C     COMBINAISON DES MATRICES DE MASSE ET DE DIFFUSION                *
C                                                                      *
C     MAK = MAK + TM1/SIGMAK                                           *
C     MAE = MAE + TM1/SIGMAE                                           *
C
C***********************************************************************
C
      CALL OM( 'M=M+CN  ' , MAK , TM1 , S , 1.D0/SIGMAK , MESH )
      CALL OM( 'M=M+CN  ' , MAE , TM1 , S , 1.D0/SIGMAE , MESH )
C
C***********************************************************************
C     DIRICHLET TYPE BOUNDARY CONDITIONS
C***********************************************************************
C
      CALL DIRICH(AK,MAK,SMK,KBOR,LIMKEP(1,1),TB,MESH,KDIR,MSK,MASKPT)
      CALL DIRICH(EP,MAE,SME,EBOR,LIMKEP(1,2),TB,MESH,KDIR,MSK,MASKPT)
C
C***********************************************************************
C     RESOLUTION DES DEUX SYSTEMES OBTENUS
C***********************************************************************
C
      CALL SOLVE(AK,MAK,SMK,TB,SLVK ,INFOKE,MESH,TM1)
      CALL SOLVE(EP,MAE,SME,TB,SLVEP,INFOKE,MESH,TM1)
C
C***********************************************************************
C     CLIPPING DES PETITES VALEURS                                     *
C***********************************************************************
C
      CALL CLIP(AK,0.D0,.TRUE.,KMAX,.FALSE.,0)
      CALL CLIP(EP,EMIN,.TRUE.,EMAX,.FALSE.,0)
C
C***********************************************************************
C
      RETURN
      END 
