C
C  ***************************************************
C  DECLARATION OF GLOBAL DATA STRUCTURE IN ARTEMIS 5.1
C  ***************************************************
C
      MODULE DECLARATIONS_ARTEMIS
C



        USE BIEF_DEF
C

C
C     TO MANAGE THE GLOBAL BUILDING OF THE BOUNDARY
        INTEGER NPTFR_TOT,NPOIN_TOT
C     TO MANAGE THE GLOBAL BUILDING OF THE BOUNDARY
        INTEGER, ALLOCATABLE :: KP1BOR_TOT(:),NBOR_TOT(:), LIHBORT(:),
     c       LIDIRT(:)
        DOUBLE PRECISION, ALLOCATABLE ::  MASK1T(:),MASK2T(:),MASK3T(:), 
     c       MASK4T(:)
        DOUBLE PRECISION, ALLOCATABLE ::XT(:),YT(:),CGT(:),KT(:),CTT(:),
     c       RPT(:),ALFAPT(:),HBT(:),TETABT(:),TETAPT(:)



C       NOTE: THIS MODULE IS ORGANISED IN 10 PARTS
C
C       1) VECTORS (WILL BE DECLARED AS BIEF_OBJ STRUCTURES)
C       2) MATRICES (WILL BE DECLARED AS BIEF_OBJ STRUCTURES) 
C       3) BLOCKS (WILL BE DECLARED AS BIEF_OBJ STRUCTURES)
C       4) INTEGERS
C       5) LOGICAL VALUES
C       6) REALS
C       7) STRINGS
C       8) SLVCFG STRUCTURES
C       9) MESH STRUCTURE
C      10) ALIASES    
C
C-----------------------------------------------------------------------
C
C       1) VECTORS
C
C-----------------------------------------------------------------------
C
C       REAL & IMAGINARY PARTS OF WAVE POTENTIAL
C





        TYPE(BIEF_OBJ), TARGET :: PHIR,PHII
C
C       WATER DEPTH AT REST
C
        TYPE(BIEF_OBJ), TARGET :: H
C
C       WAVE NUMBER
C
        TYPE(BIEF_OBJ), TARGET :: K
C
C       PHASE & GROUP CELERITIES
C
        TYPE(BIEF_OBJ), TARGET :: C,CG
C
C       WAVE HEIGHT AND WAVE PHASE (REGULAR MODE)
C
        TYPE(BIEF_OBJ), TARGET :: HHO,PHAS
C        
C       SURFACE WAVE VELOCITY COMPONENTS 
C
        TYPE(BIEF_OBJ), TARGET :: U0,V0
	
C        
C       MEAN COSINUS AND SINUS VALUES OF WAVE DIRECTION 
C
        TYPE(BIEF_OBJ), TARGET :: MCOS,MSIN
C
C       WAVE INCIDENCE (OR DIRECTION)
C
        TYPE(BIEF_OBJ), TARGET :: INCI
C
C       FREE SURFACE AND BOTTOM LEVELS
C
        TYPE(BIEF_OBJ), TARGET :: S,ZF
C
C       FRICTION FACTOR
C
        TYPE(BIEF_OBJ), TARGET :: FW 
C
C       WAVE HEIGHT (RANDOM WAVE)
C
        TYPE(BIEF_OBJ), TARGET :: HALE
C
C       WAVE PERIODS ARRAY (RANDOM MODE)
C
        TYPE(BIEF_OBJ), TARGET :: PALE           
C
C       REFLEXION COEFFICIENTS
C
        TYPE(BIEF_OBJ), TARGET :: RP,TETAP,ALFAP
C
C       INCIDENT WAVE HEIGHT AND DIRECTION AT BOUNDARY
C
        TYPE(BIEF_OBJ), TARGET :: HB,TETAB
C        
C       REAL AND IMAGINARY PARTS OF INCIDENT WAVE AT BOUNDARY
C  
        TYPE(BIEF_OBJ), TARGET :: PHIRB,PHIIB
C
C       BOUNDARY COEFFICIENTS
C
        TYPE(BIEF_OBJ), TARGET :: APHI1B,BPHI1B,CPHI1B,DPHI1B
        TYPE(BIEF_OBJ), TARGET :: APHI2B,BPHI2B,CPHI2B,DPHI2B
        TYPE(BIEF_OBJ), TARGET :: APHI3B,BPHI3B,CPHI3B,DPHI3B
        TYPE(BIEF_OBJ), TARGET :: APHI4B,BPHI4B,CPHI4B,DPHI4B
C
C       W1 WORKING ARRAYS
C
        TYPE(BIEF_OBJ), TARGET :: W1
C
C       INTEGER WORKING ARRAYS
C
        TYPE(BIEF_OBJ), TARGET :: IT1,IT2,IT3
C
C       VOID STRUCTURE
C
        TYPE(BIEF_OBJ), TARGET :: SBID
C
C       RIGHT MEMBERS OF SYSTEM TO BE SOLVED       
C
        TYPE(BIEF_OBJ), TARGET :: CV1,CV2
C
C       WAVE DISSIPATION QUANTITIES
C
        TYPE(BIEF_OBJ), TARGET :: MU,MU2
        TYPE(BIEF_OBJ), TARGET :: QB
        TYPE(BIEF_OBJ), TARGET :: HMU,HMUANC
C
C       RADIATION STRESSES QUANTITIES
C
        TYPE(BIEF_OBJ), TARGET :: SXX,SXY,SYY
        TYPE(BIEF_OBJ), TARGET :: FX,FY
C
C       MEAN WAVE PERIODS
C
        TYPE(BIEF_OBJ), TARGET :: T01,T02,TM
C  
C       WAVE DIRECTIONS AT BOUNDARY (RANDOM MODE)
C  
        TYPE(BIEF_OBJ), TARGET :: DALE
C
C       BOUNDARY CONDITIONS TYPES
C
        TYPE(BIEF_OBJ), TARGET :: LIUBOR,LIVBOR,LIHBOR,NUMLIQ
C
C
        TYPE(BIEF_OBJ), TARGET :: LIDIR , MASKEL
C
C       MASKING ARRAYS FOR BOUNDARIES
C
        TYPE(BIEF_OBJ), TARGET :: MASK1 , MASK2 , MASK3 , MASK4
	
C
C
C --> ER : DEBUT
C -->   COURANT
C
        TYPE(BIEF_OBJ), TARGET :: UC,VC
	
C -->   PULSATION RELATIVE

        TYPE(BIEF_OBJ), TARGET :: WR
        	
C -->   VECTEUR REEL INTERMEDIAIRE : VECTEUR D'ONDE ET ERREUR

        TYPE(BIEF_OBJ), TARGET :: KN1,KN2,KNANC1,KNANC2
        	
	
C --> ER : FIN      
        
	
C --> ER : FIN      

C
C-----------------------------------------------------------------------
C
C       2) MATRICES
C
C-----------------------------------------------------------------------
C
C       MATRICES FOR SYSTEM SOLVING 
C        
        TYPE(BIEF_OBJ), TARGET :: AM1,AM2,AM3
        TYPE(BIEF_OBJ), TARGET :: BM1,BM2
        TYPE(BIEF_OBJ), TARGET :: MBOR
C
C-----------------------------------------------------------------------
C
C       3) BLOCKS
C
C-----------------------------------------------------------------------
C
C
C       BLOCK OF POTENTIAL VECTORS
C
        TYPE(BIEF_OBJ), TARGET :: PHIB
C
C       BLOCK OF WORKING ARRAYS
C
        TYPE(BIEF_OBJ), TARGET :: TB,TBBD
C
C       BLOCK OF PRIVATE VECTORS
C
        TYPE(BIEF_OBJ), TARGET :: PRIVE
C
C       BLOCK OF MATRICES
C
        TYPE(BIEF_OBJ), TARGET :: MAT
C
C       BLOCK OF UNKNOWN VECTORS
C
        TYPE(BIEF_OBJ), TARGET :: UNK
C
C       BLOCK OF RIGHT HAND SIDE VECTORS IN SOLVING SYSTEM
C
        TYPE(BIEF_OBJ), TARGET :: RHS
C
C       BLOCK OF VARIABLES FOR OUTPUT
C
        TYPE(BIEF_OBJ), TARGET :: VARSOR
C
C-----------------------------------------------------------------------
C
C       4) INTEGERS
C
C-----------------------------------------------------------------------
C
C       KEY-WORDS AND PARAMETERS
C
C       MAXIMUM DE VARIABLES DE SORTIE
C
        INTEGER, PARAMETER :: MAXVAR = 100
C       MAXIMUM DE FRONTIERES LIQUIDES
        INTEGER, PARAMETER :: MAXFRO = 100
C
      INTEGER LEOPRD , LISPRD , NITMAX 
      INTEGER STDGEO , STDRES , ISOLVE(2) , LISFON , DISESP
      INTEGER NPALE  , NDALE  , OPTASS
      INTEGER IBREAK , NITDIS , LVMAC  , KFROT
      INTEGER FORMFR , REGIDO , PRODUC , NPRIV
      INTEGER PTINIG , PTINIL
      INTEGER IELM,IELM0,IELMB,IELMB0
      INTEGER MARDAT(3),MARTIM(3)
      INTEGER NFRLIQ,NFRSOL,DEBLIQ(MAXFRO),FINLIQ(MAXFRO)
      INTEGER DEBSOL(MAXFRO),FINSOL(MAXFRO)
C     ORIGIN COORDINATES
      INTEGER I_ORIG,J_ORIG
C
C-----------------------------------------------------------------------
C
C       5) LOGICAL VALUES
C
C-----------------------------------------------------------------------
C
      LOGICAL LISTIN , INFOGR  , BALAYE
      LOGICAL ALEMON , ALEMUL  , MSK            , SPHERI , DEFERL
      LOGICAL FROTTE , ENTFW   , ENTREG         , ENTRUG , LISHOU
      LOGICAL SORLEO(MAXVAR)   , SORIMP(MAXVAR) , VALID
C --> ER : DEBUT
C -->   COURANT
C
      LOGICAL COURANT

C
C-----------------------------------------------------------------------
C
C       6) REALS
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION GRAV   , HMIN   , PER    , OMEGA
      DOUBLE PRECISION TETAH
      DOUBLE PRECISION COTINI , HAUTIN , PERDEB , PERFIN , PERPAS
      DOUBLE PRECISION PERPIC , GAMMA  , TETMIN , TETMAX , EXPOS
      DOUBLE PRECISION RELAX  , FFON
      DOUBLE PRECISION EPSDIS , RELDIS , ALFABJ , GAMMAS
      DOUBLE PRECISION KDALLY , GDALLY
      DOUBLE PRECISION VISCO  , DIAM90 , DIAM50
      DOUBLE PRECISION MVSED  , MVEAU
      DOUBLE PRECISION FWCOEF , RICOEF
      DOUBLE PRECISION PMIN   , PMAX
C --> ER : DEBUT
C -->   COURANT : VALEURS EN X, Y
C
      DOUBLE PRECISION CURRENTX,CURRENTY
 
C
C-----------------------------------------------------------------------
C
C       7) STRINGS
C
C-----------------------------------------------------------------------
C
      CHARACTER*72 TITCAS , VARDES , VARIMP , CDTINI
      CHARACTER*3  BINGEO,BINRES
      CHARACTER*20 EQUA
      CHARACTER*32 VARCLA(10),TEXTE(MAXVAR),TEXTPR(MAXVAR)
C
C-----------------------------------------------------------------------
C
C       8) SLVCFG STRUCTURES
C
C-----------------------------------------------------------------------
C
      TYPE(SLVCFG) :: SLVART
C
C-----------------------------------------------------------------------
C
C       9) MESH STRUCTURE
C
C-----------------------------------------------------------------------
C
        TYPE(BIEF_MESH) :: MESH
C
C-----------------------------------------------------------------------
C
C      10) ALIASES
C
C-----------------------------------------------------------------------
C
C       DECLARATION OF POINTERS FOR ALIASES.
C       TARGETS ARE DEFINED IN POINT_ARTEMIS
C
C       ALIASES FOR WORKING VECTORS IN TB AND TBBD
C
        TYPE(BIEF_OBJ),POINTER :: T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12
        TYPE(BIEF_OBJ),POINTER :: TBD1,TBD2,TBD3,TBD4
C
C       USEFUL COMPONENTS IN STRUCTURE MESH
C
        TYPE(BIEF_OBJ), POINTER :: IKLE
        DOUBLE PRECISION, DIMENSION(:), POINTER :: X,Y 
        INTEGER, POINTER        :: NELEM
        INTEGER, POINTER        :: NELMAX
        INTEGER, POINTER        :: NPTFR
        INTEGER, POINTER        :: NPTFRX
        INTEGER, POINTER        :: DIM
        INTEGER, POINTER        :: TYPELM
        INTEGER, POINTER        :: NPOIN
        INTEGER, POINTER        :: NPMAX
        INTEGER, POINTER        :: MXPTVS
        INTEGER, POINTER        :: MXELVS
        INTEGER, POINTER        :: LV
C-----------------------------------------------------------------------
C
C      10) ART_FILES and associated
C
C-----------------------------------------------------------------------
      INTEGER :: ARTGEO,ARTCAS,ARTCLI,ARTFON,ARTRES,ARTRBI,ARTRFO,
     *           ARTREF,ARTBI1,ARTBI2,ARTFO1,ARTFO2
      TYPE(BIEF_FILE) :: ART_FILES(44)    
C
        SAVE
C
      END MODULE DECLARATIONS_ARTEMIS
C
