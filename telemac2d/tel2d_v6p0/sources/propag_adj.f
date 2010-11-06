C                       *********************
                        SUBROUTINE PROPAG_ADJ
C                       *********************
C
     *(UCONV,VCONV,CONVV,H0,C0,COTOND,PATMOS,ATMOS,
     * HPROP,UN,VN,HN,UTILD,VTILD,HTILD,DH,DU,DV,DHN,VISC,VISC_S,FU,FV,
     * SMH,MESH,ZF,AM1,AM2,AM3,BM1,BM2,CM1,CM2,TM1,A23,A32,MBOR,
     * CV1,CV2,CV3,W1,UBOR,VBOR,AUBOR,HBOR,DIRBOR,
     * TE1,TE2,TE3,TE4,TE5,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,
     * LIMPRO,MASK,GRAV,ROEAU,CF,DIFVIT,IORDRH,IORDRU,LT,AT,DT,
     * TETAH,TETAHC,TETAU,TETAD,
     * AGGLOC,AGGLOU,KDIR,INFOGR,KFROT,ICONVF,
     * PRIVE,ISOUSI,BILMAS,MASSES,YASMH,OPTBAN,CORCON,
     * OPTSUP,MSK,MASKEL,MASKPT,RO,ROVAR,
     * MAT,RHS,UNK,TB,S,BD,PRECCU,SOLSYS,CFLMAX,OPDVIT,OPTSOU,
     * NFRLIQ,SLVPRO,EQUA,VERTIC,
     * ADJO,UD,VD,HD,U,V,H,UU,VV,HH,UIT1,VIT1,HIT1,PP,QQ,RR,
     * TAM1,TAM2,TAM3,TBM1,TBM2,TCM1,
     * TCM2,MATADJ,UNKADJ,ALPHA1,ALPHA2,ALPHA3,ADJDIR,ESTIME,OPTCOST,
     * NIT,NVARRES,VARSOR,
     * NRES,NREF,ALIRE,TROUVE,MAXVAR,VARCL,VARCLA,
     * TEXTE,TEXREF,TEXRES,W,OUTINI,CHESTR,KARMAN,NDEF,
     * ITURB,LISRUG,LINDNER,SB,DP,SP,CHBORD,CFBOR,HFROT,UNSV2D)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.5    24/04/97    J-M HERVOUET (LNH) 30 87 80 18
C                                          C MOULIN   (LNH) 30 87 83 81
C			     18/09/00      A LEOPARDI (UNINA)
C           COMPLETE VERSION 13/11/00
C***********************************************************************
C
C
C      FONCTION:
C      =========
C
C      COMPUTES THE RIGHT HAND SIDE OF THE ADJOINT SYSTEM
C
C
C      ECRITURE MATRICIELLE :  Systeme adjoint
C
C            T      N-1    T     N-1    T     N-1        *
C             AM1  P     +  CM1 Q     +  CM2 R     =  CV1
C
C            T     N-1     T     N-1                     *
C             BM1 P      +  AM2 Q                  =  CV2
C
C            T     N-1                  T     N-1        *
C             BM2 P                   +  AM3 R     =  CV3
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   U ,V ,H      |<-- |  VALEURS A L' ETAPE N+1 (System direct) OR
C |   UCONV,VCONV  | -->|  CHAMP CONVECTEUR                            |
C |   CONVV        | -->|  LOGIQUES INDIQUANT LES VARIABLES QUE L'ON   |
C |                |    |  VEUT CONVECTER                              |
C |                |    |  CONVV(1):U,V CONVV(2):H                     |
C |   C0           | -->|  CELERITE DE REFERENCE
C |   COTOND       |<-- |  EXPRESSION DE CU/G DANS LA THEORIE DE L'ONDE|
C |                |    |  INCIDENTE                                   |
C |   PATMOS       | -->|  TABLEAU DE VALEURS DE LA PRESSION ATMOSPHER.|
C |   ATMOS        | -->|  LOGIQUE INDIQUANT SI PATMOS EST REMPLI.     |
C |   HPROP        | -->|  HAUTEUR DE PROPAGATION                     |
C |   UN,VN,HN     | -->|  VALEURS A L' ETAPE N.
C |   UTILD,V.,H.  | -->|  VALEURS APRES LA CONVECTION.
C |   DH,DHN       |<-- |  STOCKAGE DE LA VARIABLE DH  (DHN AU TEMPS N)
C |   DU,DV        |<-- |  STOCKAGE DES QCCROISSEMENTS EN U ET V
C |   VISC         | -->|  VISCOSITE TURBULENTE .
C |   FU,FV        |<-->|  TERMES SOURCES TRAITES EN P1
C |   SMH          | -->|  TERMES SOURCES DE L'EQUATION DE CONTINUITE
C |   ZF           | -->|  COTE DU FONT AU NOEUD DE MAILLAGE .
C |   AM1,2,3      |<-->|  MATRICES
C |   BM1,2        |<-->|  MATRICES
C |   A23,A32      |<-->|  MATRICES
C |   TM1          |<-->|  MATRICE
C |   CV1,CV2,CV3  |<-->|  SECONDS MEMBRES DU SYSTEME.
C |   W1           |<-->|  TABLEAU DE TRAVAIL.
C |   UBOR         | -->|  CONDITIONS AUX LIMITES SUR U.
C |   VBOR         | -->|  CONDITIONS AUX LIMITES SUR V.
C |   AUBOR        | -->|  CONDITIONS AUX LIMITES SUR LE FROTTEMENT.
C |   HBOR         | -->|  CONDITIONS AUX LIMITES SUR H.
C |   LIMPRO       | -->|  TYPES DE CONDITIONS AUX LIMITES
C |   MASK         | -->|  BLOC DE MASQUES POUR LES SEGMENTS :
C |                |    |    MASK(MSK1): 1. SI KDIR SUR U 0. SINON
C |                |    |    MASK(MSK2): 1. SI KDIR SUR V 0. SINON
C |                |    |    MASK(MSK3): 1. SI KDDL SUR U 0. SINON
C |                |    |    MASK(MSK4): 1. SI KDDL SUR V 0. SINON
C |                |    |    MASK(MSK6): 1. SI KNEU SUR V 0. SINON
C |                |    |    MASK(MSK7): 1. SI KOND 0. SINON
C |                |    |    MASK(MSK9): 1. SI KDIR SUR H (POINT)
C |   GRAV         | -->|  CONSTANTE DE GRAVITE .
C |   ROEAU        | -->|  MASSE VOLUMIQUE DE L'EAU.
C |   PROLIN       | -->|  INDIQUE SI LA PROPAGATION EST LINEARISEE
C |   HAULIN       | -->|  HAUTEUR DE REFERENCE.
C |   CHESTR       | -->|  COEFFICIENT DE FROTTEMENT AU FOND.
C |   DIFVIT       | -->|  INDIQUE S'IL FAUT FAIRE LA DIFFUSION DE U,V
C |   IORDRH       | -->|  ORDRE DU TIR INITIAL POUR H
C |   IORDRU       | -->|  ORDRE DU TIR INITIAL POUR U
C |   LT,AT,DT     | -->|  NUMERO D'ITERATION, TEMPS, PAS DE TEMPS
C |   NPTFR        | -->|  NOMBRE DE POINTS FRONTIERES
C |   TETAH        | -->|  IMPLICITATION SUR H DANS L'EQUATION SUR U
C |   TETAHC       | -->|  IMPLICITATION SUR H DANS LA CONTINUITE
C |   TETAU        | -->|  IMPLICITATION SUR U ET V
C |   TETAD        | -->|  IMPLICITATION SUR LA DIFFUSION (=1.)
C |   AGGLOC       | -->|  COEFFICIENT DE MASS-LUMPING SUR H
C |   AGGLOU       | -->|  COEFFICIENT DE MASS-LUMPING SUR U
C |   KDIR         | -->|  CONDITION A LA LIMITE DE TYPE DIRICHLET
C |   INFOGR       | -->|  INFORMATIONS SUR LE GRADIENT (LOGIQUE)
C |   KFROT        | -->|  LOI DE FROTTEMENT SUR LE FOND
C |   ICONVF       | -->|  FORME DE LA CONVECTION
C |                |    |  TABLEAU DE 4 VALEURS ENTIERES POUR :
C |                |    |  ICONVF(1) : U ET V
C |                |    |  ICONVF(2) : H
C |                |    |  ICONVF(3) : TRACEUR
C |                |    |  ICONVF(4) : K ET EPSILON
C |   PRIVE        | -->|  TABLEAU DE TRAVAIL DEFINI DANS PRINCI
C |   ISOUSI       | -->|  NUMERO DE LA SOUS-ITERATION DANS LE PAS
C |                |    |  DE TEMPS.
C |   BILMAS       | -->|  INDIQUE SI ON FAIT LE BILAN DE MASSE
C |   MASSES       | -->|  MASSE CREEE PAR TERME SOURCE PENDANT
C |                |    |  LE PAS DE TEMPS.
C |   YASMH        | -->|  INDIQUE SI ON PREND EN COMPTE SMH
C |   OPTBAN       | -->|  OPTION DE TRAITEMENT DES BANCS DECOUVRANTS
C |                |    |
C |                |    |
C |                |    |  NON UTILISE POUR L'INSTANT :
C |   CORCON       | -->|  CORRECTION DE CONTINUITE SUR LES POINTS A
C |                |    |  HAUTEUR IMPOSEE (ON CORRIGE LES VITESSES)
C |   COSUPG       | -->|  COEFFICIENTS DE DECENTREMENT POUR S.U.P.G
C |                |    |  TABLEAU DE 4 REELS COMPRIS ENTRE 0 ET 1
C |                |    |  COSUPG(1) : U,V          |
C |                |    |  COSUPG(2) : H            | 0 CORRESPOND AU
C |                |    |  COSUPG(3) : TRACEUR      | SCHEMA EXPLICITE
C |                |    |  COSUPG(4) : K ET EPSILON | CENTRE
C |   OPTHYB       | -->|  OPTION DU SCHEMA HYBRIDE
C |   MSK          | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |   MASKEL       | -->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |   MASKPT       | -->|  MASQUES PAR POINTS.
C |   KARMAN       | -->|  CONSTANTE DE KARMAN.
C |   RO           | -->|  MASSE VOLUMIQUE SI ELLE VARIABLE
C |   ROVAR        | -->|  OUI SI LA MASSE VOLUMIQUE EST VARIABLE.
C |   S            | -->|  STRUCTURE BIDON
C |---------------------------------------------------------------------
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMMES APPELES :
C
C**********************************************************************
C
      USE BIEF
      USE INTERFACE_TELEMAC2D, EX_PROPAG_ADJ => PROPAG_ADJ
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C    
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: LT,OPTSUP(4),KDIR,KFROT,ICONVF(4)
      INTEGER, INTENT(IN) :: IORDRH,IORDRU,ISOUSI,OPTBAN,OPTSOU,SOLSYS
      INTEGER, INTENT(IN) :: OPDVIT,NFRLIQ,LISRUG,ITURB,OPTCOST
      INTEGER, INTENT(IN)    :: NIT,NRES,NREF,MAXVAR,HFROT
      INTEGER, INTENT(INOUT) :: NVARRES,TROUVE(*),ALIRE(*)
      LOGICAL, INTENT(IN)    :: BILMAS,ATMOS,DIFVIT,INFOGR,CONVV(4),MSK
      LOGICAL, INTENT(IN)    :: YASMH,ROVAR,PRECCU,VERTIC,ADJO,CORCON
      LOGICAL, INTENT(IN)    :: OUTINI,LINDNER
      DOUBLE PRECISION, INTENT(IN)    :: TETAU,TETAD,TETAH,AGGLOC,AGGLOU
      DOUBLE PRECISION, INTENT(IN)    :: TETAHC,AT,DT,GRAV,ROEAU,CFLMAX
      DOUBLE PRECISION, INTENT(IN)    :: KARMAN,NDEF,DP,SP
      DOUBLE PRECISION, INTENT(INOUT) :: MASSES,SB
      TYPE(SLVCFG), INTENT(INOUT)     :: SLVPRO
      CHARACTER(LEN=20), INTENT(IN)   :: EQUA
      TYPE(BIEF_OBJ), INTENT(IN)      :: UCONV,VCONV,SMH,UN,VN,HN
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: RO
      TYPE(BIEF_OBJ), INTENT(IN)      :: UTILD,VTILD,PATMOS,CF,UNSV2D
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: U,V,H,CV1,CV2,CV3,PRIVE,DH,DHN
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: DU,DV,FU,FV,VISC,VISC_S,HTILD
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: UBOR,VBOR,HBOR,AUBOR
      TYPE(BIEF_OBJ), INTENT(IN)      :: MASKEL,MASKPT,ZF
      TYPE(BIEF_OBJ), INTENT(IN)      :: HPROP,H0,C0,COTOND,LIMPRO
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: T1,T2,T3,T4,T5,T6,T7,T8,T9,T10
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: T11
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: TE1,TE2,TE3,TE4,TE5
C     STRUCTURES DE MATRICES
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: TAM1,TAM2,TAM3,TBM1
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: TBM2,TCM1,TCM2
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: AM1,AM2,AM3,BM1,BM2,CM1,CM2,TM1
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: A23,A32,MBOR
C
      TYPE(BIEF_OBJ),  INTENT(INOUT)  :: MASK,MAT,RHS,UNK,TB,BD,DIRBOR
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
      TYPE(BIEF_OBJ),  INTENT(INOUT)  :: CHESTR
      TYPE (BIEF_OBJ), INTENT(INOUT)  :: HD,UD,VD,ALPHA1,ALPHA2,ALPHA3
      TYPE (BIEF_OBJ), INTENT(INOUT)  :: HH,UU,VV,UIT1,VIT1,HIT1
      TYPE (BIEF_OBJ), INTENT(INOUT)  :: PP,QQ,RR,CHBORD,CFBOR
      TYPE(BIEF_OBJ),  INTENT(INOUT)  :: W1
      TYPE(BIEF_OBJ), INTENT(IN)      :: S
      REAL,  INTENT(INOUT)            :: W(*)
      TYPE(BIEF_OBJ),  INTENT(INOUT)  :: VARSOR
      TYPE(BIEF_OBJ),  INTENT(INOUT)  :: MATADJ,UNKADJ,ADJDIR,VARCL
      CHARACTER(LEN=72), INTENT(IN)   :: ESTIME
      CHARACTER(LEN=32), INTENT(INOUT):: VARCLA(10),TEXTE(*)
      CHARACTER(LEN=32), INTENT(INOUT):: TEXREF(*),TEXRES(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER ITER,I,IELMU,IELMH
      INTEGER UDIR,UDDL,UNEU,HOND,UNONNEU,VDIR,VDDL 
C
      DOUBLE PRECISION Z(1),SL1,SL1U,C,AT1,HIST(1)
C
      LOGICAL MSKGRA
C
      CHARACTER*16 FORMULE
C
C-----------------------------------------------------------------------
C
C FH-FRDATA
      DOUBLE PRECISION, PARAMETER :: VK = 1.D-6
C FH-FRDATA
C-----------------------------------------------------------------------
C
      IELMH=HH%ELM
      IELMU=UU%ELM
C
C  ADRESSES DES TABLEAUX DANS LE BLOC DES MASQUES MASK
C
      UDIR = 1
      VDIR = 2
      UDDL = 3
      VDDL = 4
      UNEU = 5
C     VNEU = 6
      HOND = 7
      UNONNEU = 8
C
C-----------------------------------------------------------------------
C
      IF(SOLSYS.NE.1) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'TRAITEMENT DU SYSTEME LINEAIRE : ',SOLSYS
          WRITE(LU,*) 'CAS NON PREVU EN MODE ESTIMATION'
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'TREATMENT OF THE LINEAR SYSTEM : ',SOLSYS
          WRITE(LU,*) 'UNEXPECTED CASE IN ESTIMATION MODE'
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C     COMPUTATION OF MATRIX FOR ADJOINT SYSTEM
C
      CALL OM( 'M=TN    ' , TAM1, AM1, S, C, MESH )
      CALL OM( 'M=TN    ' , TAM2, AM2, S, C, MESH )
      CALL OM( 'M=TN    ' , TAM3, AM3, S, C, MESH )
      CALL OM( 'M=TN    ' , TBM1, BM1, S, C, MESH )
      CALL OM( 'M=TN    ' , TBM2, BM2, S, C, MESH )
      CALL OM( 'M=TN    ' , TCM1, CM1, S, C, MESH )
      CALL OM( 'M=TN    ' , TCM2, CM2, S, C, MESH )
C  
C=======================================================================
C
C     COMPUTATION OF RIGHT HAND SIDES FOR ADJOINT SYSTEM
C
C=======================================================================
C
C     NB: HIT1, UIT1, VIT1 ARE DIRECT VARIABLES AT TIME IT+1
C         HH  , UU  , VV   ARE DIRECT VARIABLES AT TIME IT
C         HN  , UN  , VN   ARE DIRECT VARIABLES AT TIME IT-1
C
C
C           IT    IT    IT
C  TERMS 2 W   ( X   - M   ) OR EQUIVALENT DEPENDING ON THE COST FUNCTION
C           IP    IP    IP
C
C     INITIALISES CV1, CV2 AND CV3
C     IN STEADY STATE MODE, WEIGHTS ALPHA1, ALPHA2 AND ALPHA3 ARE
C     INITIALISED BY A CALL TO "MESURES" IN HOMERE_T2D_ADJ, IN THE LOOP
C     FOR THE COMPUTATION OF THE COST FUNCTION. THEN THEY ARE CANCELLED
C     AT THE END OF THIS ROUTINE
C
      CALL COST_FUNCTION(C,OPTCOST,'RHS')
C
C-----------------------------------------------------------------------
C
C  PREPARATION OF FRICTION TERMS AND VECTOR T1 EQUAL TO 0
C
C     T10 =MASS-MATRIX LUMPED / COS(SLOPE)
      CALL SLOPES(TE3,ZF,MESH)
      CALL VECTOR(T10,'=','MASBAS          ',IELMU,1.D0,T2,
     *                T2,T2,T2,T2,T2,MESH,.TRUE.,TE3)
C     FU PUT IN T11 (AND FV=FU)
C
C     T2 VA CONTENIR CF DE L'ITERATION IT+1
C
      CALL CPSTVC(CF,T2)
C
CFH-FRDATA
C     CALL COEFRO(T2,HH,UU,VV,KARMAN,KFROT,CHESTR,GRAV,MESH,T1)
      CALL FRICTION_UNIF(MESH,HH,UU,VV,CHESTR,S,KFROT,0,LISRUG,
     &                   .FALSE.,SB,NDEF,DP,SP,VK,KARMAN,GRAV,
     &                   T1,T2,CHBORD,T2,CFBOR)
CFH-FRDATA
C                                        
      CALL FRICTI(T11,T3,T4,T5,UU,VV,HH,T2,MESH,T6,T7,VERTIC,UNSV2D,
     *            MSK,MASKEL,HFROT)
C     CALL FRICTI(T11,T3,T4,T5,UU,VV,HH,CF,MESH,T6,VERTIC)
C
C     FINAL FU OF PHD PUT IN T11
      CALL OS('X=XY    ', T11 , T10 , T10 , C )
C
C     T1 : VECTEUR NUL
      CALL CPSTVC(HH,T1)
      CALL OS('X=C     ' , T1 , T1 , T1 , 0.D0)
C
C-----------------------------------------------------------------------
C
C  COMPUTATION OF CV1 : CV1 = CV1 + T11+T12+T13+T14+T15+T16+T17
C
C-----------------------------------------------------------------------
C  TERM T11 (3 PARTS)
C-----------------------------------------------------------------------
C
C  T11_1 : M/DT  * H
C                   ADJ
C     TERME 1B DE AL
      CALL MATRIX(AM1,'M=N     ','MATMAS          ',IELMH,IELMH,
     *            1.D0/DT,T1,T1,T1,T1,T1,T1,MESH,MSK,MASKEL)
      CALL MATVEC('X=X+AY  ', CV1 , AM1 , PP , C , MESH)
C
C  T11_2 : ADVECTION 
C
C  T11_3 :  
C
C     8B DE AL
C     CALL MATRIX(BM1,'M=N     ','MATGRF         X',IELMH,IELMH,
C    *           (TETAU-1.D0),UU,T1,T1,T1,T1,T1,MESH,MSK,MASKEL)
C     CALL MATVEC('X=X+AY  ', CV1 , BM1 , PN , C , MESH)
C     9B DE AL
C     CALL MATRIX(BM2,'M=N     ','MATGRF         Y',IELMH,IELMH,
C    *           (TETAU-1.D0),VV,T1,T1,T1,T1,T1,MESH,MSK,MASKEL)
C     CALL MATVEC('X=X+AY  ', CV1 , BM2 , PN , C , MESH)
C
C     VERSION CORRIGEE JMH: NOW ICONVF(2) IS ALWAYS 5
C
!     IF(ICONVF(2).EQ.5.OR.ICONVF(2).EQ.8) THEN
        FORMULE='MATFGR          '
!     ELSE
!       FORMULE='MATGRF          '
!     ENDIF
C
      FORMULE(16:16)='X'
      CALL MATRIX(BM1,'M=N     ',FORMULE,IELMH,IELMU,
     *           (TETAU-1.D0),PP,T1,T1,T1,T1,T1,MESH,MSK,MASKEL)
      CALL MATVEC('X=X+AY  ', CV1 , BM1 , UU , C , MESH)
      FORMULE(16:16)='Y'
      CALL MATRIX(BM2,'M=N     ',FORMULE,IELMH,IELMU,
     *           (TETAU-1.D0),PP,T1,T1,T1,T1,T1,MESH,MSK,MASKEL)
      CALL MATVEC('X=X+AY  ', CV1 , BM2 , VV , C , MESH)
C
C                           H
C  T11_4 : BOUNDARY TERM TB1
C                           ADJ
C
C     JMH : MOI J'AI -1.D0 AU LIEU DE TETAU-1.D0, MAIS ENSUITE IL RETRANCHE TETAU ??
C           CA REVIENT AU MEME SAUF ONDE INCIDENTE ???
      CALL VECTOR(T2,'=','FLUBDF          ',IELBOR(IELMH,1),
     *            TETAU-1.D0,PP,T1,T1,UU,VV,T1,MESH,.TRUE.,
     *            MASK%ADR(8)%P)
      CALL OSDB('X=X+Y   ' , CV1 , T2 , T2 , C , MESH)
C     DIRICHLET ON VELOCITY :
      CALL VECTOR(T2,'=','FLUBDF          ',IELBOR(IELMH,1),
     *            -TETAU,PP,T1,T1,UU,VV,T1,MESH,
     *            .TRUE.,MASK%ADR(UDIR)%P)
      CALL OSDB('X=X+Y   ' , CV1 , T2 , T2 , C , MESH)
C     FREE FLOW ON VELOCITY :
      CALL VECTOR(T2,'=','FLUBDF          ',IELBOR(IELMH,1),
     *            -TETAU,PP,T1,T1,UU,VV,T1,MESH,
     *            .TRUE.,MASK%ADR(UDDL)%P)
      CALL OSDB('X=X+Y   ' , CV1 , T2 , T2 , C , MESH)
C
C-----------------------------------------------------------------------
C  TERM T12
C-----------------------------------------------------------------------
C
C     TERME 2B DE AL 
      CALL MATVEC('X=X+CAY ',CV1,TCM1,QQ,(TETAH-1.D0)/TETAH,MESH)
C
C-----------------------------------------------------------------------
C  TERM T13
C-----------------------------------------------------------------------
C
C     TERME 3B DE AL 
      CALL MATVEC('X=X+CAY ',CV1,TCM2,RR,(TETAH-1.D0)/TETAH,MESH)
C 
C-----------------------------------------------------------------------
C  TERM T14
C-----------------------------------------------------------------------
C 
C     TERME 1C DE AL     
C     VERSION EB+AL, NOTE JMH : PAS D'ACCORD     
C     CALL MATRIX(BM1,'M=N     ','MATGRF         X',IELMH,IELMH,
C    *            -TETAU,UN,T1,T1,T1,T1,T1,MESH,MSK,MASKEL)
C     CALL MATVEC('X=X+AY  ', CV1 , BM1 , PN , C , MESH)
C
C     VERSION JMH+CC
C                        DEBUT DE FORMULE FAIT POUR T11_3
      FORMULE(16:16)='X'
      CALL MATRIX(BM1,'M=N     ',FORMULE,IELMH,IELMU,
     *            -TETAU,PP,T1,T1,T1,T1,T1,MESH,MSK,MASKEL)
      CALL MATVEC('X=X+AY  ', CV1 , BM1 , UIT1 , C , MESH)
C
C-----------------------------------------------------------------------
C  TERM T15 + T17
C-----------------------------------------------------------------------
C
C          IT+1   IT+1    IT+1    IT+1
C     T3= U    * Q     + V     * R
C
      CALL OS('X=YZ    ', T3 , UIT1 , QQ , C )
      CALL OS('X=X+YZ  ', T3 , VIT1 , RR , C )
C
C     T4=-(4/3)/H OR -1/H
      CALL OS('X=1/Y   ', T4 , HH , HH , C ,2,0.D0,1.D-6)
      IF(KFROT.EQ.3) THEN
        CALL OS('X=CX    ', T4 , T4 , T4 , -4.D0/3.D0)
      ELSEIF(KFROT.EQ.2) THEN
        CALL OS('X=CX    ', T4 , T4 , T4 , -1.D0     )
      ELSE
        IF(LNG.EQ.1) WRITE(LU,*) 'LOI NON TRAITEE POUR L''ESTIMATION'
        IF(LNG.EQ.2) WRITE(LU,*) 'WRONG FRICTION LAW FOR ESTIMATION'
        CALL PLANTE(1)
        STOP
      ENDIF
      CALL OS('X=XY    ', T4 , T11 , T11 , C )
C     ET SI T3 EST QUASI-BULLE ?
      CALL OS('X=X+YZ  ', CV1 , T4 , T3  , C ) 
C
C-----------------------------------------------------------------------
C  TERM T16
C-----------------------------------------------------------------------
C
C     TERME 2C DE AL 
C     VERSION EB+AL, NOTE JMH : PAS D'ACCORD  
C     CALL MATRIX(BM2,'M=N     ','MATGRF         Y',IELMH,IELMH,
C    *            -TETAU,VN,T1,T1,T1,T1,T1,MESH,MSK,MASKEL)
C     CALL MATVEC('X=X+AY  ', CV1 , BM2 , PN , C , MESH)
C     VERSION JMH+CC
C                        DEBUT DE FORMULE FAIT POUR T11_3
      FORMULE(16:16)='Y'
      CALL MATRIX(BM2,'M=N     ',FORMULE,IELMH,IELMU,
     *            -TETAU,PP,T1,T1,T1,T1,T1,MESH,MSK,MASKEL)
      CALL MATVEC('X=X+AY  ', CV1 , BM2 , VIT1 , C , MESH)
C
C-----------------------------------------------------------------------
C
C  COMPUTATION OF CV2 AND CV3 :
C
C-----------------------------------------------------------------------
C  TERM  T21
C-----------------------------------------------------------------------
C
C  T21_1 : ADVECTION
C
C  T21_2 : (5B CHEZ EB+AL)
C
      FORMULE(16:16)='X'
      CALL MATRIX(BM1,'M=N     ',FORMULE,IELMH,IELMU,
     *            1.D0,HH,T1,T1,T1,T1,T1,MESH,MSK,MASKEL)
      CALL MATVEC('X=X+CTAY',CV2,BM1,PP,TETAU-1.D0,MESH)
C ANCIENNE PROGRAMMATION (TBM1 FAITE AVEC HN PAS AU BON NIVEAU DE TEMPS).
C     CALL MATVEC('X=X+CAY ',CV2,TBM1,PP,(TETAU-1.D0)/TETAU,MESH)
C
C                           U
C  T21_3 : BOUNDARY TERM TB1
C                           ADJ
C
C     JMH : MOI J'AI 1.D0 AU LIEU DE TETAU-1.D0
      CALL VECTOR(T2,'=','FLUBDF          ',IELBOR(IELMH,1),
     *            TETAU-1.D0,PP,T1,T1,HH,T1,T1,MESH,.TRUE.,
     *            MASK%ADR(8)%P)
      CALL OSDB('X=X+Y   ' , CV2 , T2 , T2 , C , MESH)
C     DIRICHLET ON VELOCITY U:
      CALL VECTOR(T2,'=','FLUBDF          ',IELBOR(IELMH,1),
     *            -TETAU,PP,T1,T1,HH,T1,T1,MESH,
     *            .TRUE.,MASK%ADR(UDIR)%P)
      CALL OSDB('X=X+Y   ' , CV2 , T2 , T2 , C , MESH)
C     FREE FLOW ON U:
      CALL VECTOR(T2,'=','FLUBDF          ',IELBOR(IELMH,1),
     *            -TETAU,PP,T1,T1,HH,T1,T1,MESH,
     *            .TRUE.,MASK%ADR(UDDL)%P)
      CALL OSDB('X=X+Y   ' , CV2 , T2 , T2 , C , MESH)
C
C-----------------------------------------------------------------------
C TERM  T31
C-----------------------------------------------------------------------
C
C  T31_1 : ADVECTION
C
C  T31_2 : (7B CHEZ EB+AL)
C
      FORMULE(16:16)='Y'
      CALL MATRIX(BM2,'M=N     ',FORMULE,IELMH,IELMU,
     *            1.D0,HH,T1,T1,T1,T1,T1,MESH,MSK,MASKEL)
      CALL MATVEC('X=X+CTAY',CV3,BM2,PP,TETAU-1.D0,MESH)
C ANCIENNE PROGRAMMATION (TBM2 FAITE AVEC HN PAS AU BON NIVEAU DE TEMPS).
C     CALL MATVEC('X=X+CAY ',CV3,TBM2,PP,(TETAU-1.D0)/TETAU,MESH)
C
C                           V
C  T31_3 : BOUNDARY TERM TB1
C                           ADJ
C
C     JMH : MOI J'AI 1.D0 AU LIEU DE TETAU-1.D0
      CALL VECTOR(T4,'=','FLUBDF          ',IELBOR(IELMH,1),
     *            TETAU-1.D0,PP,T1,T1,T1,HH,T1,MESH,.TRUE.,
     *            MASK%ADR(8)%P)
      CALL OSDB('X=X+Y   ' , CV3 , T4 , T4 , C , MESH)
C     DIRICHLET ON VELOCITY V:
      CALL VECTOR(T4,'=','FLUBDF          ',IELBOR(IELMH,1),
     *            -TETAU,PP,T1,T1,T1,HH,T1,MESH,
     *            .TRUE.,MASK%ADR(VDIR)%P)
      CALL OSDB('X=X+Y   ' , CV3 , T4 , T4 , C , MESH)
C     FREE FLOW ON V:
      CALL VECTOR(T4,'=','FLUBDF          ',IELBOR(IELMH,1),
     *            -TETAU,PP,T1,T1,T1,HH,T1,MESH,
     *            .TRUE.,MASK%ADR(VDDL)%P)
      CALL OSDB('X=X+Y   ' , CV3 , T4 , T4 , C , MESH)
C
C-----------------------------------------------------------------------
C TERM  T22 (2 PARTS) 
C-----------------------------------------------------------------------
C
C     TERME 4B DE AL       T
C     AM2 IS MASSE/DT AT TIME IT+1 
      CALL MATRIX(AM2,'M=N     ','MATMAS          ',IELMU,IELMU,
     *            1.D0/DT,T1,T1,T1,T1,T1,T1,MESH,MSK,MASKEL)
      CALL MATVEC('X=X+AY  ', CV2 , AM2 , QQ , C , MESH)
C
C     MANQUE UN TERME DE CONVECTION
C
C
C-----------------------------------------------------------------------
C TERM  T32 (2 PARTS)
C-----------------------------------------------------------------------
C
C     TERME 6B DE AL 
C     AM2 IS MASSE/DT AT TIME IT+1
      CALL MATVEC('X=X+AY  ', CV3 , AM2 , RR , C , MESH)
C
C     MANQUE UN TERME DE CONVECTION
C
C-----------------------------------------------------------------------
C TERMS  T23 AND T33
C-----------------------------------------------------------------------
C
C   ADVECTION : NOT YET IMPLEMENTED
C
C-----------------------------------------------------------------------
C TERM  T24+T25 AND T34+T35
C-----------------------------------------------------------------------
C
C     T5=U/(U^2+V^2)  C     T6=V/(U^2+V^2)
      CALL OS('X=YZ    ', T7 , UU , UU , C    )
      CALL OS('X=X+YZ  ', T7 , VV , VV , C    )
      CALL OS('X=+(Y,C)', T7 , T7 , T7 , 1.D-6)
      CALL OS('X=Y/Z   ', T5 , UU , T7 , C    )
      CALL OS('X=Y/Z   ', T6 , VV , T7 , C    )
C
C     ADD TERMS TO CV2, CV3
C
C     T3=U*Q+V*R (DEJA FAIT PLUS HAUT)        
      CALL OS('X=XY    ', T5 , T11 , T11 , C )
      CALL OS('X=XY    ', T6 , T11 , T11 , C )  
      CALL OS('X=X+YZ  ', CV2 , T5 , T3  , C )
      CALL OS('X=X+YZ  ', CV3 , T6 , T3  , C )
C
C=======================================================================
C
C     END OF COMPUTATION OF RIGHT HAND SIDE FOR ADJOINT SYSTEM
C
C=======================================================================
C
C     
C     DIRICHLET CONDITIONS FOR ADJOINT VARIABLES
      CALL OS ('X=C     ',ADJDIR,ADJDIR,ADJDIR,0.D0)           
      CALL DIRICH(UNKADJ,MATADJ,RHS,ADJDIR,LIMPRO%I,
     *            TB,MESH,KDIR,MSK,MASKPT)
C
      CALL SOLVE(UNKADJ,MATADJ,RHS,TB,SLVPRO,INFOGR,MESH,TM1)
C
C     CONTRIBUTION TO COST FUNCTION
C
      CALL COST_FUNCTION(C,OPTCOST,'GRD')
C
C     PREPARING NEXT TIME-STEP
C
      CALL OS( 'X=Y     ' , HIT1 , HH  , HH  , C )
      CALL OS( 'X=Y     ' , UIT1 , UU  , UU  , C )
      CALL OS( 'X=Y     ' , VIT1 , VV  , VV  , C )
C
      IF(     INCLU2(ESTIME,'PERMANENT')
     *    .OR.INCLU2(ESTIME,'STEADY'   )  ) THEN
C
C      STEADY STATE : NO UPDATING OF DATA AND RESULTS, 
C                     ONLY LAST TIME-STEP CONSIDERED
C
C      CALL OS( 'X=C     ' , ALPHA1 , ALPHA1 , ALPHA1 , 0.D0 )
C      CALL OS( 'X=C     ' , ALPHA2 , ALPHA2 , ALPHA2 , 0.D0 )
C      CALL OS( 'X=C     ' , ALPHA3 , ALPHA3 , ALPHA3 , 0.D0 )
C      U AND V MODIFIED BY BORD, RESET HERE (H USEFUL ??)
       CALL OS( 'X=Y     ' , H , HN  , HN  , C )
       CALL OS( 'X=Y     ' , U , UN  , UN  , C )
       CALL OS( 'X=Y     ' , V , VN  , VN  , C )
C
      ELSE
C
C      UNSTEADY STATE : UPDATING DATA AND RESULTS
C
       IF(LT.LT.NIT) THEN
C
C       HIT,.., HH,.. IN INITIAL CONDITIONS, SEE PROPIN_ADJ
        CALL OS( 'X=Y     ' , HH , HN  , HN  , C )
        CALL OS( 'X=Y     ' , UU , UN  , UN  , C )
        CALL OS( 'X=Y     ' , VV , VN  , VN  , C )
C         
C       READING OF TELEMAC2D RESULTS (RESULTS FILE - UNIT NRES)
C       SEE ALSO CONDIN_ADJ
C
        DO I=1,2*(NVARRES+1)
          BACKSPACE NRES
        ENDDO 
        CALL LITENR(VARSOR,VARCL,NRES,'STD',
     *       HIST,0,MESH%NPOIN,AT1,TEXTE,
     *       TEXRES,NVARRES,VARCLA,0,TROUVE,ALIRE,W,.FALSE.,MAXVAR)
C
C       READING THE MEASUREMENTS (REFERENCE FILE - UNIT NREF)
C
        ITER=NIT-LT
        IF(OUTINI) ITER=ITER+1
        CALL MESURES(ITER,AT-DT)
C
       ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
