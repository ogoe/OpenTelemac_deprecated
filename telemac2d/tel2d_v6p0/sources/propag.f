C                       *****************
                        SUBROUTINE PROPAG
C                       *****************
C
     *(U,V,H,UCONV,VCONV,CONVV,H0,C0,COTOND,PATMOS,ATMOS,
     * HPROP,UN,VN,HN,UTILD,VTILD,HTILD,DH,DU,DV,DHN,VISC,VISC_S,FU,FV,
     * SMH,MESH,ZF,AM1,AM2,AM3,BM1,BM2,CM1,CM2,TM1,A23,A32,MBOR,
     * CV1,CV2,CV3,W1,UBOR,VBOR,AUBOR,HBOR,DIRBOR,
     * TE1,TE2,TE3,TE4,TE5,T1,T2,T3,T4,T5,T6,T7,T8,
     * LIMPRO,MASK,GRAV,ROEAU,CF,DIFVIT,IORDRH,IORDRU,LT,AT,DT,
     * TETAH,TETAHC,TETAU,TETAD,
     * AGGLOH,AGGLOU,KDIR,INFOGR,KFROT,ICONVF,
     * PRIVE,ISOUSI,BILMAS,MASSES,YASMH,OPTBAN,CORCON,
     * OPTSUP,MSK,MASKEL,MASKPT,RO,ROVAR,
     * MAT,RHS,UNK,TB,S,BD,PRECCU,SOLSYS,CFLMAX,OPDVIT,OPTSOU,
     * NFRLIQ,SLVPRO,EQUA,VERTIC,ADJO,ZFLATS,TETAZCOMP,UDEL,VDEL,DM1,
     * ZCONV,COUPLING,FLBOR,BM1S,BM2S,CV1S,VOLU2D,V2DPAR,UNSV2D,
     * NUMDIG,NWEIRS,NPSING,HFROT)
C
C***********************************************************************
C  TELEMAC 2D VERSION 6.0  09/10/2009 J-M HERVOUET (LNHE) 01 30 87 80 18
C                                       
C  07/05/07 : MODIFICATION SUR LES SOURCES EN CAS DE MASS-LUMPING (JMH)
C  10/06/08 : CONVECTEUR VOLUMES FINIS POUR LES VITESSES
C  02/10/08 : APPEL DE CVTRVF (UN TABLEAU DE TRAVAIL EN PLUS)
C  20/07/09 : ICONVF(2) = 5 MANDATORY, ALL OTHER CASES ERASED
C  22/07/09 : EQUALITY OF FLUXES IMPOSED ON EITHER SIDE OF A WEIR
C  09/10/09 : PARAMETERIZING ADVECTION OPTIONS
C
C***********************************************************************
C
C      FONCTION:
C      =========
C
C      ETAPE DE PROPAGATION - DIFFUSION - TERMES SOURCES POUR RESOUDRE
C      LES EQUATIONS DE SAINT-VENANT.
C
C      CONDITIONS AU LIMITES :
C
C      ==>   CONDITION DE NEUMANN
C
C            * DIFFUSION   : NU DU/DN = AUBOR . U
C                            TRAITE DIRECTEMENT SUR LA MATRICE DE
C                            DIFFUSION
C
C            * PROPAGATION : LES TERMES DE BORD SONT TRAITES DANS LES
C                            SECONDS MEMBRES (IMPLICITES).
C
C      ==>   CONDITION DE DIRICHLET
C
C            * DIFFUSION, PROPAGATION :
C                            TRAITE PAR MODIFICATION DES EQUATIONS
C                            DANS " PROCLI ".
C
C      ECRITURE MATRICIELLE :
C
C                   N+1          N+1          N+1
C             AM1  H     +  BM1 U     +  BM2 V     =  CV1
C
C            T     N+1           N+1
C          -  CM1 H      +  AM2 U                  =  CV2
C
C            T     N+1                        N+1
C          -  CM2 H                   +  AM3 V     =  CV3
C
C
C
C    AVERTISSEMENT :
C
C      ==>   BM. REPRESENTENT LES MATRICES DE DIVERGENCE
C            BM1 : DERIVATION PAR RAPPORT A X
C            BM2 : DERIVATION PAR RAPPORT A Y
C
C
C      ==>   LA TRANSPOSEE DES MATRICES BM. EST EGALE A L'OPPOSE  DU
C            GRADIENT. DE CE FAIT, DES SIGNES SONT INVERSES DANS LES
C            EQUATIONS EN VITESSE.
C
C
C      ==>   LA MATRICE LAPLACIEN ( TM1) A ETE INTEGREE PAR PARTIE .
C            DE CE FAIT , LE SIGNE EST INVERSE DANS LES EQUATIONS EN
C            VITESSE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   U ,V ,H      |<-- |  VALEURS A L' ETAPE N+1.                     |
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
C |   TETAD        | -->|  IMPLICITATION SUR LA DIFFUSION
C |   AGGLOH       | -->|  COEFFICIENT DE MASS-LUMPING SUR H
C |   AGGLOU       | -->|  COEFFICIENT DE MASS-LUMPING SUR U
C |   KDIR         | -->|  CONDITION A LA LIMITE DE TYPE DIRICHLET
C |   INFOGR       | -->|  INFORMATIONS SUR LE GRADIENT (LOGIQUE)
C |   KFROT        | -->|  LOI DE FROTTEMENT SUR LE FOND
C |   ICONVF       | -->|  FORME DE LA CONVECTION
C |                |    |  TABLEAU DE 4 VALEURS ENTIERES POUR :
C |                |    |  ICONVF(1) : U ET V
C |                |    |  ICONVF(2) : H (MANDATORY VALUE = 5)
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
C |   ADJO         | -->|  IF YES : ADJOINT MODE
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
      USE DECLARATIONS_TELEMAC, ONLY : ADV_CAR,ADV_SUP,ADV_NSC,ADV_PSI,
     *   ADV_PSI_NC,ADV_NSC_NC,ADV_LPO,ADV_NSC_TF,ADV_PSI_TF,ADV_LPO_TF
C
      USE INTERFACE_TELEMAC2D, EX_PROPAG => PROPAG
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: LT,OPTSUP(4),KDIR,KFROT,ICONVF(4),NWEIRS
      INTEGER, INTENT(IN) :: IORDRH,IORDRU,ISOUSI,OPTBAN,OPTSOU,SOLSYS
      INTEGER, INTENT(IN)             :: OPDVIT,NFRLIQ,HFROT
      DOUBLE PRECISION, INTENT(IN)    :: TETAU,TETAD,TETAH,AGGLOH,AGGLOU
      DOUBLE PRECISION, INTENT(IN)    :: TETAHC,AT,DT,GRAV,ROEAU
      DOUBLE PRECISION, INTENT(IN)    :: TETAZCOMP
      DOUBLE PRECISION, INTENT(INOUT) :: CFLMAX,MASSES
      LOGICAL, INTENT(IN) :: BILMAS,ATMOS,DIFVIT,INFOGR,CONVV(4),MSK
      LOGICAL, INTENT(IN) :: YASMH,ROVAR,PRECCU,VERTIC,ADJO,CORCON
      TYPE(SLVCFG), INTENT(INOUT)     :: SLVPRO
      CHARACTER(LEN=20),  INTENT(IN)  :: EQUA
      CHARACTER(LEN=*) ,  INTENT(IN)  :: COUPLING
C                                                           NPSMAX 
      INTEGER, INTENT(IN) :: NPSING(NWEIRS),NUMDIG(2,NWEIRS,*     )
C
C  STRUCTURES DE VECTEURS
C
      TYPE(BIEF_OBJ), INTENT(IN)    :: UCONV,VCONV,SMH,UN,VN,HN
      TYPE(BIEF_OBJ), INTENT(IN)    :: VOLU2D,V2DPAR,UNSV2D
      TYPE(BIEF_OBJ), INTENT(INOUT) :: RO,UDEL,VDEL,DM1,ZCONV,FLBOR
      TYPE(BIEF_OBJ), INTENT(IN)    :: UTILD,VTILD,PATMOS,CF
      TYPE(BIEF_OBJ), INTENT(INOUT) :: U,V,H,CV1,CV2,CV3,PRIVE,DH,DHN
      TYPE(BIEF_OBJ), INTENT(INOUT) :: CV1S
      TYPE(BIEF_OBJ), INTENT(INOUT) :: DU,DV,FU,FV,VISC,VISC_S,HTILD
      TYPE(BIEF_OBJ), INTENT(INOUT) :: UBOR,VBOR,HBOR,AUBOR,COTOND
      TYPE(BIEF_OBJ), INTENT(IN)    :: MASKEL,MASKPT,ZF
      TYPE(BIEF_OBJ), INTENT(IN)    :: HPROP,H0,C0,LIMPRO
C
C     TE : PAR ELEMENT              TE4,TE5 ONLY IF OPTBAN=3
      TYPE(BIEF_OBJ), INTENT(INOUT) :: TE1,TE2,TE3,TE4,TE5,ZFLATS
C     T  : PAR POINT
      TYPE(BIEF_OBJ), INTENT(INOUT) :: T1,T2,T3,T4,T5,T6,T7,T8
      TYPE(BIEF_OBJ), INTENT(INOUT) :: W1
C     STRUCTURE BIDON
      TYPE(BIEF_OBJ), INTENT(IN)    :: S
C
C-----------------------------------------------------------------------
C
C  STRUCTURES DE MATRICES
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: AM1,AM2,AM3,BM1,BM2,CM1,CM2,TM1
      TYPE(BIEF_OBJ), INTENT(INOUT) :: A23,A32,MBOR,BM1S,BM2S
C
C-----------------------------------------------------------------------
C
C  STRUCTURES DE BLOCS
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: MASK,MAT,RHS,UNK,TB,BD,DIRBOR
C
C-----------------------------------------------------------------------
C
C  STRUCTURE DE MAILLAGE
C
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C      
      INTEGER I,IELMU,IELMH,UDIR,UDDL,UNEU,HOND,UNONNEU,IELEM,NELEM
      INTEGER I1,I2,I3,DIMLIM,N,IOPT
C
      DOUBLE PRECISION Z(1),SL1,SL1U,C
C
      DOUBLE PRECISION P_DSUM
      EXTERNAL         P_DSUM
C
      LOGICAL MSKGRA
C
      CHARACTER*16 FORMUL     
C
C-----------------------------------------------------------------------
C
      DIMLIM=LIMPRO%DIM1
C
      IELMH=H%ELM
      IELMU=U%ELM
C
C  ADRESSES DES TABLEAUX DANS LE BLOC DES MASQUES MASK
C
      UDIR = 1
C     VDIR = 2
      UDDL = 3
C     VDDL = 4
      UNEU = 5
C     VNEU = 6
      HOND = 7
      UNONNEU = 8
C
C-----------------------------------------------------------------------
C
C  CONDITIONS AUX LIMITES DIRICHLET POUR L'ACCROISSEMENT DE H
C
C  HBOR = HBOR - HN (HBOR SUR LE BORD, HN SUR LE DOMAINE)
C
      CALL OSBD( 'X=X-Y   ' , HBOR , HN , HN , C , MESH )
C
C=======================================================================
C
C   MATRICES DE TYPE GRADIENT DE L'EQUATION DE CONTINUITE :
C
C    BM1 = - TETAU ( D(HN.F1)/DX * F2)
C    BM2 = - TETAU ( D(HN.F1)/DY * F2)
C
      IF(OPTBAN.EQ.3) THEN
C       PRISE EN COMPTE DE LA POROSITE
        CALL MATRIX(BM1,'M=N     ','MATFGR         X',IELMH,IELMU,
     *              TETAU,HPROP,S,S,S,S,S,MESH,.TRUE.,TE5)
        CALL MATRIX(BM2,'M=N     ','MATFGR         Y',IELMH,IELMU,
     *              TETAU,HPROP,S,S,S,S,S,MESH,.TRUE.,TE5)
C
C       MATRICES ISSUES DE SUPG APPLIQUE AU TERME DE CONVECTION
C
        IF(OPTSUP(2).EQ.1) THEN
          CALL KSUPG(TE1,TE2,1.D0,UCONV,VCONV,MESH)
          CALL MATRIX(BM1,'M=M+N   ','MATUGH         X',IELMH,IELMU,
     *                TETAU,HPROP,S,S,TE1,TE2,S,MESH,.TRUE.,TE5)
          CALL MATRIX(BM2,'M=M+N   ','MATUGH         Y',IELMH,IELMU,
     *                TETAU,HPROP,S,S,TE1,TE2,S,MESH,.TRUE.,TE5)
        ELSEIF(OPTSUP(2).EQ.2) THEN
          CALL MATRIX(BM1,'M=M+N   ','MATUGH         X',IELMH,IELMU,
     *                0.5*DT*TETAU,HPROP,S,S,UCONV,VCONV,S,
     *                MESH,.TRUE.,TE5)
          CALL MATRIX(BM2,'M=M+N   ','MATUGH         Y',IELMH,IELMU,
     *                0.5*DT*TETAU,HPROP,S,S,UCONV,VCONV,S,
     *                MESH,.TRUE.,TE5)
        ENDIF
C
      ELSE
C
        CALL MATRIX(BM1,'M=N     ','MATFGR         X',IELMH,IELMU,
     *              TETAU,HPROP,S,S,S,S,S,MESH,MSK,MASKEL)
        CALL MATRIX(BM2,'M=N     ','MATFGR         Y',IELMH,IELMU,
     *              TETAU,HPROP,S,S,S,S,S,MESH,MSK,MASKEL)
C
C       MATRICES ISSUES DE SUPG APPLIQUE AU TERME DE CONVECTION
C
        IF(OPTSUP(2).EQ.1) THEN
          CALL KSUPG(TE1,TE2,1.D0,UCONV,VCONV,MESH)
          CALL MATRIX(BM1,'M=M+N   ','MATUGH         X',IELMH,IELMU,
     *                TETAU,HPROP,S,S,TE1,TE2,S,MESH,MSK,MASKEL)
          CALL MATRIX(BM2,'M=M+N   ','MATUGH         Y',IELMH,IELMU,
     *                TETAU,HPROP,S,S,TE1,TE2,S,MESH,MSK,MASKEL)
        ELSEIF(OPTSUP(2).EQ.2) THEN
          CALL MATRIX(BM1,'M=M+N   ','MATUGH         X',IELMH,IELMU,
     *                0.5*DT*TETAU,HPROP,S,S,UCONV,VCONV,S,
     *                MESH,MSK,MASKEL)
          CALL MATRIX(BM2,'M=M+N   ','MATUGH         Y',IELMH,IELMU,
     *                0.5*DT*TETAU,HPROP,S,S,UCONV,VCONV,S,
     *                MESH,MSK,MASKEL)
        ENDIF
C
      ENDIF
C
C=======================================================================
C
C   CONSTRUCTION DE LA MATRICE DE DIFFUSION :
C
C    TM1 =  VISC . (P03 * (DF1/DX * DF2/DX + DF1/DY * DF2/DY) )
C
      IF(DIFVIT) THEN
C
        IF(OPDVIT.EQ.2) THEN
C         SAUVEGARDE DE LA DIFFUSION
          CALL OS('X=Y     ',VISC_S,VISC,VISC,C)
C         MULTIPLICATION PAR HPROP DE LA DIFFUSION
          CALL OV_2('X=XY    ',VISC%R,1,HPROP%R,1,Z,1,C,
     *                         VISC%MAXDIM1,VISC%DIM1)
          IF(VISC%DIM2.EQ.3) THEN
            CALL OV_2('X=XY    ',VISC%R,2,HPROP%R,1,Z,1,C,
     *                           VISC%MAXDIM1,VISC%DIM1)
            CALL OV_2('X=XY    ',VISC%R,3,HPROP%R,1,Z,1,C,
     *                           VISC%MAXDIM1,VISC%DIM1)
          ENDIF
        ENDIF
C
        CALL MATRIX(TM1,'M=N     ','MATDIF          ',IELMU,IELMU,
     *              1.D0,S,S,S,VISC,S,S,MESH,MSK,MASKEL)
C
        IF(OPDVIT.EQ.2) THEN
C          MULTIPLICATION DE LA MATRICE PAR 1/HPROP
           CALL CPSTVC(HPROP,T4)
           DO I=1,HPROP%DIM1
C            CAUTION: HIDDEN PARAMETER 1.D-2, NO DIFFUSION BELOW 1 CM DEPTH
C                                             IN THIS OPTION
             IF(HPROP%R(I).GT.1.D-2) THEN
               T4%R(I)=1.D0/HPROP%R(I)
             ELSE
               T4%R(I)=0.D0
             ENDIF
           ENDDO
           IF(T4%ELM.NE.IELMU) CALL CHGDIS(T4,T4%ELM,IELMU,MESH)
           CALL OM( 'M=X(M)  ' , TM1 , TM1 , S  , C , MESH )
           CALL OM( 'M=DM    ' , TM1 , TM1 , T4 , C , MESH )
C          RESTITUTION DE LA DIFFUSION
           CALL OS('X=Y     ',VISC,VISC_S,VISC_S,C)
        ENDIF
C
C       IF ADDED ON 23/07/2002 BY JMH (MAY HAPPEN IN PARALLELISM)
        IF(MESH%NPTFR.GT.0) THEN
         CALL MATRIX(MBOR,'M=N     ','FMATMA          ',
     *               IELBOR(IELMU,1),IELBOR(IELMU,1),
     *               -1.D0,AUBOR,S,S,S,S,S,MESH,.TRUE.,MASK%ADR(UNEU)%P)
         CALL OM( 'M=M+N   ' , TM1 , MBOR , S , C , MESH )
        ENDIF
C
C       PARTIE EXPLICITE FAITE DANS LE PROCHAIN IF(DIFVIT...
C
      ENDIF
C
C=======================================================================
C
C  CALCUL DE LA SURFACE LIBRE (MISE DANS T8)
C
      CALL OS( 'X=Y+Z   ' , T8 , HN , ZF , C )
C
C  OPTION 1 DU TRAITEMENT DES BANCS DECOUVRANTS
C  ON UTILISE PARFOIS UN MASQUE POUR LES BANCS DECOUVRANTS
C  BIEN QUE CE NE SOIT PAS L'OPTION AVEC MASQUAGE.
C  CE MASQUE EST MIS DANS TE3
C
      IF(OPTBAN.EQ.1.OR.OPTBAN.EQ.3) THEN
        CALL DECVRT(TE3,T8,ZF,MESH)
      ENDIF
C
C     GRADIENT DE SURFACE LIBRE
C
      IF (OPTBAN.EQ.1.OR.OPTBAN.EQ.3) THEN
C                   SL CORRIGE PAR ELEMENT
        CALL CORRSL(ZFLATS,T8,ZF,MESH)
        CALL VECTOR(CV2,'=','GRADF          X',IELMU,
     *              -GRAV,ZFLATS,S,S,S,S,S,MESH,MSK,MASKEL)
        CALL VECTOR(CV3,'=','GRADF          Y',IELMU,
     *              -GRAV,ZFLATS,S,S,S,S,S,MESH,MSK,MASKEL)
C       CORRSL DECLARE L'ELEMENT COMME DISCONTINU, ON REMET COMME AVANT
      ELSE
        CALL VECTOR(CV2,'=','GRADF          X',IELMU,
     *              -GRAV,T8,S,S,S,S,S,MESH,MSK,MASKEL)
        CALL VECTOR(CV3,'=','GRADF          Y',IELMU,
     *              -GRAV,T8,S,S,S,S,S,MESH,MSK,MASKEL)
      ENDIF
C
C  TERMES DE GRADIENT SUPPLEMENTAIRES : PRESSION ATMOSPHERIQUE
C                                       DENSITE VARIABLE
C
C  CES TERMES MOTEURS NE DOIVENT PAS ETRE AJOUTES DANS LES BANCS
C  DECOUVRANTS.
C
      IF(ROVAR.OR.ATMOS) THEN
C
C  CAS OU IL FAUT MASQUER CES GRADIENTS : MASQUAGE DEMANDE OU
C                                         OPTBAN=1
C
      IF(MSK.OR.OPTBAN.EQ.1.OR.OPTBAN.EQ.3) THEN
        MSKGRA = .TRUE.
C       SI OPTBAN<>1 IL FAUT REMPLIR LE MASQUE QUE L'ON UTILISE ICI
        IF(OPTBAN.NE.1.AND.OPTBAN.NE.3) THEN
         CALL OV('X=Y     ',TE3%R,MASKEL%R,MASKEL%R,C,TE3%DIM1)
        ENDIF
      ELSE
        MSKGRA = .FALSE.
      ENDIF
C
C     GRADIENT DE PRESSION ATMOSPHERIQUE
C
      IF(ATMOS) THEN
        CALL VECTOR(CV2,'+','GRADF          X',IELMU,
     *              -1.D0/ROEAU,PATMOS,S,S,S,S,S,MESH,MSKGRA,TE3)
        CALL VECTOR(CV3,'+','GRADF          Y',IELMU,
     *              -1.D0/ROEAU,PATMOS,S,S,S,S,S,MESH,MSKGRA,TE3)
      ENDIF
C
C     TERMES SUPPLEMENTAIRES SI LA MASSE VOLUMIQUE EST VARIABLE
C
      IF(ROVAR) THEN
C
        CALL OS( 'X=X+C   ' , RO , S , S , -ROEAU )
C
C       PRESSION BAROCLINE
C
        CALL VECTOR(CV2,'+','GGRADF         X',IELMU,
     *              -0.5D0*GRAV/ROEAU,RO,HN,S,S,S,S,MESH,MSKGRA,TE3)
        CALL VECTOR(CV3,'+','GGRADF         Y',IELMU,
     *              -0.5D0*GRAV/ROEAU,RO,HN,S,S,S,S,MESH,MSKGRA,TE3)
        CALL OS( 'X=X+C   ' , RO , S , S , +ROEAU )
C
      ENDIF
C
      ENDIF
C
C=======================================================================
C
C   MATRICE DE MASSE / DT
C   AM1 SERA MODIFIEE PLUS TARD SI ON UTILISE SUPG SUR H
C
C    AM1 = ( 1 / DT )  * (F1 * F2)
C
      SL1 = 1.D0 / DT
      IF(OPTBAN.EQ.1.OR.OPTBAN.EQ.3) THEN
C       MATRICE DE MASSE LUMPEE LOCALEMENT SUR LES BANCS DECOUVRANTS
        FORMUL='MSLUMP          '
      ELSE
C       MATRICE DE MASSE NORMALE
        FORMUL='MATMAS          '
      ENDIF
      IF(OPTBAN.NE.3) THEN
        CALL MATRIX(AM1,'M=N     ',FORMUL,IELMH,IELMH,
     *              SL1,TE3,S,S,S,S,S,MESH,MSK,MASKEL)
      ELSE
        CALL MATRIX(AM1,'M=N     ',FORMUL,IELMH,IELMH,
     *              SL1,TE3,S,S,S,S,S,MESH,.TRUE.,TE5)
      ENDIF
C
C   MATRICE DE MASSE POUR L'EQUATION DE QTE DE MOUVEMENT
C
      IF(SOLSYS.EQ.1) THEN
C
C                           OPTBAN.NE.3 TO AVOID POROSITY IN AM2
      IF(IELMU.EQ.IELMH.AND.OPTBAN.NE.3) THEN
        CALL OM( 'M=N     ' , AM2 , AM1 , S , C , MESH )
      ELSE
        CALL MATRIX(AM2,'M=N     ',FORMUL,IELMU,IELMU,
     *              SL1,TE3,S,S,S,S,S,MESH,MSK,MASKEL)
      ENDIF
C     MASS-LUMPING DE AM2 :
      IF(AGGLOU.GT.0.001D0) THEN
        CALL LUMP(T2,AM2,MESH,AGGLOU)
        CALL OM( 'M=CN    ' , AM2 , AM2 , S  , 1.D0-AGGLOU , MESH )
        CALL OM( 'M=M+D   ' , AM2 , AM2 , T2 , C           , MESH )
      ENDIF
C
      ELSEIF(SOLSYS.EQ.2) THEN
        IF(IELMU.NE.IELMH) THEN
          CALL VECTOR(AM2%D,'=','MASBAS          ',IELMU,
     *                SL1,S,S,S,S,S,S,MESH,MSK,MASKEL)
        ELSE
          CALL OS('X=CY    ',X=AM2%D,Y=VOLU2D,C=SL1)
        ENDIF
        AM2%TYPDIA='Q'
        AM2%TYPEXT='0'
      ENDIF
C
C MASS-LUMPING DE AM1 :
C
      IF(AGGLOH.GT.0.001D0) THEN
        CALL LUMP(T1,AM1,MESH,AGGLOH)
        CALL OM( 'M=CN    ' , AM1 , AM1 , S , 1.D0-AGGLOH , MESH )
        CALL OM( 'M=M+D   ' , AM1 , AM1 , T1 , C          , MESH )
      ENDIF
C
C FIN DU MASS-LUMPING
C
C STOCKAGE PROVISOIRE DE AM2 MATRICE DE MASSE LUMPEE DANS AM3
C POUR LE CALCUL DES TERMES PROVENANT DES DERIVEES EN TEMPS
C
      IF(SOLSYS.EQ.1) THEN
        CALL OM( 'M=N     ' , AM3 , AM2 , S , C , MESH )
      ELSEIF(SOLSYS.EQ.2) THEN
        CALL OS('X=Y     ',X=AM3%D,Y=AM2%D)
        AM3%TYPDIA=AM2%TYPDIA
        AM3%TYPEXT=AM2%TYPEXT
      ENDIF
C
C=======================================================================
C
C   CALCUL DES VECTEURS SECONDS MEMBRES
C
C   ATTENTION : FAIRE TRES ATTENTION AU PARAMETRE LEGO (FALSE OU TRUE)
C               DANS LES APPELS A MATVEC ,VGRADF. SI LEGO = .FALSE. LE
C               RESULTAT DU CALCUL EST DANS W1. DANS CE CAS IL
C               FAUT TOUJOURS ENSUITE UTILISER DES OPERATIONS 'X=X+...'
C               POUR NE PAS EFFACER CE QU'IL Y A DANS W1 , ET
C               TERMINER AVEC LEGO=.TRUE. POUR L'ASSEMBLAGE FINAL.
C
C-----------------------------------------------------------------------
C
C     CV1 = CV1 + SL1U * ( BM1 * UN + BM2 * VN )
C
      SL1U   = (TETAU-1.D0)/TETAU
C                           TETAU : BM1 A ETE CONSTRUITE AVEC CE COEFF
      IF(ABS(SL1U).GT.0.0001D0) THEN
        CALL MATVEC('X=CAY   ',CV1,BM1,UN,SL1U,MESH)
        CALL MATVEC('X=X+CAY ',CV1,BM2,VN,SL1U,MESH)
      ELSE
        CALL CPSTVC(H,CV1)
        CALL OS('X=0     ',X=CV1)
      ENDIF
C
C     TERMES SOURCES DE L'EQUATION DE CONTINUITE :
C
      MASSES = 0.D0
      IF(YASMH) THEN
        IF(OPTSOU.EQ.1) THEN
C         VERSION NORMALE
C         NOTE JMH : ET LE MASS-LUMPING ?? UTILISER DT*AM1*SMH ??
C         EN TENANT COMPTE DU MASS-LUMPING
          CALL MATVEC( 'X=CAY   ',T1,AM1,SMH,DT,MESH)
C         SANS TENIR COMPTE DU MASS-LUMPING
!         CALL VECTOR(T1,'=','MASVEC          ',IELMH,
!    *                1.D0,SMH,S,S,S,S,S,MESH,MSK,MASKEL)
          CALL OS( 'X=X+Y   ' , CV1 , T1 , S , C )
          IF(BILMAS) MASSES = DT * BIEF_SUM(T1)
        ELSE
C         VERSION DIRAC
          CALL OS( 'X=X+Y   ' ,X=CV1,Y=SMH)
          IF(BILMAS) MASSES = DT * BIEF_SUM(SMH)
        ENDIF
C       LA LIGNE SUIVANTE EST MISE DANS BILAN
C       IF(NCSIZE.GT.1) MASSES=P_DSUM(MASSES)
      ENDIF
C
C-----------------------------------------------------------------------
C
C     T1 = FU * DT + UTILD
C
      IF(ICONVF(1).NE.ADV_LPO.AND.ICONVF(1).NE.ADV_LPO_TF.AND.
     *   ICONVF(1).NE.ADV_NSC.AND.ICONVF(1).NE.ADV_NSC_TF.AND.
     *   ICONVF(1).NE.ADV_PSI.AND.ICONVF(1).NE.ADV_PSI_TF     ) THEN
!     IF(ICONVF(1).NE.8) THEN
        CALL OS( 'X=CY    ' , T1 , FU , S , DT )
      ENDIF
C
      IF(ICONVF(1).EQ.ADV_CAR.OR.(.NOT.CONVV(1))) THEN
C
C       CONVECTION FORTE
C
        CALL OS( 'X=X+Y   ' , T1 , UTILD , S , C )
C
C------ SCHEMA SEMI-IMPLICITE CENTRE + S.U.P.G -----------------------
C
      ELSEIF(ICONVF(1).EQ.ADV_SUP) THEN
C
C TERME DE CONVECTION CENTRE SEMI-IMPLICITE : MATRICE
C
        CALL MATRIX(CM2,'M=N     ','MATVGR          ',IELMU,IELMU,
     *              1.D0,S,S,S,UCONV,VCONV,VCONV,MESH,MSK,MASKEL)
C
C       CONTRIBUTION DE SUPG
C
        IF(OPTSUP(1).EQ.1) THEN
C         SUPG CLASSIQUE
          CALL KSUPG(TE1,TE2,1.D0,UCONV,VCONV,MESH)
          CALL MATRIX(CM2,'M=M+N   ','MASUPG          ',IELMU,IELMU,
     *                1.D0,TE1,TE2,S,UCONV,VCONV,VCONV,
     *                MESH,MSK,MASKEL)
        ELSEIF(OPTSUP(1).EQ.2) THEN
C         SUPG MODIFIE
          CALL MATRIX(CM2,'M=M+N   ','MAUGUG          ',IELMU,IELMU,
     *                0.5D0*DT,S,S,S,UCONV,VCONV,VCONV,
     *                MESH,MSK,MASKEL)
        ENDIF
C       FIN DE LA CONTRIBUTION DE SUPG
C
C SECOND MEMBRE EXPLICITE
C
        IF(ABS(SL1U).GT.0.0001D0) THEN
          CALL MATVEC( 'X=X+CAY ',CV2,CM2,UN,TETAU-1.D0,MESH)
        ENDIF
C
C MATRICE : ON DONNE A AM2 UNE STRUCTURE NON SYMETRIQUE
C
        IF(AM2%TYPEXT.NE.'Q') THEN
          CALL OM( 'M=X(M)  ' , AM2 , AM2 , S , C , MESH )
        ENDIF
        CALL OM( 'M=M+CN  ' , AM2 , CM2 , S , TETAU , MESH )
C
        CALL OS( 'X=X+Y   ' , T1 , UN , S , C )
C
C-----------------------------------------------------------------
C
C  SCHEMA PSI
C
      ELSEIF(ICONVF(1).EQ.ADV_PSI_NC) THEN
C
C DEBUT DU CALCUL DU SECOND MEMBRE CV2
C
        CALL VGFPSI(T6,IELMU,UCONV,VCONV,UN,DT,-1.D0,CFLMAX,
     *              T4,T5,MESH,MSK,MASKEL)
        CALL OS( 'X=X+Y   ' , CV2 , T6 , S , C )
        CALL OS( 'X=X+Y   ' , T1  , UN , S , C )
C
C------ SCHEMA SEMI-IMPLICITE N --------------------------------------
C
      ELSEIF(ICONVF(1).EQ.ADV_NSC_NC) THEN
C
C TERME DE CONVECTION CENTRE SEMI-IMPLICITE : MATRICE
C
        CALL MATRIX(CM2,'M=N     ','MATVGR         N',IELMU,IELMU,
     *              1.D0,S,S,S,UCONV,VCONV,VCONV,MESH,MSK,MASKEL)
C
C SECOND MEMBRE EXPLICITE
C
        IF(ABS(SL1U).GT.0.0001D0) THEN
          CALL MATVEC( 'X=X+CAY ',CV2,CM2,UN,TETAU-1.D0,MESH)
        ENDIF
C
C MATRICE : ON DONNE A AM2 UNE STRUCTURE NON SYMETRIQUE
C
        IF(AM2%TYPEXT.NE.'Q') THEN
          CALL OM( 'M=X(M)  ' , AM2 , AM2 , S , C , MESH )
        ENDIF
        CALL OM( 'M=M+CN  ' , AM2 , CM2 , S , TETAU , MESH )
        CALL OS( 'X=X+Y   ' , T1 , UN , S , C )
C
C------ SCHEMA VOLUMES FINIS --------------------------------------
C
C     NOTE: HERE THE CONTINUITY EQUATION IS NOT SOLVED
C           BY H, HN AND UCONV,VCONV (UCONV, VCONV HAVE BEEN
C           UPDATED AND H, HN ARE STILL UNCHANGED, SO THE FINAL H
C           COMPUTED IN CVTRVF WILL NOT BE THE FINAL DEPTH OF THE
C           TIME-STEP).           
C
      ELSEIF(ICONVF(1).EQ.ADV_LPO.OR.
     *       ICONVF(1).EQ.ADV_NSC.OR.
     *       ICONVF(1).EQ.ADV_PSI     ) THEN
C
        IF(ICONVF(1).EQ.ADV_LPO.OR.ICONVF(1).EQ.ADV_NSC) IOPT=2
        IF(ICONVF(1).EQ.ADV_PSI) IOPT=3
C       ICI YASMH=.FALSE. (SOURCES PRISES EN COMPTE DANS FU)
        IF(TB%N.LT.22) THEN
          WRITE(LU,*) 'SIZE OF TB TOO SMALL IN PROPAG'
          CALL PLANTE(1)
          STOP
        ENDIF
        CALL CVTRVF(T1,UN,S,.FALSE.,.TRUE.,H,HN,
     *              HPROP,UCONV,VCONV,S,S,
     *              1,S,S,FU,S,.FALSE.,S,
     *              .FALSE.,UBOR,MASK,MESH,
     *              TB%ADR(13)%P,TB%ADR(14)%P,TB%ADR(15)%P,
     *              TB%ADR(16)%P,TB%ADR(17)%P,TB%ADR(18)%P,
     *              TB%ADR(19)%P,TB%ADR(20)%P,TB%ADR(21)%P, 
     *              TB%ADR(22)%P,    
     *              AGGLOH,TE1,DT,INFOGR,BILMAS,
     *              1,MSK,MASKEL,S,C,1,
     *              LIMPRO%I(1+DIMLIM:2*DIMLIM),
     *              KDIR,3,MESH%NPTFR,FLBOR,.FALSE.,
     *              V2DPAR,UNSV2D,IOPT,TB%ADR(12)%P,MASKPT)
        IF(IELMU.NE.11) CALL CHGDIS(T1,11,IELMU,MESH)
C
C------ SCHEMA VOLUMES FINIS AVEC BANCS DECOUVRANTS -------------------
C
C     NOTE: HERE THE CONTINUITY EQUATION IS NOT SOLVED
C           BY H, HN AND UCONV,VCONV (UCONV, VCONV HAVE BEEN
C           UPDATED AND H, HN ARE STILL UNCHANGED, SO THE FINAL H
C           COMPUTED IN CVTRVF WILL NOT BE THE FINAL DEPTH OF THE
C           TIME-STEP).
C
      ELSEIF(ICONVF(1).EQ.ADV_LPO_TF.OR.
     *       ICONVF(1).EQ.ADV_NSC_TF.OR.
     *       ICONVF(1).EQ.ADV_PSI_TF     ) THEN
C
        IF(ICONVF(1).EQ.ADV_LPO_TF) IOPT=2
        IF(ICONVF(1).EQ.ADV_NSC_TF) IOPT=2
        IF(ICONVF(1).EQ.ADV_PSI_TF) IOPT=3
C       ICI YASMH=.FALSE. (SOURCES PRISES EN COMPTE DANS FU)
        IF(TB%N.LT.22) THEN
          WRITE(LU,*) 'SIZE OF TB TOO SMALL IN PROPAG'
          CALL PLANTE(1)
          STOP
        ENDIF
        CALL CVTRVF_POS(T1,UN,S,.FALSE.,.TRUE.,H,HN,
     *                  HPROP,UCONV,VCONV,S,S,
     *                  1,S,S,FU,S,.FALSE.,S,.FALSE.,UBOR,MASK,MESH,
     *                  TB%ADR(13)%P,TB%ADR(14)%P,TB%ADR(15)%P,
     *                  TB%ADR(16)%P,TB%ADR(17)%P,TB%ADR(18)%P,
     *                  TB%ADR(19)%P,TB%ADR(20)%P,TB%ADR(21)%P, 
     *                  TB%ADR(22)%P,    
     *                  AGGLOH,TE1,DT,INFOGR,BILMAS,1,MSK,MASKEL,S,C,1,
     *                  LIMPRO%I(1+DIMLIM:2*DIMLIM),
     *                  KDIR,3,MESH%NPTFR,FLBOR,.FALSE.,
     *                  V2DPAR,UNSV2D,IOPT,TB%ADR(12)%P,MASKPT,
     *       MESH%GLOSEG%I(                 1:  MESH%GLOSEG%DIM1),
     *       MESH%GLOSEG%I(MESH%GLOSEG%DIM1+1:2*MESH%GLOSEG%DIM1),
     *       MESH%NBOR%I)
        IF(IELMU.NE.11) CALL CHGDIS(T1,11,IELMU,MESH)
C
C-----------------------------------------------------------------
C
      ELSE
C
        IF(LNG.EQ.1) WRITE(LU,2002) ICONVF(1)
        IF(LNG.EQ.2) WRITE(LU,2003) ICONVF(1)
2002    FORMAT(1X,'PROPAG : FORME DE LA CONVECTION DE U INCONNUE :',1I6)
2003    FORMAT(1X,'PROPAG : UNKNOWN TYPE OF ADVECTION FOR U: ',1I6)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C AJOUT DU TERME AM3 * T1
C ICI AM3 DOIT ETRE MASSE/DT AVEC LE MASS-LUMPING SUR U EVENTUELLEMENT
C
      IF(SOLSYS.EQ.1) THEN
        CALL MATVEC( 'X=X+AY  ',CV2,AM3,T1,C,MESH)
      ELSEIF(SOLSYS.EQ.2) THEN
        CALL OS('X=X+YZ  ',X=CV2,Y=AM3%D,Z=T1)
      ENDIF
C
C-----------------------------------------------------------------------
C
C     T1 = FV * DT + VTILD
C
      IF(ICONVF(1).NE.ADV_LPO.AND.ICONVF(1).NE.ADV_LPO_TF.AND.
     *   ICONVF(1).NE.ADV_NSC.AND.ICONVF(1).NE.ADV_NSC_TF.AND.
     *   ICONVF(1).NE.ADV_PSI.AND.ICONVF(1).NE.ADV_PSI_TF     ) THEN
        CALL OS( 'X=CY    ' , T1 , FV , S , DT )
      ENDIF
C
      IF(ICONVF(1).EQ.ADV_CAR.OR.(.NOT.CONVV(1))) THEN
C
C       CONVECTION FORTE
C
        CALL OS( 'X=X+Y   ' , T1 , VTILD , S , C )
C
C------ SCHEMA SEMI-IMPLICITE CENTRE + S.U.P.G -----------------------
C
      ELSEIF(ICONVF(1).EQ.ADV_SUP) THEN
C
C SECOND MEMBRE EXPLICITE
C
C       CM2 DEJA CALCULEE POUR U
        IF(ABS(SL1U).GT.0.0001D0) THEN
          CALL MATVEC( 'X=X+CAY ',CV3,CM2,VN,TETAU-1.D0,MESH)
        ENDIF
C
        CALL OS( 'X=X+Y   ' , T1 , VN , S , C )
C
C-----------------------------------------------------------------
C
C  SCHEMA PSI
C
      ELSEIF(ICONVF(1).EQ.ADV_PSI_NC) THEN
C
C DEBUT DU CALCUL DU SECOND MEMBRE CV3
C
        CALL VGFPSI(T6,IELMU,UCONV,VCONV,VN,DT,-1.D0,CFLMAX,
     *              T4,T5,MESH,MSK,MASKEL)
        CALL OS( 'X=X+Y   ' , CV3 , T6 , S , C )
        CALL OS( 'X=X+Y   ' , T1  , VN , S , C )
C
C------ SCHEMA SEMI-IMPLICITE N --------------------------------------
C
      ELSEIF(ICONVF(1).EQ.ADV_NSC_NC) THEN
C
C SECOND MEMBRE EXPLICITE
C
C       CM2 DEJA CALCULEE POUR U
        IF(ABS(SL1U).GT.0.0001D0) THEN
          CALL MATVEC( 'X=X+CAY ',CV3,CM2,VN,TETAU-1.D0,MESH)
        ENDIF
C
        CALL OS( 'X=X+Y   ' , T1 , VN , S , C )
C
C------ SCHEMA VOLUMES FINIS -----------------------------------------
C
      ELSEIF(ICONVF(1).EQ.ADV_LPO.OR.
     *       ICONVF(1).EQ.ADV_NSC.OR.
     *       ICONVF(1).EQ.ADV_PSI     ) THEN
C
        IF(ICONVF(1).EQ.ADV_LPO.OR.ICONVF(1).EQ.ADV_NSC) IOPT=2
        IF(ICONVF(1).EQ.ADV_PSI) IOPT=3
C
        CALL CVTRVF(T1,VN,S,.FALSE.,.TRUE.,H,HN,
     *              HPROP,UCONV,VCONV,S,S,
     *              1,S,S,FV,S,.FALSE.,S,
     *              .FALSE.,VBOR,MASK,MESH,
     *              TB%ADR(13)%P,TB%ADR(14)%P,TB%ADR(15)%P,
     *              TB%ADR(16)%P,TB%ADR(17)%P,TB%ADR(18)%P,
     *              TB%ADR(19)%P,TB%ADR(20)%P,TB%ADR(21)%P,
     *              TB%ADR(22)%P,     
     *              AGGLOH,TE1,DT,INFOGR,BILMAS,
     *              1,MSK,MASKEL,S,C,1,
     *              LIMPRO%I(1+2*DIMLIM:3*DIMLIM),
     *              KDIR,3,MESH%NPTFR,FLBOR,.FALSE.,
     *              V2DPAR,UNSV2D,IOPT,TB%ADR(12)%P,MASKPT)
        IF(IELMU.NE.11) CALL CHGDIS(T1,11,IELMU,MESH)
C
C------ SCHEMA VOLUMES FINIS AVEC BANCS DECOUVRANTS ------------------
C
      ELSEIF(ICONVF(1).EQ.ADV_LPO_TF.OR.
     *       ICONVF(1).EQ.ADV_NSC_TF.OR.
     *       ICONVF(1).EQ.ADV_PSI_TF     ) THEN
C
        IF(ICONVF(1).EQ.ADV_LPO_TF) IOPT=2
        IF(ICONVF(1).EQ.ADV_NSC_TF) IOPT=2
        IF(ICONVF(1).EQ.ADV_PSI_TF) IOPT=3
C
C       ICI YASMH=.FALSE. (SOURCES PRISES EN COMPTE DANS FU)
        IF(TB%N.LT.22) THEN
          WRITE(LU,*) 'SIZE OF TB TOO SMALL IN PROPAG'
          CALL PLANTE(1)
          STOP
        ENDIF
        CALL CVTRVF_POS(T1,VN,S,.FALSE.,.TRUE.,H,HN,
     *                  HPROP,UCONV,VCONV,S,S,1,S,S,FV,S,.FALSE.,S,
     *                  .FALSE.,VBOR,MASK,MESH,
     *                  TB%ADR(13)%P,TB%ADR(14)%P,TB%ADR(15)%P,
     *                  TB%ADR(16)%P,TB%ADR(17)%P,TB%ADR(18)%P,
     *                  TB%ADR(19)%P,TB%ADR(20)%P,TB%ADR(21)%P, 
     *                  TB%ADR(22)%P,    
     *                  AGGLOH,TE1,DT,INFOGR,BILMAS,
     *                  1,MSK,MASKEL,S,C,1,
     *                  LIMPRO%I(1+2*DIMLIM:3*DIMLIM),
     *                  KDIR,3,MESH%NPTFR,FLBOR,.FALSE.,
     *                  V2DPAR,UNSV2D,IOPT,TB%ADR(12)%P,MASKPT,
     *       MESH%GLOSEG%I(                 1:  MESH%GLOSEG%DIM1),
     *       MESH%GLOSEG%I(MESH%GLOSEG%DIM1+1:2*MESH%GLOSEG%DIM1),
     *       MESH%NBOR%I)
        IF(IELMU.NE.11) CALL CHGDIS(T1,11,IELMU,MESH)
C
C-----------------------------------------------------------------
C
      ENDIF
C
C AJOUT DU TERME AM3 * T1, PUIS ASSEMBLAGE : RESULTAT DANS T2
C ICI AM3 DOIT ETRE MASSE/DT AVEC LE MASS-LUMPING SUR U EVENTUELLEMENT
C
      IF(SOLSYS.EQ.1) THEN
        CALL MATVEC( 'X=X+AY  ',CV3,AM3,T1,C,MESH)
      ELSEIF(SOLSYS.EQ.2) THEN
        CALL OS('X=X+YZ  ',X=CV3,Y=AM3%D,Z=T1)
      ENDIF
C
C=======================================================================
C
C
C     CALCUL DES MATRICES DIAGONALES : - FU*(MASSE DES BASES)
C                                 ET : - FV*(MASSE DES BASES)
C
      IF(KFROT.NE.0.OR.VERTIC) THEN
C
C       T3,T4 : BOTTOM FRICTION      T5,T6 : VERTICAL STRUCTURES
C
        CALL FRICTI(T3,T4,T5,T6,UN,VN,HN,CF,MESH,T1,T2,VERTIC,
     *              UNSV2D,MSK,MASKEL,HFROT)
C
C       CALCUL DES MATRICES DIAGONALES : - FU*(MASSE DES BASES)
C                                   ET : - FV*(MASSE DES BASES)
C
        CALL SLOPES(TE3,ZF,MESH)
        CALL VECTOR(T2,'=','MASBAS          ',IELMU,
     *              -1.D0,S,S,S,S,S,S,MESH,.TRUE.,TE3)
C
        CALL OS( 'X=XY    ' , X=T3 , Y=T2 )
        CALL OS( 'X=XY    ' , X=T4 , Y=T2 )        
C
        IF(VERTIC) THEN
          CALL VECTOR(T2,'=','MASBAS          ',IELMU,
     *              -1.D0,S,S,S,S,S,S,MESH,.FALSE.,TE3)
          CALL OS( 'X=XY    ' , T5 , T2 , S , C )
          CALL OS( 'X=XY    ' , T6 , T2 , S , C )
          CALL OS( 'X=X+Y   ' , T3 , T5 , S , C )
          CALL OS( 'X=X+Y   ' , T4 , T6 , S , C )
        ENDIF
C
      ELSE
C
        CALL CPSTVC(U,T3)
        CALL CPSTVC(V,T4)
        CALL OS( 'X=0     ' , X=T3 )
        CALL OS( 'X=0     ' , X=T4 )
C
        IF(OPTBAN.EQ.1) THEN
C         SLOWING DOWN VELOCITIES ON TIDAL FLATS
C         HAPPENS WHEN CALLED BY TELEMAC-3D WITH OPTT2D=1
          IF(IELMU.NE.IELMH) THEN
            CALL VECTOR(T2,'=','MASBAS          ',IELMU,
     *                  1.D0,S,S,S,S,S,S,MESH,.FALSE.,TE3)
            IF(NCSIZE.GT.1) CALL PARCOM(T2,2,MESH)
            CALL OS('X=Y     ',X=T5,Y=HPROP)
            IF(T5%ELM.NE.IELMU) CALL CHGDIS(T5,T5%ELM,IELMU,MESH)
            DO I=1,U%DIM1
C             HIDDEN PARAMETER
              IF(T5%R(I).LT.1.D-3) THEN
                T3%R(I)=10.D0*T2%R(I)/DT
                T4%R(I)=T3%R(I)
              ENDIF
            ENDDO
          ELSE
            DO I=1,U%DIM1
C             HIDDEN PARAMETER
              IF(HPROP%R(I).LT.1.D-3) THEN
                T3%R(I)=10.D0*V2DPAR%R(I)/DT
                T4%R(I)=T3%R(I)
              ENDIF
            ENDDO
          ENDIF
        ENDIF
C
      ENDIF
C
C=======================================================================
C
C   CALCUL DES MATRICES
C
C    AM1: DEJA CALCULEE
C
C    AM2 = AM2 + TM1
C
      IF(DIFVIT) THEN
C
C    TEST : ESSAI D'IMPLICITATION DE LA DIAGONALE DE TM1 EN EQUATION
C           D'ONDES.
C
        IF(SOLSYS.EQ.2) THEN
C
          CALL OS('X=X+CY  ',X=AM2%D,Y=TM1%D,C=TETAD)
          CALL OS('X=CX    ',X=TM1%D,C=1.D0-TETAD)
C
          CALL MATVEC( 'X=X+CAY ',CV2,TM1,UN,-1.D0,MESH)
          CALL MATVEC( 'X=X+CAY ',CV3,TM1,VN,-1.D0,MESH)
C      
        ELSEIF(SOLSYS.EQ.1) THEN
C
          IF(TETAD.LT.0.9999D0) THEN
            IF(TETAD.GT.0.0001D0) THEN
              IF(AM2%TYPEXT.EQ.'S'.AND.TM1%TYPEXT.EQ.'Q') THEN
                CALL OM( 'M=X(M)  ' , AM2 , AM2 , S , C , MESH )
              ENDIF
              CALL OM( 'M=M+CN  ' , AM2 , TM1 , S , TETAD , MESH )
            ENDIF
C           PARTIE EXPLICITE :
            CALL MATVEC( 'X=X+CAY ',CV2,TM1,UN,TETAD-1.D0,MESH)
            CALL MATVEC( 'X=X+CAY ',CV3,TM1,VN,TETAD-1.D0,MESH)
          ELSE
C           ENTIEREMENT IMPLICITE
            IF(AM2%TYPEXT.EQ.'S'.AND.TM1%TYPEXT.EQ.'Q') THEN
              CALL OM( 'M=X(M)  ' , AM2 , AM2 , S , C , MESH )
            ENDIF
            CALL OM( 'M=M+N   ' , AM2 , TM1 , S , C , MESH )
          ENDIF
C
        ELSE
C
          IF(LNG.EQ.1) WRITE(LU,*) 'PROPAG : MAUVAIS CHOIX POUR SOLSYS'
          IF(LNG.EQ.2) WRITE(LU,*) 'PROPAG : WRONG CHOICE FOR SOLSYS'
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
      ENDIF
C
C    AM3 = AM2 (DIFFUSION HAS BEEN ADDED TO AM2, NOT TO AM3)
C
      IF(SOLSYS.EQ.1) THEN
        CALL OM( 'M=N     ' , AM3 , AM2 , S , C , MESH )
      ELSEIF(SOLSYS.EQ.2) THEN
        CALL OS('X=Y     ',X=AM3%D,Y=AM2%D)
      ENDIF
C
C=======================================================================
C
C   DEFINA METHOD CORRECTED : RIGHT HAND SIDE MODIFIED
C
C     TM1 IS DONE AS AM1, BUT WITH TE4 INSTEAD OF TE5
C  
      IF(OPTBAN.EQ.3) THEN
        SL1 = 1.D0 / DT
C       MATRICE DE MASSE LUMPEE LOCALEMENT SUR LES BANCS DECOUVRANTS
        FORMUL='MSLUMP          '
        CALL MATRIX(TM1,'M=N     ',FORMUL,IELMH,IELMH,
     *              SL1,TE3,S,S,S,S,S,MESH,.TRUE.,TE4)
C       MASS-LUMPING :
        IF(AGGLOH.GT.0.001D0) THEN
          CALL LUMP(T1,TM1,MESH,AGGLOH)
          CALL OM( 'M=CN    ' , TM1 , TM1 , S , 1.D0-AGGLOH , MESH )
          CALL OM( 'M=M+D   ' , TM1 , TM1 , T1 , C          , MESH )
        ENDIF
        CALL MATVEC('X=X+AY  ',CV1,TM1,HN,C,MESH)
      ENDIF
C
C=======================================================================
C
C  PRISE EN COMPTE DES TERMES SOURCES IMPLICITES :
C
      CALL OM( 'M=M+D   ' , AM2,AM2 , T3 , C , MESH )
      CALL OM( 'M=M+D   ' , AM3,AM3 , T4 , C , MESH )
C
C=======================================================================
C
C
C  BOUNDARY TERMS
C
C     PARALLELISM : SUBDOMAINS MAY HAVE NO BOUNDARY POINT AT ALL
C
      IF(MESH%NPTFR.GT.0) THEN
C
C  PRISE EN COMPTE DES TERMES SOMME(PSI H(N) U(N) . N ) AU BORD
C  CES TERMES NE DOIVENT PAS ETRE PRIS EN COMPTE SUR LES PAROIS
C  SOLIDES, D'OU L'UTILISATION DE MASQUE MASK(*,8) QUI VAUT 0.
C  POUR LES PAROIS DE TYPE KLOG.
C
      CALL VECTOR(FLBOR,'=','FLUBDF          ',IELBOR(IELMH,1),
     *            1.D0-TETAU,HPROP,S,S,UN,VN,S,
     *            MESH,.TRUE.,MASK%ADR(UNONNEU)%P)
C
C-----------------------------------------------------------------------
C
C  PRISE EN COMPTE DES TERMES  SOMME(PSI HN U(N+1) . N ) AU BORD
C
C  CONDITIONS DE DIRICHLET : U(N+1) = UBOR ; V(N+1) = VBOR
C
C  UDIR : 1 POUR UN SEGMENT DIRICHLET SUR U, 0 SINON
C
      CALL VECTOR(FLBOR,'+','FLUBDF          ',IELBOR(IELMH,1),
     *            TETAU,HPROP,S,S,UBOR,VBOR,S,
     *            MESH,.TRUE.,MASK%ADR(UDIR)%P)
C
C  PRISE EN COMPTE DES TERMES  SOMME(PSI HN U(N+1) . N ) AU BORD
C
C  CONDITIONS DE SORTIE LIBRE
C
C  UDDL : 1 POUR UN SEGMENT LIBRE SUR U, 0 SINON
C
      CALL VECTOR(FLBOR,'+','FLUBDF          ',IELBOR(IELMH,1),
     *            TETAU,HPROP,S,S,UN,VN,S,
     *            MESH,.TRUE.,MASK%ADR(UDDL)%P)
C
C     WITH POROSITY
C
      IF(OPTBAN.EQ.3) THEN
        DO I=1,MESH%NPTFR
          FLBOR%R(I)=FLBOR%R(I)*TE5%R(MESH%NELBOR%I(I))
        ENDDO
      ENDIF
C
C     EQUALITY OF FLUXES ON EITHER SIDE OF A WEIR
C     (WILL NOT WORK IN PARALLEL BUT WEIRS DON'T WORK IN
C      PARALLEL ANYWAY)
C
      IF(NWEIRS.GT.0) THEN
        DO N=1,NWEIRS
          DO I=1,NPSING(N)
            I1=NUMDIG(1,N,I)
            I2=NUMDIG(2,N,I)
            FLBOR%R(I1)=(FLBOR%R(I1)-FLBOR%R(I2))*0.5D0
            FLBOR%R(I2)=-FLBOR%R(I1)
          ENDDO
        ENDDO
      ENDIF
C
C     BOUNDARY TERMS IN THE RIGHT HAND SIDE
C
      CALL OSDB( 'X=X-Y   ' , CV1 , FLBOR , FLBOR , C , MESH )
C
C  PRISE EN COMPTE DES TERMES  SOMME(PSI HN U(N+1) . N ) AU BORD
C  HOND : 1 POUR LES SEGMENTS D'ONDE INCIDENTE, 0 SINON
C
      CALL INCIDE(COTOND,HN,C0,PATMOS,ATMOS,ZF,MESH,
     *            LT,AT,GRAV,ROEAU,PRIVE)
C
      CALL CPSTVC(COTOND,T5)
      CALL CPSTVC(T5,T6)
      CALL CPSTVC(T5,T2)
C  CALCUL DE CPROP (T5) ET CN (T2)
      CALL OSBD( 'X=CY    ' , T5 , HPROP , S , GRAV , MESH )
      CALL OS  ( 'X=+(Y,C)' , T5 , T5    , S , 0.D0        )
      CALL OS  ( 'X=SQR(Y)' , T5 , T5    , S , C           )
      CALL OSBD( 'X=CY    ' , T2 , HN    , S , GRAV , MESH )
      CALL OS  ( 'X=+(Y,C)' , T2 , T2    , S , 0.D0        )
      CALL OS  ( 'X=SQR(Y)' , T2 , T2    , S , C           )
C
      CALL OS  ( 'X=Y-Z   ' , T6 , T2    , C0     , C      )
      CALL OSBD( 'X=CXY   ' , T6 , HPROP , S      , -2.D0 , MESH )
      CALL OS  ( 'X=X+YZ  ' , T6 , T5    , COTOND , C      )
C
      CALL VECTOR(T2,'=','MASVEC          ',IELBOR(IELMH,1),
     *            1.D0,T6,S,S,S,S,S,MESH,.TRUE.,MASK%ADR(HOND)%P)
      IF(OPTBAN.EQ.3) THEN
        DO I=1,MESH%NPTFR
          T2%R(I)=T2%R(I)*TE5%R(MESH%NELBOR%I(I))
          T5%R(I)=T5%R(I)*TE5%R(MESH%NELBOR%I(I))
        ENDDO
      ENDIF
      CALL OSDB( 'X=X+Y   ' , CV1 , T2 , S , C , MESH )
C
C  ONDE INCIDENTE IMPLICITEE : ON AJOUTE DES TERMES A LA MATRICE AM1
C  MASK(MSK7) : 1 POUR LES SEGMENTS D'ONDE INCIDENTE, 0 SINON
C
      CALL MATRIX(MBOR,'M=N     ','FMATMA          ',
     *            IELBOR(IELMH,1),IELBOR(IELMH,1),
     *            1.D0,T5,S,S,S,S,S,
     *            MESH,.TRUE.,MASK%ADR(HOND)%P)
      CALL OM( 'M=M+N   ' , AM1 , MBOR , S , C , MESH )
C
C     END OF: IF(MESH%NPTFR.GT.0) THEN
      ENDIF
C
C END OF BOUNDARY TERMS
C
C-----------------------------------------------------------------------
C
      IF(EQUA(1:10).EQ.'BOUSSINESQ') THEN
C
C       PRISE EN COMPTE DES TERMES DE BOUSSINESQ
C
        CALL ROTNE0(MESH,CM1,
     *              AM2,A23,A32,AM3,CV2,CV3,UN,VN,H0,MSK,MASKEL,S,DT)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  MATRICES DE GRADIENT
C
C     NOT USED WITH WAVE EQUATION
      IF(SOLSYS.EQ.1) THEN
C
      CALL MATRIX(CM1,'M=N     ','MATGRA         X',IELMU,IELMH,
     *            TETAH*GRAV,S,S,S,S,S,S,MESH,MSK,MASKEL)
      CALL MATRIX(CM2,'M=N     ','MATGRA         Y',IELMU,IELMH,
     *            TETAH*GRAV,S,S,S,S,S,S,MESH,MSK,MASKEL)
C
      ENDIF
C
C=======================================================================
C
C CHOIX DU TIR INITIAL.
C
C POUR L'INSTANT U ET V NE SONT PAS CHANGES, CE QUI SIGNIFIE QUE L'ON
C A TOUJOURS POUR U ET V LE RESULTAT DU DERNIER SYSTEME RESOLU.
C
      IF(IORDRH.EQ.0) THEN
C
        CALL OS( 'X=0     ' , X=DH )
C
      ELSEIF(IORDRH.EQ.1) THEN
C
        IF(LT.EQ.1.AND.ISOUSI.EQ.1) THEN
          CALL OS( 'X=0     ' , X=DH )
        ENDIF
C
      ELSEIF(IORDRH.EQ.2) THEN
C
        IF(LT.EQ.1.AND.ISOUSI.EQ.1) THEN
          CALL OS( 'X=0     ' , X=DH )
          CALL OS( 'X=0     ' , X=DHN)
        ENDIF
        IF (LT.GT.2) CALL OS( 'X=CX    ' , DH , S , S , 2.D0 )
        CALL OS( 'X=X-Y   ' , DH , DHN , S , C )
C       STOCKAGE DE DH(N) DANS DH(N-1)
        CALL OS( 'X=X+Y   ' , DHN , DH   , S , C     )
        CALL OS( 'X=CX    ' , DHN , DHN  , S , 0.5D0 )
C
      ELSE
C
        IF(LNG.EQ.1) WRITE(LU,30) IORDRH
        IF(LNG.EQ.2) WRITE(LU,31) IORDRH
30      FORMAT(1X,'PROPAG : IORDRH=',1I6,' VALEUR NON PREVUE')
31      FORMAT(1X,'PROPAG : IORDRH=',1I6,' VALUE OUT OF RANGE')
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
      IF(IORDRU.EQ.0) THEN
C
        CALL OS( 'X=0     ' , X=U )
        CALL OS( 'X=0     ' , X=V )
C
      ELSEIF(IORDRU.EQ.1) THEN
C
C       U = UN ET V = VN DEJA FAITS
C
      ELSEIF(IORDRU.EQ.2) THEN
C
        IF(LT.EQ.1.AND.ISOUSI.EQ.1) THEN
          CALL OS( 'X=0     ' , X=DU )
          CALL OS( 'X=0     ' , X=DV )
        ENDIF
        IF(ISOUSI.EQ.1) THEN
          CALL OS( 'X=Y+Z   ' , X=U , Y=UN , Z=DU )
          CALL OS( 'X=Y+Z   ' , X=V , Y=VN , Z=DV )
        ENDIF
C
      ELSE
C
        IF(LNG.EQ.1) WRITE(LU,32) IORDRU
        IF(LNG.EQ.2) WRITE(LU,33) IORDRU
32      FORMAT(1X,'PROPAG : IORDRU=',1I6,' VALEUR NON PREVUE')
33      FORMAT(1X,'PROPAG : IORDRU=',1I6,' VALUE OUT OF RANGE')
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C=======================================================================
C
      IF(SOLSYS.EQ.2) THEN
C
        IF(NCSIZE.GT.1) THEN
          CALL PARCOM(AM2%D,2,MESH)
          CALL PARCOM(AM3%D,2,MESH)
        ENDIF
C       INVERSION DE AM2%D ET AM3%D (QUI RESSERVIRONT ENSUITE)
        CALL OS( 'X=1/Y   ' , AM2%D , AM2%D , AM2%D , C ,2,0.D0,1.D-6)
        CALL OS( 'X=1/Y   ' , AM3%D , AM3%D , AM3%D , C ,2,0.D0,1.D-6)
C
C       AJOUT DE LA MATRICE DE "DIFFUSION" A AM1
C
C       ICI AM2%D EST DEJA INVERSEE
        IF(IELMH.EQ.IELMU) THEN
!         ON VEUT DIVISER PAR (1/DT + FROT) QUI EST DANS AM2%D MAIS
!         DANS AM2%D IL EST PROJETE SUR LES BASES, DONC ON REMULTIPLIE
!         PAR LA MASSE DES BASES
          CALL OS('X=CYZ   ',X=DM1,Y=V2DPAR,Z=AM2%D,C=-TETAU*TETAH*GRAV)
        ELSE
!         ICI IL FAUT PRENDRE LA MASSE DES BASES DE U
          CALL VECTOR(T4,'=','MASBAS          ',IELMU,
     *                1.D0,S,S,S,S,S,S,MESH,.FALSE.,TE3)
          IF(NCSIZE.GT.1) CALL PARCOM(T4,2,MESH)
          CALL OS('X=CYZ   ',X=DM1,Y=T4,Z=AM2%D,C=-TETAU*TETAH*GRAV)
          CALL CHGDIS(DM1,IELMU,IELMH,MESH)
        ENDIF   
C
        CALL MATRIX(AM1,'M=M+N   ','MATDIFUV        ',IELMH,IELMH,
     *              -1.D0,S,S,S,DM1,HPROP,S,MESH,MSK,MASKEL)
C
C       NOUVEAU SECOND MEMBRE CV1
C
        IF(ABS(TETAZCOMP-1.D0).LT.1.D-6) THEN 
          CALL OS( 'X=YZ    ' , X=T2 , Y=CV2 , Z=AM2%D )
          CALL OS( 'X=YZ    ' , X=T3 , Y=CV3 , Z=AM3%D )
        ELSE
          CALL OS( 'X=Y     ' , X=T4 , Y=CV2 )
          CALL OS( 'X=Y     ' , X=T5 , Y=CV3 )
C         GRADIENT DE SURFACE LIBRE
          IF(OPTBAN.EQ.1.OR.OPTBAN.EQ.3) THEN
            CALL VECTOR(T4,'+','GRADF          X',IELMU,
     *      (1.D0-TETAZCOMP)*GRAV,ZFLATS,S,S,S,S,S,MESH,MSK,MASKEL)
            CALL VECTOR(T5,'+','GRADF          Y',IELMU,
     *      (1.D0-TETAZCOMP)*GRAV,ZFLATS,S,S,S,S,S,MESH,MSK,MASKEL)
          ELSE
            CALL VECTOR(T4,'+','GRADF          X',IELMU,
     *      (1.D0-TETAZCOMP)*GRAV,T8,S,S,S,S,S,MESH,MSK,MASKEL)
            CALL VECTOR(T5,'+','GRADF          Y',IELMU,
     *      (1.D0-TETAZCOMP)*GRAV,T8,S,S,S,S,S,MESH,MSK,MASKEL)
          ENDIF          
          CALL OS( 'X=YZ    ' , X=T2 , Y=T4 , Z=AM2%D )
          CALL OS( 'X=YZ    ' , X=T3 , Y=T5 , Z=AM3%D )
        ENDIF
C       
        IF(NCSIZE.GT.1) THEN
          CALL PARCOM(T2,2,MESH)
          CALL PARCOM(T3,2,MESH)
        ENDIF
C
C       PRISE EN COMPTE DES CONDITIONS A LA LIMITE
C       ON COMMET UNE ERREUR EN GRAD(DH)
        DO I=1,MESH%NPTFR
          IF(LIMPRO%I(I+DIMLIM).EQ.KDIR) THEN
            T2%R(MESH%NBOR%I(I)) = UBOR%R(I)
          ENDIF
          IF(LIMPRO%I(I+2*DIMLIM).EQ.KDIR) THEN
            T3%R(MESH%NBOR%I(I)) = VBOR%R(I)
          ENDIF
        ENDDO
C
C       REMEMBER THAT COEFFICIENT TETAU IS IN BM1 AND BM2
C       AND THAT SUPG UPWINDING IS ALSO IN BM1 AND BM2
C       OTHERWISE COULD BE INCLUDED IN HUGRADP BELOW
        CALL MATVEC('X=X-AY  ',CV1,BM1,T2,C,MESH)
        CALL MATVEC('X=X-AY  ',CV1,BM2,T3,C,MESH)
!
        IF(ABS(TETAZCOMP-1.D0).GT.1.D-6) THEN
          FORMUL='HUGRADP3        '
!                        3: U AND V, HERE T2 AND T3 NOT TAKEN
          CALL OS('X=CY    ',X=T1,Y=DM1,
     *                       C=(1.D0-TETAZCOMP)/TETAH/TETAU)
          IF(OPTBAN.EQ.1.OR.OPTBAN.EQ.3) THEN
            CALL VECTOR(CV1,'+',FORMUL,IELMH,TETAU,
     *                  HPROP,T1,ZFLATS,T2,T3,T3,MESH,MSK,MASKEL)
          ELSE
            CALL VECTOR(CV1,'+',FORMUL,IELMH,TETAU,
     *                  HPROP,T1,T8    ,T2,T3,T3,MESH,MSK,MASKEL)
          ENDIF
        ENDIF
C
        CALL OS('X=CY    ',X=UDEL,Y=T2,C=TETAU)
        CALL OS('X=CY    ',X=VDEL,Y=T3,C=TETAU)
        CALL OS('X=X+CY  ',X=UDEL,Y=UN,C=1.D0-TETAU)
        CALL OS('X=X+CY  ',X=VDEL,Y=VN,C=1.D0-TETAU)
C
      ENDIF
C
C=======================================================================
C
C     A CE NIVEAU, A23 ET A32 SONT NULLES
C     ELLES SONT NON NULLES SEULEMENT APRES UN PRECONDITIONNEMENT
C     DIAGONAL-BLOC OU AVEC BOUSSINESQ
C
      IF(EQUA(1:10).NE.'BOUSSINESQ') THEN
        A23%TYPDIA='0'
        A32%TYPDIA='0'
        A23%TYPEXT='0'
        A32%TYPEXT='0'
      ENDIF
C
C=======================================================================
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     IN ADJOINT MODE : THE SYSTEM IS NOT SOLVED, RETURN HERE
C
C
      IF(ADJO) RETURN      
C
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C=======================================================================
C
C   CONDITIONS AUX LIMITES DE TYPE DIRICHLET :
C
      IF(SOLSYS.EQ.1) THEN
C
        IF(CORCON.AND.NFRLIQ.GT.0) THEN
C
C         SAUVEGARDE DE L'EQUATION DE CONTINUITE
          CALL OM( 'M=N     ' , TM1  , AM1 , S , C , MESH )
          CALL OM( 'M=N     ' , BM1S , BM1 , S , C , MESH )
          CALL OM( 'M=N     ' , BM2S , BM2 , S , C , MESH )
          CALL OS( 'X=Y     ' ,X=CV1S,Y=CV1)
C
        ENDIF
C
        CALL DIRICH(UNK,MAT,RHS,DIRBOR,LIMPRO%I,TB,MESH,KDIR,MSK,MASKPT)
C
      ENDIF
C
C  RESOLUTION DU SYSTEME OBTENU :
C
C+++++++++++++++++++++++++++++++++++
C  PRECONDITIONNEMENT SPECIAL H-U  +
C+++++++++++++++++++++++++++++++++++
C
      IF(PRECCU) THEN
C
C     PREPARATION DES DIAGONALES QUI VONT PRECONDITIONNER :
C
C     HTILD : D1 (IL DOIT Y AVOIR UNE SIMPLIFICATION PAR GRAV ?
      CALL OS('X=+(Y,C)' , X=HTILD , Y=HN , C=0.D0 )
      CALL OS('X=CX    ' , X=HTILD , C=4.D0/GRAV )
      CALL OS('X=SQR(Y)' , X=HTILD , Y=HTILD )
      CALL OS('X=+(Y,C)' , X=HTILD , Y=HTILD , C=2.D0/GRAV )
C     T1 : D2 (ON NE LE CONSERVE PAS)
      CALL OS('X=1/Y   ' , X=T1 , Y=HTILD )
      CALL OS('X=CX    ' , X=T1 , C=4.D0*TETAH/TETAU )
C
C     CHANGEMENT DU SECOND MEMBRE
C
      CALL OS('X=XY    ' , X=CV1 , Y=T1 )
C
C     CHANGEMENT DE LA VARIABLE DH
C
      CALL OS('X=Y/Z   ' , X=DH , Y=DH , Z=HTILD )
C
C     PRECONDITIONNEMENT DE AM1
C
      IF(AM1%TYPEXT.EQ.'S') CALL OM( 'M=X(M)  ' , AM1,AM1 ,S,C,MESH )
      CALL OM( 'M=DM    ' , AM1 , AM1 , T1    , C , MESH )
      CALL OM( 'M=MD    ' , AM1 , AM1 , HTILD , C , MESH )
C
C     PRECONDITIONNEMENT DE BM1 ET BM2 :
C
      CALL OM( 'M=DM    ' , BM1 , BM1 , T1 , C , MESH )
      CALL OM( 'M=DM    ' , BM2 , BM2 , T1 , C , MESH )
C
C     PRECONDITIONNEMENT DE CM1 ET CM2 :
C
      CALL OM( 'M=MD    ' , CM1 , CM1 , HTILD , C , MESH )
      CALL OM( 'M=MD    ' , CM2 , CM2 , HTILD , C , MESH )
C
      ENDIF
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  FIN DU PRECONDITIONNEMENT SPECIAL H-U                            +
C  SAUF RECUPERATION DE LA VARIABLE DH (VOIR APRES APPEL A SOLV09)  +
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
      IF(SOLSYS.EQ.1) THEN
C
C       METHODE CLASSIQUE
C
C       CAS DU PRECONDITIONNEMENT BLOC-DIAGONAL
C       A23 ET A32 VONT ETRE UTILISEES, ON LES INITIALISE.
C       AVEC BOUSSINESQ COUPLE, C'EST DEJA FAIT
        IF(EQUA(1:10).NE.'BOUSSINESQ') THEN
          IF(3*(SLVPRO%PRECON/3).EQ.SLVPRO%PRECON) THEN
            A23%TYPDIA='Q'
            A32%TYPDIA='Q'
            A23%TYPEXT='Q'
            A32%TYPEXT='Q'
            CALL OM( 'M=0     ' , A23,A23 , S,C , MESH )
            CALL OM( 'M=0     ' , A32,A32 , S,C , MESH )
            IF (AM2%TYPEXT.EQ.'S') THEN
             CALL OM( 'M=X(M)  ' , AM2,AM2 , S,C , MESH )
            ENDIF
            IF (AM3%TYPEXT.EQ.'S') THEN
             CALL OM( 'M=X(M)  ' , AM3,AM3 , S,C , MESH )
            ENDIF
          ENDIF
        ENDIF
C
C       PRINT*,'DAM1=',DOTS(AM1%D,AM1%D)
C       PRINT*,'DAM2=',DOTS(AM2%D,AM2%D)
C       PRINT*,'DAM3=',DOTS(AM3%D,AM3%D)
C       PRINT*,'DBM1=',DOTS(BM1%D,BM1%D)
C       PRINT*,'DBM2=',DOTS(BM2%D,BM2%D)
C       PRINT*,'DCM1=',DOTS(CM1%D,CM1%D)
C       PRINT*,'DCM2=',DOTS(CM2%D,CM2%D)
C       PRINT*,'CV1=',DOTS(CV1,CV1)
C       PRINT*,'CV2=',DOTS(CV2,CV2)
C       PRINT*,'CV3=',DOTS(CV3,CV3)
C
        CALL SOLVE(UNK,MAT,RHS,TB,SLVPRO,INFOGR,MESH,TM1)
C
      ELSEIF(SOLSYS.EQ.2) THEN
C 
C       EQUATION D'ONDE GENERALISEE
C
C       SYSTEME EN H
C
C       SAUVEGARDE DE AM1 DANS BM1 ET DE CV1 DANS BM2%D
C
        IF(CORCON.AND.NFRLIQ.GT.0) THEN
          CALL OM('M=N     ',BM1,AM1,S,C,MESH)       
          CALL OS('X=Y     ',X=BM2%D,Y=CV1)
        ENDIF
C
        CALL DIRICH(DH,AM1,CV1,HBOR,LIMPRO%I,TB,MESH,KDIR,MSK,MASKPT)
        CALL SOLVE(DH,AM1,CV1,TB,SLVPRO,INFOGR,MESH,TM1)
C 
        NELEM=MESH%NELEM
        DO IELEM=1,NELEM
          ZCONV%R(IELEM        )=DH%R(MESH%IKLE%I(IELEM        ))
          ZCONV%R(IELEM+  NELEM)=DH%R(MESH%IKLE%I(IELEM+  NELEM))
          ZCONV%R(IELEM+2*NELEM)=DH%R(MESH%IKLE%I(IELEM+2*NELEM))
        ENDDO
        IF(ABS(1.D0-TETAZCOMP).GT.1.D-6) THEN
          C=(1.D0-TETAZCOMP)/TETAH
          IF(OPTBAN.EQ.1.OR.OPTBAN.EQ.3) THEN
!           FREE SURFACE PIECE-WISE LINEAR IN ZFLATS
            CALL OS('X=X+CY  ',X=ZCONV,Y=ZFLATS,C=C)
          ELSE
!           FREE SURFACE LINEAR          
            DO IELEM=1,NELEM
              I1=MESH%IKLE%I(IELEM        )
              I2=MESH%IKLE%I(IELEM+  NELEM)
              I3=MESH%IKLE%I(IELEM+2*NELEM)  
              ZCONV%R(IELEM        )=ZCONV%R(IELEM        )+
     *        C*(T8%R(I1))
              ZCONV%R(IELEM+  NELEM)=ZCONV%R(IELEM+  NELEM)+
     *        C*(T8%R(I2))
              ZCONV%R(IELEM+2*NELEM)=ZCONV%R(IELEM+2*NELEM)+
     *        C*(T8%R(I3))
            ENDDO
          ENDIF        
        ENDIF 
C
C       SYSTEMES EN U ET V
C  
        CALL VECTOR(CV2,'+','GRADF          X',IELMU,
     *              -GRAV*TETAH,DH,S,S,S,S,S,MESH,MSK,MASKEL)
        CALL VECTOR(CV3,'+','GRADF          Y',IELMU,
     *              -GRAV*TETAH,DH,S,S,S,S,S,MESH,MSK,MASKEL)
        IF(NCSIZE.GT.1) THEN
          CALL PARCOM(CV2,2,MESH)
          CALL PARCOM(CV3,2,MESH)
        ENDIF
C                                      AM2%D ET AM3%D DEJA INVERSEES
        CALL OS('X=YZ    ',X=U,Y=CV2,Z=AM2%D)
        CALL OS('X=YZ    ',X=V,Y=CV3,Z=AM3%D)
C
        DO I=1,MESH%NPTFR
          IF(LIMPRO%I(I+DIMLIM).EQ.KDIR) THEN
            U%R(MESH%NBOR%I(I)) = UBOR%R(I)
          ENDIF
          IF(LIMPRO%I(I+2*DIMLIM).EQ.KDIR) THEN
            V%R(MESH%NBOR%I(I)) = VBOR%R(I)
          ENDIF
        ENDDO
C
C       CORRECTION FINALE DES FLUX AU BORD POUR LES FRONTIERES
C       A HAUTEUR IMPOSEE (POUR RESOUDRE L'EQUATION DE CONTINUITE)
C
        IF(CORCON.AND.NFRLIQ.GT.0) THEN
          CALL MATVEC('X=X+CAY ',BM2%D,BM1,DH,-1.D0,MESH)
          DO I=1,MESH%NPTFR
            IF(LIMPRO%I(I).EQ.KDIR) THEN
              FLBOR%R(I)=FLBOR%R(I)+BM2%D%R(MESH%NBOR%I(I))
            ENDIF
          ENDDO
        ENDIF     
C 
      ENDIF
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PRECONDITIONNEMENT SPECIAL H-U : RECUPERATION DE LA VARIABLE DH
C
      IF(PRECCU) CALL OS('X=XY    ' , X=DH , Y=HTILD )
C
      IF(CORCON.AND.SOLSYS.EQ.1.AND.NFRLIQ.GT.0) THEN
C
C       CORRECTION FINALE DES FLUX AU BORD POUR LES FRONTIERES
C       A HAUTEUR IMPOSEE (POUR RESOUDRE L'EQUATION DE CONTINUITE)
C
        CALL MATVEC('X=X+CAY ',CV1S,TM1,DH,-1.D0,MESH)
        CALL MATVEC('X=X+CAY ',CV1S,BM1S,U,-1.D0,MESH)
        CALL MATVEC('X=X+CAY ',CV1S,BM2S,V,-1.D0,MESH)
        DO I=1,MESH%NPTFR
          IF(LIMPRO%I(I).EQ.KDIR) THEN
            FLBOR%R(I)=FLBOR%R(I)+CV1S%R(MESH%NBOR%I(I))
          ENDIF
        ENDDO     
      ENDIF
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  PASSAGE DE DH A H
C
      CALL OS( 'X=Y+Z   ' , X=H , Y=DH , Z=HN )
C
C  RECONSTRUCTION DE HBOR : HBOR = HBOR + HN
C  HBOR EST UTILISE A NOUVEAU QUAND ON FAIT DES SOUS-ITERATIONS
C
      CALL OSBD( 'X=X+Y   ' , HBOR , HN , S , C , MESH )
C
C  SAUVEGARDE DE L'ACCROISSEMENT SUR LES VITESSES
C
      IF(IORDRU.EQ.2) THEN
        CALL OS( 'X=Y-Z   ' , X=DU , Y=U , Z=UN )
        CALL OS( 'X=Y-Z   ' , X=DV , Y=V , Z=VN )
      ENDIF
C
C  COMPATIBLE VELOCITY FIELD IN CONTINUITY EQUATION
C
      IF(SOLSYS.EQ.1) THEN
        CALL OS ('X=CY    ',X=UDEL,Y=U ,C=     TETAU)
        CALL OS ('X=CY    ',X=VDEL,Y=V ,C=     TETAU)
        CALL OS ('X=X+CY  ',X=UDEL,Y=UN,C=1.D0-TETAU)
        CALL OS ('X=X+CY  ',X=VDEL,Y=VN,C=1.D0-TETAU)
      ENDIF
C
C-----------------------------------------------------------------------
C    
      RETURN
      END
