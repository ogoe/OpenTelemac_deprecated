C                       *****************                               
                        SUBROUTINE VOLFIN                               
C                       *****************                          
C                                                                       
     * (W1,AT,DT,LT,NIT,NELEM,NPTFR,                
     *  TB,ZF,CF,NPOIN,HN,H,U,V,QU,QV,G,LISTIN,                         
     *  S,MSK,MASKEL,MESH,LIMPRO,NBOR,KDIR,KNEU,KDDL,             
     *  HBOR,UBOR,VBOR,MASSES,FLUENT,FLUSOR,CFLWTD,DTVARI,KFROT,
     *  NREJET,ISCE,TSCE2,MAXSCE,MAXTRA,YASMH,SMH,
     *  NTRAC,DIMT,T,HT,TN,DLIMT,LIMTRA,
     *  TBOR,MASSOU,FLUTENT,FLUTSOR,DTHAUT,DPX,DPY,DJX,DJY,CMI,JMI,
     *  DJXT,DJYT,DIFVIT,ITURB,PROPNU,DIFT,DIFNU,
     *  DX,DY,OPTVF,
     *  HSTOK,HCSTOK,LOGFR,DSZ,FLUXT,FLUHBOR,DTN,FLUSORTN,FLUENTN,LTT,
     *  FLUXTEMP,FLUHBTEMP,HC,SMTR,AIRST,TMAX,DTT)
C                                                      
C***********************************************************************
C  TELEMAC-2D VERSION 5.8    05/08/07  J-M HERVOUET (LNH) 01 30 87 80 18 
C  COPYRIGHT EDF-1997                    INRIA                          
C-----------------------------------------------------------------------
C                                                                       
C  FONCTION  : INTERFACE POUR APPEL DU SOUS-PROGRAMME RESOLU.           
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________ 
C |      NOM       |MODE|                   ROLE                       | 
C |________________|____|______________________________________________|
C |      W1        |<-->| TABLEAU DE TRAVAIL                           | 
C |                |    | VOIR SON DIMENSIONNEMENT DANS POINT          | 
C |                |    | ICI DE TAILLE MINIMUM 9*NPOIN + 3*NPTFR      | 
C |      AT,DT,LT  | -->| TEMPS, PAS DE TEMPS, NUMERO DU PAS           | 
C |      NIT       | -->| NOMBRE TOTAL D'ITERATIONS.                   | 
C |      NELEM     | -->| NOMBRE D'ELEMENTS                            | 
C |      NPTFR     | -->| NOMBRE DE POINTS FRONTIERE                   | 
C |      TB        | -->| BLOC DE TABLEAUX DE TRAVAIL DIMENSION NPOIN  | 
C |      ZF        | -->| COTES DU FOND                                | 
C |      CF        | -->| COEFFICIENT DE FROTTEMENT                    !
C |      NPOIN     | -->| NOMBRE DE POINTS DU MAILLAGE.                | 
C |      HN,H      |<-->| HAUTEURS D'EAU AU TEMPS N ET N+1             | 
C |      U,V       |<-->| COMPOSANTES DE LA VITESSE.                   | 
C |      QU,QV     |<-->| COMPOSANTES DU DEBIT.                        | 
C |      G         | -->| GRAVITENTES DU DEBIT.                        | 
C |      LISTIN    | -->| SI OUI, MESSAGES IMPRIMES SUR LISTING.       | 
C |      S         | -->| STRUCTURE VIDE POUR APPEL A VECTOR.          | 
C |      MSK       | -->| SI OUI, MASQUAGE D'ELEMENTS.                 | 
C |      MASKEL    | -->| MASQUE DES ELEMENTS.                         | 
C |      MESH      | -->| BLOC DE TABLEAUX D'ENTIERS DU MAILLAGE.      |
C |      LIMPRO    | -->| TYPES DE CONDITIONS AUX LIMITES.             | 
C |      NBOR      | -->| NUMEROS GLOBAUX DES POINTS DE BORD.          | 
C |      KDIR      | -->| CONVENTION POUR LES POINTS DIRICHLET.        | 
C |      KNEU      | -->| CONVENTION POUR LES POINTS NEUMANN.          | 
C |      KDDL      | -->| CONVENTION POUR LES POINTS LIBRES.           | 
C |      HBOR      | -->| VALEURS IMPOSEES DE H                        | 
C |      UBOR      | -->| VALEURS IMPOSEES DE U                        | 
C |      VBOR      | -->| VALEURS IMPOSEES DE V                        | 
C |      MASSES    |<-- | MASSE AJOUTEE PAR TERME SOURCE               | 
C !  FLUENT,FLUSORT|<-- | FLUX MASSE ENTREE ET SORTIE DE TN A TN+1     |
C !  CFLWTD        ! -->! NOMBRE DE CFL                                !
C |  KFROT         | -->| LOI DE FROTTEMENT SUR LE FOND                |
C |  NREJET        | -->| NOMBRE DE SOURCES/PUITS                      |
C |  ISCE          | -->| POINTS SOURCES                               |
C |  YASMH         | -->|  INDIQUE SI ON PREND EN COMPTE SMH           |
C |  SMH           | -->|  TERMES SOURCES DE L'EQUATION DE CONTINUITE  |
C |  TRAC          | -->|  LOGIQUE INDIQUANT LA PRESENCE D'UN TRACEUR  |
C !  DIMT          ! -->!  DIMENSION DU TRACEUR                        !
C !  T             !<-- !  TRACEUR MIS A JOUR                          !
C !  HT            !<-->!  HT AU TEMPS N                               !
C !  TN            ! -->!  T  AU TEMPS N                               !
C |  TSCE          | -->|  VALEURS DU TRACEUR AUX SOURCES              |
C !  DLIMT         ! -->!  DIMENSION DU TRACEUR AU BORD                !
C |  LIMTRA        | -->|  TYPES DE CONDITIONS AUX LIMITES SUR TRACEUR |
C |  TBOR          | -->|  CONDITIONS AUX LIMITES SUR T                !
C |  MASSOU        |<-- |  MASSE DE TRACEUR AJOUTEE PAR TERME SOURCE   |
C ! FLUTENT,FLUTSOR|<-- |  FLUX TRACEUR ENTREE ET SORTIE               | 
C |  DTHAUT        ! -->!  UTILISE POUR CONDITION CFL                  !
C |  DPX, DPY      ! -->!  GRADIENTS DES FONCTIONS DE BASE             !
C |  DJX, DJY      ! -- !  TABLEAUX DE TRAVAIL                         !
C !  CMI           ! -->!  COORDONNEES DES POINTS MILIEUX D'INTERFACE  !
C !  JMI           ! -->!  NUMERO DU TRIANGLE AUQUEL APPARTIENT LE     !
C |                !    !  POINT MILIEU DE L'INTERFACE                 !
C |  DJXT,DJYT     | -- |  TABLEAUX DE TRAVAIL POUR TRACEUR            |
C |  DIFVIT        | -->|  INDIQUE S'IL FAUT FAIRE LA DIFFUSION DE U,V |
C |  ITURB         | -->|  MODELE DE TURBULENCE  1 : LAMINAIRE         |
C |  PROPNU        | -->|  COEFFICIENT DE DIFFUSION MOLECULAIRE        |
C |  DIFT          | -->|  LOGIQUE INDIQUANT S'IL Y A DIFFUSION TRACEUR|
C |  DIFNU         | -->|  COEFFICIENT DE DIFFUSION DU TRACEUR         |
C |  DX,DY         | -- |  TABLEAUX DE TRAVAIL                         |
C |  OPTVF         ! -->!  OPTION SCHEMA                               !
C |                !    !0:ROE, 1:CINETIQUE ORDRE 1,2:CINETIQUE ORDRE 2!
C |  HSTOK,HCSTOK  !<-->!  H, H CORRIGE  A STOCKER POUR TRACEUR        !
C |  LOGFR         !<-->!  REFERENCE DES NOEUDS FRONTIERE              !
C |  DSZ           !<-->!  VARIATION DE Z POUR ORDRE 2                 !
C |  DTN           !<-->!  PAS DE TEMPS   DE TN+1 A TN+2               !
C |  FLUXT,FLUHBOR !<-->!  FLUX, FLUX BORD ACCUMULES POUR TRACEUR      !
C |FLUSORTN,FLUENTN!<-->!  FLUX MASSE ENTREE ET SORTIE DE TN+1 A TN+2  !
C |  LTT           !<-->!  NOMBRE DE PAS DE TEMPS TRACEUR              !
C |  FLUXTEMP      !<-->!  FLUX POUR TRACEUR                           !
C !  FLUHBTEMP     !<-->!  FLUX BORD POUR TRACEUR                      !
C !  HC            !<-->!  H RECONSTRUIT ORDRE 2   CORRIGE             !
C |  SMTR          |<-->!  TERMES SOURCES DU TRACEUR                   !
C |  AIRST         ! -->!  AIRES DES SOUS-TRIANGLES DANS CELLULES      !
C !  TMAX          ! -->!  TEMPS DE FIN DU CALCUL                      !
C !  DTT           !<-->!  PAS DE TEMPS TRACEUR                        !
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE) 
C-----------------------------------------------------------------------
C                                                                       
C PROGRAMME APPELANT : TELMAC                                           
C PROGRAMMES APPELES : RESOLU                                           
C                                                                       
C***********************************************************************
C
      USE BIEF
      USE INTERFACE_TELEMAC2D, EX_VOLFIN => VOLFIN
C                                                                       
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NPTFR,KDIR,KNEU,KDDL,DIMT,KFROT,OPTVF
      INTEGER, INTENT(IN)    :: NELEM,NPOIN,LT,NIT,NREJET,ITURB,DLIMT
      INTEGER, INTENT(IN)    :: NTRAC,MAXSCE,MAXTRA
      INTEGER, INTENT(INOUT) :: LTT
      INTEGER, INTENT(IN)    :: LIMPRO(NPTFR,6),NBOR(NPTFR)
      INTEGER, INTENT(IN)    :: LIMTRA(DLIMT)
      INTEGER, INTENT(IN)    :: ISCE(NREJET)
      INTEGER, INTENT(INOUT) :: JMI(*),LOGFR(NPOIN)
      LOGICAL, INTENT(IN)    :: DIFVIT,DIFT,LISTIN,MSK,DTVARI,YASMH
      DOUBLE PRECISION, INTENT(IN) :: PROPNU,DIFNU
      DOUBLE PRECISION, INTENT(INOUT) :: AT,DT,MASSES,DTT
      DOUBLE PRECISION, INTENT(INOUT) :: H(NPOIN),QU(NPOIN),QV(NPOIN)           
      DOUBLE PRECISION, INTENT(INOUT) :: W1(*)
C                                              NSEG    NSEG               
      DOUBLE PRECISION, INTENT(INOUT) :: DSZ(2,*),HC(2,*)
      DOUBLE PRECISION, INTENT(INOUT) :: U(NPOIN),V(NPOIN) 
      DOUBLE PRECISION, INTENT(IN)    :: HN(NPOIN),SMH(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: CF(NPOIN),ZF(NPOIN),G 
      DOUBLE PRECISION, INTENT(INOUT) :: HSTOK(DIMT),HCSTOK(2,*)
      DOUBLE PRECISION, INTENT(IN)    :: HBOR(NPTFR),UBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: VBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: TSCE2(MAXSCE,MAXTRA)
      DOUBLE PRECISION, INTENT(INOUT) :: DX(3,*),DY(3,*)
      DOUBLE PRECISION, INTENT(IN)    :: AIRST(2,*)
      DOUBLE PRECISION, INTENT(INOUT) :: DPX(3,*),DPY(3,*)
      DOUBLE PRECISION, INTENT(INOUT) :: CMI(2,*),DJX(3,*),DJY(3,*)
      DOUBLE PRECISION, INTENT(IN)    :: CFLWTD,DTHAUT(NPOIN),TMAX 
      DOUBLE PRECISION, INTENT(INOUT) :: FLUSOR,FLUENT,DTN,MASSOU(*)
      DOUBLE PRECISION, INTENT(INOUT) :: FLUSORTN,FLUENTN
      DOUBLE PRECISION, INTENT(INOUT) :: DJXT(*),DJYT(*)
      DOUBLE PRECISION, INTENT(INOUT) :: FLUTENT(*),FLUTSOR(*)    
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: TB
      TYPE(BIEF_OBJ), INTENT(IN)      :: S,MASKEL 
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
      TYPE(BIEF_OBJ) , INTENT(IN)     :: TBOR,TN
      TYPE(BIEF_OBJ) , INTENT(INOUT)  :: T,HT,SMTR,FLUHBOR,FLUHBTEMP
      TYPE(BIEF_OBJ) , INTENT(INOUT)  :: FLUXTEMP,FLUXT             
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+                   
C                        
C     MASSE AJOUTEE PAR TERME SOURCE (NULLE POUR L'INSTANT)                
C                                                                       
      MASSES = 0.D0  
C                                                                       
C     CALCUL DES AIRES DES CELLULES                                     
C                                                                       
      CALL VECTOR(TB%ADR(1)%P,'=','MASBAS          ',11,      
     *            1.D0,S,S,S,S,S,S,MESH,MSK,MASKEL)                                             
C
      CALL RESOLU(W1,W1(1+3*NPOIN),MESH%NUBO%I,                         
     *            MESH%VNOIN%R,W1(1+9*NPOIN),AT,DT,LT,NIT,   
     *            NELEM,MESH%NSEG,NPTFR,W1(1+6*NPOIN),                       
     *            TB%ADR(1)%P%R,MESH%SURFAC%R,
     *            MESH%X%R,MESH%Y%R,MESH%IKLE%I, 
     *            ZF,CF,NPOIN,HN,H,U,V,QU,QV,G,LISTIN,                  
     *            MESH%XNEBOR%R,MESH%YNEBOR%R, 
     *            LIMPRO,NBOR,KDIR,KNEU,KDDL,HBOR,UBOR,VBOR,            
     *            FLUSOR,FLUENT,CFLWTD,DTVARI,MESH%NELMAX,KFROT,
     *            NREJET,ISCE,TSCE2,MAXSCE,MAXTRA,YASMH,SMH,MASSES,
     *            NTRAC,DIMT,T,HT,TN,DIMT,LIMTRA,
     *            TBOR,MASSOU,FLUTENT,FLUTSOR,DTHAUT,DPX,DPY,DJX,DJY,
     *            CMI,JMI,SMTR,
     *            TB%ADR(3)%P%R,TB%ADR(4)%P%R,
     *            DJXT,DJYT,DIFVIT,ITURB,PROPNU,DIFT,DIFNU,
     *            DX,DY,OPTVF,
     *            FLUSORTN,FLUENTN,DSZ,AIRST,HSTOK,HCSTOK,FLUXT,FLUHBOR,
     *            LOGFR,LTT,DTN,FLUXTEMP,FLUHBTEMP,HC,TMAX,DTT,
     *            TB%ADR(6)%P%R,TB%ADR(7)%P%R,TB%ADR(8)%P%R,
     *            TB%ADR(9)%P%R,TB%ADR(10)%P%R)
C                                                                      
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END
