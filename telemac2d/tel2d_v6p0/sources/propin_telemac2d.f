C                       ***************************
                        SUBROUTINE PROPIN_TELEMAC2D
C                       ***************************
C
     *(LIMPRO,LIMDIM,MASK,LIUBOR,LIVBOR,LIHBOR,KP1BOR,NBOR,NPTFR,
     * KENT,KENTU,KSORT,KADH,KLOG,KINC,KNEU,KDIR,KDDL,KOND,
     * CLH,CLU,CLV,IELMU,U,V,GRAV,H,LT,NPOIN,NELBOR,NELMAX,MSK,MASKEL,
     * NFRLIQ,THOMFR,DEBLIQ,FINLIQ,FRTYPE,XNEBOR,YNEBOR,ENTET,MESH)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.9  27/06/2008 J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C      FONCTIONS:
C      ==========
C
C        1) VERIFICATION DE LA COMPATIBILITE DES CONDITIONS LIMITES.
C
C        2) REMPLISSAGE DES TABLEAUX LIMPRO ET MASK.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   LIMPRO       |<-- |    TYPES DE CONDITIONS AUX LIMITES POUR LA   |
C |                |    |    PROPAGATION                               |
C |                |    |    PAR POINTS   :    .1:H  .2:U  .3:V        |
C |                |    |    PAR SEGMENTS :    .4:H  .5:U  .6:V
C |   MASK         |<-- |    MASQUES POUR LES SEGMENTS
C |                |    |    MASK(NPTFR,1) : 1. SI KDIR SUR U 0. SINON
C |                |    |    MASK(NPTFR,2) : 1. SI KDIR SUR V 0. SINON
C |                |    |    MASK(NPTFR,3) : 1. SI KDDL SUR U 0. SINON
C |                |    |    MASK(NPTFR,4) : 1. SI KDDL SUR V 0. SINON
C |                |    |    MASK(NPTFR,5) : 1. SI KNEU SUR U 0. SINON
C |                |    |    MASK(NPTFR,6) : 1. SI KNEU SUR V 0. SINON
C |                |    |    MASK(NPTFR,7) : 1. SI KOND 0. SINON
C |                |    |    MASK(NPTFR,8) : 1. - MASK( ,5)
C |                |    |    MASK(NPTFR,9) : 1. SI H DIRICHLET
C |                |    |    MASK(NPTFR,10): 1. SI H NEUMANN
C |                |    |    MASK(NPTFR,11): 1. SI H DEGRE DE LIBERTE
C |   LIUBOR       | -->|    TYPES DE CONDITIONS AUX LIMITES SUR U
C |   LIVBOR       | -->|    TYPES DE CONDITIONS AUX LIMITES SUR V
C |   LIHBOR       | -->|    TYPES DE CONDITIONS AUX LIMITES SUR H
C |   KP1BOR       | -->|    POINT SUIVANT SUR LA FRONTIERE.
C |   NBOR         | -->|    CORRESPONDANCE ENTRE NUMEROTATION DES
C |                |    |    POINTS FRONTIERES ET NUMEROTATION GLOBALE
C |   NPTFR        | -->|    DIMENSION DES TABLEAUX.
C |                |    | CONDITIONS AUX LIMITES PHYSIQUES:
C |   KENT         | -->|    INDICATEUR DE POINT D'ENTREE FLUIDE
C |   KENTU        | -->|    INDICATEUR DE VITESSE IMPOSEE.
C |   KSORT        | -->|    INDICATEUR DE POINT DE SORTIE FLUIDE
C |   KADH         | -->|    INDICATEUR DE POINT DIRICHLET
C |   KLOG         | -->|    INDICATEUR DE PAROI SOLIDE
C |   KINC         | -->|    INDICATEUR D'ONDE INCIDENTE
C |                |    | CONDITIONS AUX LIMITES TECHNIQUES:
C |   KNEU         | -->| INDICATEUR DE POINT DE NEUMANN
C |   KDIR         | -->| INDICATEUR DE POINT DE DIRICHLET
C |   KDDL         | -->| INDICATEUR DE DEGRE DE LIBERTE AU BORD
C |   KOND         | -->| INDICATEUR D'ONDE INCIDENTE
C |   CLH,CLU,CLV  |<-->| TYPES DE CONDITIONS AUX LIMITES SUR H,U,V
C |                |    | RECOPIES DE LIHBOR,LIUBOR,LIVBOR
C |   U,V, ,H      | -->| VALEURS DE U,V   ET H AU TEMPS T
C |   GRAV         | -->| PESANTEUR
C |   LT           | -->| NUMERO DE L'ITERATION COURANTE.
C |   NPOIN        | -->| NOMBRE DE NOEUD DU MAILLAGE
C |   NELBOR       | -->| NUMEROS DES ELEMENTS ADJACENTS AUX BORDS.
C |   NELMAX       | -->| NOMBRE MAXIMUM D'ELEMENTS.
C |   MSK          | -->| SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |   MASKEL       | -->| TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    | =1. : NORMAL   =0. : ELEMENT MASQUE
C |   NFRLIQ       | -->| NOMBRE DE FRONTIERES LIQUIDES
C |   THOMFR       | -->| TRAITEMENT PAR CARACTERISTIQUES DES FRONTIERES
C |                |    | LIQUIDES
C |   DEBLIQ       | -->| NUMERO DU PREMIER POINT DE LA FRONTIERE LIQUID
C |   FINLIQ       | -->| NUMERO DU DERNIER POINT DE LA FRONTIERE LIQUID
C |   FRTYPE       | -->| TYPE DE TRAITEMENT POUR LES FRONTIERES LIQUIDE
C |   ENTET        | -->| SI OUI : MESSAGES IMPRIMES
C |                |    | SAUF MESSAGES D'ERREURS QUI TOUJOURS IMPRIMES
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NPOIN,NELMAX,NPTFR,KOND,KENTU,LT,NFRLIQ
      INTEGER, INTENT(IN) :: KENT,KSORT,KADH,KLOG,KINC,KNEU,KDIR,KDDL
      INTEGER, INTENT(IN) :: LIMDIM,IELMU
      INTEGER, INTENT(IN) :: NELBOR(NPTFR),LIVBOR(NPTFR),LIHBOR(NPTFR)
      INTEGER, INTENT(IN) :: LIUBOR(NPTFR),FRTYPE(*)
      INTEGER, INTENT(IN) :: DEBLIQ(NFRLIQ),FINLIQ(NFRLIQ)
      INTEGER, INTENT(INOUT) :: LIMPRO(LIMDIM,6)
      INTEGER, INTENT(IN) :: KP1BOR(NPTFR),NBOR(NPTFR)
      INTEGER, INTENT(INOUT) :: CLH(NPTFR),CLU(NPTFR),CLV(NPTFR)
      LOGICAL, INTENT(IN) :: MSK,THOMFR,ENTET
      DOUBLE PRECISION, INTENT(IN)   :: XNEBOR(NPTFR),YNEBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)   :: U(NPOIN),V(NPOIN),H(NPOIN)
      DOUBLE PRECISION, INTENT(IN)   :: GRAV,MASKEL(NELMAX)  
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: MASK
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C    
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K,N,IFRLIQ,IELEM,KP1
C
      INTEGER, PARAMETER :: UDIR    =  1
      INTEGER, PARAMETER :: VDIR    =  2
      INTEGER, PARAMETER :: UDDL    =  3
      INTEGER, PARAMETER :: VDDL    =  4
      INTEGER, PARAMETER :: UNEU    =  5
      INTEGER, PARAMETER :: VNEU    =  6
      INTEGER, PARAMETER :: HOND    =  7
      INTEGER, PARAMETER :: UNONNEU =  8
      INTEGER, PARAMETER :: HDIR    =  9
      INTEGER, PARAMETER :: HNEU    = 10
      INTEGER, PARAMETER :: HDDL    = 11
C
      DOUBLE PRECISION YY,F2,F3   
C
      LOGICAL DEP,ALERTE1,ALERTE2
C
      INTEGER IGUILT1,IGUILT2
C
      INTRINSIC MAX
C
C-----------------------------------------------------------------------
C
      DO 1 K=1,NPTFR
        CLH(K) = LIHBOR(K)
        CLU(K) = LIUBOR(K)
        CLV(K) = LIVBOR(K)
1     CONTINUE
C
C COURT-CIRCUITAGE POUR OPTION DE TRAITEMENT PAR THOMPSON
C
      IF(NFRLIQ.NE.0.AND.THOMFR) THEN
C
      DO 2 IFRLIQ = 1 , NFRLIQ
        DEP=.FALSE.
        K = DEBLIQ(IFRLIQ)
15      IF(FRTYPE(IFRLIQ).EQ.2.AND.H(NBOR(K)).GT.1.D-3) THEN
          CLH(K) = KENT
          CLU(K) = KENTU
          CLV(K) = KENTU
        ENDIF
        IF(K.EQ.FINLIQ(IFRLIQ).AND.DEP) THEN
          GO TO 16
        ELSE
          DEP=.TRUE.
!         NE MARCHE PAS EN PARALLELE MAIS THOMSON NON PLUS
          K = KP1BOR(K)
          GO TO 15
        ENDIF
16      CONTINUE
2     CONTINUE
C
      ENDIF
C
C  CONTROLE ET MODIFICATION EVENTUELLE DES CONDITIONS POUR EVITER DES
C  CAS NON PHYSIQUES : SORTIE ENTIEREMENT LIBRE EN ECOULEMENT FLUVIAL
C                      ONDE INCIDENTE EN ECOULEMENT TORRENTIEL SORTANT
C                      ONDE INCIDENTE EN ECOULEMENT TORRENTIEL RENTRANT
C
      ALERTE1=.FALSE.
      ALERTE2=.FALSE.
C
      DO 3 K=1,NPTFR
C
        N = NBOR(K)
        F2 = (U(N)**2+V(N)**2) / GRAV / MAX(H(N),1.D-8)        
C
C       ONDE INCIDENTE EN ECOULEMENT TORRENTIEL SORTANT
C       ONDE INCIDENTE EN ECOULEMENT TORRENTIEL RENTRANT
C
        IF(CLU(K).EQ.KINC.AND.
     *     CLV(K).EQ.KINC.AND.
     *     F2.GE.1.D0) THEN
          CLU(K) = KSORT
          CLV(K) = KSORT
        ENDIF
C
C       SORTIE ENTIEREMENT LIBRE EN ECOULEMENT FLUVIAL
C
        IF(CLH(K).EQ.KSORT.AND.
     *     CLU(K).EQ.KSORT.AND.
     *     CLV(K).EQ.KSORT.AND.
     *     F2.LE.1.D0) THEN
          CLU(K) = KINC
          CLV(K) = KINC
        ENDIF
C
C       VITESSE LIBRE ENTRANTE
C
        IF(CLU(K).EQ.KSORT.AND.CLV(K).EQ.KSORT) THEN
          F3 = U(N)*XNEBOR(K)+V(N)*YNEBOR(K)
          IF(F3.LE.-1.D-2) THEN
            ALERTE1=.TRUE.
            IGUILT1=K
          ENDIF
        ENDIF
C
C       ENTREE TORRENTIELLE AVEC HAUTEUR LIBRE
C
        IF(CLH(K).EQ.KSORT.AND.F2.GE.1.D0) THEN
          F3 = U(N)*XNEBOR(K)+V(N)*YNEBOR(K)
          IF(F3.LE.-1.D-2) THEN
            ALERTE2=.TRUE.
            IGUILT2=K
          ENDIF
        ENDIF
C
3     CONTINUE
C
      IF(ALERTE1.AND.ENTET) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'PROBLEME MAL POSE, VITESSE LIBRE ENTRANTE'
          WRITE(LU,*) 'PAR EXEMPLE AU POINT DE BORD NUMERO ',IGUILT1
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'ILL-POSED PROBLEM, ENTERING FREE VELOCITY'
          WRITE(LU,*) 'FOR EXAMPLE AT BOUNDARY POINT NUMBER ',IGUILT1
        ENDIF
      ENDIF
C
      IF(ALERTE2.AND.ENTET) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'PROBLEME MAL POSE, HAUTEUR LIBRE'
          WRITE(LU,*) 'SUR FRONTIERE AVEC VITESSE ENTRANTE'
          WRITE(LU,*) 'ET ECOULEMENT TORRENTIEL'
          WRITE(LU,*) 'PAR EXEMPLE AU POINT DE BORD NUMERO ',IGUILT2
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'ILL-POSED PROBLEM, FREE DEPTH'
          WRITE(LU,*) 'ON BOUNDARY WITH ENTERING VELOCITY'
          WRITE(LU,*) 'AND SUPERCRITICAL FLOW'
          WRITE(LU,*) 'FOR EXAMPLE AT BOUNDARY POINT NUMBER ',IGUILT2
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C
C INITIALISATION DES CONDITIONS AUX LIMITES POUR LA PROPAGATION :
C
C     INITIALISATION A ZERO DE TOUS LES VECTEURS DU BLOC
C
      CALL OS('X=0     ',X=MASK)
C
      DO 4 K=1,NPTFR
C
C     SI LE SUIVANT DE K N'EST PAS DANS LE SOUS-DOMAINE EN PARALLELISME
C     ON AURA KP1=K
      KP1=KP1BOR(K)
C
C-----------------------------------------------------------------------
C
C     CONDITIONS AUX LIMITES SUR LA HAUTEUR
C
      IF(CLH(K).EQ.KENT) THEN
        LIMPRO(K,1) = KDIR
        IF(KP1.NE.K) THEN
          IF(CLH(KP1).EQ.KENT) THEN
            MASK%ADR(HDIR)%P%R(K)=1.D0
          ELSEIF(CLH(KP1).EQ.KLOG) THEN
            MASK%ADR(HNEU)%P%R(K)=1.D0
          ELSEIF(CLH(KP1).EQ.KSORT) THEN
            MASK%ADR(HDDL)%P%R(K)=1.D0
          ELSE
            IF(LNG.EQ.1) WRITE(LU,10) K
            IF(LNG.EQ.2) WRITE(LU,11) K
            CALL PLANTE(1)
            STOP
          ENDIF
        ENDIF
      ELSEIF(CLH(K).EQ.KSORT) THEN
        LIMPRO(K,1) = KDDL
        IF(KP1.NE.K) THEN
          IF(CLH(KP1).EQ.KSORT) THEN
            MASK%ADR(HDDL)%P%R(K)=1.D0
          ELSEIF(CLH(KP1).EQ.KLOG) THEN
            MASK%ADR(HNEU)%P%R(K)=1.D0
          ELSEIF(CLH(KP1).EQ.KENT) THEN
            MASK%ADR(HDDL)%P%R(K)=1.D0
          ELSE
            IF(LNG.EQ.1) WRITE(LU,10) K
            IF(LNG.EQ.2) WRITE(LU,11) K
            CALL PLANTE(1)
            STOP
          ENDIF
        ENDIF
      ELSEIF(CLH(K).EQ.KLOG ) THEN
        LIMPRO(K,1) = KNEU
        IF(KP1.NE.K) MASK%ADR(HNEU)%P%R(K)=1.D0
      ELSE
        IF(LNG.EQ.1) WRITE(LU,10) K
        IF(LNG.EQ.2) WRITE(LU,11) K
        CALL PLANTE(1)
        STOP
      ENDIF
C
C   CONDITIONS AUX LIMITES SUR U
C
      IF(CLU(K).EQ.KENT.OR.CLU(K).EQ.KENTU) THEN
        LIMPRO(K,2) = KDIR
        IF(KP1.NE.K) THEN
          IF(CLU(KP1).EQ.KENT.OR.CLU(KP1).EQ.KENTU) THEN
            MASK%ADR(UDIR)%P%R(K)=1.D0
          ELSEIF(CLU(KP1).EQ.KADH) THEN
            MASK%ADR(UDIR)%P%R(K)=1.D0
          ELSEIF(CLU(KP1).EQ.KLOG) THEN
            MASK%ADR(UNEU)%P%R(K)=1.D0
          ELSEIF(CLU(KP1).EQ.KINC) THEN
            MASK%ADR(HOND)%P%R(K)=1.D0
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSEIF(CLU(KP1).EQ.KSORT) THEN
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSE
            IF(LNG.EQ.1) WRITE(LU,20) K
            IF(LNG.EQ.2) WRITE(LU,21) K
            CALL PLANTE(1)
            STOP
          ENDIF
        ENDIF
      ELSEIF(CLU(K).EQ.KADH) THEN
        LIMPRO(K,2) = KDIR
        IF(KP1.NE.K) MASK%ADR(UNEU)%P%R(K)=1.D0
      ELSEIF(CLU(K).EQ.KSORT) THEN
        LIMPRO(K,2) = KDDL
        IF(KP1.NE.K) THEN
          IF(CLU(KP1).EQ.KSORT) THEN
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSEIF(CLU(KP1).EQ.KLOG) THEN
            MASK%ADR(UNEU)%P%R(K)=1.D0
          ELSEIF(CLU(KP1).EQ.KINC) THEN
            MASK%ADR(HOND)%P%R(K)=1.D0
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSEIF(CLU(KP1).EQ.KENT) THEN
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSEIF(CLU(KP1).EQ.KENTU) THEN
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSEIF(CLU(KP1).EQ.KADH) THEN
            MASK%ADR(UNEU)%P%R(K)=1.D0
          ELSE
            IF(LNG.EQ.1) WRITE(LU,20) K
            IF(LNG.EQ.2) WRITE(LU,21) K
            CALL PLANTE(1)
            STOP
          ENDIF
        ENDIF
      ELSEIF(CLU(K).EQ.KLOG) THEN
        LIMPRO(K,2) = KDDL
        IF(KP1.NE.K) MASK%ADR(UNEU)%P%R(K)=1.D0
      ELSEIF(CLU(K).EQ.KINC ) THEN
        LIMPRO(K,2) = KDDL
        IF(KP1.NE.K) THEN
          IF(CLU(KP1).EQ.KINC) THEN
            MASK%ADR(HOND)%P%R(K)=1.D0
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSEIF(CLU(KP1).EQ.KLOG) THEN
            MASK%ADR(UNEU)%P%R(K)=1.D0
          ELSEIF(CLU(KP1).EQ.KSORT) THEN
            MASK%ADR(HOND)%P%R(K)=1.D0
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSEIF(CLU(KP1).EQ.KENT) THEN
            MASK%ADR(HOND)%P%R(K)=1.D0
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSEIF(CLU(KP1).EQ.KENTU) THEN
            MASK%ADR(HOND)%P%R(K)=1.D0
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSEIF(CLU(KP1).EQ.KADH) THEN
            MASK%ADR(HOND)%P%R(K)=1.D0
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSE
            IF(LNG.EQ.1) WRITE(LU,20) K
            IF(LNG.EQ.2) WRITE(LU,21) K
            CALL PLANTE(1)
            STOP
          ENDIF
        ENDIF
      ELSE
        IF(LNG.EQ.1) WRITE(LU,20) K
        IF(LNG.EQ.2) WRITE(LU,21) K
        CALL PLANTE(1)
        STOP
      ENDIF
C
C   CONDITIONS AUX LIMITES SUR V
C
      IF(CLV(K).EQ.KENT.OR.CLV(K).EQ.KENTU) THEN
        LIMPRO(K,3) = KDIR
        IF(KP1.NE.K) THEN
          IF(CLV(KP1).EQ.KENT.OR.CLV(KP1).EQ.KENTU) THEN
            MASK%ADR(VDIR)%P%R(K)=1.D0
          ELSEIF(CLV(KP1).EQ.KADH) THEN
            MASK%ADR(VDIR)%P%R(K)=1.D0
          ELSEIF(CLV(KP1).EQ.KLOG) THEN
            MASK%ADR(VNEU)%P%R(K)=1.D0
          ELSEIF(CLV(KP1).EQ.KINC) THEN
            MASK%ADR(HOND)%P%R(K)=1.D0
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSEIF(CLV(KP1).EQ.KSORT) THEN
            MASK%ADR(VDDL)%P%R(K)=1.D0
          ELSE
            IF(LNG.EQ.1) WRITE(LU,30) K
            IF(LNG.EQ.2) WRITE(LU,31) K
            CALL PLANTE(1)
            STOP
          ENDIF
        ENDIF
      ELSEIF(CLV(K).EQ.KADH) THEN
        LIMPRO(K,3) = KDIR
        IF(KP1.NE.K) MASK%ADR(VNEU)%P%R(K)=1.D0
      ELSEIF(CLV(K).EQ.KSORT) THEN
        LIMPRO(K,3) = KDDL
        IF(KP1.NE.K) THEN
          IF(CLV(KP1).EQ.KSORT) THEN
            MASK%ADR(VDDL)%P%R(K)=1.D0
          ELSEIF(CLV(KP1).EQ.KLOG) THEN
            MASK%ADR(VNEU)%P%R(K)=1.D0
          ELSEIF(CLV(KP1).EQ.KINC) THEN
            MASK%ADR(HOND)%P%R(K)=1.D0
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSEIF(CLV(KP1).EQ.KENT) THEN
            MASK%ADR(VDDL)%P%R(K)=1.D0
          ELSEIF(CLV(KP1).EQ.KENTU) THEN
            MASK%ADR(VDDL)%P%R(K)=1.D0
          ELSEIF(CLV(KP1).EQ.KADH) THEN
            MASK%ADR(VDDL)%P%R(K)=1.D0
          ELSE
            IF(LNG.EQ.1) WRITE(LU,30) K
            IF(LNG.EQ.2) WRITE(LU,31) K
            CALL PLANTE(1)
            STOP
          ENDIF
        ENDIF
      ELSEIF(CLV(K).EQ.KLOG ) THEN
        LIMPRO(K,3) = KDDL
        IF(KP1.NE.K) MASK%ADR(VNEU)%P%R(K)=1.D0
      ELSEIF(CLV(K).EQ.KINC ) THEN
        LIMPRO(K,3) = KDDL
        IF(KP1.NE.K) THEN
          IF(CLV(KP1).EQ.KINC) THEN
            MASK%ADR(HOND)%P%R(K)=1.D0
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSEIF(CLV(KP1).EQ.KLOG) THEN
            MASK%ADR(VNEU)%P%R(K)=1.D0
          ELSEIF(CLV(KP1).EQ.KSORT) THEN
            MASK%ADR(HOND)%P%R(K)=1.D0
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSEIF(CLV(KP1).EQ.KENT) THEN
            MASK%ADR(HOND)%P%R(K)=1.D0
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSEIF(CLV(KP1).EQ.KENTU) THEN
            MASK%ADR(HOND)%P%R(K)=1.D0
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSEIF(CLV(KP1).EQ.KADH) THEN
            MASK%ADR(HOND)%P%R(K)=1.D0
            MASK%ADR(UDDL)%P%R(K)=1.D0
          ELSE
            IF(LNG.EQ.1) WRITE(LU,30) K
            IF(LNG.EQ.2) WRITE(LU,31) K
            CALL PLANTE(1)
            STOP
          ENDIF
        ENDIF
      ELSE
        IF(LNG.EQ.1) WRITE(LU,30) K
        IF(LNG.EQ.2) WRITE(LU,31) K
        CALL PLANTE(1)
        STOP
      ENDIF
C
4     CONTINUE
C
C-----------------------------------------------------------------------
C
C     MASQUE DES FRONTIERES LIQUIDES
C
      DO K=1,NPTFR
        KP1=KP1BOR(K)
        IF(KP1.NE.K) MASK%ADR(UNONNEU)%P%R(K)=1.D0-MASK%ADR(UNEU)%P%R(K)
      ENDDO
C
C     ON DEDUIT DES MASQUES LES TABLEAUX LIMPRO(.,4 5 ET 6)
C
      DO K=1,NPTFR
        IF(MASK%ADR(HDIR)%P%R(K).GT.0.5D0) THEN
          LIMPRO(K,4)=KDIR
        ELSEIF(MASK%ADR(HDDL)%P%R(K).GT.0.5D0) THEN
          LIMPRO(K,4)=KDDL
        ELSEIF(MASK%ADR(HNEU)%P%R(K).GT.0.5D0) THEN
          LIMPRO(K,4)=KNEU
        ELSEIF(MASK%ADR(HOND)%P%R(K).GT.0.5D0) THEN
          LIMPRO(K,4)=KOND
        ELSE
          IF(NCSIZE.GT.1) THEN
            LIMPRO(K,4)=0
          ELSE
            IF(LNG.EQ.1) WRITE(LU,10) K
            IF(LNG.EQ.2) WRITE(LU,11) K
            CALL PLANTE(1)
            STOP
          ENDIF
        ENDIF
        IF(MASK%ADR(UDIR)%P%R(K).GT.0.5D0) THEN
          LIMPRO(K,5)=KDIR
        ELSEIF(MASK%ADR(UDDL)%P%R(K).GT.0.5D0) THEN
          LIMPRO(K,5)=KDDL
        ELSEIF(MASK%ADR(UNEU)%P%R(K).GT.0.5D0) THEN
          LIMPRO(K,5)=KNEU
        ELSEIF(MASK%ADR(HOND)%P%R(K).GT.0.5D0) THEN
          LIMPRO(K,5)=KOND
        ELSE
          IF(NCSIZE.GT.1) THEN
            LIMPRO(K,5)=0
          ELSE
            IF(LNG.EQ.1) WRITE(LU,20) K
            IF(LNG.EQ.2) WRITE(LU,21) K
            CALL PLANTE(1)
            STOP
          ENDIF
        ENDIF
        IF(MASK%ADR(VDIR)%P%R(K).GT.0.5D0) THEN
          LIMPRO(K,6)=KDIR
        ELSEIF(MASK%ADR(VDDL)%P%R(K).GT.0.5D0) THEN
          LIMPRO(K,6)=KDDL
        ELSEIF(MASK%ADR(VNEU)%P%R(K).GT.0.5D0) THEN
          LIMPRO(K,6)=KNEU
        ELSEIF(MASK%ADR(HOND)%P%R(K).GT.0.5D0) THEN
          LIMPRO(K,6)=KOND
        ELSE
          IF(NCSIZE.GT.1) THEN
            LIMPRO(K,6)=0
          ELSE
            IF(LNG.EQ.1) WRITE(LU,30) K
            IF(LNG.EQ.2) WRITE(LU,31) K
            CALL PLANTE(1)
            STOP
          ENDIF
        ENDIF
      ENDDO
C
C     COMPLEMENT POUR LES VITESSES QUADRATIQUES
C     LE POINT AU MILIEU D'UN SEGMENT A LA MEME CONDITION QUE LE SEGMENT
C
      IF(IELMU.EQ.13) THEN
        DO K=1,NPTFR
          LIMPRO(K+NPTFR,2)=LIMPRO(K,5)
          LIMPRO(K+NPTFR,3)=LIMPRO(K,6)
        ENDDO
      ENDIF
C
C     MASQUAGE PAR LE MASQUE DES ELEMENTS
C
      IF(MSK) THEN
        DO K=1,NPTFR
          IELEM=NELBOR(K)
          IF(IELEM.GT.0) THEN
            YY = MASKEL(IELEM)
            MASK%ADR(UDIR   )%P%R(K) = MASK%ADR(UDIR   )%P%R(K) * YY
            MASK%ADR(VDIR   )%P%R(K) = MASK%ADR(VDIR   )%P%R(K) * YY
            MASK%ADR(UDDL   )%P%R(K) = MASK%ADR(UDDL   )%P%R(K) * YY
            MASK%ADR(VDDL   )%P%R(K) = MASK%ADR(VDDL   )%P%R(K) * YY
            MASK%ADR(UNEU   )%P%R(K) = MASK%ADR(UNEU   )%P%R(K) * YY
            MASK%ADR(VNEU   )%P%R(K) = MASK%ADR(VNEU   )%P%R(K) * YY
            MASK%ADR(HOND   )%P%R(K) = MASK%ADR(HOND   )%P%R(K) * YY
            MASK%ADR(HDIR   )%P%R(K) = MASK%ADR(HDIR   )%P%R(K) * YY
            MASK%ADR(HNEU   )%P%R(K) = MASK%ADR(HNEU   )%P%R(K) * YY
            MASK%ADR(HDDL   )%P%R(K) = MASK%ADR(HDDL   )%P%R(K) * YY
            MASK%ADR(UNONNEU)%P%R(K) = MASK%ADR(UNONNEU)%P%R(K) * YY
          ENDIF
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------
C
10    FORMAT(1X,'PROPIN : POINT DE BORD',1I5,' CAS NON PREVU SUR H')
11    FORMAT(1X,'PROPIN: BOUNDARY POINT',1I5,' UNKNOWN PARAMETER FOR H')
20    FORMAT(1X,'PROPIN : POINT DE BORD',1I5,' CAS NON PREVU SUR U')
21    FORMAT(1X,'PROPIN: BOUNDARY POINT',1I5,' UNKNOWN PARAMETER FOR U')
30    FORMAT(1X,'PROPIN : POINT DE BORD',1I5,' CAS NON PREVU SUR V')
31    FORMAT(1X,'PROPIN: BOUNDARY POINT',1I5,' UNKNOWN PARAMETER FOR V')
C
C-----------------------------------------------------------------------
C
      RETURN
      END
