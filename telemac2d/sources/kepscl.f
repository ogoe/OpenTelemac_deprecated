C                       *****************
                        SUBROUTINE KEPSCL
C                       *****************
C
     *(KBOR,EBOR,AUBOR,CF,CFBOR,DISBOR,
     * UN,VN,HN,LIMKEP,LIUBOR,LIMPRO,NBOR,NPTFR,
     * KARMAN,CMU,C2,ESTAR,SCHMIT,LISRUG,PROPNU,KMIN,EMIN,
     * KNEU,KDIR,KENT,KENTU,KADH,KLOG)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2    27/11/92    J-M HERVOUET (LNH) 30 87 80 18
C                            26/04/94    L. VAN HAREN (LNH) 30 87 84 14
C***********************************************************************
C
C     FONCTION  : CALCUL DE KBOR , EBOR ET AUBOR QUAND LE MODELE DE
C                 TURBULENCE EST LE K-EPSILON.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      KBOR      |<-- | ENERGIE TURBULENTE IMPOSEE AU BORD           |
C |      EBOR      |<-- | DISSIPATION TURBULENTE IMPOSEE AU BORD       |
C |     AUBOR      |<-- | COEFFICIENT DE FROTTEMENT SUR LES PAROIS     |
C |      CF        | -->| COEFFICIENT DE FROTTEMENT POUR K-EPSILON     |
C |      VISC      | -->| DIFFUSION TURBULENTE (FAITE DANS VISTUR)     |
C |    DISBOR      | -->| DISTANCE AU BORD DES POINTS VOISINS DU BORD  |
C |    UN , VN     | -->| COMPOSANTES DE LA VITESSE AU TEMPS N         |
C |       HN       | -->| HAUTEUR D'EAU AU TEMPS N                     |
C |     LIMKEP     | -->| CONDITIONS AUX LIMITES SUR K ET EPSILON      |
C |     LIUBOR     | -->| CONDITIONS AUX LIMITES SUR U                 |
C |     LIMPRO     | -->| CONDITIONS AUX LIMITES EN PROPAGATION        |
C |      NBOR      | -->| ADRESSES DES POINTS DE BORD                  |
C |     NELBOR     | -->| NUMEROS DES ELEMENTS ADJACENTS AUX BORDS     |
C |     KP1BOR     | -->| NUMERO DU POINT SUIVANT SUR LE BORD          |
C |     NPTFR      | -->| NOMBRE DE POINTS FRONTIERES                  |
C |     KARMAN     | -->| CONSTANTE DE KARMAN                          |
C |     CMU        | -->| CONSTANTE DU MODELE K-EPSILON                |
C |     C1,C2      | -->| CONSTANTES DU MODELE K-EPSILON               |
C |     SIGMAK     | -->| CONSTANTE DU MODELE K-EPSILON                |
C |     SIGMAE     | -->| CONSTANTE DU MODELE K-EPSILON                |
C |     ESTAR      | -->| CONSTANTE DU MODELE K-EPSILON                |
C |     SCHMIT     | -->| CONSTANTE DU MODELE K-EPSILON                |
C |     LISRUG     | -->| REGIME DE TURBULENCE 1: LISSE 2: RUGUEUX     |
C |     GRAV       | -->| ACCELERATION DE LA PESANTEUR                 |
C |     PROPNU     | -->| COEFFICIENT DE DIFFUSION MOLECULAIRE         |
C |     KMIN,KMAX  | -->| K MINIMUM ET MAXIMUM EN CAS DE CLIPPING      |
C |     EMIN,EMAX  | -->| EPSILON MINIMUM ET MAXIMUM EN CAS DE CLIPPING|
C |     KNEU       | -->| CONDITION A LA LIMITE DE TYPE NEUMANN        |
C |     KDIR       | -->| CONDITION A LA LIMITE DE TYPE DIRICHLET      |
C |     KDDL       | -->| CONDITION A LA LIMITE DE DEGRE DE LIBERTE    |
C |     KOND       | -->| CONDITION A LA LIMITE DE TYPE ONDE INCIDENTE |
C |     KENT       | -->| CONVENTION POUR UNE ENTREE LIQUIDE           |
C |     KENTU      | -->| CONVENTION POUR DES VITESSES IMPOSEES        |
C |     KSORT      | -->| CONVENTION POUR UNE SORTIE LIQUIDE           |
C |     KADH       | -->| CONVENTION POUR UNE PAROI AVEC ADHERENCE     |
C |     KLOG       | -->| CONVENTION POUR UNE PAROI LOGARITHMIQUE      |
C |     KINC       | -->| CONVENTION POUR UNE ONDE INCIDENTE           |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : TELMAC
C
C  SOUS-PROGRAMME APPELES : OV
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C      
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NPTFR,LISRUG
      INTEGER, INTENT(IN) :: KNEU,KDIR,KENT,KADH,KLOG,KENTU
      INTEGER, INTENT(IN) :: LIMPRO(NPTFR,6),NBOR(NPTFR)
      INTEGER, INTENT(IN) :: LIMKEP(NPTFR,2),LIUBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: KMIN,EMIN
      DOUBLE PRECISION, INTENT(IN)    :: CF(*),CFBOR(*)
      DOUBLE PRECISION, INTENT(IN)    :: UN(*),VN(*),HN(*),DISBOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: KBOR(*),EBOR(*),AUBOR(*)
      DOUBLE PRECISION, INTENT(IN) :: KARMAN,CMU,C2,ESTAR,SCHMIT,PROPNU
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C      
      INTEGER N,K,IT
C
      DOUBLE PRECISION KFOND,EFOND,TIERS,UETUTA,CEPS,USTAR
      DOUBLE PRECISION SSQCMU,UTANG,DIST,DENOM,EBORD,KBORD 
C
C-----------------------------------------------------------------------
C
      INTRINSIC SQRT,MAX,LOG
C
C-----------------------------------------------------------------------
C
      TIERS  = 1.D0/3.D0
      SSQCMU = 1.D0/ SQRT(CMU)
C
C=======================================================================
C
C  BOUCLE SUR LA FRONTIERE
C
C  CALCUL DE LA VITESSE DE FROTTEMENT SUR LA PAROI
C
C  CALCUL DE KBOR,EBOR, ET AUBOR
C
C=======================================================================
C
      DO 1000 K=1,NPTFR
C
         KBOR(K) = 0.D0
         EBOR(K) = 0.D0
         N     = NBOR(K)
         UTANG = SQRT( UN(N)**2 + VN(N)**2 )
C        ATTENTION : MODIF PAR RAPPORT A LA NOTE DE PRINCIPE
C        DIST  = DISBOR(K)*0.1D0
         DIST  = DISBOR(K)*0.33D0
C
C        CALCUL DE UETOIL POUR LES PAROIS SOLIDES
C        ----------------------------------------
C
C                                   UETOIL
C        UETUTA REPRESENTE PARTOUT  ------
C                                   UTANG
C
C        UETUTA A L'AVANTAGE DE GARDER UN SENS MEME SI UTANG=0
C
C        ********************
         IF(LISRUG.EQ.1) THEN
C        ********************
C
C           TIR INITIAL ; PUIS 5 ITERATIONS
            UETUTA = 6.D-2
            DO 60 IT=1,5
C
             IF(DIST*UETUTA*UTANG/PROPNU .LT. 30.D0) THEN
               UETUTA = 7.25D-2
             ELSE
               UETUTA=1.D0/(5.5D0+LOG(DIST*UETUTA*UTANG/PROPNU)/KARMAN)
             ENDIF
C
60          CONTINUE
C
C        ************************
         ELSEIF(LISRUG.EQ.2) THEN
C        ************************
C
            UETUTA = SQRT( 0.5D0 * CFBOR(K) )
C
C        ****
         ELSE
C        ****
C
            IF(LNG.EQ.1) WRITE(LU,400) LISRUG
            IF(LNG.EQ.2) WRITE(LU,401) LISRUG
400         FORMAT(1X,'KEPSCL : REGIME DE TURBULENCE INCONNU : ',1I6)
401         FORMAT(1X,'KEPSCL : UNKNOWN TURBULENCE MODEL : ',1I6)
            CALL PLANTE(1)
            STOP
C
C        *****
         ENDIF
C        *****
C
C                                                        DIRICHLET SUR K
C                                                        ---------------
C
C
         IF(LIMKEP(K,1).EQ.KDIR) THEN
C        ----------------------------
C
C           ************************************************
            IF(LIUBOR(K).EQ.KENT.OR.LIUBOR(K).EQ.KENTU) THEN
C           ************************************************
C
C              ENTREE DE DOMAINE : TURBULENCE DUE AU FOND
C
               CEPS    = C2*SQRT(CMU)/SQRT(ESTAR*SCHMIT) /
     *                   (0.5D0*CF(N))**0.75D0
               DENOM   = CEPS * 0.5D0*CF(N)
               USTAR = SQRT( 0.5D0 * CF(N) * ( UN(N)**2 + VN(N)**2 ) )
               KBOR(K) = C2 * USTAR**2 / MAX(DENOM,1.D-10)
C
C           ***************************************************
            ELSEIF(LIUBOR(K).EQ.KLOG.OR.LIUBOR(K).EQ.KADH) THEN
C           ***************************************************
C
C              PAROI
C
               CEPS    = C2*SQRT(CMU)/SQRT(ESTAR*SCHMIT) /
     *                   (0.5D0*CF(N))**0.75D0
               DENOM   = CEPS * 0.5D0*CF(N)
               USTAR = SQRT( 0.5D0 * CF(N) * ( UN(N)**2 + VN(N)**2 ) )
               KFOND   = C2 * USTAR**2 / MAX(DENOM,1.D-10)
               KBORD   = SSQCMU*(UETUTA*UTANG)**2
               KBOR(K) = KBORD + KFOND
C
C           ****
            ELSE
C           ****
C
               IF(LNG.EQ.1) WRITE(LU,500) K,LIUBOR(K)
               IF(LNG.EQ.2) WRITE(LU,501) K,LIUBOR(K)
500            FORMAT(1X,'KEPSCL: POINT DE BORD ',1I6,
     *                   'CAS NON PREVU POUR KBOR',1X,'LIUBOR=',1I6)
501            FORMAT(1X,'KEPSCL: BOUNDARY POINT ',1I6,
     *                   'UNKNOWN CASE FOR KBOR',1X,'LIUBOR=',1I6)
               CALL PLANTE(1)
               STOP
C
C           *****
            ENDIF
C           *****
C
         ENDIF
C        -----
C
C                                                  DIRICHLET SUR EPSILON
C                                                  ---------------------
C
         IF(LIMKEP(K,2).EQ.KDIR) THEN
C        ----------------------------
C
C           ************************************************
            IF(LIUBOR(K).EQ.KENT.OR.LIUBOR(K).EQ.KENTU) THEN
C           ************************************************
C
C              ENTREE DE DOMAINE : TURBULENCE DUE AU FOND
C
               DENOM   = SQRT(0.5D0*CF(N)) * HN(N)
               USTAR   = SQRT(0.5D0*CF(N) * ( UN(N)**2 + VN(N)**2 ) )
               EFOND   = USTAR**3 / MAX(DENOM,1.D-10)
               EBOR(K) = MAX( EFOND , EMIN )
C
C           ***************************************************
            ELSEIF(LIUBOR(K).EQ.KLOG.OR.LIUBOR(K).EQ.KADH) THEN
C           ***************************************************
C
C              PAROI
C
               DENOM   = SQRT(0.5D0*CF(N)) * HN(N)
               USTAR   = SQRT(0.5D0*CF(N) * ( UN(N)**2 + VN(N)**2 ) )
               EFOND   = USTAR**3 / MAX(DENOM,1.D-10)
               EBORD   = (UETUTA*UTANG)**3 / ( KARMAN*DIST )
               EBOR(K) = MAX( EBORD + EFOND, EMIN )
C
C           ****
            ELSE
C           ****
C
C              AUTRE
C
               IF(LNG.EQ.1) WRITE(LU,600) K,LIUBOR(K)
               IF(LNG.EQ.2) WRITE(LU,601) K,LIUBOR(K)
600            FORMAT(1X,'KEPSCL: POINT DE BORD ',1I6,
     *                   'CAS NON PREVU POUR EBOR',1X,'LIUBOR=',1I6)
601            FORMAT(1X,'KEPSCL: BOUNDARY POINT ',1I6,
     *                   'UNKNOWN CASE FOR EBOR',1X,'LIUBOR=',1I6)
               CALL PLANTE(1)
               STOP
C
C           *****
            ENDIF
C           *****
C
         ENDIF
C        -----
C
C                                                        CALCUL DE AUBOR
C                                                        ---------------
C
C
C  AUBOR COMPTE POUR LE SEGMENT ENTRE K ET KP1BOR(K)
C
C  LOI UTILISEE  : NUT * DU/DN = UETOIL**2 = -AUBOR*U(N+1)
C  TRANSFORME EN : NUT * DU/DN = UETOIL**2  *  U(N+1) / U(N)
C                              = UETOIL * (UETOIL/UTANG) * U(N+1)
C
         IF (LIMPRO(K,5).EQ.KNEU) THEN
            AUBOR(K) = - UTANG * UETUTA**2
         ELSE
            AUBOR(K) = 0.D0
         ENDIF
C
1000  CONTINUE
C
C=======================================================================
C
C                        /* FIN BOUCLE SUR LA FRONTIERE */
C
C=======================================================================
C
      RETURN
      END 
