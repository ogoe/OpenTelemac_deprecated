C                       *****************
                        SUBROUTINE CARACT
C                       *****************
C
     * ( U , UTILD , UCONV , VCONV , WCONV , X , Y , ZSTAR ,
     *   T1 , T2 , ZCONV , DX , DY , DZ , Z , SHP , SHZ , SURDET ,
     *   DT , IKLE , IFABOR , ELT , ETA , ITRAV1 , ITRAV2 , IELM ,IELMU,
     *   NELEM , NELMAX , NOMB , NPOIN , NPOIN2 , NDP , NPLAN , 
     *   LV , MSK , MASKEL , MESH , FAC , TEST , STEST,INITLOC)
C
C***********************************************************************
C BIEF VERSION 5.9        27/08/08    J-M HERVOUET (LNHE) 01 30 87 80 18
C                                        ALGIANE FROEHLY (MATMECA)
C***********************************************************************
C
C     FONCTION:
C
C     RESOUT LES EQUATIONS DE CONVECTION PAR LA METHODE DES
C     CARACTERISTIQUES, POUR UN ENSEMBLE DE FONCTIONS.
C
C     ATTENTION LA COMPATIBILITE AVEC 3.0 N'EST PLUS ASSUREE A CAUSE
C     DE L'APPEL A PARCOM
C
C     EN REVANCHE U ET UTILD PEUVENT MAINTENANT ETRE DES VECTEURS,
C     DANS CE CAS NOMB SERA CONSIDERE EGAL A 1
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   U            | -->| VARIABLES A L'ETAPE N .                      |
C |   UTILD        |<-- | VARIABLES APRES LA CONVECTION .              |
C |   UCONV,VCONV..| -->| COMPOSANTES DES VITESSES DU CONVECTEUR.      |
C |   X,Y,ZSTAR    | -->| COORDONNEES DU MAILLAGE .                    |
C |   XCONV,YCONV..| -- | COORDONNEES AU PIED DES CARACTERISTIQUES.    |
C |   DX,DY,DZ     | -- | STOCKAGE DES SOUS-PAS .                      |
C |   Z            | -->| COTE DANS LE MAILLAGE REEL (POUR TEL3D) .    |
C |   SHP          | -- | COORDONNEES BARYCENTRIQUES 2D AU PIED DES    |
C |                |    | COURBES CARACTERISTIQUES.                    |
C |   SHZ          | -- | COORDONNEES BARYCENTRIQUES SUIVANT Z AU PIED |
C |                |    | DES COURBES CARACTERISTIQUES (POUR TEL3D)    |
C |   SURDET       | -->| 1/DETERMINANT POUR LES ELEMENTS 2D.          |
C |   DT           | -->| PAS DE TEMPS                                 |
C |   IKLE         | -->| NUMEROS GLOBAUX DES POINTS DES ELEMENTS 2D.  |
C |   IFABOR       | -->| NUMEROS DES ELEMENTS VOISINS (ATTENTION, POUR|
C |                |    | TEL3D, IFABOR EST LE TABLEAU IBOR DE MITRID).|
C |   ELT          | -- | NUMEROS DES ELEMENTS 2D AU PIED DES COURBES  |
C |                |    | CARACTERISTIQUES.                            |
C |   ETA          | -- | NUMEROS DES ETAGES AU PIED DES COURBES       |
C |                |    | CARACTERISTIQUES (POUR TEL3D).               |
C |   ITRAV1       | -- | TABLEAU DE TRAVAIL ENTIER.                   |
C |   ITRAV2       | -- | TABLEAU DE TRAVAIL ENTIER.                   |
C |   IELM         | -->| TYPE D'ELEMENT : 11 : TRIANGLE P1            |
C |                |    |                  21 : QUADRANGLE P1          |
C |                |    |                  41 : PRISME DE TEL3D        |
C |   NELEM        | -->| NOMBRE TOTAL D'ELEMENTS DANS LE MAILLAGE 2D. |
C |   NELMAX       | -->| NOMBRE MAXIMAL D'ELEMENTS DANS LE MAILLAGE 2D|
C |   NOMB         | -->| NOMBRE DE VARIABLES A CONVECTER.             |
C |   NPOIN        | -->| NOMBRE TOTAL DE POINTS DU MAILLAGE.          |
C |   NPOIN2       | -->| NOMBRE DE POINTS DU MAILLAGE 2D (POUR TEL3D).|
C |   NDP          | -->| NOMBRE DE POINTS PAR ELEMENT 2D.             |
C |   NPLAN        | -->| NOMBRE DE PLAN SUIVANT Z (POUR TEL3D).       |
C |   NPLINT       | -->| PLAN DE REFERENCE INTERMEDIAIRE (POUR TEL3D).|
C |   LV           | -->| LONGUEUR DU VECTEUR POUR LA VECTORISATION.   |
C |   MSK          | -->| SI OUI, PRESENCE D'ELEMENTS MASQUES.         |
C |   MASKEL       | -->| TABLEAU DE MASQUAGE DES ELEMENTS             |
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE.        |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC , MITRID
C
C SOUS-PROGRAMMES APPELES : CHAR11 , CHAR41 ,
C                           GTSH11 , GTSH41
C
C***********************************************************************
C
      USE BIEF, EX_CARACT => CARACT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NELEM,NELMAX,NOMB,NPOIN,NPOIN2
      INTEGER, INTENT(IN)             :: NDP,NPLAN,IELM,IELMU,LV
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: U
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: UTILD
      DOUBLE PRECISION, INTENT(IN)    :: UCONV(*),VCONV(*),WCONV(*)
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: Y(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: Z(NPOIN2,NPLAN),ZSTAR(NPLAN)
      DOUBLE PRECISION, INTENT(OUT)   :: ZCONV(NPOIN2,NPLAN)
      DOUBLE PRECISION, INTENT(OUT)   :: DX(*),DY(*),DZ(NPOIN)
      DOUBLE PRECISION, INTENT(OUT)   :: SHP(NDP,*),SHZ(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: MASKEL(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: SURDET(NELEM)
      DOUBLE PRECISION, INTENT(IN)    :: DT
      INTEGER, INTENT(IN)             :: IKLE(NELMAX,*)
      INTEGER, INTENT(IN)             :: IFABOR(NELMAX,*)
      INTEGER, INTENT(OUT)            :: ELT(*),ETA(*)
      INTEGER, INTENT(OUT)            :: ITRAV1(*),ITRAV2(*)
      LOGICAL, INTENT(IN)             :: MSK,INITLOC
      DOUBLE PRECISION, INTENT(OUT)   :: TEST(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: FAC(NPOIN)
      TYPE(BIEF_MESH) , INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: STEST
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: T1,T2
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NRK,I,NPOINT,NPOINT2,IPOIN,IELMI
      DOUBLE PRECISION C
      LOGICAL QUAD
C
C***********************************************************************
C
C NOMBRE DE SOUS-PAS DE RUNGE-KUTTA PAR ELEMENT TRAVERSE
C
      QUAD = .FALSE.
      NRK = 3
      CALL OV( 'X=Y     ' , T1%R , X , Z , C , NPOIN )
      CALL OV( 'X=Y     ' , T2%R , Y , Z , C , NPOIN )
C
C-----------------------------------------------------------------------
C
      IF(NCSIZE.GT.1) CALL OV('X=C     ',TEST,Y,Z,1.D0,NPOIN)
C
C-----------------------------------------------------------------------
C
      IF(IELM.EQ.11) THEN
C
C ON RECHERCHE SI L'UN DES ELEMENTS A CONVECTER EST DISCRETISE P2
C      
        DO I=1,U%N
          IF(U%ADR(I)%P%ELM.EQ.13) THEN
            QUAD = .TRUE.
          ENDIF
        ENDDO
C      
        IF(QUAD) THEN       
          CALL CHGDIS(T1,IELM,13,MESH)
          CALL CHGDIS(T2,IELM,13,MESH)
        ENDIF
C
C-----------------------------------------------------------------------
C
        IF(.NOT.QUAD) THEN
C
C    TRIANGLES P1
C    ============
C
C       REMPLISSAGE DES SHP ET DES ELT OPTIMISE
C
        NPOINT  = NPOIN
        NPOINT2 = NPOIN
        IF(INITLOC) THEN
          CALL GTSH11(UCONV,VCONV,X,Y,SHP,ELT,IKLE,
     *                ITRAV1,ITRAV2,NPOINT,NELEM,NELMAX,LV,MSK,MASKEL)
        ENDIF
C
C       APPEL DU SOUS-PROGRAMME DE REMONTEE DES COURBES CARATERISTIQUES
C 
        IF(NCSIZE.GT.1) CALL OV('X=C     ',TEST,Y,Z,1.D0,NPOINT)    
        CALL CHAR11(UCONV,VCONV,DT,NRK,X,Y,IKLE,
     *              IFABOR,T1%R,T2%R,DX,DY,SHP,ELT,ITRAV1,
     *              NPOINT,NPOINT2,NELEM,NELMAX,SURDET,-1,TEST)
C
C-----------------------------------------------------------------------
C
        ELSEIF(QUAD) THEN
C
C         TRIANGLES P2 POUR L'UNE DES VARIABLES CONVECTEE
C         ===============================================
C
C         REMPLISSAGE DES SHP ET DES ELT OPTIMISE
C
          NPOINT  = NPOIN+MESH%NSEG 
          NPOINT2 = NPOIN+MESH%NSEG 
C
C         CAS NON PREVU D'UN TRACEUR QUADRATIQUE ET D'UNE VITESSE LINEAIRE
C      
          IF(IELMU.NE.13)THEN
            IF(LNG.EQ.1) WRITE(LU,21) 
            IF(LNG.EQ.2) WRITE(LU,22) 
            CALL PLANTE(1)
            STOP
          ENDIF
C 
          IF(INITLOC) THEN      
            CALL GTSH13( UCONV , VCONV , X , Y, SHP , ELT , IKLE,
     *                   ITRAV1  , ITRAV2  , NPOINT2, NELEM , 
     *                   NELMAX  , LV ,MSK , MASKEL )
          ENDIF
C
C         APPEL DU SOUS-PROGRAMME DE REMONTEE DES COURBES CARACTERISTIQUES
C
C         CALL CHAR13( UCONV , VCONV , DT    , NRK , X , Y ,
          IF(NCSIZE.GT.1) CALL OV('X=C     ',TEST,Y,Z,1.D0,NPOINT)  
          CALL CHAR11( UCONV , VCONV , DT    , NRK , X , Y ,
     *                 IKLE    , IFABOR  ,
     *                 T1%R    , T2%R    , DX    , DY ,SHP,ELT ,ITRAV1,
     *                 NPOINT  , NPOINT2 , NELEM ,
     *                 NELMAX  , SURDET  , -1    , TEST)
C  
        ENDIF
C
C-----------------------------------------------------------------------
C   
      ELSEIF(IELM.EQ.41) THEN
C
C    PRISMES DE TELEMAC-3D
C    =====================
C
        NPOINT = NPOIN
        NPOINT2 = NPOIN2
        DO I = 1,NPLAN
          CALL OV('X=C     ' ,ZCONV(1,I),Y,Z,ZSTAR(I),NPOIN2)
        ENDDO
C
C      REMPLISSAGE DES SHP ET DES ELT OPTIMISE
C
        IF(INITLOC) THEN
          CALL GTSH41( UCONV , VCONV , WCONV , X , Y , SHP , SHZ ,
     *                 ELT , ETA , IKLE , ITRAV1 , ITRAV2 , NPOINT2 ,
     *                 NELEM , NPLAN , LV , MSK , MASKEL )
        ENDIF
C
C      APPEL DU SOUS-PROGRAMME DE REMONTEE DES COURBES CARATERISTIQUES
C
         IF(NCSIZE.GT.1) CALL OV('X=C     ',TEST,Y,Z,1.D0,NPOINT)
         CALL CHAR41( UCONV , VCONV , WCONV , DT , NRK ,
     *                X , Y , ZSTAR , Z , IKLE , IFABOR , T1%R , T2%R , 
     *                ZCONV , DX , DY , DZ , SHP , SHZ , ELT , ETA , 
     *                ITRAV1 , NPOINT ,
     *                NPOINT2 , NELEM , NPLAN , SURDET , -1 ,
     *                ITRAV2 , TEST )
C
C-----------------------------------------------------------------------
C
      ELSE
C
        IF(LNG.EQ.1) WRITE(LU,11) IELM
        IF(LNG.EQ.2) WRITE(LU,12) IELM
        CALL PLANTE(1)
        STOP
C
      ENDIF

C
C  PROVISOIRE
C  TEST = NOMBRE DE SOUS-DOMAINES QUI ONT TRAITE UN POINT
C
C-----------------------------------------------------------------------
C
C  AJUSTEMENT DES SHP EN FONCTION DE CE QUE L'ON A TROUVE
C  DANS LES AUTRES SOUS-DOMAINES (DONNE PAR TEST)
C
      IF(NCSIZE.GT.1) THEN
C
        IF(QUAD) THEN
          IF(LNG.EQ.1) WRITE(LU,19) 
          IF(LNG.EQ.2) WRITE(LU,20) 
          CALL PLANTE(1)
          STOP
        ENDIF
C
        DO IPOIN = 1,NPOIN
C         A CHARACTERISTIC WHICH DID NOT STOP IN THE SUB-DOMAIN
C         WILL GIVE 0
          IF(TEST(IPOIN).LT.0.5D0) THEN
            SHP(1,IPOIN) = 0.D0
            SHP(2,IPOIN) = 0.D0
            SHP(3,IPOIN) = 0.D0
          ENDIF
        ENDDO
C
C       THIS IS JUST FOR PRINTING A WARNING
C
        CALL PARCOM(STEST,2,MESH)
        DO IPOIN = 1,NPOIN
          IF(TEST(IPOIN).LT.0.5D0) THEN
            IF(LNG.EQ.1) WRITE(LU,13) IPOIN
            IF(LNG.EQ.2) WRITE(LU,14) IPOIN
          ENDIF
        ENDDO
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  INTERPOLATION AU PIED DES CARACTERISTIQUES SI DEMANDEE
C
      IF(NOMB.GT.0) THEN
C
        IF(U%TYPE.EQ.2.AND.UTILD%TYPE.EQ.2) THEN
C
C         U ET UTILD VECTEURS (NOMB VAUT ALORS 1)
C
          CALL INTERP(U%R,UTILD%R,SHP,NDP,SHZ,ETA,ELT,
     *                NPOINT,NPOINT2,NPLAN,IELM,IKLE,NELMAX)
C
          IF(NCSIZE.GT.1) THEN
            IF(QUAD) THEN
              IF(LNG.EQ.1) WRITE(LU,19) 
              IF(LNG.EQ.2) WRITE(LU,20) 
              CALL PLANTE(1)
              STOP
            ENDIF
C           CHOOSING THE RESULT WITH MAXIMUM
C           ABSOLUTE VALUE
            CALL PARCOM(UTILD,1,MESH)             
          ENDIF 
C
        ELSEIF(U%TYPE.EQ.4.AND.UTILD%TYPE.EQ.4) THEN
C
C     U ET UTILD BLOCS DE VECTEURS
C
         IF(U%N.LT.NOMB.OR.UTILD%N.LT.NOMB) THEN
            IF(LNG.EQ.1) WRITE(LU,15) U%N,UTILD%N
            IF(LNG.EQ.2) WRITE(LU,16) U%N,UTILD%N
            CALL PLANTE(1)
            STOP
          ENDIF
C
          DO 61 I=1,NOMB
C
C INTERPOLATION DES VITESSES EN PRESENCE D'UNE VARIABLE QUADRATIQUE
C          
          IF(QUAD) THEN
            IELMI = U%ADR(I)%P%ELM 
            NPOINT2=U%ADR(I)%P%DIM1
            CALL INTERP(U%ADR(I)%P%R,UTILD%ADR(I)%P%R,SHP,NDP,SHZ,
     *           ETA,ELT,NPOINT2,NPOINT2,NPLAN,IELMI,IKLE,NELMAX)
C
C INTERPOLATION DES VITESSES DANS LES AUTRES CAS
C         
          ELSE 
            CALL INTERP(U%ADR(I)%P%R,UTILD%ADR(I)%P%R,SHP,NDP,SHZ,
     *            ETA,ELT,NPOINT,NPOINT2,NPLAN,IELM,IKLE,NELMAX)
          ENDIF
C
          IF(NCSIZE.GT.1) THEN
            IF(QUAD)THEN
              IF(LNG.EQ.1) WRITE(LU,19) 
              IF(LNG.EQ.2) WRITE(LU,20) 
              CALL PLANTE(1)
              STOP
            ELSE
C             CHOOSING THE RESULT WITH MAXIMUM
C             ABSOLUTE VALUE
              CALL PARCOM(UTILD%ADR(I)%P,1,MESH)
            ENDIF
          ENDIF
C
61        CONTINUE
C
        ELSE
C
          IF(LNG.EQ.1) WRITE(LU,17) U%TYPE,UTILD%TYPE
          IF(LNG.EQ.2) WRITE(LU,18) U%TYPE,UTILD%TYPE
          CALL PLANTE(1)
C
        ENDIF
C
      ENDIF
C
C-----------------------------------------------------------------------
C
11    FORMAT(1X,'CARACT : TYPE D''ELEMENT INCONNU : ',I6)
12    FORMAT(1X,'CARACT: UNKNOWN TYPE OF ELEMENT : ',I6)
C
13    FORMAT(1X,'CARACT : (PARALLELE) REMONTEE INCOMPLETE POUR : ',I6)
14    FORMAT(1X,'CARACT: (PARALLEL) INCOMPLETE PATH LINE FOR : ',I6)
C
15    FORMAT(1X,'CARACT : MAUVAIS BLOC DES VARIABLES : ',2I6)
16    FORMAT(1X,'CARACT: WRONG BLOCK OF VARIABLES : ',2I6)
C
17    FORMAT(1X,'CARACT : TYPE D''OBJET INCONNU : ',2I6)
18    FORMAT(1X,'CARACT: UNKNOWN TYPE OF OBJECT: ',2I6)
C
19    FORMAT(1X,'CARACT : PARALLELISME NON PREVU EN QUADRATIQUE')
20    FORMAT(1X,'CARACT : PARALLELISM NOT TREATED WITH QUADRATIC ')
C
21    FORMAT(1X,'CARACT : VITESSES LINEAIRES ET TRACEUR QUADRATIQUE')
22    FORMAT(1X,'CARACT : LINEAR VELOCITY AND QUADRATIC TRACER')
C
C-----------------------------------------------------------------------
C
      RETURN
      END
