C                       *****************
                        SUBROUTINE CHARAC
C                       *****************
C
     *( FN  , FTILD  , NOMB   , UCONV  , VCONV , WCONV  , ZSTAR ,
     *  DT  , IFAMAS , IELM   , NPOIN2 , NPLAN , NPLINT ,
     *  MSK , MASKEL , SHP,SHZ , TB    , IT1,IT2,IT3,IT4,MESH ,
     *  NELEM2,NELMAX2,IKLE2,SURDET2   , INILOC)
C
C***********************************************************************
C BIEF VERSION 6.0      12/02/2010    J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C  FONCTION : APPEL DE LA METHODE DES CARACTERISTIQUES
C             (SOUS-PROGRAMME CARACT)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   FN           | -->| VARIABLES A L'ETAPE N .
C |   FTILD        |<-- | VARIABLES APRES LA CONVECTION .
C |   NOMB         | -->| NOMBRE DE VARIABLES A CONVECTER.
C |   UCONV,VCONV..| -->| COMPOSANTES DES VITESSES DU CONVECTEUR.
C |   ZSTAR        | -->| COORDONNEES VERTICALES EN 3D.
C |   DT           | -->| PAS DE TEMPS
C |   IFAMAS       | -->| IFABOR MODIFIE QUAND DES ELEMENTS SONT MASQUES
C |   IELM         | -->| TYPE D'ELEMENT : 11 : TRIANGLE P1
C |                |    |                  41 : PRISME DE TEL3D
C |   NPOIN2       | -->| NOMBRE DE POINTS DU MAILLAGE 2D (POUR TEL3D).
C |   NPLAN        | -->| NOMBRE DE PLAN SUIVANT Z (POUR TEL3D).
C |   NPLINT       | -->| PLAN DE REFERENCE INTERMEDIAIRE (POUR TEL3D).
C |   MSK          | -->| SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |   MASKEL       | -->| TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE.
C |   MAT          | -->| MATRICE DE TRAVAIL
C |   TB           | -->| BLOC DE TABLEAUX DE TRAVAIL (AU MOINS 8)
C |   MESH         | -->| BLOC DES ENTIERS DU MAILLAGE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : CARACT
C
C**********************************************************************
C
      USE BIEF, EX_CHARAC => CHARAC
      USE STREAMLINE, ONLY : SCARACT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(IN)         :: NOMB
      INTEGER         , INTENT(IN)         :: NPLAN,NPLINT,NELEM2
      INTEGER         , INTENT(IN)         :: NPOIN2,IELM,NELMAX2
      INTEGER         , INTENT(INOUT)      :: IT1(*),IT2(*)
      INTEGER         , INTENT(INOUT)      :: IT3(*),IT4(*)
      TYPE(BIEF_OBJ)  , INTENT(IN)         :: FN,UCONV,VCONV,WCONV
      TYPE(BIEF_OBJ)  , INTENT(IN)         :: ZSTAR,MASKEL,IKLE2,SURDET2
      TYPE(BIEF_OBJ)  , INTENT(INOUT)      :: FTILD,TB,SHP,SHZ
      LOGICAL         , INTENT(IN)         :: MSK
      DOUBLE PRECISION, INTENT(IN)         :: DT
      TYPE(BIEF_MESH) , INTENT(INOUT)      :: MESH
      TYPE(BIEF_OBJ)  , INTENT(IN), TARGET :: IFAMAS
      LOGICAL, OPTIONAL, INTENT(IN)        :: INILOC
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NPOIN,IELMU
      LOGICAL INITLOC
C
C-----------------------------------------------------------------------
C
      TYPE(BIEF_OBJ), POINTER :: T1,T2,T3,T4,T5,T6,T7
      INTEGER, DIMENSION(:), POINTER :: IFA
      INTEGER I,J,K
C
C-----------------------------------------------------------------------
C  TABLEAUX DE TRAVAIL PRIS DANS LE BLOC TB
C-----------------------------------------------------------------------
C
      T1 =>TB%ADR( 1)%P
      T2 =>TB%ADR( 2)%P
      T3 =>TB%ADR( 3)%P
      T4 =>TB%ADR( 4)%P
      T5 =>TB%ADR( 5)%P
      T6 =>TB%ADR( 6)%P
      T7 =>TB%ADR( 7)%P
C
C-----------------------------------------------------------------------
C  INITIALISING THE LOCATION OF POINTS OR NOT
C-----------------------------------------------------------------------
C
      IF(PRESENT(INILOC)) THEN
        INITLOC=INILOC
      ELSE
        INITLOC=.TRUE.
      ENDIF
C
C-----------------------------------------------------------------------
C  DEPLOIEMENT DE LA STRUCTURE DE MAILLAGE
C-----------------------------------------------------------------------
C
      NPOIN = MESH%NPOIN
      IELMU = UCONV%ELM
C
C-----------------------------------------------------------------------
C     CHECKING SHP SIZE (ONCE A BUG...)
C-----------------------------------------------------------------------
C
      IF(3*NPOIN.GT.SHP%MAXDIM1*SHP%MAXDIM2) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'TAILLE DE SHP:',SHP%MAXDIM1*SHP%MAXDIM2
          WRITE(LU,*) 'TROP PETITE DANS CHARAC, ',3*NPOIN
          WRITE(LU,*) 'REQUISE'
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'SIZE OF SHP:',SHP%MAXDIM1*SHP%MAXDIM2
          WRITE(LU,*) 'TOO SMALL IN CHARAC, ',3*NPOIN
          WRITE(LU,*) 'REQUESTED'
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF   
C
C-----------------------------------------------------------------------
C  APPEL DE CARACT
C-----------------------------------------------------------------------
C
      IF(MSK) THEN
C       APPEL AVEC IFAMAS
        IFA=>IFAMAS%I
      ELSE
C       APPEL AVEC IFABOR
        IFA=>MESH%IFABOR%I
      ENDIF
!
      CALL CPSTVC(MESH%X,T1)
      CALL CPSTVC(MESH%Y,T2)
!
      IF(NCSIZE.EQ.0) THEN
!
        CALL CARACT( FN , FTILD , UCONV%R , VCONV%R , WCONV%R ,
     *               MESH%X%R,MESH%Y%R,ZSTAR%R,
     *               T1,T2,T3%R,T4%R,T5%R,T6%R,
     *               MESH%Z%R,SHP%R,SHZ%R,
     *               SURDET2%R,DT,IKLE2%I,IFA,
     *               IT1,IT2,IT3,IT4,
     *               IELM,IELMU,NELEM2,NELMAX2,NOMB,NPOIN,NPOIN2,
     *               3,NPLAN,MESH%LV,
     *               MSK,MASKEL%R,MESH,MESH%FAC%R,T7%R,T7,INITLOC)
!
      ELSEIF(NCSIZE.GE.1) THEN
!
        CALL SCARACT( FN , FTILD , UCONV%R , VCONV%R , WCONV%R ,
     &                MESH%X%R,MESH%Y%R,ZSTAR%R,
     &                T1,T2,T3%R,T4%R,T5%R,T6%R,
     &                MESH%Z%R,SHP%R,SHZ%R,
     &                SURDET2%R,DT,IKLE2%I,IFA,
     &                IT1,IT2,IT3,IT4,
     &                IELM,IELMU,NELEM2,NELMAX2,NOMB,NPOIN,NPOIN2,
     &                3,NPLAN,MESH%LV,MSK,MASKEL%R,
     &                MESH,MESH%FAC%R,T7%R,T7,INITLOC)
!
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
