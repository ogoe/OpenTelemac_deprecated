C                       *****************
                        SUBROUTINE RESCUE
C                       *****************
C
     *(U,V,H,S,ZF,T,TRAC0,NTRAC,ITURB,NPOIN,AKEP,TROUVE)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.8   31/08/07   J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C   FONCTION  : RECONSTITUTION DE DONNEES MANQUANTES EN CAS DE SUITE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   U , V        |<-- | COMPOSANTES DES VECTEURS VITESSE.            |
C |   C            |<-- | CELERITE.                                    |
C |   H            |<-- | HAUTEUR.                                     |
C |   S            |<-- | SURFACE LIBRE.                               |
C |   ZF           |<-- | COTE DES POINTS DU FOND.                     |
C |   F            |<-- | NOMBRE DE FROUDE.                            |
C |   Q            |<-- | DEBIT SCALAIRE.                              |
C |   T            |<-- | TRACEUR.                                     |
C |   AK           |<-- | ENERGIE TURBULENTE.                          |
C |   EP           |<-- | DISSIPATION.                                 |
C |   VISC         |<-- | VISCOSITE TURBULENTE.                        |
C |   TRAC0        | -->| VALEUR INITIALE DU TRACEUR.                  |
C |   TRAC         | -->| INDIQUE LA PRESENCE D'UN TRACEUR.            |
C |   ITURB        | -->| MODELE DE TURBULENCE.                        |
C |   NPOIN        | -->| NOMBRE DE POINTS DANS LE MAILLAGE            |
C |   GRAV         | -->| ACCELERATION DE LA PESANTEUR.                |
C |   AKEP         | -->| LOGIQUE QUI INDIQUE S'IL FAUT INITIALISER K  |
C |                |    | ET EPSILON.                                  |
C |   TROUVE       | -->| TABLEAU INDIQUANT LES VARIABLES TROUVEES DANS|
C |                |    | LE FICHIER DES RESULTATS DU CALCUL PRECEDENT |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMME APPELE : OV
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
      INTEGER, INTENT(IN)             :: TROUVE(*),ITURB,NPOIN,NTRAC
      LOGICAL, INTENT(INOUT)          :: AKEP
      DOUBLE PRECISION, INTENT(INOUT) :: U(NPOIN),V(NPOIN),H(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: S(NPOIN),ZF(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: TRAC0(NTRAC)
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: T
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C  
      INTEGER ITRAC    
      DOUBLE PRECISION BID
C
C-----------------------------------------------------------------------
C
C  VITESSE U
C
      IF(TROUVE(1).NE.1 )  THEN
          IF(LNG.EQ.1) WRITE(LU,190)
          IF(LNG.EQ.2) WRITE(LU,191)
190       FORMAT(1X,'RESCUE : FICHIER DE RESULTATS DU CALCUL PRECEDENT',
     *         /,1X,'         SANS LA VITESSE U, ON LA PREND NULLE')
191       FORMAT(1X,'RESCUE : PREVIOUS COMPUTATION RESULTS FILE',
     *         /,1X,'         WITHOUT VELOCITY U, WE FIX IT TO ZERO')
          CALL OV( 'X=C     ' , U , U , U , 0.D0 , NPOIN )
      ENDIF
C
C-----------------------------------------------------------------------
C
C  VITESSE V
C
      IF(TROUVE(2).NE.1 )  THEN
          IF(LNG.EQ.1) WRITE(LU,200)
          IF(LNG.EQ.2) WRITE(LU,201)
200       FORMAT(1X,'RESCUE : FICHIER DE RESULTATS DU CALCUL PRECEDENT',
     *         /,1X,'         SANS LA VITESSE V, ON LA PREND NULLE')
201       FORMAT(1X,'RESCUE : PREVIOUS COMPUTATION RESULTS FILE',
     *         /,1X,'         WITHOUT VELOCITY V, WE FIX IT TO ZERO')
          CALL OV( 'X=C     ' , V , V , V , 0.D0 , NPOIN )
      ENDIF
C
C-----------------------------------------------------------------------
C
C  HAUTEUR
C
      IF(TROUVE(4).NE.1) THEN
        IF(TROUVE(5).EQ.1) THEN
          IF(LNG.EQ.1) WRITE(LU,400)
          IF(LNG.EQ.2) WRITE(LU,401)
400       FORMAT(1X,'RESCUE : HAUTEUR D''EAU CALCULEE AVEC LE FOND',
     *         /,1X,'         ET LA SURFACE LIBRE')
401       FORMAT(1X,'RESCUE : WATER DEPTH COMPUTED WITH BATHYMETRY',
     *         /,1X,'         AND SURFACE ELEVATION')
          CALL OV( 'X=Y-Z   ' , H , S , ZF , BID , NPOIN )
        ELSE
          IF(LNG.EQ.1) WRITE(LU,420)
          IF(LNG.EQ.2) WRITE(LU,421)
420       FORMAT(1X,'RESCUE : IMPOSSIBLE DE CALCULER LA HAUTEUR D''EAU')
421       FORMAT(1X,'RESCUE : WATER DEPTH CANNOT BE COMPUTED')
          CALL PLANTE(1)
          STOP
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C
C  TRACEUR
C
      IF(NTRAC.GT.0) THEN
        DO ITRAC=1,NTRAC
          IF(TROUVE(31+ITRAC).EQ.0) THEN
            IF(LNG.EQ.1) WRITE(LU,900)
            IF(LNG.EQ.2) WRITE(LU,901)
900         FORMAT(1X,'RESCUE : CALCUL PRECEDENT SANS TRACEUR,',
     *           /,1X,'         ON PREND TRAC0')
901         FORMAT(1X,'RESCUE : PREVIOUS CALCULATION WITHOUT TRACER',
     *           /,1X,'         WE FIX IT TO TRAC0')
            CALL OS( 'X=C     ' , X=T%ADR(ITRAC)%P,C=TRAC0(ITRAC))
          ENDIF
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------
C
C  K ET EPSILON
C
      IF(ITURB.EQ.3.AND.TROUVE(10).EQ.1.AND.TROUVE(11).EQ.1) THEN
        AKEP=.FALSE.
      ENDIF
      IF(ITURB.EQ.3.AND.(TROUVE(10).EQ.0.OR.TROUVE(11).EQ.0)) THEN
        IF(LNG.EQ.1) WRITE(LU,950)
        IF(LNG.EQ.2) WRITE(LU,951)
950     FORMAT(1X,'RESCUE : K ET EPSILON SERONT REINITIALISES')
951     FORMAT(1X,'RESCUE : K ET EPSILON WILL BE SET AGAIN')
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
