C                       *****************
                        SUBROUTINE FONSTR
C                       *****************
C
     *(H,ZF,Z,CHESTR,NGEO,NFON,NOMFON,MESH,FFON,LISTIN)
C
C***********************************************************************
C  BIEF VERSION 5.6           17/08/94    J-M HERVOUET (LNH) 30 71 80 18
C                                       
C***********************************************************************
C
C  FONCTION  :  RECHERCHE DES FONDS DANS LE FICHIER DE GEOMETRIE.
C               RECHERCHE DES COEFFICIENTS DE FROTTEMENT.
C
C
C  ATTENTION :  LES NOMS DES VARIABLES LUES ONT ETE ICI ECRITS
C               DIRECTEMENT ET NE SONT PLUS PRIS DANS 'TEXTE'
C
C               CECI PERMET DE REPRENDRE UNE GEOMETRIE FAITE DANS UNE
C               AUTRE LANGUE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   H            |<-- | HAUTEUR D'EAU
C |   ZF           |<-- | FOND
C |   Z            |<-- | COTE DE LA SURFACE LIBRE
C |   CHESTR       |<-- | COEFFICIENT DE FROTTEMENT.
C |   NGEO         | -->| NUMERO DU CANAL DU FICHIER DE GEOMETRIE
C |   NFON         | -->| NUMERO DU CANAL DU FICHIER DES FONDS
C |   NOMFON       | -->| NOM DU FICHIER DES FONDS
C |   MESH         | -->| 
C |   FFON         | -->| 
C |   LISTIN       | -->| 
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT :  
C PROGRAMMES APPELES : LIT , OV
C
C***********************************************************************
C
C    ON SUPPOSE QUE LA PARTIE ENTETE DU FICHIER A ETE LUE ET QUE L'ON
C    VA COMMENCER A LIRE LES ENREGISTREMENTS DES RESULTATS
C
C-----------------------------------------------------------------------
C
      USE BIEF, EX_FONSTR => FONSTR
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: H,ZF,Z,CHESTR
      CHARACTER(LEN=72), INTENT(IN) :: NOMFON
      TYPE(BIEF_MESH), INTENT(IN)   :: MESH
      DOUBLE PRECISION, INTENT(IN)  :: FFON
      LOGICAL, INTENT(IN)           :: LISTIN
      INTEGER, INTENT(IN)           :: NGEO,NFON
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER ERR
C
      DOUBLE PRECISION BID
      REAL, ALLOCATABLE :: W(:)
C
      LOGICAL CALFON,CALFRO,OK,LUZF,LUH,LUZ
C
C-----------------------------------------------------------------------
C
      ALLOCATE(W(MESH%NPOIN),STAT=ERR)
      IF(ERR.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'FONSTR : MAUVAISE ALLOCATION DE W'
        IF(LNG.EQ.2) WRITE(LU,*) 'FONSTR: WRONG ALLOCATION OF W'
        STOP
      ENDIF
C
C   INITIALISATION
C
      LUH  =  .FALSE.
      LUZ  =  .FALSE.
      LUZF =  .FALSE.
      CALFRO = .TRUE.
C
C-----------------------------------------------------------------------
C
C     RECHERCHE DU COEFFICIENT DE FROTTEMENT
C
      IF(LNG.EQ.1) CALL FIND_IN_SEL(CHESTR,'FROTTEMENT      ',NGEO,W,OK,
     *                              TIME=BID)
      IF(LNG.EQ.2) CALL FIND_IN_SEL(CHESTR,'BOTTOM FRICTION ',NGEO,W,OK,
     *                              TIME=BID)
C     CAS D'UNE GEOMETRIE AVEC UNE LANGUE ETRANGERE
      IF(.NOT.OK.AND.LNG.EQ.1) THEN  
        CALL FIND_IN_SEL(CHESTR,'BOTTOM FRICTION ',NGEO,W,OK,TIME=BID)
      ENDIF
      IF(.NOT.OK.AND.LNG.EQ.2) THEN
        CALL FIND_IN_SEL(CHESTR,'FROTTEMENT      ',NGEO,W,OK,TIME=BID)
      ENDIF
      IF(OK) THEN
        CALFRO = .FALSE.
        IF(LNG.EQ.1) WRITE(LU,5)
        IF(LNG.EQ.2) WRITE(LU,6)
5       FORMAT(1X,'FONSTR : COEFFICIENTS DE FROTTEMENT LUS DANS',/,
     *         1X,'         LE FICHIER DE GEOMETRIE')
6       FORMAT(1X,'FONSTR : FRICTION COEFFICIENTS READ IN THE',/,
     *         1X,'         GEOMETRY FILE')
      ENDIF
C
C     RECHERCHE DU FOND
C
      IF(LNG.EQ.1) CALL FIND_IN_SEL(ZF,'FOND            ',NGEO,W,OK,
     *                              TIME=BID)
      IF(LNG.EQ.2) CALL FIND_IN_SEL(ZF,'BOTTOM          ',NGEO,W,OK,
     *                              TIME=BID)
      IF(.NOT.OK.AND.LNG.EQ.1) THEN  
        CALL FIND_IN_SEL(ZF,'BOTTOM          ',NGEO,W,OK,TIME=BID)
      ENDIF
      IF(.NOT.OK.AND.LNG.EQ.2) THEN
        CALL FIND_IN_SEL(ZF,'FOND            ',NGEO,W,OK,TIME=BID)
      ENDIF
C     MESHES FROM BALMAT ?
      IF(.NOT.OK) CALL FIND_IN_SEL(ZF,'altimetrie      ',NGEO,W,OK,
     *                             TIME=BID)
C     TOMAWAC FRANCAIS ?
      IF(.NOT.OK) CALL FIND_IN_SEL(ZF,'COTE_DU_FOND    ',NGEO,W,OK,
     *                             TIME=BID)
C     TOMAWAC ANGLAIS ?
      IF(.NOT.OK) CALL FIND_IN_SEL(ZF,'BOTTOM_LEVEL    ',NGEO,W,OK,
     *                             TIME=BID)
      LUZF = OK 
C
      IF(.NOT.LUZF) THEN
C       RECHERCHE DE LA HAUTEUR ET DE LA SURFACE LIBRE
        IF(LNG.EQ.1) CALL FIND_IN_SEL(H,'HAUTEUR D''EAU   ',NGEO,W,OK,
     *                                TIME=BID)
        IF(LNG.EQ.2) CALL FIND_IN_SEL(H,'WATER DEPTH     ',NGEO,W,OK,
     *                                TIME=BID)
        IF(.NOT.OK.AND.LNG.EQ.1) THEN  
          CALL FIND_IN_SEL(H,'WATER DEPTH     ',NGEO,W,OK,TIME=BID)
        ENDIF
        IF(.NOT.OK.AND.LNG.EQ.2) THEN
          CALL FIND_IN_SEL(H,'HAUTEUR D''EAU   ',NGEO,W,OK,TIME=BID)
        ENDIF
        LUH = OK
        IF(LNG.EQ.1) CALL FIND_IN_SEL(Z,'SURFACE LIBRE   ',NGEO,W,OK,
     *                                TIME=BID)
        IF(LNG.EQ.2) CALL FIND_IN_SEL(Z,'FREE SURFACE    ',NGEO,W,OK,
     *                                TIME=BID)
        IF(.NOT.OK.AND.LNG.EQ.1) THEN  
          CALL FIND_IN_SEL(Z,'FREE SURFACE    ',NGEO,W,OK,TIME=BID)
        ENDIF
        IF(.NOT.OK.AND.LNG.EQ.2) THEN
          CALL FIND_IN_SEL(Z,'SURFACE LIBRE   ',NGEO,W,OK,TIME=BID)
        ENDIF
        LUZ = OK
      ENDIF
C
C     INITIALISATION DU FOND
C
      IF(LUZF) THEN
C
         CALFON = .FALSE.
C
      ELSE
C
         IF (LUZ.AND.LUH) THEN
C
            CALL OS( 'X=Y-Z   ' , ZF , Z , H , BID )
            IF(LNG.EQ.1) WRITE(LU,24)
            IF(LNG.EQ.2) WRITE(LU,25)
24          FORMAT(1X,'FONSTR (BIEF) : ATTENTION, FOND CALCULE AVEC',/,
     *                '                PROFONDEUR ET SURFACE LIBRE',/,
     *                '                DU FICHIER DE GEOMETRIE')
25          FORMAT(1X,'FONSTR (BIEF): ATTENTION, THE BOTTOM RESULTS',/,
     *                '               FROM DEPTH AND SURFACE ELEVATION',
     *              /,'               FOUND IN THE GEOMETRY FILE')
            CALFON = .FALSE.
C
         ELSE
C
            CALFON = .TRUE.
C
         ENDIF
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C CONSTRUCTION DU FOND SI ON NE L'A PAS TROUVE DANS LA GEOMETRIE
C
      IF(NOMFON(1:1).NE.' ') THEN
C       UN FICHIER DES FONDS A ETE DONNE, ON (RE)CALCULE DONC LES FONDS
        IF(LISTIN) THEN
          IF(LNG.EQ.1) WRITE(LU,2223) NOMFON
          IF(LNG.EQ.2) WRITE(LU,2224) NOMFON
          IF(.NOT.CALFON) THEN
            IF(LNG.EQ.1) WRITE(LU,2225)
            IF(LNG.EQ.2) WRITE(LU,2226)
          ENDIF
        ENDIF
2223    FORMAT(/,1X,'FONSTR (BIEF) : FOND DANS LE FICHIER : ',A72)
2224    FORMAT(/,1X,'FONSTR (BIEF): BATHYMETRY GIVEN IN FILE : ',A72)
2225    FORMAT(  1X,'                LE FOND TROUVE DANS LE FICHIER',/,
     *           1X,'                DE GEOMETRIE EST IGNORE',/)
2226    FORMAT(  1X,'               BATHYMETRY FOUND IN THE',/,
     *           1X,'               GEOMETRY FILE IS IGNORED',/)
C
        CALL FOND(ZF%R,MESH%X%R,MESH%Y%R,MESH%NPOIN,NFON,
     *            MESH%NBOR%I,MESH%KP1BOR%I,MESH%NPTFR)
C
      ELSEIF(CALFON) THEN
        IF(LISTIN) THEN
          IF(LNG.EQ.1) WRITE(LU,2227)
          IF(LNG.EQ.2) WRITE(LU,2228)
        ENDIF
2227    FORMAT(/,1X,'FONSTR (BIEF) : PAS DE FOND DANS LE FICHIER DE',
     *         /,1X,'                GEOMETRIE ET PAS DE FICHIER DES',
     *         /,1X,'                FONDS. LE FOND EST INITIALISE A'
     *         /,1X,'                ZERO MAIS PEUT ENCORE ETRE MODIFIE'
     *         /,1X,'                DANS CORFON.',
     *         /,1X)
2228    FORMAT(/,1X,'FONSTR (BIEF): NO BATHYMETRY IN THE GEOMETRY FILE',
     *         /,1X,'               AND NO BATHYMETRY FILE. THE BOTTOM',
     *         /,1X,'               LEVEL IS FIXED TO ZERO BUT STILL',
     *         /,1X,'               CAN BE MODIFIED IN CORFON.',
     *         /,1X)
        CALL OS( 'X=C     ' , ZF , ZF , ZF , 0.D0 )
      ENDIF
C
C-----------------------------------------------------------------------
C
C CALCUL DU COEFFICIENT DE FROTTEMENT SUR LE FOND
C
      IF(CALFRO) THEN
        CALL OS( 'X=C     ' , CHESTR , CHESTR , CHESTR , FFON )
      ENDIF
      CALL STRCHE
C
C-----------------------------------------------------------------------
C
      DEALLOCATE(W)
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
