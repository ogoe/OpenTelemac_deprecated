C                       ****************
                        SUBROUTINE BILAN
C                       ****************
C
     *(MESH,H,WORK,MASK,AT,DT,LT,NIT,INFO,MASSES,MSK,MASKEL,EQUA,POROS,
     * OPTBAN,NPTFR,FLBOR,FLUX_BOUNDARIES,NUMLIQ,NFRLIQ)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.9  27/03/08   J-M HERVOUET (LNHE) 01 30 87 80 18
C
C                          14/01/05 : COMPATIBLE COMPUTATION OF FLUXES
C                                     AT EXITS
C
C                          27/03/08 : PRINTING FLUXES PER BOUNDARY
C                                     INSTEAD OF FREE AND IMPOSED FLUX
C***********************************************************************
C
C  FONCTION:  EFFECTUE LE BILAN DE MASSE D'EAU
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   WORK         | -->|  TABLEAU DE TRAVAIL.                         |
C |   AT           | -->|  TEMPS                                       |
C |   DT           | -->|  PAS DE TEMPS                                |
C |   LT,NIT       | -->|  NUMERO DU PAS DE TEMPS, NOMBRE TOTAL DE PAS.|
C |   INFO         | -->|  LOGIQUE INDIQUANT SI ON FAIT LES IMPRESSIONS|
C |   MASSES       | -->|  MASSE APPORTEE PAR TERME SOURCE.
C |   MSK          | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |   MASKEL       | -->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMMES APPELES : VECTOR , BIEF_SUM
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
C     SIZE OF NUMLIQ AND FLUX_BOUNDARIES IS NFRLIQ BUT NFRLIQ
C     CAN BE 0.
C
      INTEGER, INTENT(IN)            :: LT,NIT,OPTBAN,NPTFR,NFRLIQ
      INTEGER, INTENT(IN)            :: NUMLIQ(*)
      CHARACTER(LEN=20), INTENT(IN)  :: EQUA
      LOGICAL, INTENT(IN)            :: INFO,MSK
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: WORK,FLBOR
      TYPE(BIEF_OBJ), INTENT(IN)     :: H,MASKEL,POROS,MASK
      DOUBLE PRECISION, INTENT(IN)   :: AT,DT
      DOUBLE PRECISION, INTENT(INOUT):: MASSES,FLUX_BOUNDARIES(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I      
C
      DOUBLE PRECISION P_DSUM
      EXTERNAL         P_DSUM
C
      DOUBLE PRECISION ERREUR,FLUX1,PERDUE,DENOM
      DOUBLE PRECISION MASSE0,MASSE1,MASSE2,MASENT,RELATI,MASSET      
C
      INTRINSIC ABS      
C
      SAVE MASSE0,MASSE1,MASSE2,MASENT,MASSET
C
C-----------------------------------------------------------------------
C
C  CALCUL COMPATIBLE DE LA MASSE D'EAU
C
      IF(LT.NE.0) MASSE1 = MASSE2
C
      IF(EQUA(1:15).EQ.'SAINT-VENANT VF') THEN
        CALL VECTOR(WORK,'=','MASBAS          ',H%ELM,
     *              1.D0,H,H,H,H,H,H,MESH,MSK,MASKEL)
        CALL OS( 'X=XY    ' , X=WORK , Y=H )
      ELSEIF(OPTBAN.EQ.3) THEN
        CALL VECTOR(WORK,'=','MASVEC          ',H%ELM,
     *              1.D0,H,H,H,H,H,H,MESH,.TRUE.,POROS)
      ELSE
        CALL VECTOR(WORK,'=','MASVEC          ',H%ELM,
     *              1.D0,H,H,H,H,H,H,MESH,MSK,MASKEL)
      ENDIF
      MASSE2 = BIEF_SUM(WORK)
C
      IF(NCSIZE.GT.1) MASSE2 = P_DSUM(MASSE2)
C
      IF(LT.EQ.0) THEN
        MASSE0 = MASSE2
        MASSE1 = MASSE2
        MASENT = 0.D0
        MASSET = 0.D0
C
C       FOR THE FIRST CALL, RETURN HERE
C
        CALL OS('X=0     ',X=FLBOR)
        IF(NFRLIQ.GT.0) THEN
          DO I=1,NFRLIQ
            FLUX_BOUNDARIES(I)=0.D0
          ENDDO
        ENDIF
        RETURN
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C   MASSE AJOUTEE PAR TERME SOURCE
C
      IF(NCSIZE.GT.1) MASSES = P_DSUM(MASSES)
      MASSET = MASSET + MASSES
C
C=======================================================================
C
C   CALCUL DES FLUX AUX FRONTIERES LIQUIDES
C
      IF(NFRLIQ.GT.0) THEN
        DO I=1,NFRLIQ
          FLUX_BOUNDARIES(I)=0.D0
        ENDDO
        IF(NPTFR.GT.0) THEN
          DO I=1,NPTFR
!           NOTE: ON POURRAIT DEFINIR FLUX_BOUNDARIES ENTRE 0 ET NFRLIQ
            IF(NUMLIQ(I).GT.0) THEN
              FLUX_BOUNDARIES(NUMLIQ(I))=
     *        FLUX_BOUNDARIES(NUMLIQ(I))+FLBOR%R(I)
            ENDIF
          ENDDO
        ENDIF
        IF(NCSIZE.GT.1) THEN
          DO I=1,NFRLIQ
            FLUX_BOUNDARIES(I)=P_DSUM(FLUX_BOUNDARIES(I))
          ENDDO
        ENDIF
      ENDIF
C
C=======================================================================
C
C   FLUX TOTAL AUX FRONTIERES LIQUIDES
C
      FLUX1=0.D0
      IF(NFRLIQ.GT.0) THEN
        DO I=1,NFRLIQ
          FLUX1=FLUX1+FLUX_BOUNDARIES(I)
        ENDDO        
      ENDIF
C
C=======================================================================
C
      MASENT = MASENT - FLUX1*DT
C
C=======================================================================
C
C   CALCUL DE L'ERREUR SUR LA MASSE POUR CE PAS
C
      ERREUR = MASSE1 + MASSES - MASSE2 - DT*FLUX1
C
C=======================================================================
C
C   IMPRESSIONS :
C
      IF(INFO) THEN
C
C-----------------------------------------------------------------------
C
C     IMPRESSIONS POUR LA MASSE D'EAU
C
        IF(LT.EQ.0) THEN
C
          CALL ENTETE(7,AT,LT)
          IF(LNG.EQ.1) WRITE(LU,1000) MASSE0
          IF(LNG.EQ.2) WRITE(LU,2000) MASSE0
C
        ELSE
C
          CALL ENTETE(7,AT,LT) 
          IF(LNG.EQ.1) THEN
            WRITE(LU,1010) MASSE2
            IF(NFRLIQ.GT.0) THEN
              DO I=1,NFRLIQ
                WRITE(LU,3020) I,-FLUX_BOUNDARIES(I)
              ENDDO
            ENDIF
          ENDIF 
          IF(LNG.EQ.2) THEN
            WRITE(LU,2010) MASSE2
            IF(NFRLIQ.GT.0) THEN
              DO I=1,NFRLIQ
                WRITE(LU,4020) I,-FLUX_BOUNDARIES(I)
              ENDDO
            ENDIF
          ENDIF         
          IF(ABS(MASSES).GT.1.D-6) THEN
            IF(LNG.EQ.1) WRITE(LU,1031) MASSES
            IF(LNG.EQ.2) WRITE(LU,2031) MASSES
          ENDIF
C         CALCUL DE L'ERREUR RELATIVE OU ABSOLUE
          DENOM = MAX(MASSE2,ABS(FLUX1*DT))
          IF(DENOM.GT.1.D-8) THEN
            ERREUR = ERREUR / DENOM
            IF(LNG.EQ.1) WRITE(LU,1040) AT,ERREUR
            IF(LNG.EQ.2) WRITE(LU,2040) AT,ERREUR
          ELSE
            IF(LNG.EQ.1) WRITE(LU,1050) AT,ERREUR
            IF(LNG.EQ.2) WRITE(LU,2050) AT,ERREUR
          ENDIF
C
        ENDIF
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  BILAN DE MASSE FINAL
C
      IF(LT.EQ.NIT.AND.INFO) THEN
C
        CALL ENTETE(8,AT,LT)
C       PERDUE = MASSE0+MASSET+MASENT+MASAJT-MASSE2
        PERDUE = MASSE0+MASSET+MASENT-MASSE2
        DENOM = MAX( MASSE0 , MASSE2 , ABS(MASENT) )
        IF(DENOM.GT.1.D-8) THEN
          RELATI = PERDUE / DENOM
          IF(LNG.EQ.1) WRITE(LU,1060) RELATI
          IF(LNG.EQ.2) WRITE(LU,2060) RELATI
        ELSE
          RELATI = PERDUE
          IF(LNG.EQ.1) WRITE(LU,1070) RELATI
          IF(LNG.EQ.2) WRITE(LU,2070) RELATI
        ENDIF
        IF(LNG.EQ.1) THEN
          WRITE(LU,1080) MASSE0,MASSE2
          IF(ABS(MASENT).GT.1.D-8) WRITE(LU,1081) MASENT
          IF(ABS(MASSET).GT.1.D-8) WRITE(LU,1082) MASSET
C         IF(ABS(MASAJT).GT.1.D-8) WRITE(LU,1083) MASAJT
          WRITE(LU,1084) PERDUE
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,2080) MASSE0,MASSE2
          IF(ABS(MASENT).GT.1.D-8) WRITE(LU,2081) MASENT
          IF(ABS(MASSET).GT.1.D-8) WRITE(LU,2082) MASSET
C         IF(ABS(MASAJT).GT.1.D-8) WRITE(LU,2083) MASAJT
          WRITE(LU,2084) PERDUE
        ENDIF
C
      ENDIF
C
C  FIN DES IMPRESSIONS :
C
C=======================================================================
C
C  FORMATS D'IMPRESSION :
C
1000  FORMAT(5X,'VOLUME D''EAU INITIAL DANS LE DOMAINE: ',G16.7,' M3')
2000  FORMAT(5X,'INITIAL WATER VOLUME IN THE DOMAIN: ',G16.7,' M3')
C
1010  FORMAT(5X,'VOLUME DANS LE DOMAINE :',G16.7,' M3')
2010  FORMAT(5X,'VOLUME IN THE DOMAIN :',G16.7,' M3')
C
1031  FORMAT(5X,'VOLUME AJOUTE PAR TERME SOURCE : ',G16.7,' M3')
2031  FORMAT(5X,'ADDITIONAL VOLUME DUE TO SOURCE TERMS: ',G16.7,' M3')
C
1040  FORMAT(5X,'ERREUR RELATIVE EN VOLUME A T = ',G16.4,' S : ',G16.7)
2040  FORMAT(5X,'RELATIVE ERROR IN VOLUME AT T = ',G16.4,' S : ',G16.7)
C
1050  FORMAT(5X,'ERREUR ABSOLUE EN VOLUME A T = ',G16.4,' S: ',G16.7)
2050  FORMAT(5X,'ABSOLUTE ERROR IN VOLUME AT T = ',G16.4,'S: ',G16.7)
C
1060  FORMAT(/,5X,'ERREUR RELATIVE CUMULEE SUR LE VOLUME : ',G16.7)
2060  FORMAT(/,5X,'RELATIVE ERROR CUMULATED ON VOLUME: ',G16.7)
C
1070  FORMAT(/,5X,'ERREUR ABSOLUE CUMULEE SUR LE VOLUME : ',G16.7)
2070  FORMAT(/,5X,'ABSOLUTE ERROR CUMULATED ON VOLUME: ',G16.7)
C
1080  FORMAT(/,5X,'VOLUME INITIAL              : ',G16.7,' M3',
     *       /,5X,'VOLUME FINAL                : ',G16.7,' M3')
1081  FORMAT(  5X,'VOLUME ENTRE AUX FRONTIERES : ',G16.7,' M3',
     *            '  ( SI <0 VOLUME SORTI )')
1082  FORMAT(  5X,'VOLUME AJOUTE ( SOURCES   ) : ',G16.7,' M3')
C1083  FORMAT(  5X,'VOLUME AJOUTE ( CDT. LIM. ) : ',G16.7,' M3')
1084  FORMAT(  5X,'VOLUME TOTAL PERDU          : ',G16.7,' M3')
2080  FORMAT(/,5X,'INITIAL VOLUME              : ',G16.7,' M3',
     *       /,5X,'FINAL VOLUME                : ',G16.7,' M3')
2081  FORMAT(  5X,'VOLUME THAT ENTERED THE DOMAIN: ',G16.7,' M3',
     *            '  ( IF <0 EXIT )')
2082  FORMAT(  5X,'VOLUME ADDED BY SOURCE TERM   : ',G16.7,' M3')
2084  FORMAT(  5X,'TOTAL VOLUME LOST             : ',G16.7,' M3')
3020  FORMAT(5X,'FLUX FRONTIERE ',I4,' : ', G16.7 ,' M3/S',
     *          '  ( >0 : ENTRANT  <0 : SORTANT )')
4020  FORMAT(5X,'FLUX BOUNDARY ',I4,': ', G16.7 ,' M3/S',
     *          '  ( >0 : ENTERING  <0 : EXITING )')
C
C=======================================================================
C
      RETURN
      END
