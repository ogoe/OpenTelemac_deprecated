C                       *****************
                        SUBROUTINE BILANT
C                       *****************
C
     *(H,WORK2,WORK3,DT,LT,NIT,INFO,
     * T,AGGLOT,MASSOU,MASTR0,MASTR2,MASTEN,
     * MASTOU,MSK,MASKEL,MESH,
     * FLBOR,NUMLIQ,NFRLIQ,NPTFR,NAMETRAC,FLBORTRA)
C
C***********************************************************************
C  BIEF VERSION 5.9     10/06/08      J-M HERVOUET (LNHE) 01 30 87 80 18
C                                        C MOULIN   (LNH) 01 30 87 83 81
C***********************************************************************
C
C  FONCTION:  EFFECTUE LE BILAN DE MASSE DE TRACEUR
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   H            | -->|  VALEURS DE H A L' ETAPE N+1.
C |   WORK2,3      | -->|  TABLEAUX DE TRAVAIL.                        |                                      |
C |   DT           | -->|  PAS DE TEMPS                                |
C |   LT,NIT       | -->|  NUMERO DU PAS DE TEMPS, NOMBRE TOTAL DE PAS.|
C |   INFO         | -->|  LOGIQUE INDIQUANT SI ON FAIT LES IMPRESSIONS|
C |   T            | -->|  TRACEUR AU TEMPS T(N+1) 
C |   TETAT        | -->|  SEMI-IMPLICITATION DU TRACEUR.              |
C |   MASSOU       | -->|  QUANTITE DE TRACEUR APPORTEE PAR LE TERME   |
C |                |    |  SOURCE
C |   MSK          | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |   MASKEL       | -->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELEMAC2D, CVDFTR2
C
C SOUS-PROGRAMMES APPELES : 
C
C***********************************************************************
C
      USE BIEF, EX_BILANT => BILANT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)            :: LT,NIT,NFRLIQ,NPTFR
      INTEGER, INTENT(IN)            :: NUMLIQ(NFRLIQ)
      DOUBLE PRECISION, INTENT(IN)   :: DT,MASSOU,AGGLOT
      LOGICAL, INTENT(IN)            :: INFO,MSK
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: WORK2,WORK3
      TYPE(BIEF_OBJ), INTENT(IN)     :: H,T,MASKEL,FLBOR
      TYPE(BIEF_OBJ), INTENT(IN)     :: FLBORTRA
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      DOUBLE PRECISION, INTENT(INOUT):: MASTR0,MASTR2,MASTEN,MASTOU
      CHARACTER(LEN=32), INTENT(IN)  :: NAMETRAC
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION P_DSUM
      EXTERNAL         P_DSUM
C
      INTEGER I,IFRLIQ,IELMT,IELMH
C
      DOUBLE PRECISION ERREUT,PERDUE,FLUXT,MASBOR,RELATI,DENOM,MASTR1
C     300 EST ICI MAXFRO, LE NOMBRE MAXIMUM DE FRONTIERES LIQUIDES
      DOUBLE PRECISION FLT_BOUND(300)
C
      INTRINSIC ABS,MAX
C
C-----------------------------------------------------------------------
C
      IELMT = T%ELM
      IELMH = H%ELM
C
C-----------------------------------------------------------------------
C
C  CALCUL COMPATIBLE DE LA QUANTITE DE TRACEUR AU TEMPS N+1:
C  ON TIENT COMPTE DU MASS-LUMPING MAIS IL FAUT QUE AGGLOC=AGGLOT
C
      IF(LT.NE.0) MASTR1 = MASTR2
C
      CALL VECTOR(WORK2,'=','MASVEC          ',IELMT,
     *            1.D0-AGGLOT,T,T,T,T,T,T,MESH,MSK,MASKEL)
C     H EST MIS ICI POUR UNE STRUCTURE BIDON
      CALL VECTOR(WORK3,'=','MASBAS          ',IELMT,
     *                 AGGLOT,H,H,H,H,H,H,MESH,MSK,MASKEL)
C
      CALL OS('X=X+YZ  ',X=WORK2,Y=WORK3,Z=T)
C
      MASTR2 = DOTS(WORK2,H)
      IF(NCSIZE.GT.1) MASTR2=P_DSUM(MASTR2) 
C
      IF(LT.EQ.0) THEN
        MASTR0 = MASTR2
        MASTR1 = MASTR2
        MASTEN = 0.D0
        MASTOU = 0.D0
      ENDIF
C
C=======================================================================
C
C   CALCUL DES FLUX (IL MANQUE LE FLUX DIFFUSIF,...A VOIR)
C
C=======================================================================
C
      FLUXT=0.D0
C
      IF(LT.GT.0.AND.NFRLIQ.GT.0) THEN
        DO IFRLIQ=1,NFRLIQ
          FLT_BOUND(IFRLIQ)=0.D0
        ENDDO
        IF(NPTFR.GT.0) THEN
          DO I=1,NPTFR
!           NOTE: ON POURRAIT DEFINIR FLUX_BOUNDARIES ENTRE 0 ET NFRLIQ
            IFRLIQ=NUMLIQ(I)
            IF(IFRLIQ.GT.0) THEN
!             FLBORTRA MUST NOT BE ASSEMBLED IN //
              FLT_BOUND(IFRLIQ)=FLT_BOUND(IFRLIQ)+FLBORTRA%R(I)
            ENDIF
          ENDDO
        ENDIF
        IF(NCSIZE.GT.1) THEN
          DO IFRLIQ=1,NFRLIQ
            FLT_BOUND(IFRLIQ)=P_DSUM(FLT_BOUND(IFRLIQ))
          ENDDO
        ENDIF
        DO IFRLIQ=1,NFRLIQ
          FLUXT=FLUXT+FLT_BOUND(IFRLIQ)
        ENDDO
      ENDIF
C
C=======================================================================
C
C     CALCUL DES FLUX AUX FRONTIERES LIQUIDES
C
      MASTEN = MASTEN - FLUXT * DT
      MASTOU = MASTOU + MASSOU
C
C=======================================================================
C
C     CALCUL DU FLUX DE TRACEUR A TRAVERS LES PAROIS, PAR LOI DE FLUX
C
C     PROVISOIRE, A PROGRAMMER
      MASBOR = 0.D0
C
C=======================================================================
C
C     CALCUL DE L'ERREUR SUR LA MASSE POUR CE PAS
C
      ERREUT = MASTR1 + MASSOU - MASTR2 - DT*FLUXT
C
C=======================================================================
C
C     IMPRESSIONS :
C
      IF(INFO) THEN
C
C-----------------------------------------------------------------------
C
C     IMPRESSIONS POUR LE TRACEUR
C
        WRITE(LU,*)
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) '                      BILAN DE QUANTITE DE ',
     *    TRIM(NAMETRAC(1:16)),' (UNITE : ',TRIM(NAMETRAC(17:32)),')'
        ENDIF                 
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) '                           BALANCE OF ',
     *    TRIM(NAMETRAC(1:16)),' (UNIT: ',TRIM(NAMETRAC(17:32)),')'
        ENDIF 
C
        IF(LT.EQ.0) THEN
C
          IF(LNG.EQ.1) WRITE(LU,1090) MASTR0
          IF(LNG.EQ.2) WRITE(LU,2090) MASTR0
C
        ELSE
C
          IF(LNG.EQ.1) WRITE(LU,1100) MASTR2
          IF(LNG.EQ.2) WRITE(LU,2100) MASTR2
          IF(NFRLIQ.GT.0) THEN
            DO IFRLIQ=1,NFRLIQ
              IF(LNG.EQ.1) WRITE(LU,1110) IFRLIQ,-FLT_BOUND(IFRLIQ)
              IF(LNG.EQ.2) WRITE(LU,2110) IFRLIQ,-FLT_BOUND(IFRLIQ)
            ENDDO
          ENDIF
          IF(ABS(MASSOU).GT.1.D-8) THEN
            IF(LNG.EQ.1) WRITE(LU,1113) MASSOU
            IF(LNG.EQ.2) WRITE(LU,2113) MASSOU
          ENDIF
          DENOM = MAX(MASTR2,ABS(FLUXT*DT))
          IF(DENOM.GT.1.D-8) THEN
            ERREUT = ERREUT / DENOM
            IF(LNG.EQ.1) WRITE(LU,1120) ERREUT
            IF(LNG.EQ.2) WRITE(LU,2120) ERREUT
          ELSE
            IF(LNG.EQ.1) WRITE(LU,1130) ERREUT
            IF(LNG.EQ.2) WRITE(LU,2130) ERREUT
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
      IF(LT.EQ.NIT) THEN
C
        WRITE(LU,*)
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) '                BILAN FINAL DE QUANTITE DE ',
     *    TRIM(NAMETRAC(1:16)),' (UNITE : ',TRIM(NAMETRAC(17:32)),')'
        ENDIF                 
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) '                     FINAL BALANCE OF ',
     *    TRIM(NAMETRAC(1:16)),' (UNIT: ',TRIM(NAMETRAC(17:32)),')'
        ENDIF 
C
          PERDUE = MASTR0+MASTEN+
     *             MASBOR+MASTOU-MASTR2
          DENOM = MAX(MASTR0,MASTR2,ABS(MASTEN),ABS(MASTOU))
          IF(DENOM.GT.1.D-8) THEN
            RELATI = PERDUE / DENOM
            IF(LNG.EQ.1) WRITE(LU,1140) RELATI
            IF(LNG.EQ.2) WRITE(LU,2140) RELATI
          ELSE
            RELATI = PERDUE
            IF(LNG.EQ.1) WRITE(LU,1150) RELATI
            IF(LNG.EQ.2) WRITE(LU,2150) RELATI
          ENDIF
          IF(LNG.EQ.1) THEN
            WRITE(LU,1160) MASTR0,MASTR2
            IF(ABS(MASTEN).GT.1.D-8) WRITE(LU,1161) MASTEN
            IF(ABS(MASTOU).GT.1.D-8) WRITE(LU,1164) MASTOU
            WRITE(LU,1165) PERDUE
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,2160) MASTR0,MASTR2
            IF(ABS(MASTEN).GT.1.D-8) WRITE(LU,2161) MASTEN
            IF(ABS(MASTOU).GT.1.D-8) WRITE(LU,2164) MASTOU
            WRITE(LU,2165) PERDUE
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
1090  FORMAT(5X,'QUANTITE INITIALE DE TRACEUR :',G16.7)
2090  FORMAT(5X,'INITIAL QUANTITY OF TRACER:',G16.7)
1100  FORMAT(/,5X,'QUANTITE DE TRACEUR :',G16.7)
2100  FORMAT(/,5X,'QUANTITY OF TRACER:',G16.7)
1110  FORMAT(5X,'FRONTIERE ',1I3,' FLUX :           ',G16.7,
     *          ' ( >0 : ENTRANT  <0 : SORTANT )')
1113  FORMAT(5X,'QUANTITE CREEE PAR TERME SOURCE :  ' , G16.7 )
2110  FORMAT(5X,'BOUNDARY ',1I3,' FLUX:         ',G16.7,
     *          ' ( >0 : ENTERING  <0 : EXITING )')
2113  FORMAT(5X,'QUANTITY CREATED BY SOURCE TERM:   ' , G16.7 )
1120  FORMAT(5X,'ERREUR RELATIVE : ',G16.7)
2120  FORMAT(5X,'RELATIVE ERROR: ',G16.7)
1130  FORMAT(5X,'ERREUR ABSOLUE : ',G16.7)
2130  FORMAT(5X,'ABSOLUTE ERROR: ',G16.7)
1140  FORMAT(/,5X,'ERREUR RELATIVE CUMULEE : ',G16.7)
2140  FORMAT(/,5X,'RELATIVE ERROR CUMULATED: ',G16.7)
1150  FORMAT(/,5X,'ERREUR ABSOLUE  CUMULEE: ',G16.7)
2150  FORMAT(/,5X,'ABSOLUTE CUMULATED ERROR: ',G16.7)
1160  FORMAT(/,5X,'QUANTITE INITIALE                 : ',G16.7,
     *       /,5X,'QUANTITE FINALE                   : ',G16.7)
1161  FORMAT(  5X,'QUANTITE ENTREE AUX FRONT. LIQUID.: ',G16.7,
     *            '  ( SI <0 QUANTITE SORTIE )')
1164  FORMAT(  5X,'QUANTITE CREEE PAR TERME SOURCE   : ',G16.7)
1165  FORMAT(  5X,'QUANTITE TOTALE PERDUE            : ',G16.7)
2160  FORMAT(/,5X,'INITIAL QUANTITY                  : ',G16.7,
     *       /,5X,'FINAL QUANTITY                    : ',G16.7)
2161  FORMAT(  5X,'QUANTITY ENTERED THROUGH LIQ. BND.: ',G16.7,
     *            '  ( IF <0 EXIT )')
2164  FORMAT(  5X,'QUANTITY CREATED BY SOURCE TERM   : ',G16.7)
2165  FORMAT(  5X,'TOTAL QUANTITY LOST               : ',G16.7)
C
C=======================================================================
C
      RETURN
      END
