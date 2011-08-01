C                       ******************
                        SUBROUTINE BILANT1
C                       ******************
C
     *(H,UCONV,VCONV,HPROP,WORK1,WORK2,WORK3,WORK4,WORK5,DT,LT,NIT,INFO,
     * MASKTR,T,TN,TETAT,MASSOU,MSK,MASKEL,MESH,FLUSOR,FLUENT,EQUA,LTT,
     * ITRAC)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.8    05/11/07  J-M HERVOUET (LNH) 01 30 87 80 18
C                                          
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
C |   UCONV,VCONV  | -->|  CHAMP CONVECTEUR                            |
C |   HPROP        | -->|  HAUTEUR DE PROPAGATION
C |   WORK1,2,6    | -->|  TABLEAUX DE TRAVAIL.                        |
C |   AT           | -->|  TEMPS                                       |
C |   DT           | -->|  PAS DE TEMPS                                |
C |   LT,NIT       | -->|  NUMERO DU PAS DE TEMPS, NOMBRE TOTAL DE PAS.|
C |   LTT          | -->|  NUMERO DU PAS DE TEMPS TRACEUR              |
C |   INFO         | -->|  LOGIQUE INDIQUANT SI ON FAIT LES IMPRESSIONS|
C |   T,TN         | -->|  TRACEUR AU TEMPS T(N+1) ET T(N)             |
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
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMMES APPELES : CHVARC , CLIP , DIVF , FLUBOR , LISSAG ,
C                           MASVEC , OV
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
      INTEGER, INTENT(IN)            :: LT,NIT,LTT,ITRAC
      DOUBLE PRECISION, INTENT(IN)   :: DT,TETAT,MASSOU,FLUSOR,FLUENT
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: WORK1,WORK2,WORK3,WORK4,WORK5
      TYPE(BIEF_OBJ), INTENT(IN)     :: HPROP,UCONV,VCONV,H,T,TN,MASKEL
      TYPE(BIEF_OBJ), INTENT(IN)     :: MASKTR
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      LOGICAL, INTENT(IN)            :: MSK,INFO
      CHARACTER(LEN=20), INTENT(IN)  :: EQUA
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION P_DSUM
      EXTERNAL         P_DSUM
C
      INTEGER DIR,DDL,OND,IELMT,IELMH      
C
      DOUBLE PRECISION ERREUT,PERDUE
      DOUBLE PRECISION FLUXT,MASBOR
      DOUBLE PRECISION FLTDIR,FLTDDL,FLTOND
      DOUBLE PRECISION C,RELATI,DENOM
C
      DOUBLE PRECISION MASTR0(100),MASTR1(100),MASTR2(100),MASTEN(100)
      DOUBLE PRECISION MASTOU(100),DIRTOT(100)
C
      INTRINSIC ABS,MAX
C
C-----------------------------------------------------------------------
C
      SAVE MASTR0,MASTR1,MASTR2,MASTEN,MASTOU,DIRTOT
C
C-----------------------------------------------------------------------
C
      IELMT = T%ELM
      IELMH = H%ELM
C
C-----------------------------------------------------------------------
C
C PROVISOIRE : H ET HPROP SONT REMPLACES PAR WORK4 ET WORK5 PARTOUT
C
      CALL OS ('X=Y     ' , WORK4 , H     , H , C )
      CALL OS ('X=Y     ' , WORK5 , HPROP , H , C )
C
      IF(IELMT.NE.IELMH) THEN
        CALL CHGDIS(WORK4,IELMH,IELMT,MESH)
        CALL CHGDIS(WORK5,IELMH,IELMT,MESH)
      ENDIF
C
C FIN DU PROVISOIRE, SAUF REMPLACEMENT DE H ET HPROP PLUS LOIN
C
C-----------------------------------------------------------------------
C
C  CALCUL COMPATIBLE DE LA QUANTITE DE TRACEUR :
C
      IF(LT.NE.0) MASTR1(ITRAC) = MASTR2(ITRAC)
C     H EST MIS ICI POUR UNE STRUCTURE BIDON
C
      IF(EQUA(1:15).EQ.'SAINT-VENANT VF') THEN
        CALL VECTOR(WORK2,'=','MASBAS          ',IELMT,
     *              1.D0,H,H,H,H,H,H,MESH,MSK,MASKEL)
        CALL OS( 'X=XY    ' , WORK2 , T , H , C )
      ELSE
        CALL VECTOR(WORK2,'=','MASVEC          ',IELMT,
     *              1.D0,T,H,H,H,H,H,MESH,MSK,MASKEL)
      ENDIF
      MASTR2(ITRAC) = DOTS(WORK2,WORK4)
      IF(NCSIZE.GT.1) MASTR2(ITRAC)=P_DSUM(MASTR2(ITRAC))
C
      IF(LT.EQ.0) THEN
        MASTR0(ITRAC) = MASTR2(ITRAC)
        MASTR1(ITRAC) = MASTR2(ITRAC)
        MASTEN(ITRAC) = 0.D0
        MASTOU(ITRAC) = 0.D0
        DIRTOT(ITRAC) = 0.D0
      ENDIF
C
C=======================================================================
C
C   CALCUL DES FLUX (IL MANQUE LE FLUX DIFFUSIF,...A VOIR)
C
C=======================================================================
C
      DIR=1
      DDL=2
      OND=4
C
C=======================================================================
C   CALCUL DES FLUX IMPOSES (DEBIT IMPOSE OU VITESSES IMPOSEES)
C
      IF(EQUA(1:15).EQ.'SAINT-VENANT VF') THEN
        FLTDIR = FLUENT
      ELSE
        CALL VECTOR(WORK2,'=','FLUBDF          ',IELBOR(IELMT,1),
     *              1.D0,WORK5,H,H,UCONV,VCONV,VCONV,
     *              MESH,.TRUE.,MASKTR%ADR(DIR)%P)
C
        CALL CPSTVC(WORK2,WORK3)
        CALL OSBD( 'X=CY    ' , WORK3 , T  ,  T , TETAT      , MESH )
        CALL OSBD( 'X=X+CY  ' , WORK3 , TN ,  T , 1.D0-TETAT , MESH )
        FLTDIR=DOTS(WORK2,WORK3)
        IF(NCSIZE.GT.1) FLTDIR=P_DSUM(FLTDIR)
      ENDIF
C
C=======================================================================
C
C   CALCUL DES FLUX LIBRES
C
      IF(EQUA(1:15).EQ.'SAINT-VENANT VF') THEN
        FLTDDL = FLUSOR
      ELSE
        CALL VECTOR(WORK2,'=','FLUBDF          ',IELBOR(IELMT,1),
     *              1.D0,WORK5,H,H,UCONV,VCONV,VCONV,
     *              MESH,.TRUE.,MASKTR%ADR(DDL)%P)
        CALL CPSTVC(WORK2,WORK3)
        CALL OSBD( 'X=CY    ' , WORK3 , T  ,  T , TETAT      , MESH )
        CALL OSBD( 'X=X+CY  ' , WORK3 , TN ,  T , 1.D0-TETAT , MESH )
        FLTDDL=DOTS(WORK2,WORK3)
        IF(NCSIZE.GT.1) FLTDDL=P_DSUM(FLTDDL)
      ENDIF
C
C=======================================================================
C
C   CALCUL DES FLUX PAR ONDE INCIDENTE
C
      CALL VECTOR(WORK2,'=','FLUBDF          ',IELBOR(IELMT,1),
     *            1.D0,WORK5,H,H,UCONV,VCONV,VCONV,
     *            MESH,.TRUE.,MASKTR%ADR(OND)%P)
      CALL CPSTVC(WORK2,WORK3)
      CALL OSBD( 'X=CY    ' , WORK3 , T  ,  T , TETAT      , MESH )
      CALL OSBD( 'X=X+CY  ' , WORK3 , TN ,  T , 1.D0-TETAT , MESH )
      FLTOND=DOTS(WORK2,WORK3)
      IF(NCSIZE.GT.1) FLTOND=P_DSUM(FLTOND)
C
C=======================================================================
C
C   CALCUL DES FLUX AUX FRONTIERES LIQUIDES
C
      FLUXT = FLTDIR + FLTDDL + FLTOND
      MASTEN(ITRAC) = MASTEN(ITRAC) - FLUXT
      MASTOU(ITRAC) = MASTOU(ITRAC) + MASSOU
      DIRTOT(ITRAC) = DIRTOT(ITRAC) - FLTDIR
C
C=======================================================================
C
C   CALCUL DU FLUX DE TRACEUR A TRAVERS LES PAROIS, PAR LOI DE FLUX
C
C     PROVISOIRE, A PROGRAMMER
      MASBOR = 0.D0
C
C=======================================================================
C
C   CALCUL DE L'ERREUR SUR LA MASSE POUR CE PAS
C   est nulle pour tous les pas de temps hydro ou il n'y a pas 
C   de mise a jour traceur
C
      ERREUT = MASTR1(ITRAC) + MASSOU - MASTR2(ITRAC) - FLUXT
C
C=======================================================================
C
C  IMPRESSIONS :
C
      IF(INFO) THEN
C
C-----------------------------------------------------------------------
C
C     IMPRESSIONS POUR LE TRACEUR
C
        IF(LNG.EQ.1) WRITE(LU,500) ITRAC
        IF(LNG.EQ.2) WRITE(LU,501) ITRAC
C
        IF(LT.EQ.0) THEN
C
          IF(LNG.EQ.1) WRITE(LU,1090) MASTR0(ITRAC)
          IF(LNG.EQ.2) WRITE(LU,2090) MASTR0(ITRAC)
C
        ELSE
C
          IF(LNG.EQ.1) THEN
            WRITE(LU,1100) MASTR2(ITRAC)
            IF(ABS(FLTDIR).GT.1.D-8) WRITE(LU,1110) -FLTDIR
            IF(ABS(FLTDDL).GT.1.D-8) WRITE(LU,1111) -FLTDDL
            IF(ABS(FLTOND).GT.1.D-8) WRITE(LU,1112) -FLTOND
            IF(ABS(MASSOU).GT.1.D-8) WRITE(LU,1113) MASSOU
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,2100) MASTR2(ITRAC)
            IF(ABS(FLTDIR).GT.1.D-8) WRITE(LU,2110) -FLTDIR
            IF(ABS(FLTDDL).GT.1.D-8) WRITE(LU,2111) -FLTDDL
            IF(ABS(FLTOND).GT.1.D-8) WRITE(LU,2112) -FLTOND
            IF(ABS(MASSOU).GT.1.D-8) WRITE(LU,2113) MASSOU
          ENDIF
C
          PERDUE = MASTR0(ITRAC)+MASTEN(ITRAC)+
     *             MASBOR+MASTOU(ITRAC)-MASTR2(ITRAC)
          DENOM = MAX(MASTR0(ITRAC),MASTR2(ITRAC),ABS(DIRTOT(ITRAC)))
          IF(DENOM.GT.1.D-8) THEN
            RELATI = PERDUE / DENOM
            IF(LNG.EQ.1) WRITE(LU,1140) RELATI
            IF(LNG.EQ.2) WRITE(LU,2140) RELATI
          ELSE
            RELATI = PERDUE
            IF(LNG.EQ.1) WRITE(LU,1150) RELATI
            IF(LNG.EQ.2) WRITE(LU,2150) RELATI
          ENDIF
c
          IF(LNG.EQ.1) THEN
            IF(ABS(MASTEN(ITRAC)).GT.1.D-8) WRITE(LU,1161) MASTEN(ITRAC)
            IF(ABS(MASTOU(ITRAC)).GT.1.D-8) WRITE(LU,1164) MASTOU(ITRAC)
          ENDIF
          IF(LNG.EQ.2) THEN
            IF(ABS(MASTEN(ITRAC)).GT.1.D-8) WRITE(LU,2161) MASTEN(ITRAC)
            IF(ABS(MASTOU(ITRAC)).GT.1.D-8) WRITE(LU,2164) MASTOU(ITRAC)
          ENDIF
C
        ENDIF
C
      ENDIF
C
C-----------------------------------------------------------------------
C  BILAN DE MASSE FINAL
C
       IF(LT.EQ.NIT) THEN
C
          IF(LNG.EQ.1) WRITE(LU,600) ITRAC
          IF(LNG.EQ.2) WRITE(LU,601) ITRAC
C
          PERDUE = MASTR0(ITRAC)+MASTEN(ITRAC)+
     *             MASBOR+MASTOU(ITRAC)-MASTR2(ITRAC)
C
          IF(LNG.EQ.1) THEN
            WRITE(LU,1160) MASTR0(ITRAC),MASTR2(ITRAC)
            WRITE(LU,1165) PERDUE
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,2160) MASTR0(ITRAC),MASTR2(ITRAC)
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
500   FORMAT(80(' '),/,22X,'BILAN DE QUANTITE DU TRACEUR ',1I2)
501   FORMAT(80(' '),/,21X,'BALANCE OF TRACER ',1I2)
600   FORMAT(80('-'),/,20X,'BILAN FINAL DE QUANTITE DU TRACEUR ',1I2)
601   FORMAT(80('-'),/,19X,'FINAL BALANCE OF TRACER ',1I2)
1090  FORMAT(5X,'QUANTITE INITIALE DE TRACEUR :',G16.7,' UNITE M3')
2090  FORMAT(5X,'INITIAL QUANTITY OF TRACER:',G16.7,' TRACER UNIT M3')
1100  FORMAT(/,5X,'QUANTITE DE TRACEUR :',G16.7,' UNITE M3')
2100  FORMAT(/,5X,'QUANTITY OF TRACER:',G16.7,' TRACER UNIT M3')
1110  FORMAT(5X,'FLUX IMPOSE DE TRACEUR :           ' , G16.7 ,
     *          '  ( >0 : ENTRANT  <0 : SORTANT )')
1111  FORMAT(5X,'FLUX LIBRE DE TRACEUR :            ' , G16.7 ,
     *          '  ( >0 : ENTRANT  <0 : SORTANT )')
1112  FORMAT(5X,'FLUX INCIDENT DE TRACEUR :         ' , G16.7 ,
     *          '  ( >0 : ENTRANT  <0 : SORTANT )')
1113  FORMAT(5X,'QUANTITE CREEE PAR TERME SOURCE :  ' , G16.7 )
C1114  FORMAT(5X,'QUANTITE AJOUTEE ( CDT. LIM. )  :  ' , G16.7 )
2110  FORMAT(5X,'PRESCRIBED FLUX OF TRACER:         ' , G16.7 ,
     *          '  ( >0 : ENTERING  <0 : EXITING )')
2111  FORMAT(5X,'FREE FLUX OF TRACER:               ' , G16.7 ,
     *          '  ( >0 : ENTERING  <0 : EXITING )')
2112  FORMAT(5X,'INCIDENT FLUX OF TRACER:           ' , G16.7 ,
     *          '  ( >0 : ENTERING  <0 : EXITING )')
2113  FORMAT(5X,'QUANTITY CREATED BY SOURCE TERM:   ' , G16.7 )
1120  FORMAT(5X,'ERREUR RELATIVE SUR LE TRACEUR : ',G16.7)
2120  FORMAT(5X,'RELATIVE ERROR ON TRACER : ',G16.7)
1130  FORMAT(5X,'ERREUR ABSOLUE SUR LE TRACEUR : ',G16.7)
2130  FORMAT(5X,'ABSOLUTE ERROR ON TRACER : ',G16.7)
1140  FORMAT(/,5X,'ERREUR RELATIVE CUMULEE SUR LE TRACEUR : ',G16.7)
2140  FORMAT(/,5X,'RELATIVE ERROR CUMULATED ON TRACER: ',G16.7)
1150  FORMAT(/,5X,'ERREUR ABSOLUE  CUMULEE SUR LE TRACEUR : ',G16.7)
2150  FORMAT(/,5X,'ABSOLUTE ERROR CUMULATED ON TRACER: ',G16.7)
1160  FORMAT(/,5X,'QUANTITE INITIALE DU TRACEUR      : ',G16.7,
     *       /,5X,'QUANTITE FINALE                   : ',G16.7)
1161  FORMAT(  5X,'QUANTITE ENTREE AUX FRONT. LIQUID.: ',G16.7,
     *            '  ( SI <0 QUANTITE SORTIE )')
1164  FORMAT(  5X,'QUANTITE CREEE PAR TERME SOURCE   : ',G16.7)
1165  FORMAT(  5X,'QUANTITE TOTALE PERDUE            : ',G16.7)
2160  FORMAT(/,5X,'INITIAL QUANTITY OF TRACER        : ',G16.7,
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
