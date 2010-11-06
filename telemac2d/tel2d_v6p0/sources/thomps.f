C                       *****************
                        SUBROUTINE THOMPS
C                       *****************
C
     *(HBOR,UBOR,VBOR,TBOR,U,V,H,T,ZF,X,Y,NBOR,FRTYPE,UNA,C,
     * UCONV,VCONV,T6,FU,FV,LIHBOR,LIUBOR,LIVBOR,LITBOR,LISPFR,T8,W1,
     * ITRAV2,
     * W1R,W2R,W3R,W4R,HBTIL,UBTIL,VBTIL,TBTIL,ZBTIL,SURDET,IKLE,
     * CF,SMH,IFABOR,NULONE,NELEM,MESH,
     * KP1BOR,XNEBOR,YNEBOR,NPOIN,NPTFR,LT,NIT,TEMPS,DT,GRAV,
     * DEBLIQ,FINLIQ,NTRAC,NFRLIQ,KSORT,LV,MSK,MASKEL,MASKPT,
     * NELBOR,NELMAX,IELM,NORD,FAIR,WINDX,WINDY,
     * VENT,HWIND,CORIOL,FCOR,SPHERI,
     * OPTPRO,MAREE,MARDAT,MARTIM,PHI0,OPTSOU,ISCE,DSCE,USCE,VSCE,T5,
     * COUROU,NPTH,VARCL,NVARCL,VARCLA,NUMLIQ,SHP,UNSV2D,HFROT)
C
C***********************************************************************
C  TELEMAC 2D VERSION 6.0    05/09/08    E DAVID (LHF) 04 76 33 42 36
C
C JMH 01/09/2008 : POINTS GROUPED REGARDLESS OF THEIR BOUNDARY NUMBER
C                  THIS IS TO HAVE AN ALGORITHM THAT WORK ALSO IN
C                  PARALLEL (BUT THE GROUPS WILL BE DIFFERENT) 
C                  CALLING GTSH11 ONCE AT THE BEGINNING AND NOT IN
C                  CARAFR (NOW GTSH11 IS INDEPENDENT OF THE VELOCITY)
C
C                  OTHER DIFFICULTIES FORBID PARALLELISM SO FAR
C
C***********************************************************************
C
C      FONCTION:    TRAITEMENT DES FRONTIERES LIQUIDES PAR LA
C                   METHODE DE THOMPSON - RESOLUTION PAR
C                   REMONTEE DES CARACTERISTIQUES
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   HBOR         |<-- |  HAUTEUR IMPOSEE.                            |
C |   UBOR         |<-- |  VITESSE U IMPOSEE.                          |
C |   VBOR         |<-- |  VITESSE V IMPOSEE.                          |
C |   TBOR         |<-- |  TRACEUR IMPOSE AU BORD                      |
C |    U,V         | -->|  COMPOSANTES DE LA VITESSE AU TEMPS N        |
C |    H           | -->|  HAUTEUR AU TEMPS N                          |
C |    T           | -->|  TRACEUR AU TEMPS N                          |
C |    ZF          | -->|  FOND                                        |
C |    X,Y         | -->|  COORDONNEES DES POINTS DU MAILLAGE          |
C |    NBOR        | -->|  ADRESSES DES POINTS DE BORD                 |
C |    FRTYPE      | -->|  TYPE DE FRONTIERES LIQUIDES                 |
C |    UNA         | -->|  TABLEAU DE TRAVAIL                          |
C |    C           | -->|  TABLEAU DE TRAVAIL : CELERITE DES ONDES     |
C |    UCONV,VCONV | -->|  TABLEAU DE TRAVAIL : CHAMPS DE VITESSE      |
C |                |    |  CONVECTEUR DES INVARIANTS DE RIEMANN        |
C |    FU,FV       | -->|  TABLEAU DE TRAVAIL : TERMES SOURCES         |
C |   LIHBOR       | -->|  CONDITIONS AUX LIMITES SUR H                |
C | LIUBOR,LIVBOR  | -->|  CONDITIONS AUX LIMITES SUR U ET V           |
C |   LITBOR       | -->|  CONDITIONS AUX LIMITES SUR LE TRACEUR       |
C |   LISPFR       | -->|  LISTE DES POINTS FRONTIERES CONTIGUS TRAITES|
C |                |    |  ENSEMBLES PAR LES CARACTERISTIQUES          |
C |   W1R,..,W4R   | -->|  INVARIANTS DE RIEMANN TRANSPORTES           |
C |   HBTIL..TBTIL | -->|  VALEURS DE H,..,T AU PIEDS DES              |
C |                |    |  CARACTERISTIQUES                            |
C |   KP1BOR       | -->|  NUMERO DU POINT FRONTIERE SUIVANT           |
C | XNEBOR,YNEBOR  | -->|  NORMALES EXTERIEURES AUX POINTS.
C |   NPOIN        | -->|  NOMBRE DE POINTS DU MAILLAGE.               |
C |   NPTFR        | -->|  NOMBRE DE POINTS FRONTIERE.                 |
C |   LT           | -->|  NUMERO DE L'ITERATION EN COURS              |
C |   TEMPS        | -->|  TEMPS                                       |
C |   DT           | -->|  PAS DE TEMPS                                |
C |   GRAV         | -->|  GRAVITE                                     |
C |   TRAC         | -->|  LOGIQUE INDIQUANT LA PRESENCE D'UN TRACEUR  |
C |   DEBLIQ       | -->|  TABLEAU D'INDICES DE DEBUT DE FRONTIERE LIQ.|
C |   FINLIQ       | -->|  TABLEAU D'INDICES DE FIN DE FRONTIERE LIQUI.|
C |   NFRLIQ       | -->|  NOMBRE DE FRONTIERES LIQUIDES
C |   KENT,KENTU,  | -->|  CONVENTION POUR LES TYPES DE CONDITIONS AUX |
C |   KSORT,       |    |  LIMITES PHYSIQUES                           |
C |   KINC         |    |  KENT:VALEURS IMPOSEES (SAUF U ET V)         |
C |                |    |  KENTU:U ET V IMPOSES                        |
C |                |    |  KSORT:VALEURS LIBRES                        |
C |                |    |  KINC:ONDE INCIDENTE                         |
C |   MSK          | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.        |
C |   MASKEL       | -->|  TABLEAU DE MASQUAGE DES ELEMENTS            |
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE         |
C |   NELBOR       | -->|  NUMEROS DES ELEMENTS ADJACENTS AUX BORDS    |
C |   NELMAX       | -->|  NOMBRE MAXIMUM D'ELEMENTS                   |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMME APPELE : FRICTI,PROSOU,CARAFR
C
C***********************************************************************
C
      USE BIEF
      USE INTERFACE_TELEMAC2D, EX_THOMPS => THOMPS
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NPTFR,LT,NIT,NPOIN,NELEM,NELMAX,NFRLIQ,LV
      INTEGER, INTENT(IN) :: NVARCL,NPTH,KSORT,IELM,NTRAC,HFROT
      INTEGER, INTENT(IN) :: OPTPRO,MARDAT(3),MARTIM(3),OPTSOU,ISCE(*)
      INTEGER, INTENT(IN) :: DEBLIQ(NFRLIQ),FINLIQ(NFRLIQ)
      INTEGER, INTENT(IN) :: NBOR(NPTFR),KP1BOR(NPTFR,2),NELBOR(NPTFR)
      INTEGER, INTENT(IN) :: IKLE(*),IFABOR(*),NULONE(*)
      INTEGER, INTENT(IN) :: LIHBOR(NPTFR),LIUBOR(NPTFR),LIVBOR(NPTFR)
      INTEGER, INTENT(IN) :: FRTYPE(NFRLIQ),NUMLIQ(NFRLIQ)
      INTEGER, INTENT(INOUT) :: LISPFR(NPTFR)  
C     ITRAV2 : TAILLE NPOIN
      INTEGER, INTENT(INOUT) :: ITRAV2(*)
      LOGICAL, INTENT(IN) :: VENT,MAREE,CORIOL,SPHERI,MSK,COUROU
      DOUBLE PRECISION, INTENT(IN) :: HWIND
      DOUBLE PRECISION, INTENT(INOUT) :: HBOR(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: UBOR(NPTFR),VBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: XNEBOR(NPTFR),YNEBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: SURDET(*),DSCE(*)
      DOUBLE PRECISION, INTENT(IN)    :: USCE(*),VSCE(*)
      DOUBLE PRECISION, INTENT(IN)  :: TEMPS,GRAV,DT,FAIR,FCOR,NORD,PHI0
      DOUBLE PRECISION, INTENT(INOUT) :: W1R(NPTFR),W2R(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: W3R(NPTFR),W4R(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: HBTIL(NPTFR),UBTIL(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: VBTIL(NPTFR),ZBTIL(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: T5(NPOIN),SHP(*)   
      TYPE(BIEF_OBJ), INTENT(IN)      :: WINDX,WINDY,MASKEL,MASKPT
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: W1,VARCL
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: FU,FV,T8,UNA,UCONV,VCONV,C,U,V
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: H,T,SMH,TBOR,TBTIL,T6   
      TYPE(BIEF_OBJ), INTENT(IN)      :: ZF,CF,LITBOR,UNSV2D 
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
      CHARACTER(LEN=32), INTENT(IN)   :: VARCLA(NVARCL)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K,NDEB,NFIN,IFRLIQ,NPT,KP,J,ITRAC,N                 
C                      
      DOUBLE PRECISION EPSIL,HMIN,HHBOR
C
      LOGICAL TSI
C
      DATA EPSIL /1.D-5/
      DATA TSI   /.FALSE./
      DATA HMIN  /2.D-2/
C
      INTRINSIC ABS
C
C-----------------------------------------------------------------------
C
      IF(NCSIZE.GT.1) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'THOMPSON NE MARCHE PAS EN PARALLELE'
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'THOMPSON NOT YET IMPLEMENTED IN PARALLEL'
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
C
C     SEULEMENT SI IL Y A DES FRONTIERES LIQUIDES
C
      IF(NFRLIQ.NE.0) THEN
C                                                  
C
C CALCUL DU TERME DE FROTTEMENT DANS UCONV
C "C,C" AJOUTE PAR JMH LE 08/08/2000 (STRUCTURES VERTICALES)
C MAIS NON PRIS EN COMPTE ICI (DERNIER ARGUMENT FALSE)
C
        CALL FRICTI(FU,FV,C,C,U,V,H,CF,MESH,T8,T6,.FALSE.,UNSV2D,
     *              MSK,MASKEL,HFROT)
C
C CALCUL DU TERME DT*UCONV*U
C
        CALL OS('X=CYZ   ', UCONV , FU , U , DT )
        CALL OS('X=CYZ   ', VCONV , FV , V , DT )
C
C CALCUL DES TERMES SOURCES DANS FU
C
        CALL PROSOU(FU,FV,SMH,U,V,H,GRAV,NORD,
     *              FAIR,WINDX,WINDY,VENT,HWIND,CORIOL,FCOR,
     *              SPHERI,TSI,MESH%COSLAT,MESH%SINLAT,
     *              TEMPS,LT,0,0,DSCE,ISCE,UNA,MESH,MSK,MASKEL,
     *              MAREE,MARDAT,MARTIM,PHI0,OPTSOU,COUROU,NPTH,
     *              VARCL,NVARCL,VARCLA,UNSV2D)
C
C ASSEMBLAGE DANS FU
C
        CALL OS('X=Y+CZ  ', FU , UCONV , FU , DT )
        CALL OS('X=Y+CZ  ', FV , VCONV , FV , DT )
C
C CALCUL DE LA CELERITE
C
        CALL OS('X=CY    ' , C , H , H , GRAV )
        CALL CLIP(C,0.D0,.TRUE.,1.D6,.FALSE.,0)
        CALL OS('X=SQR(Y)',X=C,Y=C )
C
C CORRECTION POUR LES BANCS DECOUVRANTS
C
        DO 9 K=1,NPOIN
          IF(H%R(K).LT.HMIN) THEN
            FU%R(K)=0.D0
            FV%R(K)=0.D0
          ENDIF
9       CONTINUE
C
C AVANCEMENT TEMPOREL (SPLITTING DU/DT=FU)
C
        CALL OS('X=X+Y   ',X=U,Y=FU)
        CALL OS('X=X+Y   ',X=V,Y=FV)
C
C REGROUPEMENT DES POINTS POSSEDANT LA MEME NORMALE
C A LA PRECISION 'EPSIL' PRES
C NPT : NOMBRE DE POINT CONTINUS
C LISPFR : LISTE DE CES POINTS DANS LA NUMEROTATION DES POINTS FRONTIERE
C
      NDEB=0
C
19    CONTINUE
      K=NDEB
20    CONTINUE
C
      K=K+1
      IF(K.GT.NPTFR) GO TO 1000
      IF(NUMLIQ(K).EQ.0) GO TO 20
      IF(FRTYPE(NUMLIQ(K)).EQ.2) THEN
C       PREMIER POINT THOMSON DE LA LISTE TROUVE
        NPT=1
        LISPFR(NPT)=K
      ELSE
        GO TO 20
      ENDIF
      NDEB = K
      KP   = K
30    CONTINUE
      KP=KP+1
      IF(KP.GT.NPTFR) GO TO 999
      IF(NUMLIQ(KP).EQ.0) GO TO 999
      IF(FRTYPE(NUMLIQ(KP)).EQ.2.AND.
     *   ABS(XNEBOR(KP)-XNEBOR(K)).LT.EPSIL.AND.
     *   ABS(YNEBOR(KP)-YNEBOR(K)).LT.EPSIL     ) THEN
        NPT=NPT+1
        LISPFR(NPT)=KP
        GO TO 30
      ENDIF
999   CONTINUE
      NDEB=LISPFR(NPT)
C
C MISE A JOUR DES VALEURS AUX BORDS SI SORTIE LIBRE
C
      DO J=1,NPT
        K=LISPFR(J)
        N=NBOR(K)
        IF(LIHBOR(K).EQ.KSORT) THEN
          HBOR(K)=H%R(N)
        ENDIF
        IF(LIUBOR(K).EQ.KSORT) THEN
          UBOR(K)=U%R(N)
        ENDIF
        IF(LIVBOR(K).EQ.KSORT) THEN
          VBOR(K)=V%R(N)
        ENDIF
      ENDDO
      IF(NTRAC.GT.0) THEN
        DO J=1,NPT
          K=LISPFR(J)
          N=NBOR(K)
          DO ITRAC=1,NTRAC
            IF(LITBOR%ADR(ITRAC)%P%I(K).EQ.KSORT) THEN
              TBOR%ADR(ITRAC)%P%R(K)=T%ADR(ITRAC)%P%R(N)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
C CALCUL DU CHAMP CONVECTEUR U SELON LA DIRECTION NORMALE
C A LA FRONTIERE
C
      CALL OS('X=CY    ',UNA  ,U  ,U  ,XNEBOR(LISPFR(1)))
      CALL OS('X=X+CY  ',UNA  ,V  ,V  ,YNEBOR(LISPFR(1)))
      CALL OS('X=CY    ',UCONV,UNA,UNA,XNEBOR(LISPFR(1)))
      CALL OS('X=CY    ',VCONV,UNA,UNA,YNEBOR(LISPFR(1)))
C
C CARACTERISTIQUES POUR LES POINTS GROUPES , CHAMPS CONVECTEUR U
C
      CALL GTSH11(UCONV%R,VCONV%R,X,Y,SHP,ITRAV2,
C                      INDIC  NLOC   (NE SERVENT PLUS)
     *            IKLE,ITRAV2,ITRAV2,NPOIN,NELEM,NELMAX,1,MSK,MASKEL%R)
      CALL CARAFR
     * ( U%R,V%R,H%R,T,UCONV%R,VCONV%R,X,Y,SHP, 
     *   SURDET , DT , IKLE , IFABOR , ITRAV2 ,
     *   NBOR , NELBOR , NULONE , IELM , NELEM , NELMAX , 
     *   NPOIN , 3 , NPTFR , 
     *   MSK , MASKEL%R , MASKPT%R ,  NPT , LISPFR , NTRAC ,
     *   HBTIL , UBTIL , VBTIL , TBTIL , ZBTIL , ZF%R,T5)
C
C CALCUL DES INVARIANTS DE RIEMANN W1 ET W4 (DEUXIEME DIMENSION DE TBOR) 
C TRANSPORTES PAR CE CHAMP
C
      DO J=1,NPT
       K=LISPFR(J)
       IF(UNA%R(NBOR(K)).GE.0.D0) THEN
        W1R(K)=-HBTIL(K)*(XNEBOR(LISPFR(1))*(VBTIL(K)-V%R(NBOR(K)))-
     *                    YNEBOR(LISPFR(1))*(UBTIL(K)-U%R(NBOR(K))))
        IF(NTRAC.GT.0) THEN
          DO ITRAC=1,NTRAC
            TBOR%ADR(ITRAC)%P%R(K+NPTFR)=HBTIL(K)*
     *      (TBTIL%ADR(ITRAC)%P%R(K)-T%ADR(ITRAC)%P%R(NBOR(K)))
          ENDDO
        ENDIF
       ELSE
        W1R(K)=-HBOR(K)*(XNEBOR(LISPFR(1))*(VBOR(K)-V%R(NBOR(K)))-
     *                   YNEBOR(LISPFR(1))*(UBOR(K)-U%R(NBOR(K))) )
        IF(NTRAC.GT.0) THEN
          DO ITRAC=1,NTRAC
            TBOR%ADR(ITRAC)%P%R(K+NPTFR)=
     *      HBOR(K)*(TBOR%ADR(ITRAC)%P%R(K)-T%ADR(ITRAC)%P%R(NBOR(K)))
          ENDDO
        ENDIF
       ENDIF
      ENDDO
C
C CALCUL DU CHAMP CONVECTEUR U+C SELON LA DIRECTION NORMALE
C A LA FRONTIERE
C
      CALL OS('X=X+CY  ',UNA  , C   , C   ,             1.D0 )
      CALL OS('X=CY    ',UCONV, UNA , UNA , XNEBOR(LISPFR(1)) )
      CALL OS('X=CY    ',VCONV, UNA , UNA , YNEBOR(LISPFR(1)) )
C
C CARACTERISTIQUES POUR LES POINTS GROUPES , CHAMPS U+C
C
      CALL GTSH11(UCONV%R,VCONV%R,X,Y,SHP,ITRAV2,
C                      INDIC  NLOC   (NE SERVENT PLUS)
     *            IKLE,ITRAV2,ITRAV2,NPOIN,NELEM,NELMAX,1,MSK,MASKEL%R)
      CALL CARAFR
     * ( U%R,V%R,H%R,T,UCONV%R,VCONV%R,X,Y,SHP, 
     *   SURDET,DT,IKLE,IFABOR,ITRAV2,
     *   NBOR,NELBOR,NULONE,IELM,NELEM,NELMAX, 
     *   NPOIN,3,NPTFR, 
     *   MSK,MASKEL%R,MASKPT%R,NPT,LISPFR,NTRAC,
     *   HBTIL,UBTIL,VBTIL,TBTIL,ZBTIL,ZF%R,T5)
C
C CALCUL DES INVARIANTS DE RIEMANN W2 TRANSPORTE PAR CE CHAMP CONVECTEUR
C
      DO 50 J=1,NPT
       K=LISPFR(J)
       IF (UNA%R(NBOR(K)).GE.0.D0) THEN
        W2R(K)=(-ZF%R(NBOR(K))+HBTIL(K)+ZBTIL(K))*C%R(NBOR(K))+
     *         HBTIL(K)*(XNEBOR(LISPFR(1))*(UBTIL(K)-U%R(NBOR(K)))+
     *                   YNEBOR(LISPFR(1))*(VBTIL(K)-V%R(NBOR(K))) )
       ELSE
        W2R(K)=HBOR(K)*(C%R(NBOR(K))+
     *                 (XNEBOR(LISPFR(1))*(UBOR(K)-U%R(NBOR(K)))+
     *                  YNEBOR(LISPFR(1))*(VBOR(K)-V%R(NBOR(K)))) )
       ENDIF
50    CONTINUE
C
C CALCUL DU CHAMP CONVECTEUR U-C SELON LA DIRECTION NORMALE
C A LA FRONTIERE
C
      CALL OS('X=X+CY  ',X=UNA, Y=C , C=-2.D0 )
      CALL OS('X=CY    ',X=UCONV,Y=UNA , C=XNEBOR(LISPFR(1)) )
      CALL OS('X=CY    ',X=VCONV,Y=UNA , C=YNEBOR(LISPFR(1)) )
C
C CARACTERISTIQUES POUR LES POINTS GROUPES , CHAMPS U+C
C
      CALL GTSH11(UCONV%R,VCONV%R,X,Y,SHP,ITRAV2,
C                      INDIC  NLOC   (NE SERVENT PLUS)
     *            IKLE,ITRAV2,ITRAV2,NPOIN,NELEM,NELMAX,1,MSK,MASKEL%R)
      CALL CARAFR
     * ( U%R,V%R,H%R,T,UCONV%R,VCONV%R,X,Y,SHP,
     *   SURDET,DT,IKLE,IFABOR,ITRAV2,
     *   NBOR,NELBOR,NULONE,IELM,NELEM,NELMAX,NPOIN,3,NPTFR, 
     *   MSK,MASKEL%R,MASKPT%R,NPT,LISPFR,NTRAC,
     *   HBTIL,UBTIL,VBTIL,TBTIL,ZBTIL,ZF%R,T5)
C
C CALCUL DES INVARIANTS DE RIEMANN W3 TRANSPORTE PAR CE CHAMP CONVECTEUR
C
      DO 60 J=1,NPT
       K=LISPFR(J)
       IF(UNA%R(NBOR(K)).GE.0.D0) THEN
        W3R(K)=(-ZF%R(NBOR(K))+HBTIL(K)+ZBTIL(K))*C%R(NBOR(K))-
     *         HBTIL(K)*(XNEBOR(LISPFR(1))*(UBTIL(K)-U%R(NBOR(K)))+
     *                   YNEBOR(LISPFR(1))*(VBTIL(K)-V%R(NBOR(K))) )
       ELSE
        W3R(K)=HBOR(K)*(C%R(NBOR(K))-
     *                 (XNEBOR(LISPFR(1))*(UBOR(K)-U%R(NBOR(K)))+
     *                  YNEBOR(LISPFR(1))*(VBOR(K)-V%R(NBOR(K)))) )
       ENDIF
60    CONTINUE
C
C RECONSTRUCTION DES VARIABLES DE TELEMAC-2D
C
C POUR LES BANCS DECOUVRANTS (ICI H<HMIN) IL FAUT LAISSER FAIRE CE QUE
C L'ANCIENNE VERSION AVAIT PREVU
C
      DO 70 J=1,NPT
C
        K=LISPFR(J)
        IF(C%R(NBOR(K))**2.GT.GRAV*HMIN) THEN
          HBOR(K)=(W2R(K)+W3R(K))/(2*C%R(NBOR(K)))
          IF(HBOR(K).GT.HMIN) THEN
C           BEWARE TIDAL FLATS, AND HIDDEN PARAMETER 0.1
            HHBOR=MAX(0.1D0,HBOR(K))
            UBOR(K)=(YNEBOR(LISPFR(1))*W1R(K)+XNEBOR(LISPFR(1))*W2R(K)-
     *HBOR(K)*C%R(NBOR(K))*XNEBOR(LISPFR(1)))/HHBOR+U%R(NBOR(K))
            VBOR(K)=(YNEBOR(LISPFR(1))*W2R(K)-XNEBOR(LISPFR(1))*W1R(K)-
     *HBOR(K)*C%R(NBOR(K))*YNEBOR(LISPFR(1)))/HHBOR+V%R(NBOR(K))
            IF(NTRAC.GT.0) THEN
              DO ITRAC=1,NTRAC
                TBOR%ADR(ITRAC)%P%R(K)=
     *          TBOR%ADR(ITRAC)%P%R(K+NPTFR)/HHBOR+
     *          T%ADR(ITRAC)%P%R(NBOR(K))
              ENDDO
            ENDIF
          ELSE
C           ON DEVIENT DECOUVERT
            HBOR(K)=MAX(0.D0,HBOR(K))
            UBOR(K)=0.D0
            VBOR(K)=0.D0
          ENDIF
        ELSE
C         ON ETAIT DECOUVERT, H EST DONNE PAR BORD
          UBOR(K)=0.D0
          VBOR(K)=0.D0
        ENDIF
C
70    CONTINUE
C
      IF(NDEB.LE.NPTFR) GO TO 19
C
C     TEST IF(NFRLIQ.GT.0)...
      ENDIF
C
C
1000  CONTINUE
C
C
C RECUPERATION DE LA VALEUR EXACTE DE U ET V
C CAR BESOIN POUR LES SOUS ITERATIONS
C
      CALL OS('X=X-Y   ' , X=U , Y=FU )
      CALL OS('X=X-Y   ' , X=V , Y=FV )
C
C-----------------------------------------------------------------------
C
      RETURN
      END
