C                       *************
                        SUBROUTINE OM
C                       *************
C
     *( OP , M , N , D , C , MESH )
C
C***********************************************************************
C BIEF VERSION 6.0      05/02/2010    J-M HERVOUET (LNHE) 01 30 87 80 18
C                                      F  LEPEINTRE (LNH) 30 87 78 54
C
C  ALGIANE FROEHLY LE 13/02/2008 : AJOUT OM1113 ET OM1311
C  JMH 05/02/2010 : CALL TO OMSEGBOR MODIFIED, OMSEGPAR SUPPRESSED
C***********************************************************************
C
C FONCTION : OPERATIONS SUR DES MATRICES
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES MATRICES M ET N, LA MATRICE DIAGONALE D ET LA
C   CONSTANTE C.
C
C   LE RESULTAT EST LA MATRICE M.
C
C      OP = 'M=N     '  : COPIE DE N DANS M
C      OP = 'M=CN    '  : PRODUIT DE N PAR LA CONSTANTE C
C      OP = 'M=M+CN  '  : ON AJOUTE CN A M
C      OP = 'M=MD    '  : PRODUIT DE M PAR D A DROITE
C      OP = 'M=DM    '  : PRODUIT DE M PAR D A GAUCHE
C      OP = 'M=DMD   '  : PRODUIT DE M A DROITE ET A GAUCHE PAR D
C      OP = 'M=0     '  : ANNULATION DE M (A VERIFIER)
C      OP = 'M=X(M)  '  : PASSAGE A UNE FORME NON SYMETRIQUE
C      OP = 'M=TN    '  : TRANSPOSEE DE N MISE DANS M
C      OP = 'M=MSK(M)'  : MASQUAGE DES TERMES EXTRADIAGONAUX
C
C   ATTENTION | SI ON AJOUTE UNE NOUVELLE OPERATION
C   SI OP CONTIENT N CELA SIGNIFIE QUE LA MATRICE N EST UTILISEE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C |    OP          | -->| OPERATION A EFFECTUER
C |    DM          | -->| DIAGONALE DE M
C |    TYPDIM      | -->| TYPE DE LA DIAGONALE DE M ('Q','I','0')
C |    XM          | -->| TERMES EXTRA-DIAGONAUX DE M
C |    TYPEXM      | -->| TYPE DE TERMES EXTRADIAGONAUX ('Q','S','0')
C |    DN          | -->| DIAGONALE DE N
C |    TYPDIN      | -->| TYPE DE LA DIAGONALE DE N ('Q','I','0')
C |    XN          | -->| TERMES EXTRA-DIAGONAUX DE N
C |    TYPEXN      | -->| TYPE DE TERMES EXTRADIAGONAUX
C |    D           | -->| MATRICE DIAGONALE : PEUT ETRE UNE STRUCTURE
C |                |    | OU UN TABLEAU PROVISOIREMENT (VOIR TEST SUR D)
C |    C           | -->| CONSTANTE DONNEE
C |    IKLE        | -->| CORRESPONDANCE NUMEROTATIONS LOCALE ET GLOBALE
C |    NELEM       | -->| NOMBRE D'ELEMENTS DU MAILLAGE
C |    NELMAX      | -->| NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    | (CAS D'UN MAILLAGE ADAPTATIF)
C |    IELM1       | -->| TYPE D'ELEMENT
C |    IELM2       | -->| TYPE D'ELEMENT
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   : OM1111 , OM2121 , OM4141 , PLANTE
C
C***********************************************************************
C
      USE BIEF, EX_OM => OM
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(IN)    :: C
      CHARACTER(LEN=8), INTENT(IN)    :: OP
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: M
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: N,D
      TYPE(BIEF_MESH) , INTENT(IN)    :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELM1,IELM2,IELN1,IELN2,NELEM,NELMAX
      INTEGER NDIAGM,NDIAGN,NDIAGX,NSEG1,NSEG2
      INTEGER STOM,STON,NPOIN,NPTFR,NPTFX,MDIAGX
      INTEGER SIZXN,SZMXN,NETAGE
C
      CHARACTER*1 TYPDIM,TYPEXM,TYPDIN,TYPEXN
C
      INTEGER, DIMENSION(:), POINTER :: IKLE
C
C-----------------------------------------------------------------------
C
C  CAS OU LA STRUCTURE DE M DEVIENT CELLE DE N
C

      IF(OP(3:8).EQ.'N     '.OR.OP(3:8).EQ.'CN    ') THEN
        CALL CPSTMT(N,M)
      ELSEIF(OP(3:8).EQ.'TN    ') THEN
        CALL CPSTMT(N,M,TRANS=.TRUE.)
      ENDIF
C
C  EXTRACTION DES CARACTERISTIQUES DE LA MATRICE M
C
      TYPDIM = M%TYPDIA
      TYPEXM = M%TYPEXT
      STOM = M%STO
      NDIAGM = M%D%DIM1
      MDIAGX = M%D%MAXDIM1
      IELM1 = M%ELMLIN
      IELM2 = M%ELMCOL
C
      IF(OP(3:8).EQ.'X(M)  ') THEN
        IF(M%ELMLIN.NE.M%ELMCOL) THEN
          IF(LNG.EQ.1) WRITE(LU,900) M%NAME
          IF(LNG.EQ.2) WRITE(LU,901) M%NAME
900       FORMAT(1X,'OM (BIEF) : M (NOM REEL : ',A6,') NON CARREE')
901       FORMAT(1X,'OM (BIEF) : M (REAL NAME: ',A6,') NOT SQUARE')
          IF(LNG.EQ.1) WRITE(LU,700)  
          IF(LNG.EQ.2) WRITE(LU,701)  
700       FORMAT(1X,'            EST DEJA NON SYMETRIQUE')
701       FORMAT(1X,'            IS ALREADY NON SYMMETRICAL')
          CALL PLANTE(1)
          STOP
        ENDIF
        IF(M%X%MAXDIM2*M%X%MAXDIM1.LT.
     *      DIM1_EXT(IELM1,IELM2,STOM,'Q')*
     *      DIM2_EXT(IELM1,IELM2,STOM,'Q')    ) THEN
            IF(LNG.EQ.1) WRITE(LU,400) M%NAME
            IF(LNG.EQ.2) WRITE(LU,401) M%NAME
            IF(LNG.EQ.1) WRITE(LU,800)  
            IF(LNG.EQ.2) WRITE(LU,801)  
800         FORMAT(1X,'            POUR DEVENIR NON SYMETRIQUE')
801         FORMAT(1X,'            TO BECOME NON SYMMETRICAL')
            CALL PLANTE(1)
            STOP
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C
C  EXTRACTION EVENTUELLE DES CARACTERISTIQUES DE LA MATRICE N
      IF(INCLUS(OP,'N')) THEN
        TYPDIN = N%TYPDIA
        TYPEXN = N%TYPEXT
        STON   = N%STO
        NDIAGN = N%D%DIM1
        NDIAGX = N%D%MAXDIM1
        
C
        SIZXN = N%X%DIM1
        SZMXN = N%X%MAXDIM1
C
C 07/02/03 : DIVISION BY NDIAGN=0 AVOIDED (SUBDOMAIN WITHOUT BOUNDARY POINTS
C            IN PARALLEL). COURTESY OLIVER GOETHEL (HANNOVER UNIVERSITY)
C
        IF(NDIAGN.GT.0) THEN
          NETAGE = NDIAGM/NDIAGN - 1
        ELSE
          NETAGE = 0
        ENDIF
C
        IELN1 = N%ELMLIN
        IELN2 = N%ELMCOL
        IF(NDIAGN.GT.MDIAGX) THEN
         IF(LNG.EQ.1) WRITE(LU,400) M%NAME
         IF(LNG.EQ.2) WRITE(LU,401) M%NAME
400      FORMAT(1X,'OM (BIEF) : M (NOM REEL : ',A6,') TROP PETITE')
401      FORMAT(1X,'OM (BIEF) : M (REAL NAME: ',A6,') TOO SMALL')
         STOP
        ENDIF
      ELSE
        IELN1 = IELM1
        IELN2 = IELM2
        STON = STOM
      ENDIF      
C
C-----------------------------------------------------------------------
C
C  DEPLOIEMENT DE LA STRUCTURE DE MAILLAGE
C
C     MATRICE NORMALE
      IF(DIMENS(IELM1).EQ.MESH%DIM) THEN
        IKLE=>MESH%IKLE%I
        NELEM = MESH%NELEM
        NELMAX= MESH%NELMAX
      ELSE
C     MATRICE DE BORD
        IKLE=>MESH%IKLBOR%I
        NELEM  = MESH%NELEB
        NELMAX = MESH%NELEBX
      ENDIF
C
      NPOIN= MESH%NPOIN
      NPTFR= MESH%NPTFR
      NPTFX= MESH%NPTFRX
C  
C-----------------------------------------------------------------------
C
C  STOCKAGE EBE CLASSIQUE :
C
      IF(STOM.EQ.1.AND.STON.EQ.1) THEN
C
      IF(IELM1.EQ.1.AND.IELM2.EQ.1) THEN
C
C     ELEMENTS A 2 POINTS
C
      
      CALL OM0101(OP , M%D%R,TYPDIM,M%X%R,TYPEXM,
     *                 N%D%R,TYPDIN,N%X%R,TYPEXN, D%R,C,
     *                 IKLE,NELEM,NELMAX,NDIAGM)
    
       
C
C     ELEMENTS A 3 POINTS
C
      ELSEIF( (IELM1.EQ.2 .AND.IELM2.EQ.2 ) .OR.  
     *        (IELM1.EQ.11.AND.IELM2.EQ.11) .OR.  
     *        (IELM1.EQ.61.AND.IELM2.EQ.61) .OR.
     *        (IELM1.EQ.81.AND.IELM2.EQ.81)       ) THEN
C
        IF( (IELN1.EQ.2 .AND.IELN2.EQ.2 ) .OR.  
     *      (IELN1.EQ.11.AND.IELN2.EQ.11) .OR.
     *      (IELN1.EQ.61.AND.IELN2.EQ.61) .OR.       
     *      (IELN1.EQ.81.AND.IELN2.EQ.81)         ) THEN
C
          CALL OM1111(OP , M%D%R,TYPDIM,M%X%R,TYPEXM ,
     *                     N%D%R,TYPDIN,N%X%R,TYPEXN,D%R,C,
     *                     IKLE,NELEM,NELMAX,NDIAGM)
C
        ELSEIF(IELN1.EQ.1.AND.IELN2.EQ.1) THEN
C
          CALL OM1101(OP , M%D%R,TYPDIM,M%X%R,TYPEXM,
     *                     N%D%R,TYPDIN,N%X%R,TYPEXN, C,
     *                     MESH%NULONE%I,MESH%NELBOR%I,
     *                     MESH%NBOR%I,
     *                     NELMAX,NDIAGM,NDIAGN,NDIAGX)
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,100) M%NAME
          IF (LNG.EQ.2) WRITE(LU,101) M%NAME
          IF (LNG.EQ.1) WRITE(LU,200) IELM1,IELM2
          IF (LNG.EQ.2) WRITE(LU,201) IELM1,IELM2
          IF (LNG.EQ.1) WRITE(LU,150) N%NAME
          IF (LNG.EQ.2) WRITE(LU,151) N%NAME
          IF (LNG.EQ.1) WRITE(LU,250) IELN1,IELN2
          IF (LNG.EQ.2) WRITE(LU,251) IELN1,IELN2
          IF (LNG.EQ.1) WRITE(LU,300)
          IF (LNG.EQ.2) WRITE(LU,301)
          CALL PLANTE(1)
          STOP
        ENDIF
C
C     ELEMENTS A 4 POINTS
C
      ELSEIF( (IELM1.EQ.21.AND.IELM2.EQ.21) .OR.
     *        (IELM1.EQ.71.AND.IELM2.EQ.71) .OR. 
     *        (IELM1.EQ.31.AND.IELM2.EQ.31) .OR. 
     *        (IELM1.EQ.51.AND.IELM2.EQ.51) .OR. 
     *        (IELM1.EQ.12.AND.IELM2.EQ.12)      ) THEN
C
        IF(   (IELN1.EQ.21.AND.IELN2.EQ.21) .OR.
     *        (IELN1.EQ.71.AND.IELN2.EQ.71) .OR. 
     *        (IELN1.EQ.31.AND.IELN2.EQ.31) .OR. 
     *        (IELN1.EQ.51.AND.IELN2.EQ.51) .OR. 
     *        (IELN1.EQ.12.AND.IELN2.EQ.12)      ) THEN
C
          CALL OM2121(OP , M%D%R,TYPDIM,M%X%R,TYPEXM ,
     *                N%D%R,TYPDIN,N%X%R,TYPEXN, D%R,C,
     *                IKLE,NELEM,NELMAX,NDIAGM)
        ELSEIF(  (IELM1.EQ.12.AND.IELM2.EQ.12) .AND.
     *           (IELN1.EQ.1 .AND.IELN2.EQ.1 )   ) THEN
C      
          CALL OM1201(OP , M%D%R,TYPDIM,M%X%R,TYPEXM,
     *                     N%D%R,TYPDIN,N%X%R,TYPEXN, C,
     *                     MESH%NULONE%I,MESH%NELBOR%I,
     *                     MESH%NBOR%I,
     *                     NELMAX,NDIAGM,NDIAGN,NDIAGX)
        ELSEIF(( (IELM1.EQ.51.AND.IELM2.EQ.51) .AND.
     *           (IELN1.EQ.61.AND.IELN2.EQ.61))   ) THEN    
C         PRISMES DECOUPES EN TETRAEDRES M MATRICE INTERIEURE
C                                        N MATRICE DE BORD LATERAL
          CALL OM5161(OP,M%D%R,TYPDIM,M%X%R,TYPEXM,
     *                N%D%R,TYPDIN,N%X%R,TYPEXN,C,
     *                MESH%NULONE%I,MESH%NELBOR%I,MESH%NBOR%I,
     *                NELMAX,NDIAGN,SIZXN,SZMXN)
C
        ELSEIF(( (IELM1.EQ.31.AND.IELM2.EQ.31) .AND.
     *           (IELN1.EQ.81.AND.IELN2.EQ.81))   ) THEN    
C         TETRAEDRES NON STRUCTURES      M MATRICE INTERIEURE
C                                        N MATRICE DE BORD LATERAL
          CALL OM3181(OP,M%D%R,TYPDIM,M%X%R,TYPEXM,
     *                N%D%R,TYPDIN,N%X%R,TYPEXN,C,
     *                MESH%NULONE%I,MESH%NELBOR%I,MESH%NBOR%I,
     *                NELMAX,NDIAGN,MESH%NELEB,SZMXN)
C
        ELSEIF(  (IELM1.EQ.51.AND.IELM2.EQ.51) .AND.
     *           (IELN1.EQ.11.AND.IELN2.EQ.11)   ) THEN
C         PRISMES DECOUPES EN TETRAEDRES M MATRICE INTERIEURE
C                                        N MATRICE DE BORD FOND OU SURFACE
C                                        OPERATIONS M+NF ET M+NS
          CALL OM5111(OP , M%D%R,TYPDIM,M%X%R,TYPEXM,
     *                N%D%R,TYPDIN,N%X%R,TYPEXN, C,
     *                NDIAGN,NDIAGX,SIZXN,SZMXN,NETAGE,NELMAX)
C     
        ELSE
          IF (LNG.EQ.1) WRITE(LU,100) M%NAME
          IF (LNG.EQ.2) WRITE(LU,101) M%NAME
          IF (LNG.EQ.1) WRITE(LU,200) IELM1,IELM2
          IF (LNG.EQ.2) WRITE(LU,201) IELM1,IELM2
          IF (LNG.EQ.1) WRITE(LU,150) N%NAME
          IF (LNG.EQ.2) WRITE(LU,151) N%NAME
          IF (LNG.EQ.1) WRITE(LU,250) IELN1,IELN2
          IF (LNG.EQ.2) WRITE(LU,251) IELN1,IELN2
          IF (LNG.EQ.1) WRITE(LU,300)
          IF (LNG.EQ.2) WRITE(LU,301)
          CALL PLANTE(1)
          STOP
        ENDIF
C
C     MATRICES RECTANGULAIRES 3 ET 4 POINTS
C
      ELSEIF(IELM1.EQ.11.AND.IELM2.EQ.12) THEN
C
        IF((IELN1.EQ.11.AND.IELN2.EQ.12).OR.
     *     (IELN1.EQ.12.AND.IELN2.EQ.11.AND.OP(1:4).EQ.'M=TN')) THEN
          CALL OM1112(OP , M%D%R,TYPDIM,M%X%R,TYPEXM ,
     *                     N%D%R,TYPDIN,N%X%R,TYPEXN, D%R,C,
     *                     IKLE,NELEM,NELMAX,NDIAGM)
        ELSE
          IF (LNG.EQ.1) WRITE(LU,100) M%NAME
          IF (LNG.EQ.2) WRITE(LU,101) M%NAME
          IF (LNG.EQ.1) WRITE(LU,200) IELM1,IELM2
          IF (LNG.EQ.2) WRITE(LU,201) IELM1,IELM2
          IF (LNG.EQ.1) WRITE(LU,150) N%NAME
          IF (LNG.EQ.2) WRITE(LU,151) N%NAME
          IF (LNG.EQ.1) WRITE(LU,250) IELN1,IELN2
          IF (LNG.EQ.2) WRITE(LU,251) IELN1,IELN2
          IF (LNG.EQ.1) WRITE(LU,300)
          IF (LNG.EQ.2) WRITE(LU,301)
          CALL PLANTE(1)
          STOP
        ENDIF
C
C     MATRICES RECTANGULAIRES 4 ET 3 POINTS
C
      ELSEIF(IELM1.EQ.12.AND.IELM2.EQ.11) THEN
C
        IF((IELN1.EQ.12.AND.IELN2.EQ.11).OR.
     *     (IELN1.EQ.11.AND.IELN2.EQ.12.AND.OP(1:4).EQ.'M=TN')) THEN
          CALL OM1211(OP , M%D%R,TYPDIM,M%X%R,TYPEXM ,
     *                N%D%R,TYPDIN,N%X%R,TYPEXN, D%R,C,
     *                IKLE,NELEM,NELMAX,NDIAGM)
        ELSE
          IF (LNG.EQ.1) WRITE(LU,100) M%NAME
          IF (LNG.EQ.2) WRITE(LU,101) M%NAME
          IF (LNG.EQ.1) WRITE(LU,200) IELM1,IELM2
          IF (LNG.EQ.2) WRITE(LU,201) IELM1,IELM2
          IF (LNG.EQ.1) WRITE(LU,150) N%NAME
          IF (LNG.EQ.2) WRITE(LU,151) N%NAME
          IF (LNG.EQ.1) WRITE(LU,250) IELN1,IELN2
          IF (LNG.EQ.2) WRITE(LU,251) IELN1,IELN2
          IF (LNG.EQ.1) WRITE(LU,300)
          IF (LNG.EQ.2) WRITE(LU,301)
          CALL PLANTE(1)
          STOP
        ENDIF
C
C     ELEMENTS A 6 POINTS
C
      ELSEIF( (IELM1.EQ.41.AND.IELM2.EQ.41).OR.
     *        (IELM1.EQ.13.AND.IELM2.EQ.13)      ) THEN
C
        IF( (IELN1.EQ.41.AND.IELN2.EQ.41).OR.
     *      (IELN1.EQ.13.AND.IELN2.EQ.13)        ) THEN
C
          CALL OM4141(OP , M%D%R,TYPDIM,M%X%R,TYPEXM ,
     *                N%D%R,TYPDIN,N%X%R,TYPEXN, D%R,C,
     *                IKLE,NELEM,NELMAX,NDIAGM)
C
        ELSEIF(IELN1.EQ.71.AND.IELN2.EQ.71) THEN
C
          CALL OM4121(OP,M%D%R,TYPDIM,M%X%R,TYPEXM,
     *                   N%D%R,TYPDIN,N%X%R,TYPEXN,C,
     *                   MESH%NULONE%I,MESH%NELBOR%I,MESH%NBOR%I,
     *                   NELMAX,NDIAGN,SIZXN,SZMXN)
C
        ELSEIF(IELN1.EQ.11.AND.IELN2.EQ.11) THEN
C
          CALL OM4111(OP , M%D%R,TYPDIM,M%X%R,TYPEXM,
     *                     N%D%R,TYPDIN,N%X%R,TYPEXN, C,
     *                     NDIAGN,NDIAGX,SIZXN,SZMXN,NETAGE,NELMAX)
C
        ELSEIF(  (IELM1.EQ.13.AND.IELM2.EQ.13) .AND.
     *           (IELN1.EQ.2 .AND.IELN2.EQ.2 )   ) THEN
C      
          CALL OM1302(OP , M%D%R,TYPDIM,M%X%R,TYPEXM,
     *                     N%D%R,TYPDIN,N%X%R,TYPEXN, C,
     *                     MESH%NULONE%I,MESH%NELBOR%I,
     *                     MESH%NBOR%I,
     *                     NELMAX,NDIAGM,NPTFR,NPTFX)       
C       
        ELSE
          IF (LNG.EQ.1) WRITE(LU,100) M%NAME
          IF (LNG.EQ.2) WRITE(LU,101) M%NAME
          IF (LNG.EQ.1) WRITE(LU,200) IELM1,IELM2
          IF (LNG.EQ.2) WRITE(LU,201) IELM1,IELM2
          IF (LNG.EQ.1) WRITE(LU,150) N%NAME
          IF (LNG.EQ.2) WRITE(LU,151) N%NAME
          IF (LNG.EQ.1) WRITE(LU,250) IELN1,IELN2
          IF (LNG.EQ.2) WRITE(LU,251) IELN1,IELN2
          IF (LNG.EQ.1) WRITE(LU,300)
          IF (LNG.EQ.2) WRITE(LU,301)
          CALL PLANTE(1)
          STOP
        ENDIF
C
C     MATRICES RECTANGULAIRES 3 ET 6 POINTS
C
      ELSEIF(IELM1.EQ.11.AND.IELM2.EQ.13) THEN
C
        IF((IELN1.EQ.11.AND.IELN2.EQ.13).OR.
     *     (IELN1.EQ.13.AND.IELN2.EQ.11.AND.OP(1:4).EQ.'M=TN')) THEN
          CALL OM1113(OP , M%D%R,TYPDIM,M%X%R,TYPEXM ,
     *                     N%D%R,TYPDIN,N%X%R,TYPEXN, D%R,C,
     *                     IKLE,NELEM,NELMAX,NDIAGM)
        ELSE
          IF (LNG.EQ.1) WRITE(LU,100) M%NAME
          IF (LNG.EQ.2) WRITE(LU,101) M%NAME
          IF (LNG.EQ.1) WRITE(LU,200) IELM1,IELM2
          IF (LNG.EQ.2) WRITE(LU,201) IELM1,IELM2
          IF (LNG.EQ.1) WRITE(LU,150) N%NAME
          IF (LNG.EQ.2) WRITE(LU,151) N%NAME
          IF (LNG.EQ.1) WRITE(LU,250) IELN1,IELN2
          IF (LNG.EQ.2) WRITE(LU,251) IELN1,IELN2
          IF (LNG.EQ.1) WRITE(LU,300)
          IF (LNG.EQ.2) WRITE(LU,301)
          CALL PLANTE(1)
          STOP
        ENDIF
C
C     MATRICES RECTANGULAIRES A 6 ET 3 POINTS
C
      ELSEIF(IELM1.EQ.13.AND.IELM2.EQ.11) THEN
C
        IF((IELN1.EQ.13.AND.IELN2.EQ.11).OR.
     *     (IELN1.EQ.11.AND.IELN2.EQ.13.AND.OP(1:4).EQ.'M=TN')) THEN
          CALL OM1311(OP , M%D%R,TYPDIM,M%X%R,TYPEXM ,
     *                N%D%R,TYPDIN,N%X%R,TYPEXN, D%R,C,
     *                IKLE,NELEM,NELMAX,NDIAGM)
        ELSE
          IF (LNG.EQ.1) WRITE(LU,100) M%NAME
          IF (LNG.EQ.2) WRITE(LU,101) M%NAME
          IF (LNG.EQ.1) WRITE(LU,200) IELM1,IELM2
          IF (LNG.EQ.2) WRITE(LU,201) IELM1,IELM2
          IF (LNG.EQ.1) WRITE(LU,150) N%NAME
          IF (LNG.EQ.2) WRITE(LU,151) N%NAME
          IF (LNG.EQ.1) WRITE(LU,250) IELN1,IELN2
          IF (LNG.EQ.2) WRITE(LU,251) IELN1,IELN2
          IF (LNG.EQ.1) WRITE(LU,300)
          IF (LNG.EQ.2) WRITE(LU,301)
          CALL PLANTE(1)
          STOP
        ENDIF
C
C  COMBINAISON DE IELM1 ET IELM2 NON PREVUE : ERREUR
C
      ELSE
         IF (LNG.EQ.1) WRITE(LU,100) M%NAME
         IF (LNG.EQ.2) WRITE(LU,101) M%NAME
         IF (LNG.EQ.1) WRITE(LU,200) IELM1,IELM2
         IF (LNG.EQ.2) WRITE(LU,201) IELM1,IELM2
         IF (LNG.EQ.1) WRITE(LU,410) STOM,STON
         IF (LNG.EQ.2) WRITE(LU,411) STOM,STON
         IF (LNG.EQ.1) WRITE(LU,300)
         IF (LNG.EQ.2) WRITE(LU,301)
         CALL PLANTE(1)
         STOP
      ENDIF
C
      ELSEIF(STOM.EQ.3.AND.STON.EQ.3) THEN
C
C  STOCKAGE PAR SEGMENT
C
        IF(M%ELMCOL.NE.N%ELMCOL.OR.M%ELMLIN.NE.N%ELMLIN) THEN
          WRITE(LU,*) 'M ET N DE STRUCTURES DIFFERENTES'
          CALL PLANTE(1)
          STOP
        ENDIF
C
        NSEG1 = NBSEG(M%ELMLIN)
        NSEG2 = NBSEG(M%ELMCOL)
C
C       IN LINEAR-QUADRATIC RECTANGULAR MATRICES, PURELY QUADRATIC
C       SEGMENTS ARE NOT CONSIDERED (NUMBER 13,14 AND 15, SO 3 PER ELEMENT)
C
        IF(IELM1.EQ.11.AND.IELM2.EQ.13) THEN
          NSEG2=NSEG2-3*NELEM
        ELSEIF(IELM1.EQ.13.AND.IELM2.EQ.11) THEN
          NSEG1=NSEG1-3*NELEM
        ENDIF
C
C       IN LINEAR-QUADRATIC RECTANGULAR MATRICES, PURELY QUADRATIC
C       SEGMENTS ARE NOT CONSIDERED (NUMBER 13,14 AND 15, SO 3 PER ELEMENT)
C
        CALL OMSEG(OP , M%D%R,TYPDIM,M%X%R,TYPEXM,
     *                  N%D%R,TYPDIN,N%X%R,TYPEXN,D%R,C,
     *                  NDIAGM,NSEG1,NSEG2,MESH%GLOSEG%I,
     *                  MESH%GLOSEG%MAXDIM1)
C
      ELSEIF(STOM.EQ.3.AND.STON.EQ.1) THEN
C
C       EDGE-BASED STORAGE FOR M AND EBE FOR N
C       THIS CAN HAPPEN ONLY WHEN N IS A BOUNDARY MATRIX   
C
        IF(  (M%ELMLIN.EQ.11.AND.M%ELMCOL.EQ.11.AND.
     *        N%ELMLIN.EQ.1 .AND.N%ELMCOL.EQ.1) .OR.
     *       (M%ELMLIN.EQ.12.AND.M%ELMCOL.EQ.12.AND.
     *        N%ELMLIN.EQ.1 .AND.N%ELMCOL.EQ.1) .OR.
     *       (M%ELMLIN.EQ.13.AND.M%ELMCOL.EQ.13.AND.
     *        N%ELMLIN.EQ.1 .AND.N%ELMCOL.EQ.1) .OR.
     *       (M%ELMLIN.EQ.13.AND.M%ELMCOL.EQ.13.AND.
     *        N%ELMLIN.EQ.2 .AND.N%ELMCOL.EQ.2)      ) THEN
C
          NSEG1 = NBSEG(M%ELMLIN)
          NSEG2 = NBSEG(M%ELMCOL)
C
          CALL OMSEGBOR(OP , M%D%R,TYPDIM,M%X%R,TYPEXM,
     *                       N%D%R,TYPDIN,N%X%R,TYPEXN,D%R,C,
     *                       NDIAGM,NSEG1,NSEG2,MESH%NBOR%I,
     *                       MESH%KP1BOR%I,NPTFR,
     *                       M%ELMLIN,N%ELMLIN)
C
        ELSE
          WRITE(LU,*) 'OM : UNEXPECTED CASE IN SEGMENT STORAGE'
          WRITE(LU,*) '     M%ELMLIN=',M%ELMLIN
          WRITE(LU,*) '     M%ELMCOL=',M%ELMCOL
          WRITE(LU,*) '     M%NAME=',M%NAME
          WRITE(LU,*) '     N%ELMLIN=',N%ELMLIN
          WRITE(LU,*) '     N%ELMCOL=',N%ELMCOL
          WRITE(LU,*) '     N%NAME=',N%NAME
          WRITE(LU,*) '     IMPLEMENTATION MISSING'
          CALL PLANTE(1)
          STOP 
        ENDIF
C
C  COMBINAISON DE IELM1 ET IELM2 NON PREVUE : ERREUR
C
C     ELSE
C        IF (LNG.EQ.1) WRITE(LU,100) M%NAME
C        IF (LNG.EQ.2) WRITE(LU,101) M%NAME
C        IF (LNG.EQ.1) WRITE(LU,200) IELM1,IELM2
C        IF (LNG.EQ.2) WRITE(LU,201) IELM1,IELM2
C        IF (LNG.EQ.1) WRITE(LU,410) STOM,STON
C        IF (LNG.EQ.2) WRITE(LU,411) STOM,STON
C        IF (LNG.EQ.1) WRITE(LU,300)
C        IF (LNG.EQ.2) WRITE(LU,301)
C        CALL PLANTE(1)
C        STOP
C     ENDIF
C
      ELSE
C
C  COMBINAISON DE STOCKAGES NON PREVUE
C
         IF (LNG.EQ.1) WRITE(LU,100) M%NAME
         IF (LNG.EQ.2) WRITE(LU,101) M%NAME
         IF (LNG.EQ.1) WRITE(LU,410) STOM,STON
         IF (LNG.EQ.2) WRITE(LU,411) STOM,STON
         IF (LNG.EQ.1) WRITE(LU,300)
         IF (LNG.EQ.2) WRITE(LU,301)
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  REENCODAGE DU NOUVEAU TYPE
C
      M%TYPDIA = TYPDIM
      M%TYPEXT = TYPEXM
      IF(OP(3:8).EQ.'X(M)  ') THEN
        M%X%DIM1=DIM1_EXT(IELM1,IELM2,STOM,'Q')
        M%X%DIM2=DIM2_EXT(IELM1,IELM2,STOM,'Q')
      ENDIF
C
C-----------------------------------------------------------------------
C
100      FORMAT(1X,'OM (BIEF) : MATRICE M (NOM REEL : ',A6,')')
150      FORMAT(1X,'OM (BIEF) : MATRICE N (NOM REEL : ',A6,')')
200      FORMAT(1X,'            IELM1 = ',1I6,' IELM2 = ',1I6)
250      FORMAT(1X,'            IELN1 = ',1I6,' IELN2 = ',1I6)
300      FORMAT(1X,'            CAS NON PREVU')
410      FORMAT(1X,'ET STOCKAGES   M  : ',1I6,'    N  : ',1I6)
C
101      FORMAT(1X,'OM (BIEF) : MATRIX  M (REAL NAME:',A6,')')
151      FORMAT(1X,'OM (BIEF) : MATRIX  N (REAL NAME:',A6,')')
201      FORMAT(1X,'            IELM1 = ',1I6,' IELM2 = ',1I6)
251      FORMAT(1X,'            IELN1 = ',1I6,' IELN2 = ',1I6)
301      FORMAT(1X,'            THIS CASE IS NOT IMPLEMENTED')
411      FORMAT(1X,'AND STORAGES   M  : ',1I6,' STON  : ',1I6)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
