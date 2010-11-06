C                       *****************
                        SUBROUTINE INBIEF
C                       *****************
C
     *(LIHBOR,KLOG,IT1,IT2,IT3,LVMAC,IELMX,
     * LAMBD0,SPHERI,MESH,T1,T2,OPTASS,PRODUC,EQUA)
C
C***********************************************************************
C BIEF VERSION 6.0      05/02/2010    J-M HERVOUET (LNHE) 01 30 87 80 18
C COPYRIGHT EDF 2008                  REGINA NEBAUER
C                                     LAM MINH PHUONG
C                                     EMILE RAZAFINDRAKOTO
C
C 22/01/2008 : ALLOCATION DYNAMIQUE DE IKLESTR SUPPRIMEE
C 29/02/2008 : NORMAB MODIFIEE, XSEG, YSEG NE SONT PLUS UTILISES
C 20/03/2008 : ADAPTATION DE NBOR, IKLBOR AUX TRIANGLES QUADRATIQUES
C 14/08/2008 : MODIFICATION DE PARINI
C 05/02/2010 : COMP_SEG UPDATED FOR QUADRATIC TRIANGLES
C***********************************************************************
C
C    FONCTION : PREPARATION DE LA STRUCTURE DE DONNEES POUR BIEF
C               LES TABLEAUX D'ENTIERS ET DE REELS DECRIVANT LE
C               MAILLAGE SONT CONSTRUITS ET RANGES DANS MESH.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    LIHBOR      | -->| TYPES DE CONDITIONS AUX LIMITES SUR H
C |    KLOG        | -->| CONVENTION POUR LA CONDITION LIMITE DE PAROI
C |    W1          | -->| TABLEAU DE TRAVAIL
C |    LVMAC       | -->| LONGUEUR DU VECTEUR SI MACHINE VECTORIELLE
C |    IELMX       | -->| TYPE DE L'ELEMENT LE PLUS COMPLEXE UTILISE
C |    LAMBD0      | -->| LATITUDE ORIGINE (COORDONNEES SPHERIQUES)
C |    SPHERI      | -->| LOGIQUE, SI OUI : COORDONNEES SPHERIQUES
C |    MESH        | -->| BLOCS DES TABLEAUX DU MAILLAGE
C |    T1,2,3      | -->| TABLEAUX DE TRAVAIL.
C |    OPTASS      | -->| OPTION D'ASSEMBLAGE.
C |    PRODUC      | -->| OPTION DU PRODUIT MATRICE X VECTEUR
C |    EQUA        | -->| IDENTIFICATION DU CODE OU DES EQUATIONS
C |________________|____|______________________________________________|
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PREDAT
C
C SOUS-PROGRAMME APPELE :
C
C***********************************************************************
C     
      USE BIEF, EX_INBIEF => INBIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)            :: IELMX,OPTASS,PRODUC,KLOG,LVMAC
      INTEGER, INTENT(IN)            :: LIHBOR(*)
      DOUBLE PRECISION, INTENT(IN)   :: LAMBD0
      LOGICAL, INTENT(IN)            :: SPHERI
      CHARACTER(LEN=20)              :: EQUA
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: T1,T2,IT1,IT2,IT3
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,IELEM,NELEM,NELMAX,NPTFR,NPOIN,IELM
      INTEGER MXPTVS,MXELVS,NPLAN
      INTEGER LV,NDP,IDP,I1,I2,I3,NPOIN2
      INTEGER NPTFR2,NELEM2,NELMAX2,NELEB2,NELEB    
C      
      DOUBLE PRECISION C,Z(1),X2,X3,Y2,Y3
C
C-----------------------------------------------------------------------
C     POUR APPEL A VOISIN31
      INTEGER IKLESTR(1,3)
C
C     DEPLOIEMENT DE LA STRUCTURE DE DONNEES
C
      NELEM = MESH%NELEM
      NELMAX= MESH%NELMAX
      NPOIN = MESH%NPOIN
      IELM  = MESH%X%ELM    
      NDP   = NBPEL(IELM)
      NPTFR = MESH%NPTFR
      NELEB = MESH%NELEB
C
C     WITH PRISMS, DIFFERENT FROM 2D VALUES, OTHERWISE
C
      IF(IELM.EQ.41.OR.IELM.EQ.51) THEN
        NPOIN2  =NBPTS(11)
        NELEM2  =NBPTS(10)
        NELMAX2 =NBMPTS(10)
        NPTFR2  =NBPTS(1)
        NPLAN   =NPOIN/NPOIN2
      ELSEIF(IELM.EQ.11.OR.IELM.EQ.31) THEN
        NPOIN2  =NPOIN
        NELEM2  =NELEM
        NELMAX2 =NELMAX
        NPTFR2  =NPTFR
        NELEB2  =NELEB
        NPLAN   =1
      ELSE
        WRITE(LU,*) 'UNEXPECTED ELEMENT IN INBIEF:',IELM
        CALL PLANTE(1)
        STOP
      ENDIF
C
C  PARALLELISME : 1) CONSTRUCTION DE KNOLG ET PASSAGE AUX NUMEROS
C                    GLOBAUX DU SOUS-DOMAINE.
C                 2) INITIALISATION DES TABLEAUX NHP,NHM
C                    INDPU,FAC, ETC.
C
C
      IF(NCSIZE.GT.1) THEN
C       1)
        CALL PARAGL(MESH%KNOGL%I,MESH%KNOGL%DIM1,MESH%KNOLG%I,
     *              MESH%NBOR%I ,MESH%NACHB%I,NPTFR2,NPOIN2)
C
C       2)
        CALL PARINI(MESH%NHP%I,MESH%NHM%I,MESH%INDPU%I,MESH%FAC,
     *              NPOIN2,MESH%NACHB%I,NPLAN,MESH,
     *              MESH%NB_NEIGHB,MESH%NB_NEIGHB_SEG,
     *              NELEM2,MESH%IFAPAR%I)
C
C       PRISMES : COMPLEMENT DE FAC
        IF(IELM.EQ.41.OR.IELM.EQ.51) THEN
          DO I = 2,NPLAN                        
            CALL OV_2('X=Y     ',MESH%FAC%R,I,MESH%FAC%R,1,
     *                           MESH%FAC%R,1,0.D0,NPOIN2,NPOIN2)
          ENDDO
        ENDIF       
C
      ELSE
C       THESE STUCTURES ARE ALLOCATED IN PARINI
        CALL ALLVEC(2,MESH%NB_NEIGHB_PT,'NBNGPT',0,1,0)
        CALL ALLVEC(2,MESH%LIST_SEND   ,'LSSEND',0,1,0)                
        CALL ALLVEC(2,MESH%NH_COM      ,'NH_COM',0,1,0)
        CALL ALLVEC(2,MESH%NB_NEIGHB_PT_SEG,'NBNGSG',0,1,0)
        CALL ALLVEC(2,MESH%LIST_SEND_SEG,'LSSESG',0,1,0)                
        CALL ALLVEC(2,MESH%NH_COM_SEG  ,'NH_CSG',0,1,0)
        CALL ALLVEC(1,MESH%BUF_SEND    ,'BUSEND',0,1,0)
        CALL ALLVEC(1,MESH%BUF_RECV    ,'BURECV',0,1,0)        
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C     CALCUL DES VOISINS DES FACES DE BORD POUR LE MAILLAGE DE TRIANGLES
C
C     NOTE : VOIR CPIKLE2 ET CPIKLE3, EN 3D IKLE PEUT ICI ETRE 3D CAR
C            LE DEBUT DE IKLE EN 3D EST LE MEME QUE EN 2D, CAR LES 3
C            PREMIERS POINTS DES PRISMES OU DES TETRAEDRES CORRESPONDENT 
C            AUX TROIS POINTS DES TRIANGLES DU FOND.
C
C
      IF(IELM.EQ.11.OR.IELM.EQ.41.OR.IELM.EQ.51) THEN
        CALL VOISIN(MESH%IFABOR%I,NELEM2,NELMAX2,IELM,MESH%IKLE%I,
     *              MESH%IKLE%DIM1,
     *              NPOIN2,MESH%NACHB%I,MESH%NBOR%I,NPTFR2,IT1%I,IT2%I)
C    
      ELSEIF (IELM.EQ.31) THEN
        CALL VOISIN31(MESH%IFABOR%I,NELEM2,NELMAX2,IELM,MESH%IKLE%I,
     *                MESH%IKLE%DIM1,
     *                NPOIN2,MESH%NACHB%I,MESH%NBOR%I,NPTFR2,
     *                LIHBOR,KLOG,IKLESTR,1,NELEB2)
      ELSE
        WRITE(LU,*) 'UNEXPECTED ELEMENT IN INBIEF:',IELM
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C  
      IF(IELM.EQ.11.OR.IELM.EQ.41.OR.IELM.EQ.51) THEN
C
      MXPTVS = MESH%MXPTVS
      MXELVS = MESH%MXELVS
      CALL ELEBD(MESH%NELBOR%I,MESH%NULONE%I,MESH%KP1BOR%I,
     *           MESH%IFABOR%I,MESH%NBOR%I,MESH%IKLE%I,MESH%IKLE%DIM1,
     *           MESH%IKLBOR%I,NELEM2,NELMAX2,
     *           NPOIN2,NPTFR2,IELM,
     *           LIHBOR,KLOG,MESH%IFANUM%I,
     *           OPTASS,MESH%ISEG%I,
     *           IT1%I,IT2%I,IT3%I)
C     
      ELSEIF(IELM.EQ.31) THEN
      CALL ELEBD31(MESH%NELBOR%I,MESH%NULONE%I,MESH%IKLBOR%I,
     &             MESH%IFABOR%I,MESH%NBOR%I,MESH%IKLE%I,
     &             NELEM2,NELEB2,NELMAX2,NPOIN2,NPTFR2,IELM)
           
      ELSE 
        WRITE(LU,*) 'UNEXPECTED ELEMENT IN INBIEF:',IELM
        CALL PLANTE(1)
        STOP
      ENDIF
C
C     COMPLETION OF ARRAYS FOR PRISMS
C
      IF(IELM.EQ.41) THEN
        CALL ELEB3D(MESH%IKLE%I,MESH%NBOR%I,MESH%KP1BOR%I,                         
     *              MESH%NELBOR%I,MESH%IKLBOR%I,           
     *              MESH%NULONE%I,NELEM2,NPOIN2,NPLAN,NPLAN-1,NPTFR2)
C
C     COMPLETION OF ARRAYS FOR TETRAHEDRONS
C
      ELSEIF(IELM.EQ.51) THEN
        CALL ELEB3DT(MESH%IKLE%I,MESH%NBOR%I,MESH%KP1BOR%I,                         
     *               MESH%NELBOR%I,MESH%IKLBOR%I,           
     *               MESH%NULONE%I,NELEM2,NELMAX2,
     *               NPOIN2,NPLAN,NPLAN-1,NPTFR2)
     
      ELSEIF(IELM.NE.11.AND.IELM.NE.31) THEN
        WRITE(LU,*) 'INBIEF UNEXPECTED ELEMENT: ',IELM
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C RECHERCHE DES POSSIBILITES DE VECTORISATION
C
      IF(IELM.EQ.11) THEN
C
      IF(LVMAC.NE.1) THEN
        IF(LNG.EQ.1) WRITE(LU,200) LVMAC
        IF(LNG.EQ.2) WRITE(LU,201) LVMAC
200     FORMAT(1X,'INBIEF (BIEF) : MACHINE VECTORIELLE',/,1X,
     *  'AVEC LONGUEUR DE VECTEUR :',1I6,
     *  ' (SELON VOS DONNEES OU DANS LE DICTIONNAIRE DES MOTS-CLES)')
201     FORMAT(1X,'INBIEF (BIEF): VECTOR MACHINE',/,1X,
     *  'WITH VECTOR LENGTH :',1I6,
     *  ' (ACCORDING TO YOUR DATA OR IN THE DICTIONNARY OF KEY-WORDS)')
        CALL VECLEN(LV,NDP,MESH%IKLE%I,NELEM,NELMAX,NPOIN,T1%R)
        IF(LV.LT.LVMAC) THEN
          IF(LNG.EQ.1) WRITE(LU,300) LV
          IF(LNG.EQ.2) WRITE(LU,301) LV
300       FORMAT(1X,'LONGUEUR LIMITEE A ',1I4,' PAR LA NUMEROTATION DES
     *ELEMENTS (VOIR LA DOCUMENTATION DE STBTEL)')
301       FORMAT(1X,'THIS LENGTH IS REDUCED TO ',1I4,' BY THE NUMBERING
     *OF THE ELEMENTS (SEE STBTEL DOCUMENTATION)')
        ENDIF
      ELSE
        LV = 1
        IF(LNG.EQ.1) WRITE(LU,400)
        IF(LNG.EQ.2) WRITE(LU,401)
400     FORMAT(1X,'INBIEF (BIEF) : MACHINE NON VECTORIELLE',
     *                                           ' (SELON VOS DONNEES)')
401     FORMAT(1X,'INBIEF (BIEF): NOT A VECTOR MACHINE',
     *                                      ' (ACCORDING TO YOUR DATA)')
      ENDIF
C
      MESH%LV = LV
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(SPHERI.AND.IELM.EQ.11) THEN
C
        CALL LATITU(MESH%COSLAT%R,MESH%SINLAT%R,LAMBD0,MESH%Y%R,NPOIN)
        CALL CORLAT
C
C VERSION 5.5
        CALL CPSTVC(MESH%X,T1)
        CALL CPSTVC(MESH%Y,T2)
        DO I=1,NPOIN
          T1%R(I)=MESH%X%R(I)*MESH%COSLAT%R(I)
          T2%R(I)=MESH%Y%R(I)*MESH%COSLAT%R(I)
        ENDDO
C
C VERSION 5.4 (IN PARAMETER ESTIMATION COSLAT MAY HAVE BEEN CHANGED
C              INTO QUASI-BUBBLE BY A PREVIOUS RUN, AND OS STOPS)
C       CALL OS( 'X=YZ    ' , T1 , MESH%X , MESH%COSLAT , C )
C       CALL OS( 'X=YZ    ' , T2 , MESH%Y , MESH%COSLAT , C )
C
      ELSE
C
        CALL OS( 'X=Y     ' , X=T1 , Y=MESH%X )
        CALL OS( 'X=Y     ' , X=T2 , Y=MESH%Y )
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  PASSAGE AUX COORDONNEES PAR ELEMENTS (SEULEMENT EN TRIANGLES)
C
      IF(IELM.EQ.11) THEN
C
        CALL PTTOEL(MESH%XEL,T1,MESH)
        CALL PTTOEL(MESH%YEL,T2,MESH)
C
C  PASSAGE A UN REPERE LOCAL EN X ET Y
C
        DO 10 IDP=2,NDP
        CALL OV_2('X=X-Y   ',MESH%XEL%R,IDP,
     *                       MESH%XEL%R,1  ,
     *                       MESH%XEL%R,1  , C , NELMAX , NELEM )
        CALL OV_2('X=X-Y   ',MESH%YEL%R,IDP,
     *                       MESH%YEL%R,1  ,
     *                       MESH%YEL%R,1  , C , NELMAX , NELEM )
10      CONTINUE
C
        CALL OV('X=C     ', MESH%XEL%R , Z , Z , 0.D0 , NELEM )
        CALL OV('X=C     ', MESH%YEL%R , Z , Z , 0.D0 , NELEM )
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C CALCUL DE COEFFICIENTS GEOMETRIQUES POUR CHAQUE ELEMENT
C
      IF(IELM.EQ.11) THEN
C
        CALL GEOELT(MESH%SURDET%R,MESH%SURFAC%R,
     *              MESH%XEL%R   ,MESH%YEL%R   ,NELEM,NELMAX,IELM)
C
C POUR L'INSTANT, SURDET N'EST UTILISE QUE PAR CARACT, QUI NE TRAVAILLE
C PAS SUR LE MAILLAGE TRANSFORME EN COORDONNEES SPHERIQUES.
C ICI ON EFFACE DONC SURDET CALCULE PAR GEOELE AVEC XEL ET YEL.
C
        IF(SPHERI) THEN
C
         DO IELEM = 1 , NELEM
C
         I1 = MESH%IKLE%I(IELEM)
         I2 = MESH%IKLE%I(IELEM+NELMAX)
         I3 = MESH%IKLE%I(IELEM+2*NELMAX)
         X2 = - MESH%X%R(I1) + MESH%X%R(I2)
         X3 = - MESH%X%R(I1) + MESH%X%R(I3)
         Y2 = - MESH%Y%R(I1) + MESH%Y%R(I2)
         Y3 = - MESH%Y%R(I1) + MESH%Y%R(I3)
C
         MESH%SURDET%R(IELEM) = 1.D0 / (X2*Y3 - X3*Y2)
C
         ENDDO
C
        ENDIF
C
      ELSEIF(IELM.EQ.41.OR.IELM.EQ.51.OR.IELM.EQ.31) THEN
C
C        POUR LES PRISMES, SURFAC EST LA SURFACE DES TRIANGLES
C
         DO IELEM = 1 , NELEM
C
         I1 = MESH%IKLE%I(IELEM)
         I2 = MESH%IKLE%I(IELEM+NELMAX)
         I3 = MESH%IKLE%I(IELEM+2*NELMAX)
         X2 = - MESH%X%R(I1) + MESH%X%R(I2)
         X3 = - MESH%X%R(I1) + MESH%X%R(I3)
         Y2 = - MESH%Y%R(I1) + MESH%Y%R(I2)
         Y3 = - MESH%Y%R(I1) + MESH%Y%R(I3)
C
         MESH%SURFAC%R(IELEM) = 0.5D0 * (X2*Y3 - X3*Y2)
C
         ENDDO
      ELSE
        WRITE(LU,*) 'UNEXPECTED ELEMENT IN INBIEF:',IELM
      ENDIF
C
C-----------------------------------------------------------------------
C
C DEFINITION DES NORMALES EXTERIEURES SUR LES BORDS.
C            ET DES DISTANCES AU BORD.
C
      IF(IELM.EQ.11) THEN
C
      CALL NORMAB(MESH%XNEBOR%R,MESH%YNEBOR%R,
     *            MESH%XSGBOR%R,MESH%YSGBOR%R,
     *            MESH%DISBOR%R,MESH%SURFAC%R,NELMAX,
     *            MESH%NBOR%I,MESH%KP1BOR%I,MESH%NELBOR%I,
     *            MESH%LGSEG%R,NPTFR,
     *            MESH%X%R,MESH%Y%R,MESH,T1)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  DATA STRUCTURE FOR EDGE-BASED STORAGE (FROM 5.9 ON ALWAYS DONE IN 2D)
C  SEE CALL TO COMP_SEG BELOW FOR COMPLETING THE STRUCTURE
C
      IF(IELM.EQ.11) THEN
C
      CALL STOSEG(MESH%IFABOR%I,NELEM,NELMAX,NELMAX,IELMX,MESH%IKLE%I,
     *            MESH%NBOR%I,NPTFR,
     *            MESH%GLOSEG%I,MESH%GLOSEG%MAXDIM1,
     *            MESH%ELTSEG%I,MESH%ORISEG%I,MESH%NSEG,
     *            MESH%KP1BOR%I,MESH%NELBOR%I,MESH%NULONE%I,
     *            MESH%KNOLG%I)
C
      ELSEIF(IELM.EQ.41.AND.OPTASS.EQ.3) THEN
C
      CALL STOSEG41(MESH%IFABOR%I,NELEM,NELMAX,IELMX,MESH%IKLE%I,
     *              MESH%NBOR%I,NPTFR,
     *              MESH%GLOSEG%I,MESH%GLOSEG%MAXDIM1,
     *              MESH%ELTSEG%I,MESH%ORISEG%I,MESH%NSEG,
     *              MESH%KP1BOR%I,MESH%NELBOR%I,MESH%NULONE%I,
     *              NELMAX2,NELEM2,NPTFR2,NPOIN2,NPLAN,MESH%KNOLG%I)    
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(NCSIZE.GT.1.AND.IELM.EQ.11) THEN
C
C       COMPLETING NH_COM_SEG WITH SEGMENT NUMBERS ONCE ELTSEG IS KNOWN
C
        CALL COMP_NH_COM_SEG(MESH%ELTSEG%I,NELEM,MESH%NH_COM_SEG%I,
     *                       MESH%NH_COM_SEG%DIM1,MESH%NB_NEIGHB_SEG,
     *                       MESH%NB_NEIGHB_PT_SEG%I,
     *                       MESH%GLOSEG%I,MESH%GLOSEG%DIM1,
     *                       MESH%KNOLG%I,NPOIN)
C
C       COMPLETING FAC ONCE IFABOR AND ELTSEG ARE KNOWN
C
        IF(IELM.EQ.11.AND.IELMX.EQ.13) THEN
          CALL COMP_FAC(MESH%ELTSEG%I,MESH%IFABOR%I,NELEM,
     *                  NPOIN,MESH%FAC)
        ENDIF 
C       
      ENDIF 
C
C-----------------------------------------------------------------------
C
C  DATA STRUCTURE FOR EDGE-BASED STORAGE
C
      IF(IELM.EQ.11.AND.PRODUC.EQ.2) THEN
C
      CALL FROPRO(MESH%NBOR%I,MESH%IKLE%I,
     *            NELEM,NELMAX,NPOIN,MESH%NPMAX,NPTFR,IELM,
     *            MESH%IKLEM1%I,MESH%LIMVOI%I,OPTASS,PRODUC,MXPTVS,
     *            IT1%I,MESH%GLOSEG%I,MESH%GLOSEG%DIM1,MESH%NSEG)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  COMPLEMENT DE IKLE AU DELA DES ELEMENTS LINEAIRES
C
      IF(IELM.EQ.11.AND.IELM.NE.IELMX) THEN
        IF(MESH%IKLE%DIM2.NE.NBPEL(IELMX)) THEN
          IF(LNG.EQ.1) WRITE(LU,100) IELMX
          IF(LNG.EQ.2) WRITE(LU,101) IELMX
100       FORMAT(1X,'INBIEF (BIEF) : IKLE MAL DIMENSIONNE',/,1X,
     *              'POUR UN ELEMENT DE TYPE :',1I6)
101       FORMAT(1X,'INBIEF (BIEF): WRONG DIMENSION OF IKLE',/,1X,
     *              'FOR AN ELEMENT WITH TYPE :',1I6)
          CALL PLANTE(1)
          STOP
        ENDIF
        CALL COMP_IKLE(MESH%IKLE%I,MESH%IKLBOR%I,
     *                 MESH%ELTSEG%I,MESH%NBOR%I,
     *                 IELMX,NELEM,NELMAX,NPOIN,NPTFR)
      ENDIF
C
C-----------------------------------------------------------------------
C
C COMPLEMENT DE LA STRUCTURE DE SEGMENT AU DELA DES ELEMENTS LINEAIRES
C
      IF(IELM.NE.IELMX) THEN
        CALL COMP_SEG(NELEM,NELMAX,IELMX,MESH%IKLE%I,MESH%GLOSEG%I,
     *                MESH%GLOSEG%MAXDIM1,MESH%ELTSEG%I,MESH%ORISEG%I,
     *                MESH%NSEG)
      ENDIF
C
C-----------------------------------------------------------------------
C
C COMPLEMENT DE LA STRUCTURE DE DONNEES POUR LES VOLUMES FINIS
C
      IF(EQUA(1:15).EQ.'SAINT-VENANT VF') THEN
C
          CALL INFCEL(MESH%X%R,MESH%Y%R,MESH%IKLE%I,
     *                MESH%NUBO%I,MESH%VNOIN%R,
     *                NPOIN,MXPTVS,NELEM,NELMAX,MESH%NSEG,MESH%CMI%R,
     *                MESH%JMI%I,MESH%AIRST%R)
C                                                                      
C         CALCUL DES AIRES DES CELLULES                                     
C                                                                       
          CALL VECTOR(T1,'=','MASBAS          ',11,      
     *                1.D0,T2,T2,T2,T2,T2,T2,MESH,.FALSE.,T2)               
C                                                                    
C         CALCUL DU PAS D'ESPACE LOCAL PAR CELLULE
C
          CALL HLOC(NPOIN,MESH%NSEG,MESH%NPTFR,MESH%NUBO%I,
     *              MESH%NBOR%I,MESH%VNOIN%R,
     *              MESH%XNEBOR%R,MESH%YNEBOR%R,T1%R,MESH%DTHAUT%R)
C
C         CALCUL DES GRADIENTS DES FONCTIONS DE BASE
C
          CALL GRADP(NPOIN,MESH%NELMAX,MESH%IKLE%I,MESH%SURFAC%R,
     *               MESH%X%R,MESH%Y%R,MESH%DPX%R,MESH%DPY%R)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN      
      END
