C                       *****************
                        SUBROUTINE ALMESH
C                       *****************
C
     *(MESH,NOM,IELM,SPHERI,CFG,NFIC,EQUA,NPLAN,NPMAX,NPTFRX,NELMAX,
     * I3,I4)
C
C***********************************************************************
C BIEF VERSION 6.0      05/02/2010    J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C  FONCTION  : ALLOCATION D'UNE STRUCTURE DE MAILLAGE BIEF_MESH.
C
C
C 05/02/2010 : EDGE-BASED STUCTURES ALWAYS ALLOCATED.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   MESH         | -->| STRUCTURE A ALLOUER
C |   NOM          | -->| NAME OF THE MESH
C |   IELM         | -->| ELEMENT QUI CONTIENT LE PLUS GRAND NOMBRE
C |                |    | DE POINTS QUI SERA UTILISE.
C |   SPHERI       | -->| LOGIQUE, SI OUI : COORDONNEES SPHERIQUES
C |   CFG(1)       | -->| OPTION DE STOCKAGE  1 : EBE CLASSIQUE
C |                |    |                     2 : EBE ASSEMBLE
C |   CFG(2)       | -->| MATRICE X VECTEUR   1 : EBE CLASSIQUE
C |                |    |                     2 : FRONTAL
C |   MXPTVS       | -->| NOMBRE MAXIMUM DE POINTS VOISINS D'UN POINT.
C |   MXELVS       | -->| NOMBRE MAXIMUM D'ELEMENTS VOISINS D'UN POINT.
C |   EQUA         | -->| NAME IN 20 CHARACTERS TO ENABLE DIFFERENT
C |                |    | OPTIONS. OPTIONS ARE:
C |                |    |     "SAINT-VENANT EF"
C |                |    |     "SAINT-VENANT VF"
C |                |    |     "BOUSSINESQ"
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_ALMESH => ALMESH
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_MESH)  , INTENT(INOUT)           :: MESH
      INTEGER          , INTENT(IN)              :: IELM
      INTEGER          , INTENT(IN)              :: NFIC
      LOGICAL          , INTENT(IN)              :: SPHERI
      CHARACTER(LEN=6) , INTENT(IN)              :: NOM
      CHARACTER(LEN=20), INTENT(IN)              :: EQUA
      INTEGER          , INTENT(INOUT)           :: CFG(2)
      INTEGER          , INTENT(IN),    OPTIONAL :: NPLAN
      INTEGER          , INTENT(IN),    OPTIONAL :: NPMAX
      INTEGER          , INTENT(IN),    OPTIONAL :: NPTFRX
      INTEGER          , INTENT(IN),    OPTIONAL :: NELMAX
      INTEGER          , INTENT(INOUT), OPTIONAL :: I3,I4
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER D,IELM0,IELM1,STOCFG,IELB0,IELB1,NSEG,NNPMAX
      INTEGER NNPTFRX,NNELEB,ERR,NNELMAX,NNPLAN,NPOIN,NPTFR,NELEM
      INTEGER MXPTVS,MXELVS,NDP,IB(10),IELEM,NSEGBOR
C
C     Temporary connectivity table 
C
      INTEGER, ALLOCATABLE :: IKLES(:)
C
C     Temporary table for numbering the boundary nodes
C
      INTEGER, ALLOCATABLE :: IPOBO(:)
!
      INTEGER IELB0V,IELB1V
C
      INTEGER I
C
C     FH-jaj
C     For size of KNOGL
      INTEGER :: NPOIN_MAX
      INTEGER, EXTERNAL :: P_ISUM
C
C-----------------------------------------------------------------------
C
C     FIRST READING OF THE GEOMETRY FILE TO GET NPOIN,.. IB
C
      MESH%NAME = NOM
C
CCCCCCCCCCCCCCCCCCCCC
C 3D : en plus, on a NNELEB en sortie.
CCCCCCCCCCCCCCCCCCCCC
C
C     IN PARALLEL, THIS IS WHERE NPTIR IS READ
C
      CALL READGEO1(NPOIN,NELEM,NPTFR,NDP,IB,NFIC,NNELEB)
C
      IF(PRESENT(I3)) I3=IB(3)
      IF(PRESENT(I4)) I4=IB(4)
C
C Allocation of the temporary connectivity table
C
      ALLOCATE(IKLES(NELEM*NDP),STAT=ERR)
C
      IF(ERR.NE.0) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'ALMESH : ALLOCATION DE IKLES DEFECTUEUSE'
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'ALMESH : WRONG ALLOCATION OF IKLES'
        ENDIF
        STOP
      ENDIF
C
C Allocation of the temporary table for the boundary nodes.
C
      ALLOCATE(IPOBO(NPOIN),STAT=ERR)
C
      IF(ERR.NE.0) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'ALMESH : ALLOCATION DE IPOBO DEFECTUEUSE'
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'ALMESH : WRONG ALLOCATION OF IPOBO'
        ENDIF
        STOP
      ENDIF
C
CCCCCCCCCCCCCCCCCCCCCCC
C READGEO2 prendra en plus l'argument IPOBO, tableau lu dans cette
C procedure. Actuellement il est utilise uniquement en interne de 
C cette procedure, mais la, en on a besoin plus tard.
C comme ca fait partie du fichier .geo, autant pouvoir l'utiliser.
C Les arguments MXPTVS et MXELVS seront supprimes et transferes dans 
C l'appel MXPTEL qui suit. Le calcul de ces nombres n'est pas en rapport
C avec la lecture de donnees dans le fichier .geo.
CCCCCCCCCCCCCCCCCCCCCCC
C
C Reading IPOBO
C
      CALL READGEO2(NPOIN,NELEM,NPTFR,NDP,IKLES,IPOBO,IB,NFIC)
C
CCCCCCCCCCCCCCCCCCCCCCC
C ICI on va mettre l'appel a mxptel. (au lieu de l'appel dans readgeo2.
C On calculera MXPTVS et MXELVS en fonction du type de l'element.
CCCCCCCCCCCCCCCCCCCCCCC
C
C Calculating the maximal number of elements around a node MXELVS
C and the maximum number of surrounding nodes, MXPTVS 
C
      CALL MXPTEL(MXPTVS,MXELVS,IKLES,IELM,
     *              NPOIN,NELEM,NDP,IPOBO,.TRUE.)
C     
      DEALLOCATE(IPOBO)
C
C
C-----------------------------------------------------------------------
C      
C
C     INITIALISATION OF COMMONS DIMS AND NODES
C
      IF(PRESENT(NPMAX)) THEN
        NNPMAX = NPMAX
      ELSE
        NNPMAX = NPOIN
      ENDIF
      IF(PRESENT(NPTFRX)) THEN
        NNPTFRX = NPTFRX
      ELSE
        NNPTFRX = NPTFR
      ENDIF 
      IF(PRESENT(NELMAX)) THEN
        NNELMAX = NELMAX
      ELSE
        NNELMAX = NELEM
      ENDIF 
      IF(PRESENT(NPLAN)) THEN
        NNPLAN = NPLAN
      ELSE
        NNPLAN = 1
      ENDIF
C
      IF(NCSIZE.GT.1) THEN
C
C       IN // NSEGBOR IS NOT NPTFR, NSEGBOR COMPUTED IN SEGBOR      
C
        CALL SEGBOR(NSEGBOR,IKLES,NELEM,NELMAX,NPOIN)
      ELSE
        NSEGBOR=NPTFR 
      ENDIF
C
C     IN CALL ININDS, ALL VALUES ARE 2D VALUES, ONLY NNPLAN TELLS
C     IF IT IS 2D OR 3D.
C
CC MAILLAGE 3D TETRAEDRES
CC ADD the optional argument NNELEB 
C      
      CALL ININDS(NPOIN,NPTFR,NELEM,NNPMAX,NNPTFRX,NNELMAX,NNPLAN,
     *            NSEGBOR,NNELEB)
C
C     TYPE OF ELEMENTS P0 AND P1
C
      IELM0  = 10*(IELM/10)
      IELM1  = IELM0 + 1
      D      = DIMENS(IELM0)
C
C BOUNDARY ELEMENTS (at the surface and at the bottom FOR PRISMS)
C
      IELB0  = IELBOR(IELM0,1)
      IELB1  = IELBOR(IELM1,1)
C
C LATERAL BOUNDARY ELEMENTS (DIFFERENT FOR PRISMS)
C
      IELB0V = IELBOR(IELM0,2)
      IELB1V = IELBOR(IELM1,2)
C
C-----------------------------------------------------------------------
C
C  ALLOCATION DES TABLEAUX DE REELS
C
C     COORDONNEES PAR ELEMENTS : XEL,YEL,ZEL
C
      ALLOCATE(MESH%XEL)
      ALLOCATE(MESH%YEL)
      ALLOCATE(MESH%ZEL)
C
      CALL ALLVEC(1,MESH%XEL,'XEL   ',IELM0,NBPEL(IELM1),1)
      CALL ALLVEC(1,MESH%YEL,'YEL   ',IELM0,NBPEL(IELM1),1)
C
      IF(D.GE.3) THEN
        CALL ALLVEC(1,MESH%ZEL,'ZEL   ',IELM0,NBPEL(IELM1),1)
      ELSE
        CALL ALLVEC(1,MESH%ZEL,'ZEL   ',    0,           1,0)
      ENDIF
C
C     SURFACES DES ELEMENTS : SURFAC
C     jaj can be used for element volumes...
C
      ALLOCATE(MESH%SURFAC)
      CALL ALLVEC(1,MESH%SURFAC,'SURFAC',IELM0,1,1)
C
C     1/DET : SURDET ! not used in 3D, why?
C
      ALLOCATE(MESH%SURDET)
      CALL ALLVEC(1,MESH%SURDET,'SURDET',IELM0,1,1)
C
C     LONGUEURS DES SEGMENTS : LGSEG
! can be used (in theory) for lateral surfaces in 3D, 
! but then ielb0v instead of ielb0! (2D case not affected)
C
      ALLOCATE(MESH%LGSEG)
      CALL ALLVEC(1,MESH%LGSEG,'LGSEG ',ielb0v,1,1)
C
C     NORMALES AUX SEGMENTS : XSGBOR,YSGBOR,ZSGBOR
! can be (in theory) used for "non-sigma" mesh for lateral normal vectors 
! per lateral boundary element, but then ielb0v instead of ielb0!
! 2D case not affected
C
      ALLOCATE(MESH%XSGBOR)
      ALLOCATE(MESH%YSGBOR)
      ALLOCATE(MESH%ZSGBOR)
C     SEE NORMAB FOR MEANING OF 4 DIMENSIONS
      CALL ALLVEC(1,MESH%XSGBOR,'XSGBOR',IELB0V,4,1)
      CALL ALLVEC(1,MESH%YSGBOR,'YSGBOR',IELB0V,4,1)
      IF(D.GE.3) THEN
        CALL ALLVEC(1,MESH%ZSGBOR,'ZSGBOR',IELB0V,4,1)
      ELSE
        CALL ALLVEC(1,MESH%ZSGBOR,'ZSGBOR',    0,4,0)
      ENDIF
C
C     NORMALES AUX NOEUDS : XNEBOR,YNEBOR,ZNEBOR
!
! in 3D they are normal vectors at the bottom
! so that ielb1 remains
!
      ALLOCATE(MESH%XNEBOR)
      ALLOCATE(MESH%YNEBOR)
      ALLOCATE(MESH%ZNEBOR)
!
      CALL ALLVEC(1,MESH%XNEBOR,'XNEBOR',IELB1,2,1)
      CALL ALLVEC(1,MESH%YNEBOR,'YNEBOR',IELB1,2,1)
!
      IF(D.GE.3) THEN !jaj not used, actually 
        CALL ALLVEC(1,MESH%ZNEBOR,'ZNEBOR',IELB1,2,1)
      ELSE
        CALL ALLVEC(1,MESH%ZNEBOR,'ZNEBOR',    0,2,0)
      ENDIF
C
C     COORDONNEES PAR POINTS : X, Y ET Z
C
      ALLOCATE(MESH%X)
      ALLOCATE(MESH%Y)
      ALLOCATE(MESH%Z)
      CALL ALLVEC(1,MESH%X,'X     ',IELM1,1,1)
      CALL ALLVEC(1,MESH%Y,'Y     ',IELM1,1,1)
      IF(D.GE.3) THEN
        CALL ALLVEC(1,MESH%Z,'Z     ',IELM1,1,1)
      ELSE
        CALL ALLVEC(1,MESH%Z,'Z     ',    0,1,0)
      ENDIF
C
C     COS ET SIN DE LA LATITUDE
C     ON PREND IELM (EXEMPLE : VITESSE DANS CORIOLIS)
C     COSLAT ET SINLAT SONT DES TABLEAUX DE TRAVAIL
C     AUXQUELS ON DONNE AU DEPART LA STRUCTURE DE X
C     ON PEUT LES ETENDRE PLUS TARD A L'ELEMENT IELM.
C
      ALLOCATE(MESH%COSLAT)
      ALLOCATE(MESH%SINLAT)
      IF(SPHERI) THEN ! different compared to the v2.3
        CALL ALLVEC(1,MESH%COSLAT,'COSLAT',IELM,1,2)
        CALL ALLVEC(1,MESH%SINLAT,'SINLAT',IELM,1,2)
        CALL CPSTVC(MESH%X,MESH%COSLAT)
        CALL CPSTVC(MESH%X,MESH%SINLAT)
      ELSE
        CALL ALLVEC(1,MESH%COSLAT,'COSLAT',   0,1,0)
        CALL ALLVEC(1,MESH%SINLAT,'SINLAT',   0,1,0)
      ENDIF
C
C     DISTANCES TO BOUNDARIES : DISBOR
C
      ALLOCATE(MESH%DISBOR)
      CALL ALLVEC(1,MESH%DISBOR,'DISBOR',IELB0,1,1)
C
C     MATRICE DE TRAVAIL INTERNE A BIEF, AVEC STOCKAGE CLASSIQUE
C
      STOCFG = CFG(1)
      CFG(1) = 1
      ALLOCATE(MESH%M)
      CALL ALLMAT(MESH%M,'M     ',IELM,IELM,CFG,'Q','Q')
      CFG(1) = STOCFG
C
C     MATRICE DE TRAVAIL PAR SEGMENT
C
      ALLOCATE(MESH%MSEG)
C     FROM 5.9 ON ALWAYS DONE IN 2D (IELM=11,12,13 OR 14)
      IF(CFG(1).EQ.3.OR.10*(IELM/10).EQ.10) THEN
        CALL ALLMAT(MESH%MSEG,'MSEG  ',IELM,IELM,CFG,'Q','Q')
      ELSE
        CALL ALLMAT(MESH%MSEG,'MSEG  ',IELM,IELM,CFG,'0','0')
      ENDIF
C
C     TABLEAU DE TRAVAIL POUR UN VECTEUR NON ASSEMBLE
C
      ALLOCATE(MESH%W)
      CALL ALLVEC(1,MESH%W,'W     ',IELM0,NBPEL(IELM),2)
C
C     TABLEAU DE TRAVAIL POUR UN VECTEUR NORMAL.
C
      ALLOCATE(MESH%T)
      CALL ALLVEC(1,MESH%T,'T     ',IELM,1,2)
C
C     VNOIN : TABLEAU DES NORMALES VNOIN POUR VOLUMES FINIS.
C     CMI : COORDONNEES MILIEU DE SEGMENTS POUR SCHEMAS CINETIQUES.
C     DTHAUT :
C     DPX,DPY : GRADIENTS DES FONCTIONS DE BASE. 
C
      ALLOCATE(MESH%VNOIN)
      ALLOCATE(MESH%CMI)
      ALLOCATE(MESH%AIRST)
      ALLOCATE(MESH%DTHAUT)
      ALLOCATE(MESH%DPX)
      ALLOCATE(MESH%DPY)
      IF(EQUA(1:15).EQ.'SAINT-VENANT VF') THEN
        NSEG=NBSEG(IELM1)
        CALL ALLVEC(1,MESH%VNOIN ,'VNOIN ',3*NSEG,1,0)
        CALL ALLVEC(1,MESH%CMI   ,'CMI   ',2*NSEG,1,0)
        CALL ALLVEC(1,MESH%AIRST ,'AIRST ',2*NSEG,1,0)
        CALL ALLVEC(1,MESH%DTHAUT,'DTHAUT',IELM1,1,2)
        CALL ALLVEC(1,MESH%DPX   ,'DPX   ',IELM0,3,2)
        CALL ALLVEC(1,MESH%DPY   ,'DPY   ',IELM0,3,2)
      ELSE
        CALL ALLVEC(1,MESH%VNOIN ,'VNOIN ',     0,1,0)
        CALL ALLVEC(1,MESH%CMI   ,'CMI   ',     0,1,0)
        CALL ALLVEC(1,MESH%AIRST ,'AIRST ',     0,1,0)
        CALL ALLVEC(1,MESH%DTHAUT,'DTHAUT',     0,1,0)
        CALL ALLVEC(1,MESH%DPX   ,'DPX   ',     0,1,0)
        CALL ALLVEC(1,MESH%DPY   ,'DPY   ',     0,1,0)
      ENDIF
C
C     FOR PARALLELISM
C
      ALLOCATE(MESH%XSEG)
      ALLOCATE(MESH%YSEG)
      ALLOCATE(MESH%FAC)
C     THEIR ALLVEC IS IN PARINI
      ALLOCATE(MESH%BUF_SEND)
      ALLOCATE(MESH%BUF_RECV)
C
      IF(NCSIZE.GT.1) THEN
C
C       XSEG
        CALL ALLVEC(1,MESH%XSEG,'XSEG  ',ielbor(ielm1,2),1,2)
C       YSEG
        CALL ALLVEC(1,MESH%YSEG,'YSEG  ',ielbor(ielm1,2),1,2)
C       FAC
        CALL ALLVEC(1,MESH%FAC,'FAC   ',IELM ,1,2)
C
      ELSE
        CALL ALLVEC(1,MESH%XSEG  ,'XSEG  ',0,1,0)
        CALL ALLVEC(1,MESH%YSEG  ,'YSEG  ',0,1,0)
        CALL ALLVEC(1,MESH%FAC   ,'FAC   ',0,1,0)
      ENDIF
C
C-----------------------------------------------------------------------
C
C     1) INTEGER VALUES (ALLOCATE BECAUSE THEY ARE POINTERS)
C
      ALLOCATE(MESH%NELEM)
      MESH%NELEM  = NBPTS(IELM0)
      ALLOCATE(MESH%NELMAX)
      MESH%NELMAX = NBMPTS(IELM0)
!
!
! I do use MESH%NPTFR for the number of lateral boundary NODES
! ielbor(ielm0,1) changed to ielbor(ielm1,1)
! the problem is, that for 3D (ielm=41):
! nbpts(ielbor(ielm0,1)) contains the number of HORIZONTAL boundary elements
! nbpts(ielbor(ielm0,2)) contains the number of VERTICAL boundary elements
! nbpts(ielbor(ielm1,1)) contains the number of HORIZONTAL boundary nodes
! nbpts(ielbor(ielm1,2)) contains the number of VERTICAL boundary nodes
!
! Funny, but the 2D case is NOT affected, because the number of boundary
! segments is equal to the number of boundary nodes. 
!
      ALLOCATE(MESH%NPTFR)
      MESH%NPTFR  = nbpts(ielbor(ielm1,2))
      ALLOCATE(MESH%NPTFRX)
      MESH%NPTFRX = nbmpts(ielbor(ielm1,2))
!
! number of lateral boundary elements
!
      ALLOCATE(MESH%NELEB)
      ALLOCATE(MESH%NELEBX)
CC MAILLAGE-3D
      IF(IELM.EQ.31) THEN
        MESH%NELEB   = NNELEB
        MESH%NELEBX  = NNELEB
      ELSEIF(IELM.EQ.11.OR.IELM.EQ.12.OR.IELM.EQ.13.OR.IELM.EQ.14) THEN
        MESH%NELEB   = NPTFR
        MESH%NELEBX  = NNPTFRX
      ELSEIF(IELM.EQ.41) THEN
        MESH%NELEB   = NPTFR*(NNPLAN-1)
        MESH%NELEBX  = NNPTFRX*(NNPLAN-1)
      ELSEIF(IELM.EQ.51) THEN
        MESH%NELEB   = 2*NPTFR*(NNPLAN-1)
        MESH%NELEBX  = 2*NNPTFRX*(NNPLAN-1)
      ELSE
        WRITE(LU,*) 'ALMESH, UNEXPECTED ELEMENT FOR NELEB:',IELM
        CALL PLANTE(1)
        STOP
      ENDIF
!
      ALLOCATE(MESH%DIM)
      MESH%DIM    = DIMENS(IELM0)
      ALLOCATE(MESH%TYPELM)
      MESH%TYPELM = IELM0
      ALLOCATE(MESH%NPOIN)
      MESH%NPOIN  = NBPTS(IELM1)
      ALLOCATE(MESH%NPMAX)
      MESH%NPMAX  = NBMPTS(IELM1)
      ALLOCATE(MESH%MXPTVS)
      MESH%MXPTVS = MXPTVS
      ALLOCATE(MESH%MXELVS)
      MESH%MXELVS = MXELVS
C     LV WILL BE RECOMPUTED LATER
      ALLOCATE(MESH%LV)
      MESH%LV     = 1
      ALLOCATE(MESH%NSEG)
      MESH%NSEG = NBSEG(IELM1)
C
C     2) INTEGER ARRAYS
C
C     ALLOCATION DE IKLE ET KLEI (MEME TAILLE, DEUX DIMENSIONS INVERSEES)
C
      ALLOCATE(MESH%IKLE)
      CALL ALLVEC(2,MESH%IKLE,'IKLE  ',IELM0,NBPEL(IELM),1)
      ALLOCATE(MESH%KLEI)
      CALL ALLVEC(2,MESH%KLEI,'KLEI  ',IELM0,NBPEL(IELM),1)
C
C     IFABOR
C
      ALLOCATE(MESH%IFABOR)
      CALL ALLVEC(2,MESH%IFABOR,'IFABOR',IELM0,NBFEL(IELM),1)
C
C     NELBOR
C
! nelbor & nulone
! it is now changed to VERTICAL boundary element... 
! 2D not affected, 3D usage nelbo3(nptfr,netage) - no. of lat. bd. elements
!
      ALLOCATE(MESH%NELBOR)
      CALL ALLVEC(2,MESH%NELBOR,'NELBOR',ielbor(ielm0,2),1,1)
!
!     NULONE
!
!
! extraordinarily strange geometrically
! in 2D number of boundary nodes is equal to the number of boundary
! elements... in 3D it is NOT the case! 
! nulone 3D is used internally as: nulone(nptfr,netage,4)
! "associe la numerotation locale de bord a la numerotation locale 3D"
!
      ALLOCATE(MESH%NULONE)
!    
      CALL ALLVEC(2,mesh%nulone,'NULONE',
     *            ielbor(ielm0,2),nbpel(ielbor(ielm,2)),1)
!      
!!! NOTE : for the tetraedrons, this is no longer the case. We read the 
!          boundary connectivity table BEFORE initializing the number of
!          boundary nodes and elements. This is well defined now.
!
! in 2D it is call ALLVEC(2, mesh%nulone, 'NULONE', 0,  2, 1) 
! which, in 2D only is equivalent to 
!             call ALLVEC(2, mesh%nulone, 'NULONE', 1,  2, 1)
! in 3D it is call ALLVEC(2, mesh%nulone, 'NULONE', 20, 4, 1)
!
!
C     KP1BOR
C
      ALLOCATE(MESH%KP1BOR)
      CALL ALLVEC(2,MESH%KP1BOR,'KP1BOR', IELBOR(IELM,1),2,1)
C
C     NBOR : NUMEROS GLOBAUX DES NOEUDS DE BORD
C     ALLOCATION OF NBOR, IT WILL BE READ IN LECLIM
C
      ALLOCATE(MESH%NBOR)
      CALL ALLVEC(2, MESH%NBOR,'NBOR  ', IELBOR(IELM,2) , 1, 1)
!
!     IKLBOR : IKLE DES SEGMENTS OU FACES DE BORD
!
!
! in alme3d:
! in 2D: call ALLVEC(2, mesh%iklbor,'IKLBOR',  1, 2, 1) ! why 1, 2?
! in 3D: call ALLVEC(2, mesh%iklbor,'IKLBOR', 20, 4, 1) ! proper 
! usage iklbor(nptfr,netage,4)
!
      ALLOCATE(MESH%IKLBOR)
!                    
      IF(IELM.EQ.41.OR.IELM.EQ.51) THEN
!       FACES DE BORD LATERALES
!       VOIR ININDS POUR LES ELEMENTS 20 ET 21
        CALL ALLVEC(2,MESH%IKLBOR,'IKLBOR',
     *              IELBOR(IELM0,2),NBPEL(IELBOR(IELM1,2)),1)
      ELSEIF(IELM.EQ.11.OR.IELM.EQ.12.OR.IELM.EQ.13
     *   .OR.IELM.EQ.31) THEN
        CALL ALLVEC(2,MESH%IKLBOR,'IKLBOR',
     *              IELBOR(IELM0,1),NBPEL(IELBOR(IELM ,2)),1)
      ELSE
        WRITE(LU,*) 'ALMESH : UNKNOWN ELEMENT FOR IKLBOR:',IELM
        CALL PLANTE(1)
        STOP
      ENDIF
C
C     IFANUM : NUMERO DE LA FACE DANS L'ELEMENT ADJACENT
C
      ALLOCATE(MESH%IFANUM)
      IF(CFG(1).EQ.2) THEN
        CALL ALLVEC(2,MESH%IFANUM,'IFANUM',
     &              3*NBMPTS(IELM0)+NBMPTS(01) ,1,0 )
      ELSEIF(CFG(1).NE.1.AND.CFG(1).NE.3) THEN
        IF(LNG.EQ.1) WRITE(LU,98) CFG(1)
        IF(LNG.EQ.2) WRITE(LU,99) CFG(1)
98      FORMAT(1X,'ALMESH : STOCKAGE INCONNU :',1I6)
99      FORMAT(1X,'ALMESH : UNKNOWN STORAGE:',1I6)
        CALL PLANTE(1)
        STOP
      ELSE
        CALL ALLVEC(2,MESH%IFANUM,'IFANUM', 0 , 1,0 )
      ENDIF
C
C     IKLEM1 : TABLES DE CONNECTIVITE INVERSE POUR PRODUIT FRONTAL
C     LIMVOI : NUMERO LIMITE D'UN NOMBRE DE VOISINS DONNE
C
      ALLOCATE(MESH%IKLEM1)
      ALLOCATE(MESH%LIMVOI)
      IF(CFG(2).EQ.2) THEN
C       CALL ALLVEC(2,MESH%IKLEM1,'IKLEM1',NBMPTS(IELM1)*MXPTVS,4,0)
C       POUR OPTASS=3: SYM et NON SYM DIFFERENTS
        CALL ALLVEC(2,MESH%IKLEM1,'IKLEM1',NBMPTS(IELM1)*MXPTVS,8,0)
        CALL ALLVEC(2,MESH%LIMVOI,'LIMVOI',MXPTVS,2,0)
      ELSEIF(CFG(2).NE.1) THEN
        IF(LNG.EQ.1) WRITE(LU,96) CFG(2)
        IF(LNG.EQ.2) WRITE(LU,97) CFG(2)
96      FORMAT(1X,'ALMESH : PRODUIT MATRICE-VECTEUR INCONNU :',1I6)
97      FORMAT(1X,'ALMESH: UNKNOWN MATRIX-VECTOR PRODUCT:',1I6)
        STOP
      ELSE
        CALL ALLVEC(2,MESH%IKLEM1,'IKLEM1',0,4,0)
        CALL ALLVEC(2,MESH%LIMVOI,'LIMVOI',0,2,0)
      ENDIF
C
C     INTEGER ARRAYS FOR SEGMENT-BASED STORAGE
C
      ALLOCATE(MESH%GLOSEG)
      ALLOCATE(MESH%ELTSEG)
      ALLOCATE(MESH%ORISEG)
C     FROM 5.9 ON ALWAYS DONE IN 2D (IELM=11,12,13 OR 14)
C     FROM 6.0 ON ALWAYS DONE
!     IF(CFG(1).EQ.3.OR.10*(IELM/10).EQ.10) THEN
        CALL ALLVEC(2,MESH%GLOSEG,'GLOSEG',NBSEG(IELM)  ,2,0)        
        CALL ALLVEC(2,MESH%ELTSEG,'ELTSEG',NBMPTS(IELM0),
     *                                     NBSEGEL(IELM),0)  
        CALL ALLVEC(2,MESH%ORISEG,'ORISEG',NBMPTS(IELM0),
     *                                     NBSEGEL(IELM),0)
!     ELSE
!       CALL ALLVEC(2,MESH%GLOSEG,'GLOSEG',0,2,0)        
!       CALL ALLVEC(2,MESH%ELTSEG,'ELTSEG',0,3,0)        
!       CALL ALLVEC(2,MESH%ORISEG,'ORISEG',0,3,0)
!     ENDIF
C      
C     INTEGER ARRAYS FOR PARALLELISM
C
C     KNOLG
C     NACHB
C     ISEG
C     KNOGL
C     INDPU
C     NHP
C     NHM
C
      ALLOCATE(MESH%KNOLG)
      ALLOCATE(MESH%NACHB)
      ALLOCATE(MESH%ISEG)
      ALLOCATE(MESH%KNOGL)
      ALLOCATE(MESH%INDPU)
      ALLOCATE(MESH%NHP)
      ALLOCATE(MESH%NHM)
      ALLOCATE(MESH%IFAPAR)
      ALLOCATE(MESH%NB_NEIGHB)
      ALLOCATE(MESH%NB_NEIGHB_SEG)
C     THEIR ALLVEC IS IN PARINI
      ALLOCATE(MESH%NB_NEIGHB_PT)
      ALLOCATE(MESH%LIST_SEND)
      ALLOCATE(MESH%NH_COM)      
      ALLOCATE(MESH%NB_NEIGHB_PT_SEG)
      ALLOCATE(MESH%LIST_SEND_SEG)
      ALLOCATE(MESH%NH_COM_SEG)
C
      IF(NCSIZE.GT.1) THEN
C
        CALL ALLVEC(2,MESH%KNOLG,'KNOLG ',IELM1            ,1,1)
        CALL ALLVEC(2,MESH%NACHB,'NACHB ',NBMAXNSHARE*NPTIR,1,0)
        CALL ALLVEC(2,MESH%ISEG ,'ISEG  ',IELBOR(IELM1,1)  ,1,1)
! FH-jaj
! For the size of KNOGL, we need the number of nodes
! We can't have this value, then we take the sum of the nodes
! of each sub-mesh (value bigger than the number of nodes in the mesh)
! The size of KNOGL is thus ever bigger than the number of nodes in the mesh
        NPOIN_MAX = P_ISUM(MESH%NPOIN)
        CALL ALLVEC(2,MESH%KNOGL,'KNOGL ',NPOIN_MAX,1,0)
! FH-jaj
        CALL ALLVEC(2,MESH%INDPU,'INDPU ',IELM1              ,1,1)
        CALL ALLVEC(2,MESH%NHP  ,'NHP   ',NBMAXDSHARE*2*NPTIR,1,0)
        CALL ALLVEC(2,MESH%NHM  ,'NHM   ',NBMAXDSHARE*2*NPTIR,1,0)
        CALL ALLVEC(2,MESH%IFAPAR,'IFAPAR',10,6,1)
        DO I=1,6*NBPTS(10)
          MESH%IFAPAR%I(I)=0
        ENDDO        
C
      ELSE
C
        CALL ALLVEC(2,MESH%KNOLG ,'KNOLG ',0,1,0)
        CALL ALLVEC(2,MESH%NACHB ,'NACHB ',0,1,0)
        CALL ALLVEC(2,MESH%ISEG  ,'ISEG  ',0,1,0)
        CALL ALLVEC(2,MESH%KNOGL ,'KNOGL ',0,1,0)
        CALL ALLVEC(2,MESH%INDPU ,'INDPU ',0,1,0)
        CALL ALLVEC(2,MESH%NHP   ,'NHP   ',0,1,0)
        CALL ALLVEC(2,MESH%NHM   ,'NHM   ',0,1,0)
        CALL ALLVEC(2,MESH%IFAPAR,'IFAPAR',0,1,0)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  VOLUMES FINIS
C
      ALLOCATE(MESH%NUBO)
      IF(EQUA(1:15).EQ.'SAINT-VENANT VF') THEN
        CALL ALLVEC(2,MESH%NUBO,'NUBO  ',2*NSEG,1,0)
      ELSE
        CALL ALLVEC(2,MESH%NUBO,'NUBO  ',     0,1,0)
      ENDIF
C
      ALLOCATE(MESH%JMI)
      IF(EQUA(1:15).EQ.'SAINT-VENANT VF') THEN
        CALL ALLVEC(2,MESH%JMI,'JMI   ',NSEG,1,0)
      ELSE
        CALL ALLVEC(2,MESH%JMI,'JMI   ',   0,1,0)
      ENDIF
C
C-----------------------------------------------------------------------
C
C     FILLING ARRAYS IKLE, X AND Y (AND Z)
C
C     PRISMS CUT INTO TETRAHEDRONS
C
      IF(IELM.EQ.51) THEN 
C
        CALL CPIKLE3(MESH%IKLE%I,IKLES,NELEM,NNELMAX,NPOIN,NNPLAN)
C
C       NOTE : NO Z HERE, AS IELM.EQ.41, SEE NOTE BELOW       
        CALL READGEO3(MESH%KNOLG%I,MESH%X%R,MESH%Y%R,NPOIN,NFIC,IB)
C
C     PRISMS
C
      ELSEIF(IELM.EQ.41) THEN
C 
        CALL CPIKLE2(MESH%IKLE%I,MESH%KLEI%I,IKLES,
     *               NELEM,NNELMAX,NPOIN,NNPLAN)
C       NOTE : WITH PRISMS Z IS COMPUTED WITH ZF AND H, OR
C              READ IN THE PREVIOUS COMPUTATION FILE, HENCE NO Z HERE      
        CALL READGEO3(MESH%KNOLG%I,MESH%X%R,MESH%Y%R,NPOIN,NFIC,IB)
C
C     TRIANGLES OR TETRAHEDRONS
C
      ELSEIF(IELM.EQ.11.OR.IELM.EQ.12.OR.IELM.EQ.13.OR.IELM.EQ.14
     *   .OR.IELM.EQ.31) THEN
C
C       IKLES(NDP,NELEM) COPIED INTO IKLE(NELMAX,NDP) AND KLEI(NDP,NELMAX)
        DO I = 1,NDP
          DO IELEM  = 1,NELEM
            MESH%IKLE%I((I-1)*NNELMAX+IELEM) = IKLES((IELEM-1)*NDP+I)
            MESH%KLEI%I((IELEM-1)*NDP+I)     = IKLES((IELEM-1)*NDP+I)
          ENDDO
        ENDDO
        IF(IELM.EQ.11.OR.IELM.EQ.12.OR.IELM.EQ.13.OR.IELM.EQ.14) THEN
          CALL READGEO3(MESH%KNOLG%I,MESH%X%R,MESH%Y%R,NPOIN,NFIC,IB)
        ELSEIF(IELM.EQ.31) THEN
C         TETRAHEDRONS: READING THE Z COORDINATE AFTER X AND Y
          CALL READGEO3(MESH%KNOLG%I,MESH%X%R,MESH%Y%R,NPOIN,NFIC,IB,
     *                  MESH%Z%R)
        ENDIF
C
      ELSE
C
C OTHER ELEMENT TYPES
C
        WRITE(LU,*) 'ALMESH : UNKNOWN ELEMENT:',IELM
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C     COMPLEMENT DES TABLEAUX X , Y  POUR LES PRISMES ET TETRAEDRES
C
      IF(IELM.EQ.41.OR.IELM.EQ.51) THEN
        DO I = 2,NNPLAN
          CALL OV_2( 'X=Y     ' , MESH%X%R,I, MESH%X%R,1,
     *                            MESH%X%R,1, 0.D0, NNPMAX,NPOIN)
          CALL OV_2( 'X=Y     ' , MESH%Y%R,I, MESH%Y%R,1,
     *                            MESH%Y%R,1, 0.D0, NNPMAX,NPOIN)
        END DO
      ENDIF
!
!  WATCH OUT - D Y N A M I T E 
!
! for 3D matrices computation, x,y,z (per node) are used instead
! of xel,yel,zel (per element)...  -- requires a thorough checking 
! of BIEF5, where xel,yel,zel are used in calls for 3D cases
! I set simply:
!                      PROVISOIRE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC DEBUT MAILAGE-3D
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(IELM.EQ.31.OR.IELM.EQ.41.OR.IELM.EQ.51) THEN
        mesh%xel => mesh%x
        mesh%yel => mesh%y
        mesh%zel => mesh%z
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC FIN MAILAGE-3D
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
C
C-----------------------------------------------------------------------
C
C DEALLOCATE TEMPORARY TABLES
C
      DEALLOCATE(IKLES)
    
C-----------------------------------------------------------------------
C
      IF(LNG.EQ.1) WRITE(LU,*) 'MAILLAGE : ',NOM,' ALLOUE'
      IF(LNG.EQ.2) WRITE(LU,*) 'MESH: ',NOM,' ALLOCATED'
C
C-----------------------------------------------------------------------
C
      RETURN
      END
