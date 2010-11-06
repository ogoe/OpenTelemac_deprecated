C                       ************************
                        SUBROUTINE POINT_SISYPHE
C                       ************************
C
C
C***********************************************************************
C SISYPHE VERSION 6.0                            C. LENORMANT   11/09/95
C                                                J.-M. HERVOUET
C                                                C. MACHET      10/06/02
C
C
C  16/06/2008 JMH : BOUNDARY_COLOUR AJOUTE
C  16/09/2009 JMH : AVAIL(NPOIN,10,NSICLA)
C  18/09/2009 JMH : SEE AVAI AND LAYTHI
C
C***********************************************************************
C
C
C                                                
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT   
C***********************************************************************
C
C     FONCTION  : ALLOCATIONS DES STRUCTURES
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |                |    |  
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C PROGRAMME APPELANT : HOMERE
C PROGRAMMES APPELES : ININDS,ALLVEC,ALLBLO,ALLMAT,NBMPTS,NBFEL,ALMESH
C***********************************************************************
C
      ! 1/ MODULES
      ! ----------
      USE BIEF
      USE DECLARATIONS_SISYPHE
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU


      ! 2/ LOCAL VARIABLES
      ! ------------------
      INTEGER :: I,K,NTR,IELM0,IELM1,IELBT,IELM0_SUB
      INTEGER :: CFG(2),CFGBOR(2)

C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------

      IF (LNG == 1) WRITE(LU,11)
      IF (LNG == 2) WRITE(LU,12)

      ! ************************************** !
      ! I - DISCRETISATION ET TYPE DE STOCKAGE !
      ! ************************************** !
      ! IELMT, IELMH_SIS et IELMU_SIS codee en dur dans lecdon
      IELM0     = 10
      IELM1     = 11
      IELBT     = IELBOR(IELMT,1)
      IELM0_SUB = 10*(IELMT/10)

      CFG(1)    = OPTASS                                                  
      CFG(2)    = PRODUC
      CFGBOR(1) = 1 ! CFG impose pour les matrices de bord
      CFGBOR(2) = 1 ! CFG impose pour les matrices de bord

      IF(VF) EQUA(1:15)='SAINT-VENANT VF'

      ! ******************************************* !
      ! II - ALLOCATION DE LA STRUCTURE DU MAILLAGE ! 
      ! ******************************************* !
      CALL ALMESH(MESH,'MESH_S',IELMT,SPHERI,CFG,
     *            SIS_FILES(SISGEO)%LU,EQUA)

      IKLE  => MESH%IKLE
      X     => MESH%X%R
      Y     => MESH%Y%R
      NELEM => MESH%NELEM
      NELMAX=> MESH%NELMAX
      NPTFR => MESH%NPTFR
      NPTFRX=> MESH%NPTFRX
      DIM   => MESH%DIM
      TYPELM=> MESH%TYPELM
      NPOIN => MESH%NPOIN
      NPMAX => MESH%NPMAX
      MXPTVS=> MESH%MXPTVS
      MXELVS=> MESH%MXELVS
      LV    => MESH%LV


      ! ******************** !
      ! III - TABLEAUX REELS !
      ! ******************** !
      CALL ALLVEC(1, S     , 'S     ', 0    , 1, 1) ! STRUCTURE VIDE
!
      CALL ALLVEC(1, E     , 'E     ', IELMT, 1, 2) ! RESULTAT
      CALL ALLVEC(1, Z     , 'Z     ', IELMT, 1, 2) ! RESULTAT
      CALL ALLVEC(1, DEL_Z , 'DEL_Z ', IELMT, 1, 2) ! INCREMENT OF Z IF HYDRO
      CALL ALLVEC(1, ZF_C  , 'ZF_C  ', IELMT, 1, 2) ! VARIABLES E SOMMEES
      CALL ALLVEC(1, ZF_S  , 'ZF_S  ', IELMT, 1, 2) ! VARIABLES E SOMMEES
      CALL ALLVEC(1, ESOMT , 'ESOMT ', IELMT, 1, 2) ! VARIABLES E SOMMEES
      CALL ALLVEC(1, EMAX  , 'EMAX  ', IELMT, 1, 2) ! VARIABLES E SOMMEES
      CALL ALLVEC(1, Q     , 'Q     ', IELMT, 1, 2) ! FLOWRATE
      CALL ALLVEC(1, QU    , 'QU    ', IELMT, 1, 2) ! X FLOWRATE
      CALL ALLVEC(1, QV    , 'QV    ', IELMT, 1, 2) ! Y FLOWRATE
      CALL ALLVEC(1, DEL_QU, 'DEL_QU', IELMT, 1, 2) ! INCREMENT OF QU IF HYDRO
      CALL ALLVEC(1, DEL_QV, 'DEL_QV', IELMT, 1, 2) ! INCREMENT OF QV IF HYDRO
      CALL ALLVEC(1, U2D   , 'U2D   ', IELMT, 1, 2) ! X VELOCITY
      CALL ALLVEC(1, V2D   , 'V2D   ', IELMT, 1, 2) ! Y VELOCITY
      CALL ALLVEC(1, QS    , 'QS    ', IELMT, 1, 2) ! TRANSPORT RATE
      CALL ALLVEC(1, QSX   , 'QSX   ', IELMT, 1, 2) ! X TRANSPORT RATE
      CALL ALLVEC(1, QSY   , 'QSY   ', IELMT, 1, 2) ! Y TRANSPORT RATE
      CALL ALLVEC(1, QS_C  , 'QS_C  ', IELMT, 1, 2) ! BEDLOAD RATE
      CALL ALLVEC(1, QSXC  , 'QSXC  ', IELMT, 1, 2) ! X BEDLOAD RATE
      CALL ALLVEC(1, QSYC  , 'QSYC  ', IELMT, 1, 2) ! Y BEDLOAD RATE
      CALL ALLVEC(1, QS_S  , 'QS_S  ', IELMT, 1, 2) ! SUSPENSION RATE
      CALL ALLVEC(1, QSXS  , 'QSXS  ', IELMT, 1, 2) ! X SUSPENSION RATE
      CALL ALLVEC(1, QSYS  , 'QSYS  ', IELMT, 1, 2) ! Y SUSPENSION RATE
      CALL ALLVEC(1, HIDING, 'HIDING', IELMT, 1, 2) ! HIDING FACTOR
      CALL ALLVEC(1, ZF    , 'ZF    ', IELMT, 1, 2) ! COTES DU FOND
      CALL ALLVEC(1, ZR    , 'ZR    ', IELMT, 1, 2) ! COTES DU FOND NON ERODABLE
      CALL ALLVEC(1, ZREF  , 'ZREF  ', IELMT, 1, 2) ! REFERENCE ELEVATION
      CALL ALLVEC(1, CHESTR, 'CHESTR', IELMT, 1, 2) ! COEFFICIENT DE FROTTEMENT
      CALL ALLVEC(1, COEFPN, 'COEFPN', IELMT, 1, 2) ! EFFET DE PENTE
      CALL ALLVEC(1, CALFA , 'CALFA ', IELMT, 1, 2)
      CALL ALLVEC(1, SALFA , 'SALFA ', IELMT, 1, 2)
      CALL ALLVEC(1, CF    , 'CF    ', IELMT, 1, 2) ! FROTTEMENt ADIMENSIONNEL
      CALL ALLVEC(1, TOB   , 'TOB   ', IELMT, 1, 2) ! FROTTEMENT TOTAL 
      CALL ALLVEC(1, TOBW  , 'TOBW  ', IELMT, 1, 2) ! varaibles houle
      CALL ALLVEC(1, MU    , 'MU    ', IELMT, 1, 2) ! Frott. de peau
      CALL ALLVEC(1, KSP   , 'KSP   ', IELMT, 1, 2) ! RUGOSITE de peau
      CALL ALLVEC(1, KS    , 'KS    ', IELMT, 1, 2) !  RUGOSITE totale
      CALL ALLVEC(1, KSR   , 'KSR   ', IELMT, 1, 2) !  RUGOSITE due aux rides
      CALL ALLVEC(1, THETAW, 'THETAW', IELMT, 1, 2) !variables houle
      CALL ALLVEC(1, FW    , 'FW    ', IELMT, 1, 2) ! variables houle
      CALL ALLVEC(1, UW    , 'UW    ', IELMT, 1, 2) !variables houle
      CALL ALLVEC(1, HW    , 'HW    ', IELMT, 1, 2)
      CALL ALLVEC(1, TW    , 'TW    ', IELMT, 1, 2)
      CALL ALLVEC(1, DZF_GF, 'DZF_GF', IELMT, 1, 2) ! BED LEVEL CHANGE FOR GRAIN-FEEDING
      CALL ALLVEC(1, ACLADM, 'ACLADM', IELMT, 1, 2) ! DIAMETRE MOYEN DE LA COUCHE ACTIVE
      CALL ALLVEC(1, UNLADM, 'UNLADM', IELMT, 1, 2) ! DIAMETRE MOYEN DE LA DEUXIEME COUCHE
      CALL ALLVEC(1, HCPL  , 'HCPL  ', IELMT, 1, 2) ! WATER DEPTH SAVED FOR CONSTANT FLOW DISCHARGE
      CALL ALLVEC(1, ECPL  , 'ECPL  ', IELMT, 1, 2) ! EVOLUTION SAVED FOR CONSTANT FLOW DISCHARGE
      CALL ALLVEC(1, ELAY  , 'ELAY  ', IELMT, 1, 2) ! EPAISSEUR DE LA COUCHE ACTIVE
      CALL ALLVEC(1, ESTRAT, 'ESTRAT', IELMT, 1, 2) ! EPAISSEUR DE LA DEUXIEME COUCHE
      CALL ALLVEC(1, KX    , 'KX    ', IELMT, 1, 1)
      CALL ALLVEC(1, KY    , 'KY    ', IELMT, 1, 1)
      CALL ALLVEC(1, KZ    , 'KZ    ', IELMT, 1, 1)
      CALL ALLVEC(1, UCONV , 'UCONV ', IELMT, 1, 1)
      CALL ALLVEC(1, VCONV , 'VCONV ', IELMT, 1, 1) 
      CALL ALLVEC(1, UNORM , 'UNORM ', IELMT, 1, 2)
      CALL ALLVEC(1, DISP  , 'DISP  ', IELMT, 3, 1)
      CALL ALLVEC(1, DISP_C, 'DISP_C', IELMT, 3, 1)
      CALL ALLVEC(1, MASKB , 'MASKB ', IELM0, 1, 2)
      CALL ALLVEC(1, MASK  , 'MASK  ', IELBT, 1, 2)
      CALL ALLVEC(1, AFBOR , 'AFBOR ', IELBT, 1, 1)
      CALL ALLVEC(1, BFBOR , 'BFBOR ', IELBT, 1, 1)
      CALL ALLVEC(1, FLBOR , 'FLBOR ', IELBT, 1, 1)
!     FLUX AUX BORDS POUR APPEL A CVDFTR
      CALL ALLVEC(1, FLBOR_SIS , 'FLBORS', IELBT, 1, 1)
      CALL ALLVEC(1, FLBORTRA  , 'FLBTRA', IELBT, 1, 1)
      CALL ALLVEC(1, CSTAEQ, 'CSTAEQ', IELMT, 1, 2)
      CALL ALLVEC(1, HN    , 'HN    ', IELMH_SIS, 1, 2) ! WATER DEPTH
      CALL ALLVEC(1, HCLIP , 'HCLIP ', IELMH_SIS, 1, 2) ! CLIPPING WATER DEPTH
      CALL ALLVEC(1, HPROP , 'HPROP ', IELMH_SIS, 1, 1)
      CALL ALLVEC(1, VOLU2D, 'VOLU2D', IELMH_SIS, 1, 1)
      CALL ALLVEC(1, V2DPAR, 'V2DPAR', IELMH_SIS, 1, 1)
      CALL ALLVEC(1, UNSV2D, 'UNSV2D', IELMH_SIS, 1, 1)
!
      IF(MSK) THEN
        CALL ALLVEC(1,MASKEL,'MASKEL', IELM0 , 1 , 2 )
        CALL ALLVEC(1,MSKTMP,'MSKTMP', IELM0 , 1 , 2 )
        CALL ALLVEC(1,MASKPT,'MASKPT', IELMT , 1 , 2 )
      ELSE
        CALL ALLVEC(1,MASKEL,'MASKEL', 0 , 1 , 0 )
        CALL ALLVEC(1,MSKTMP,'MSKTMP', 0 , 1 , 0 )
        CALL ALLVEC(1,MASKPT,'MASKPT', 0 , 1 , 0 )
      ENDIF
!
!     POUR LES SEDIMENTS MIXTES
!
      IF(MIXTE.OR.TASS) THEN
        CALL ALLVEC(1,FLUER_VASE,'FRMIXT',IELMT,1,2)
        CALL ALLVEC(1,TOCE_MIXTE ,'TCMIXT',IELMT,10,2)
        CALL ALLVEC(1,MS_SABLE   ,'MSSABL',IELMT,10,2)
        CALL ALLVEC(1,MS_VASE    ,'MSVASE',IELMT,10,2)
      ELSE
        CALL ALLVEC(1,FLUER_VASE,'FRMIXT',0,1,0)
        CALL ALLVEC(1,TOCE_MIXTE ,'TCMIXT',0,1,0)
        CALL ALLVEC(1,MS_SABLE   ,'MSSABL',0,1,0)
        CALL ALLVEC(1,MS_VASE    ,'MSVASE',0,1,0)
      ENDIF
!
      ! *********************** !
      ! IV - TABLEAUX D'ENTIERS ! (_IMP_)
      ! *********************** !
      CALL ALLVEC(2, LIEBOR, 'LIEBOR', IELBOR(IELM1,1), 1, 1)
      CALL ALLVEC(2, LIQBOR, 'LIQBOR', IELBOR(IELM1,1), 1, 1)
      CALL ALLVEC(2, LIMTEC, 'LIMTEC', IELBOR(IELM1,1), 1, 1)
      CALL ALLVEC(2, NUMLIQ, 'NUMLIQ', IELBOR(IELM1,1), 1, 1)
      CALL ALLVEC(2, CLT   , 'CLT   ', IELBOR(IELMT,1), 1, 1)   
      CALL ALLVEC(2, CLU   , 'CLU   ', IELBOR(IELMT,1), 1, 1)   
      CALL ALLVEC(2, CLV   , 'CLV   ', IELBOR(IELMT,1), 1, 1)
      CALL ALLVEC(2, LIMDIF, 'LIMDIF', IELBOR(IELMT,1), 1, 1)
      CALL ALLVEC(2, LICBOR, 'LICBOR', IELBOR(IELMT,1), 1, 1)
      CALL ALLVEC(2, LIHBOR, 'LIHBOR', IELBOR(IELMT,1), 1, 1)
      CALL ALLVEC(2, BOUNDARY_COLOUR,
     *                       'BNDCOL', IELBOR(IELMT,1), 1, 1)
      CALL ALLVEC(2, LIMPRO, 'LIMPRO', IELBOR(IELMT,1), 6, 1)
      CALL ALLVEC(2, INDIC , 'INDIC ', IELM1          , 1, 1)
      CALL ALLVEC(2, IT1   , 'IT1   ', IELM1          , 1, 2)
      CALL ALLVEC(2, IT2   , 'IT2   ', IELM1          , 1, 2)
      CALL ALLVEC(2, IT3   , 'IT3   ', IELM1          , 1, 2)
      CALL ALLVEC(2, IT4   , 'IT4   ', IELM1          , 1, 2)
      CALL ALLVEC(2, NLAYER, 'NLAYE ', IELMT          , 1, 2) ! NOMBRE DE COUCHES

      IF(VF) THEN
        CALL ALLVEC(2,BREACH,'BREACH',IELM1,1,2)
      ELSE
        CALL ALLVEC(2,BREACH,'BREACH',0,1,0)
      ENDIF

      IF(MSK) THEN
        CALL ALLVEC(2,IFAMAS,'IFAMAS',IELM0,NBFEL(IELM0),1)
      ELSE
        CALL ALLVEC(2,IFAMAS,'IFAMAS',0,1,0)
      ENDIF

      ! ******************* !
      ! V - BLOCK OF ARRAYS !
      ! ******************* !
      ALLOCATE(AVAIL(NPOIN,10,NSICLA)) ! fraction of each class for each layer
      ALLOCATE(ES(NPOIN,10))           ! thickness of each class ???

      !================================================================!
      CALL ALLBLO(MASKTR, 'MASKTR') ! mask of the boundary conditions
      CALL ALLBLO(EBOR  , 'EBOR  ') ! boundary conditions
      CALL ALLBLO(QBOR  , 'QBOR  ') ! boundary conditions
      CALL ALLBLO(AVAI  , 'AVAI  ') ! fraction of each class for the two first layers
      CALL ALLBLO(LAYTHI, 'LAYTHI') ! layer thicknesses
      !================================================================!
      CALL ALLBLO(QSCL  , 'QSCL  ') ! transport rate for each class
      CALL ALLBLO(QSCLX , 'QSCLX ') ! transport rate for each class along X
      CALL ALLBLO(QSCLY , 'QSCLY ') ! transport rate for each class along Y
      CALL ALLBLO(QSCL_C, 'QSCL_C') ! bedload transport rate for each class
      CALL ALLBLO(QSCLXC, 'QSCLXC') ! bedload transport rate for each class along X
      CALL ALLBLO(QSCLYC, 'QSCLYC') ! bedload transport rate for each class along Y
      CALL ALLBLO(ZFCL  , 'ZFCL  ') ! evolution for each class
      CALL ALLBLO(ZFCL_C, 'ZFCL_C') ! evolution for each class thanks bedload transport
      !================================================================!
      CALL ALLBLO(CBOR  , 'CBOR  ') ! boundary conditions
      CALL ALLBLO(QSCL_S, 'QSCL_S') ! suspension transport rate for each class
      CALL ALLBLO(QSCLXS, 'QSCLXS') ! suspension transport rate for each class along X
      CALL ALLBLO(QSCLYS, 'QSCLYS') ! suspension transport rate for each class along Y
      CALL ALLBLO(ZFCL_S, 'ZFCL_S') ! evolution for each class thanks suspension transport
      CALL ALLBLO(FLUDP , 'FLUDP ') ! deposition flux
      CALL ALLBLO(FLUDPT, 'FLUDPT') ! deposition flux for implicitation
      CALL ALLBLO(FLUER , 'FLUER ') ! erosion flux
      CALL ALLBLO(FLUERT, 'FLUERT') ! erosion flux for implicitation
      CALL ALLBLO(CS    , 'CS    ') ! concentration at time n
      CALL ALLBLO(CTILD , 'CTILD ') ! concentration at time n+1/2 (=> advection step)
      CALL ALLBLO(CST   , 'CST   ') ! concentration at time n+1   (=> result)
      !================================================================!
      !================================================================!
      CALL ALLVEC_IN_BLOCK(MASKTR, 4       , 1, 'MSKTR ', IELBT, 1, 2)
      CALL ALLVEC_IN_BLOCK(EBOR  , NSICLA  , 1, 'EBOR  ', IELBT, 1, 2)
      CALL ALLVEC_IN_BLOCK(QBOR  , NSICLA  , 1, 'QBOR  ', IELBT, 1, 2)
!
!     JMH 18/09/09 AVAI ALLOCATED WITH SIZE 0 AND POINTING TO
!                  RELEVANT SECTIONS OF AVAIL
!     CALL ALLVEC_IN_BLOCK(AVAI,NOMBLAY*NSICLA,1,'AVAI  ',IELMT,1,2)
      CALL ALLVEC_IN_BLOCK(AVAI,NOMBLAY*NSICLA,1,'AVAI  ',    0,1,0)
      DO I=1,NSICLA
        DO K=1,NOMBLAY
          AVAI%ADR(K+(I-1)*NOMBLAY)%P%R=>AVAIL(1:NPOIN,K,I)
          AVAI%ADR(K+(I-1)*NOMBLAY)%P%MAXDIM1=NPOIN 
          AVAI%ADR(K+(I-1)*NOMBLAY)%P%DIM1=NPOIN  
        ENDDO
      ENDDO
!     LAYTHI ALLOCATED WITH SIZE 0 AND POINTING TO RELEVANT SECTIONS OF ES
!     CALL ALLVEC_IN_BLOCK(LAYTHI,NOMBLAY, 1, 'LAYTHI', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(LAYTHI,NOMBLAY, 1, 'LAYTHI',     0, 1, 0)
      DO K=1,NOMBLAY
        LAYTHI%ADR(K)%P%R=>ES(1:NPOIN,K)
        LAYTHI%ADR(K)%P%MAXDIM1=NPOIN 
        LAYTHI%ADR(K)%P%DIM1=NPOIN 
      ENDDO
!
      !================================================================!
      CALL ALLVEC_IN_BLOCK(QSCL  , NSICLA  , 1, 'QSCL  ', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(QSCLX , NSICLA  , 1, 'QSCLX ', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(QSCLY , NSICLA  , 1, 'QSCLY ', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(QSCL_C, NSICLA  , 1, 'QSCL_C', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(QSCLXC, NSICLA  , 1, 'QSCLXC', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(QSCLYC, NSICLA  , 1, 'QSCLYC', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(ZFCL  , NSICLA  , 1, 'ZFCL  ', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(ZFCL_C, NSICLA  , 1, 'ZFCL_C', IELMT, 1, 2)
      !================================================================!
      CALL ALLVEC_IN_BLOCK(CBOR  , NSICLA  , 1, 'CBOR  ', IELBT, 1, 2)
      CALL ALLVEC_IN_BLOCK(QSCL_S, NSICLA  , 1, 'QSCL_S', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(QSCLXS, NSICLA  , 1, 'QSCLXS', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(QSCLYS, NSICLA  , 1, 'QSCLYS', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(ZFCL_S, NSICLA  , 1, 'ZFCL_S', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(FLUDP , NSICLA  , 1, 'FLUDP ', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(FLUDPT, NSICLA  , 1, 'FLUDPT', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(FLUER , NSICLA  , 1, 'FLUER ', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(FLUERT, NSICLA  , 1, 'FLUERT', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(CS    , NSICLA  , 1, 'CS    ', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(CTILD , NSICLA  , 1, 'CTILD ', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(CST   , NSICLA  , 1, 'CST   ', IELMT, 1, 2)
      !================================================================!


      ! ************* !
      ! VI - MATRICES !
      ! ************* !


      !================================================================!
      CALL ALLMAT(AM1_S, 'AM1_S ', IELMT, IELMT, CFG   , 'Q', 'Q') ! suspension work matrix
      CALL ALLMAT(AM2_S, 'AM2_S ', IELMT, IELMT, CFG   , 'Q', 'Q') ! suspension work matrix
      CALL ALLMAT(MBOR , 'MBOR  ', IELBT, IELBT, CFGBOR, 'Q', 'Q') ! suspension boundray matrix
      !================================================================!


      ! ****************** !
      ! VII - OTHER ARRAYS !
      ! ****************** !
!
!     NTR SHOULD BE AT LEAST THE NUMBER OF VARIABLES IN VARSOR THAT WILL BE READ IN 
!     VALIDA. HERE UP TO THE LAYER THICKNESSES 
      NTR   = 26+(NOMBLAY+4)*NSICLA+NOMBLAY+NPRIV   
      IF(SLVSED%SLV == 7) NTR = MAX(NTR,2+2*SLVSED%KRYLOV)
      IF(SLVTRA%SLV == 7) NTR = MAX(NTR,2+2*SLVTRA%KRYLOV)
      IF(3*(SLVSED%PRECON/3) == SLVSED%PRECON) NTR = NTR + 2 ! si precond. BLOC-DIAG (+2 diag)
      IF(3*(SLVTRA%PRECON/3) == SLVTRA%PRECON) NTR = NTR + 2 ! si precond. BLOC-DIAG (+2 diag)
!
!     W1 NO LONGER USED (IS SENT TO CVDFTR BUT CVDFTR DOES NOTHING WITH IT)
      CALL ALLVEC(1, W1 , 'W1    ', IELM0    , 1    , 1) ! work array 
      CALL ALLVEC(1, TE1, 'TE1   ', IELM0_SUB, 1    , 1) ! work array by element
      CALL ALLVEC(1, TE2, 'TE2   ', IELM0_SUB, 1    , 1) ! work array by element
      CALL ALLVEC(1, TE3, 'TE3   ', IELM0_SUB, 1    , 1) ! work array by element
!
      CALL ALLBLO(VARCL, 'VARCL ') ! variables clandestines
      CALL ALLBLO(PRIVE, 'PRIVE ') ! user's array
      CALL ALLBLO(TB   , 'TB    ') ! tableau de travail
!
      CALL ALLVEC_IN_BLOCK(TB   , NTR   , 1, 'T     ', IELMT, 1, 2)
      CALL ALLVEC_IN_BLOCK(VARCL, NVARCL, 1, 'CL    ', IELMT, 1, 2)
      IF(NPRIV.GT.0) THEN
        CALL ALLVEC_IN_BLOCK(PRIVE,MAX(NPRIV,4),1,'PRIV  ',IELMT,1, 2)
      ELSE
        CALL ALLVEC_IN_BLOCK(PRIVE,4           ,1,'PRIV  ',    0,1, 0)
      ENDIF
!     POUR EVITER UNE ECRITURE SUR FICHIER DE TABLEAUX NON INITIALISES
      CALL OS('X=0     ',X=PRIVE)
!
      ! ************ !
      ! VIII - ALIAS !
      ! ************ !
!
      T1   => TB%ADR( 1)%P ! work array
      T2   => TB%ADR( 2)%P ! work array
      T3   => TB%ADR( 3)%P ! work array
      T4   => TB%ADR( 4)%P ! work array
      T5   => TB%ADR( 5)%P ! work array
      T6   => TB%ADR( 6)%P ! work array
      T7   => TB%ADR( 7)%P ! work array
      T8   => TB%ADR( 8)%P ! work array
      T9   => TB%ADR( 9)%P ! work array
      T10  => TB%ADR(10)%P ! work array
      T11  => TB%ADR(11)%P ! work array
      T12  => TB%ADR(12)%P ! work array
      T13  => TB%ADR(13)%P ! work array
      T14  => TB%ADR(14)%P ! work array
!
      ! ****************************************************************** !
      ! IX - ALLOCATION D'UN BLOC RELIANT UN NOM DE VARIABLE A SON TABLEAU ! 
      ! ****************************************************************** !
!
      CALL ALLBLO(VARSOR, 'VARSOR')
      CALL ADDBLO(VARSOR, U2D    )            ! 01
      CALL ADDBLO(VARSOR, V2D    )            ! 02
      CALL ADDBLO(VARSOR, HN    )             ! 03
      CALL ADDBLO(VARSOR, Z     )             ! 04
      CALL ADDBLO(VARSOR, ZF    )             ! 05
      CALL ADDBLO(VARSOR, Q     )             ! 06
      CALL ADDBLO(VARSOR, QU    )             ! 07
      CALL ADDBLO(VARSOR, QV    )             ! 08
      CALL ADDBLO(VARSOR, ZR    )             ! 09
      CALL ADDBLO(VARSOR, CHESTR)             ! 10
      CALL ADDBLO(VARSOR, TOB   )             ! 11
      CALL ADDBLO(VARSOR, HW    )             ! 12
      CALL ADDBLO(VARSOR, TW    )             ! 13
      CALL ADDBLO(VARSOR, THETAW)             ! 14
      CALL ADDBLO(VARSOR, QS    )             ! 15
      CALL ADDBLO(VARSOR, QSX   )             ! 16
      CALL ADDBLO(VARSOR, QSY   )             ! 17
      CALL ADDBLO(VARSOR, ESOMT )             ! 18
      CALL ADDBLO(VARSOR, KS)                 ! 19
      CALL ADDBLO(VARSOR, MU)                 ! 20
C
C     AVAI: FROM 21 TO 20+NOMBLAY*NSICLA
C
      DO I = 1,NOMBLAY*NSICLA
        CALL ADDBLO(VARSOR, AVAI%ADR(I)%P)   
      ENDDO
C
C     QSCL: FROM 21+NOMBLAY*NSICLA TO 20+(NOMBLAY+1)*NSICLA
C
      DO I = 1, NSICLA
        CALL ADDBLO(VARSOR, QSCL%ADR(I)%P)  
      ENDDO
C
C     CS: FROM 21+(NOMBLAY+1)*NSICLA TO 20+(NOMBLAY+2)*NSICLA
C
      DO I=1,NSICLA
        CALL ADDBLO(VARSOR, CS%ADR(I)%P)    
      ENDDO
      CALL ADDBLO(VARSOR,QS_C)               ! 21+(NOMBLAY+2)*NSICLA
      CALL ADDBLO(VARSOR,QSXC)               ! 22+(NOMBLAY+2)*NSICLA
      CALL ADDBLO(VARSOR,QSYC)               ! 23+(NOMBLAY+2)*NSICLA
      CALL ADDBLO(VARSOR,QS_S)               ! 24+(NOMBLAY+2)*NSICLA
      CALL ADDBLO(VARSOR,QSXS)               ! 25+(NOMBLAY+2)*NSICLA
      CALL ADDBLO(VARSOR,QSYS)               ! 26+(NOMBLAY+2)*NSICLA
C
C     QSCL_C: FROM 27+(NOMBLAY+2)*NSICLA TO 26+(NOMBLAY+3)*NSICLA
C
      DO I=1,NSICLA
        CALL ADDBLO(VARSOR,QSCL_C%ADR(I)%P)  
      ENDDO
C
C     QSCL_S: FROM 27+(NOMBLAY+3)*NSICLA TO 26+(NOMBLAY+4)*NSICLA
C
      DO I=1,NSICLA
        CALL ADDBLO(VARSOR,QSCL_S%ADR(I)%P)  
      ENDDO
C
C     LAYTHI: FROM 27+(NOMBLAY+4)*NSICLA TO 26+(NOMBLAY+4)*NSICLA+NOMBLAY
C
      DO I=1,NOMBLAY   
        CALL ADDBLO(VARSOR,LAYTHI%ADR(I)%P) ! 26+(NOMBLAY+4)*NSICLA+NOMBLAY
      ENDDO
C
C     PRIVE: FROM 27+(NOMBLAY+4)*NSICLA+NOMBLAY TO
C                 26+(NOMBLAY+4)*NSICLA+MAX(4,NPRIV)+NOMBLAY
C
      DO I=1,MAX(4,NPRIV)
        CALL ADDBLO(VARSOR,PRIVE%ADR(I)%P)
      ENDDO
!
      IF(VARCL%N.GT.0) THEN
        DO I=1,VARCL%N
          CALL ADDBLO(VARSOR,VARCL%ADR(I)%P)
          SORLEO(26+MAX(4,NPRIV)+NSICLA*(NOMBLAY+4)+NOMBLAY+I)=.TRUE.
        ENDDO
      ENDIF
!
!
!-----------------------------------------------------------------------                  
! !jaj #### if required, in this place we read the input sections file                         
!      and modify NCP and CTRLSC(1:NCP) accordingly in read_sections                      
!                                                                                         
      IF(TRIM(SIS_FILES(SISSEC)%NAME).NE.'') THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*)
     &   'POINT_SISYPHE: SECTIONS DEFINIES PAR FICHIER'
        ELSEIF(LNG.EQ.2) THEN
          WRITE(LU,*)
     &   'POINT_SISYPHE: SECTIONS DEFINED IN THE SECTIONS INPUT FILE'
        ENDIF
        CALL READ_SECTIONS_SISYPHE 
      ELSE ! the previously existing way of doing things 
        IF(NCP.NE.0) THEN 
          IF(LNG.EQ.1) THEN                                 
            WRITE(LU,*)
     &      'POINT_SISYPHE: SECTIONS DEFINED IN THE PARAMETER FILE'
          ELSEIF(LNG.EQ.2) THEN
            IF(NCP.NE.0) WRITE(LU,*)
     &      'POINT_SISYPHE: SECTIONS DEFINED IN THE PARAMETER FILE'
          ENDIF
        ENDIF
      ENDIF
!
      IF(LNG == 1) WRITE(LU,21)
      IF(LNG == 2) WRITE(LU,22)
!
11    FORMAT(1X,///,21X,'*******************************',/,
     *21X,              '* ALLOCATION DE LA MEMOIRE    *',/,
     *21X,              '*******************************',/)
21    FORMAT(1X,///,21X,'****************************************',/,
     *21X,              '* FIN DE L''ALLOCATION DE LA MEMOIRE  : *',/,
     *21X,              '****************************************',/)

12    FORMAT(1X,///,21X,'*******************************',/,
     *21X,              '*     MEMORY ORGANISATION     *',/,
     *21X,              '*******************************',/)
22    FORMAT(1X,///,21X,'*************************************',/,
     *21X,              '*    END OF MEMORY ORGANIZATION:    *',/,
     *21X,              '*************************************',/)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
