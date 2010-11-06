C                       **************************
                        SUBROUTINE POINT_TELEMAC2D
C                       **************************
C
C***********************************************************************
C TELEMAC 2D VERSION 6.0  24/03/2010  J-M HERVOUET (LNHE) 01 30 87 80 18
C
C
C 11/07/2008 : DIMENSION DE LIMPRO
C 02/10/2008 : NTR=22
C 02/04/2009 : NGEO REPLACED BY T2D_FILES(T2DGEO)%LU
C 26/11/2009 : ADVECTION SPECIFIC IF EQUA='SAINT-VENANT VF', SO NO
C              OTHER METHOD USED (CHARACTERISTICS, ETC.)
C
C***********************************************************************
C
C     FONCTION  : ALLOCATIONS DES STRUCTURES
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |                |    |  
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : HOMERE
C
C  SOUS-PROGRAMME APPELE :
C
C***********************************************************************
C              
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER MEMW1,NTR,NTRT,NTRKE,I,J,I3,I4,ITRAC
      INTEGER IELMX,IELMC1,IELMC2,IELMUT,IELMHT
      INTEGER IELBU,IELBH,IELBT,IELBK,IELBE,IELB1
      INTEGER IELBX,CFG(2),CFGBOR(2),ERR
C
      CHARACTER*1 TYP
C
      INTRINSIC MAX
C
C-----------------------------------------------------------------------
C
      IF(LISTIN) THEN
         IF(LNG.EQ.1) WRITE(LU,20)
         IF(LNG.EQ.2) WRITE(LU,21)
      ENDIF
20    FORMAT(1X,///,26X,'*****************************',/,
     *26X,              '* ALLOCATION DE LA MEMOIRE  *',/,
     *26X,              '*****************************',/)
21    FORMAT(1X,///,26X,'*****************************',/,
     *26X,              '*    MEMORY ORGANIZATION    *',/,
     *26X,              '*****************************',/)
C
C-----------------------------------------------------------------------
C
C     TYPES DE DISCRETISATIONS
C
      IELM0 = 10*(IELMH/10)
      IELM1 = IELM0 + 1
C
      IELB1 = IELBOR(IELM1,1)
      IELBU = IELBOR(IELMU,1)
      IELBH = IELBOR(IELMH,1)
      IELBT = IELBOR(IELMT,1)
      IELBK = IELBOR(IELMK,1)
      IELBE = IELBOR(IELME,1)
C
      IELMX=MAX(IELMU,IELMH,IELMT,IELMK,IELME)
C
C TYPE DE STOCKAGE ET PRODUIT MATRICE X VECTEUR
C
      CFG(1) = OPTASS
      CFG(2) = PRODUC
C     CFG IMPOSE POUR LES MATRICES DE BORD
      CFGBOR(1) = 1
      CFGBOR(2) = 1
C
C=======================================================================
C
C     ALLOCATION DE LA STRUCTURE DU MAILLAGE
C
      CALL ALMESH(MESH,'MESH  ',IELMX,SPHERI,CFG,T2D_FILES(T2DGEO)%LU,
     *            EQUA,I3=I3,I4=I4)
C     IF ORIGIN COORDINATES IN GEOMETRY FILE AND NOT IN PARAMETER FILE,
C     VALUES OF GEOMETRY FILE ARE TAKEN
      IF(I3.NE.0.AND.I_ORIG.EQ.0) I_ORIG=I3
      IF(I4.NE.0.AND.J_ORIG.EQ.0) J_ORIG=I4
C
C     ALIAS POUR CERTAINES COMPOSANTES DU MAILLAGE MESH
C
      IKLE  => MESH%IKLE
      X     => MESH%X%R
      Y     => MESH%Y%R
C
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
C
C=======================================================================
C
C                     **********************
C                     *   TABLEAUX REELS   *
C                     **********************
C
C-----------------------------------------------------------------------
C
      ALLOCATE(W(NPOIN),STAT=ERR)
      IF(ERR.NE.0) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'POINT_TELEMAC2D : MAUVAISE ALLOCATION DE W'
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'POINT_TELEMAC2D: WRONG ALLOCATION OF W'
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C                       ******************
C                       *   STRUCTURES   *
C                       ******************
C
C-----------------------------------------------------------------------
C
C  ALLOCATION D'UNE STRUCTURE VIDE
C
      CALL ALLVEC(1,S,'S     ',0,1,1)
C
C  TABLEAUX CONTENANT LES VARIABLES QUI FIGURERONT DANS LE RESULTAT:
C
      CALL ALLVEC(1,U,'U     ',IELMU,1,1)
      CALL ALLVEC(1,V,'V     ',IELMU,1,1)
      CALL ALLVEC(1,H,'H     ',IELMH,1,1)
C
C  TABLEAUX CONTENANT LES VARIABLES U,V T, K ET EPSILON CONVECTEES.
C
      CALL ALLVEC(1,UTILD,'UTILD ',IELMU,1,2)
      CALL ALLVEC(1,VTILD,'VTILD ',IELMU,1,2)
      CALL ALLVEC(1,HTILD,'HTILD ',IELMH,1,2)
C
C  TABLEAUX CONTENANT LES VARIABLES U,V,H SAUVEGARDEES AU TEMPS N
C
      CALL ALLVEC(1,UN,'UN    ', IELMU,1,1 )
      CALL ALLVEC(1,VN,'VN    ', IELMU,1,1 )
      CALL ALLVEC(1,HN,'HN    ', IELMH,1,1 )
C
C  TABLEAUX POUR LA SAUVEGARDE DES ACCROISSEMENTS
C
      CALL ALLVEC(1,DH  ,'DH    ' , IELMH ,1,2 )
      IF(IORDRU.EQ.2) THEN
        CALL ALLVEC(1,DU  ,'DU    ' , IELMU , 1,2 )
        CALL ALLVEC(1,DV  ,'DV    ' , IELMU , 1,2 )
      ELSE
        CALL ALLVEC(1,DU  ,'DU    ' , 0     , 1,0 )
        CALL ALLVEC(1,DV  ,'DV    ' , 0     , 1,0 )
      ENDIF
      IF(IORDRH.EQ.2) THEN
        CALL ALLVEC(1,DHN ,'DHN   ' , IELMH , 1,2 )
      ELSE
        CALL ALLVEC(1,DHN ,'DHN   ' , 0     , 1,0 )
      ENDIF
C
C  BLOC DES INCONNUES DANS PROPAG
C
      CALL ALLBLO(UNK,'UNK   ')
      CALL ADDBLO(UNK,DH)
      CALL ADDBLO(UNK, U)
      CALL ADDBLO(UNK, V)
C
C  TABLEAUX DES CONDITIONS AUX LIMITES (TABLEAUX DE BORD)
C  POUR UBOR ET VBOR, DIMENSION 2 POUR SAUVEGARDE DANS LE
C  CAS DE VITESSE OU DEBIT IMPOSES PAR FONCTION
C
      CALL ALLVEC(1,UBOR    ,'UBOR  ',IELBU,2,1)
      CALL ALLVEC(1,VBOR    ,'VBOR  ',IELBU,2,1)
      CALL ALLVEC(1,HBOR    ,'HBOR  ',IELBH,1,1)
      CALL ALLVEC(1,AUBOR   ,'AUBOR ',IELBU,1,1)
      CALL ALLVEC(1,FLBOR   ,'FLBOR ',IELBH,1,1)
      CALL ALLVEC(1,FLBORTRA,'FLBTRA',IELBT,1,1)
C
C  CONDITIONS DE DIRICHLET MISES DANS UN BLOC POUR APPEL DE DIRICH
C
      CALL ALLBLO(DIRBOR,'DIRBOR')
      CALL ADDBLO(DIRBOR,HBOR)
      CALL ADDBLO(DIRBOR,UBOR)
      CALL ADDBLO(DIRBOR,VBOR)
C
C TABLEAU DES COTES DU FOND :
C
      CALL ALLVEC(1,ZF,'ZF    ',IELMH,1,1)
C
C TABLEAU DES COTES DU FOND PAR ELEMENT POUR BANCS DECOUVRANTS
C
      IF(MSK) THEN
        CALL ALLVEC(1,ZFE,'ZFE   ',IELM0,1,1)
      ELSE
        CALL ALLVEC(1,ZFE,'ZFE   ',    0,1,0)
      ENDIF
C
C VISCOSITES : POUR L'INSTANT MIS EN P1
C              MAIS A DEUX DIMENSIONS POUR LE MODELE DE ELDER
C
      IF(ITURB.EQ.2) THEN
        CALL ALLVEC(1,VISC ,'VISC  ',IELM1,3,1)
      ELSE
        CALL ALLVEC(1,VISC ,'VISC  ',IELM1,1,1)
      ENDIF
C     TABLEAU DE SAUVEGARDE DE LA VISCOSITE
      IF(OPDVIT.EQ.2.OR.(NTRAC.GT.0.AND.OPDTRA.EQ.2)) THEN
        IF(ITURB.EQ.2) THEN
          CALL ALLVEC(1,VISC_S,'VISC_S',IELM1,3,1)
        ELSE
          CALL ALLVEC(1,VISC_S,'VISC_S',IELM1,1,1)
        ENDIF
      ENDIF
C
C  COEFFICIENT DE FROTTEMENT
C
      CALL ALLVEC(1,CHESTR,'CHESTR',IELMU,1,1)
C
C  TABLEAUX POUR LES CONDITIONS ATMOSPHERIQUES ET D'ONDE INCIDENTE
C
      CALL ALLVEC(1,C0    ,'C0    ',IELBH,1,1)
      CALL ALLVEC(1,COTOND,'COTOND',IELBH,1,1)
      CALL ALLVEC(1,PATMOS,'PATMOS',IELMH,1,1)
      IF(ROVAR) THEN
        CALL ALLVEC(1,RO,'RO    ',IELMH,1,1)
      ELSE
        CALL ALLVEC(1,RO,'RO    ',    0,1,0)
      ENDIF
C     VENT DONNE EN P1
      IF(VENT) THEN
        CALL ALLVEC(1,WINDX,'WINDX ',IELM1,1,1)
        CALL ALLVEC(1,WINDY,'WINDY ',IELM1,1,1)
      ELSE
        CALL ALLVEC(1,WINDX,'WINDX ',    0,1,0)
        CALL ALLVEC(1,WINDY,'WINDY ',    0,1,0)
      ENDIF
C
C  TABLEAUX DES TERMES SOURCES
C
      CALL ALLVEC(1,FU,'FU    ',IELMU,1,2)
      CALL ALLVEC(1,FV,'FV    ',IELMU,1,2)
C
C  POINTEURS POUR LES MATRICES
C
      IELMHT = IELMH
C     AM1 UTILISEE POUR LES TRACEURS
      IF(NTRAC.GT.0) IELMHT = MAX(IELMHT,IELMT)
      CALL ALLMAT(AM1,'AM1   ',IELMHT,IELMHT,CFG,'Q','Q')
C
      TYP='Q'
      IF(ICONVF(1).NE.ADV_SUP    .AND.
     *   ICONVF(1).NE.ADV_NSC_NC .AND.
     *   3*(SLVPRO%PRECON/3).NE.SLVPRO%PRECON) TYP = 'S'
C
      IF(OPDVIT.EQ.2) TYP='Q'
C
      IELMUT = IELMU
C     AM2 ET AM3 UTILISEES POUR LES TRACEURS
      IF(NTRAC.GT.0) THEN
        IELMUT = MAX(IELMU,IELMT)
        TYP='Q'
      ENDIF
C     AM2 ET AM3 MODIFIEES POUR BOUSSINESQ
      IF(EQUA(1:10).EQ.'BOUSSINESQ') THEN
        TYP='Q'
      ENDIF
      CALL ALLMAT(AM2,'AM2   ',IELMUT,IELMUT,CFG,'Q',TYP)
      CALL ALLMAT(AM3,'AM3   ',IELMUT,IELMUT,CFG,'Q',TYP)
C
C  BM1 ET BM2 :
C
      CALL ALLMAT(BM1,'BM1   ',IELMH,IELMU,CFG,'Q','Q')
      CALL ALLMAT(BM2,'BM2   ',IELMH,IELMU,CFG,'Q','Q')
C
C  SAUVEGARDES DE CV1,BM1 ET BM2 POUR CORRECTION DE CONTINUITE
C
      IF(CORCON.AND.SOLSYS.EQ.1) THEN
        CALL ALLMAT(BM1S,'BM1S  ',IELMH,IELMU,CFG,'Q','Q')
        CALL ALLMAT(BM2S,'BM2S  ',IELMH,IELMU,CFG,'Q','Q')
        CALL ALLVEC(1,CV1S,'CV1S  ',IELMX,1,2)
      ELSE
        CALL ALLMAT(BM1S,'BM1S  ',IELMH,IELMU,CFG,'0','0')
        CALL ALLMAT(BM2S,'BM2S  ',IELMH,IELMU,CFG,'0','0')
        CALL ALLVEC(1,CV1S,'CV1S  ',0,1,0)
      ENDIF
C
C  CM1 ET CM2 :
C
      IELMC1 = IELMH
      IELMC2 = IELMU
C     CM2 UTILISEE POUR U DANS CERTAINS CAS
      IF(ICONVF(1).EQ.ADV_SUP.OR.ICONVF(1).EQ.ADV_NSC_NC) THEN
        IELMC1 = MAX(IELMC1,IELMU)
      ENDIF
      IF(EQUA(1:10).EQ.'BOUSSINESQ') IELMC1 = MAX(IELMC1,IELMU)
C
      CALL ALLMAT(CM1,'CM1   ',IELMC1,IELMC2,CFG,'Q','Q')
      CALL ALLMAT(CM2,'CM2   ',IELMC1,IELMC2,CFG,'Q','Q')
      CALL ALLMAT(TM1,'TM1   ',IELMU ,IELMU ,CFG,'Q','Q')
C
C  MATRICE DE BORD
C
      IELBX = MAX(IELBU,IELBH,IELBT,IELBK,IELBE)
      CALL ALLMAT(MBOR,'MBOR  ',IELBX,IELBX,CFGBOR,'Q','Q')
C
C  MATRICES A23 ET A32 UTILISEES PAR LE PRECONDITIONNEMENT DIAGONAL-BLOC
C  OU POUR LES EQUATIONS DE BOUSSINESQ
C
      TYP = '0'
      IF(3*(SLVPRO%PRECON/3).EQ.SLVPRO%PRECON) TYP = 'Q'
      IF(EQUA(1:10).EQ.'BOUSSINESQ') TYP = 'Q'
      CALL ALLMAT(A23,'A23   ',IELMU,IELMU,CFG,TYP,TYP)
      CALL ALLMAT(A32,'A32   ',IELMU,IELMU,CFG,TYP,TYP)
C
C BLOC DES MATRICES DANS PROPAG
C
      CALL ALLBLO(MAT,'MAT   ')
      CALL ADDBLO(MAT,AM1)
      CALL ADDBLO(MAT,BM1)
      CALL ADDBLO(MAT,BM2)
      CALL ADDBLO(MAT,CM1)
      CALL ADDBLO(MAT,AM2)
      CALL ADDBLO(MAT,A23)
      CALL ADDBLO(MAT,CM2)
      CALL ADDBLO(MAT,A32)
      CALL ADDBLO(MAT,AM3)
C
C TABLEAU DE TRAVAIL W1 (DIMENSIONNEMENT A VERIFIER)
C
C     MEMOIRE NECESSAIRE POUR W1 DANS VALIDA
      MEMW1 = 9*NPOIN
C     VOLUMES FINIS
      IF(EQUA(1:15).EQ.'SAINT-VENANT VF') THEN
        MEMW1 = MAX(MEMW1,9*NPOIN+3*NPTFR,2*MXPTVS*NPOIN)
      ENDIF
C     CETTE PLACE MEMOIRE EST RESERVEE SOUS FORME D'UN
C     TABLEAU P0 AVEC UNE DEUXIEME DIMENSION
      MEMW1 = MAX(3,1+MEMW1/NBMPTS(IELM0))
      CALL ALLVEC(1,W1,'W1    ',IELM0,MEMW1,1)
C
C_______________________________________________________________________
C
C  POINTEURS POUR LES SECONDS MEMBRES DE L'ETAPE DE PROPAGATION
C
C_______________________________________________________________________
C
      CALL ALLVEC(1,CV1,'CV1   ',IELMX,1,2)
      CALL ALLVEC(1,CV2,'CV2   ',IELMU,1,2)
      CALL ALLVEC(1,CV3,'CV3   ',IELMU,1,2)
C
C  BLOC DES SECONDS MEMBRES DANS PROPAG
C
      CALL ALLBLO(RHS,'RHS   ')
      CALL ADDBLO(RHS,CV1)
      CALL ADDBLO(RHS,CV2)
      CALL ADDBLO(RHS,CV3)
C_______________________________________________________________________
C
C  POINTEURS POUR LES TERMES SOURCES DE L'ETAPE DE PROPAGATION
C
C_______________________________________________________________________
C
      CALL ALLVEC(1,SMH,'SMH   ',IELMX,1,2)
C_______________________________________________________________________
C
C  POINTEURS POUR CONVECTEUR ET PROPAGATEUR
C_______________________________________________________________________
C
      CALL ALLVEC(1,UCONV,'UCONV ',IELMU,1,1)
      CALL ALLVEC(1,VCONV,'VCONV ',IELMU,1,1)
      CALL ALLVEC(1,HPROP,'HPROP ',IELMH,1,1)
C_______________________________________________________________________
C
C  POINTEURS POUR INTEGRALE DES BASES, EN PARALLELE, ET INVERSE
C_______________________________________________________________________
C
      CALL ALLVEC(1,VOLU2D,'VOLU2D',IELMH,1,1)
      CALL ALLVEC(1,V2DPAR,'V2DPAR',IELMH,1,1)
      CALL ALLVEC(1,UNSV2D,'UNSV2D',IELMH,1,1)
C_______________________________________________________________________
C
C  POINTEURS UTILISES POUR LES DERIVES LAGRANGIENNES
C_______________________________________________________________________
C
      CALL ALLVEC(1,XLAG  ,'XLAG  ',NPOIN*NLAG             ,1,0)
      CALL ALLVEC(1,YLAG  ,'YLAG  ',NPOIN*NLAG             ,1,0)
      CALL ALLVEC(1,SHPLAG,'SHPLAG',NPOIN*NBPEL(IELM1)*NLAG,1,0)
C
C-----------------------------------------------------------------------
C
C  POINTEURS DES TABLEAUX DE TRAVAIL:
C
C-----------------------------------------------------------------------
C
C  NOMBRE DE TABLEAUX A ALLOUER : NTR
C           21 : POUR CGSTAB =3 X 7, 22 POUR CVDFTR (APPEL DE CVTRVF)
      NTR = 22
      IF(SLVPRO%SLV.EQ.7) NTR = MAX(NTR,6+6*SLVPRO%KRYLOV)
C     6 DIAGONALES DE PLUS A STOCKER EN PRECONDITIONNEMENT BLOC-DIAGONAL
      IF(3*(SLVPRO%PRECON/3).EQ.SLVPRO%PRECON) NTR = NTR + 6
C
C  TAILLE MAXIMUM UTILE
C
      NTRT=0
      IF(NTRAC.GT.0) THEN
C       NTRT = 7
C       A CAUSE DE LA POSITION DES TRACEURS DANS VARSOR QUI SERA
C       LA MEME DANS TB (UTILISE PAR VALIDA)
        NTRT = 31+NTRAC
        IF(SLVTRA%SLV.EQ.7) NTRT = MAX(2+2*SLVTRA%KRYLOV,NTRT)
        NTR = MAX(NTR,NTRT)
      ENDIF
      NTRKE=0
      IF(ITURB.EQ.3) THEN
        NTRKE=7
        IF(SLVK%SLV.EQ.7) NTRKE = MAX(NTRKE,2+2*SLVK%KRYLOV)
        NTR  = MAX(NTR,NTRKE)
      ENDIF
C
C  ALLOCATIONS : TABLEAUX DE TRAVAIL DE DIMENSION LE NOMBRE MAXIMUM
C                DE DEGRES DE LIBERTE, AU NOMBRE DE NTR
C
C     TB CONTIENDRA LES TABLEAUX T1,T2,...
C
      CALL ALLBLO(TB ,'TB    ')
C
      CALL ALLVEC_IN_BLOCK(TB,NTR,1,'TB    ',IELMX,1,2)
C
C     ALIAS POUR LES 12 PREMIERS TABLEAUX DE TRAVAIL DU BLOC TB
C
      T1 =>TB%ADR( 1)%P
      T2 =>TB%ADR( 2)%P
      T3 =>TB%ADR( 3)%P
      T4 =>TB%ADR( 4)%P
      T5 =>TB%ADR( 5)%P
      T6 =>TB%ADR( 6)%P
      T7 =>TB%ADR( 7)%P
      T8 =>TB%ADR( 8)%P
      T9 =>TB%ADR( 9)%P
      T10=>TB%ADR(10)%P
      T11=>TB%ADR(11)%P
      T12=>TB%ADR(12)%P
C
C  ALLOCATIONS : TABLEAUX DE TRAVAIL DE DIMENSION LE NOMBRE MAXIMUM
C                D'ELEMENTS.
C
      CALL ALLVEC(1,TE1,'TE1   ',IELM0,1,1)
      CALL ALLVEC(1,TE2,'TE2   ',IELM0,1,1)
      CALL ALLVEC(1,TE3,'TE3   ',IELM0,1,1)
      IF(OPTBAN.EQ.1.OR.OPTBAN.EQ.3) THEN
C       PIECE-WISE LINEAR FREE SURFACE
        CALL ALLVEC(1,ZFLATS, 'ZFLATS',IELM0,3,1)
      ELSE
        CALL ALLVEC(1,ZFLATS, 'ZFLATS',    0,1,1)
      ENDIF
      IF(OPTBAN.EQ.3) THEN
        CALL ALLVEC(1,TE4,'TE4   ',IELM0,1,1)
        CALL ALLVEC(1,TE5,'TE5   ',IELM0,1,1)
      ELSE
        CALL ALLVEC(1,TE4,'TE4   ',    0,1,0)
        CALL ALLVEC(1,TE5,'TE5   ',    0,1,0)
      ENDIF
C
C-----------------------------------------------------------------------
C !jaj #### if required, in this place we read the input sections file 
C      and modify NCP and CTRLSC(1:NCP) accordingly in read_sections 
C
      IF (TRIM(T2D_FILES(T2DSEC)%NAME)/='') THEN 
        WRITE(LU,*) 
     &   'POINT_TELEMAC2D: SECTIONS DEFINED IN THE SECTIONS INPUT FILE'
        CALL READ_SECTIONS_TELEMAC2D
      ELSE ! the previously existing way of doing things 
        IF (NCP.NE.0) WRITE(LU,*) 
     &   'POINT_TELEMAC2D: SECTIONS DEFINED IN THE PARAMETER FILE'
      ENDIF
C
C     BLOCK OF MASKS FOR THE COMPUTATION OF FLUXES ACCROSS SECTIONS
C     ONLY WITH COMPATIBLE FLUXES
C
      CALL ALLBLO(MSKSEC,'MSKSEC')
      IF(NCP.GT.1.AND.COMFLU) THEN
        CALL ALLVEC_IN_BLOCK(MSKSEC,NCP/2,1,'MSKS  ',IELM0,1,1)
      ENDIF
C
C-----------------------------------------------------------------------
C
C POINTEURS DES MASQUES
C
C     BLOC DES MASQUES DES CONDITIONS AUX LIMITES
C     POUR LA PROPAGATION.
C
      CALL ALLBLO(MASK,'MASK  ')
      CALL ALLVEC_IN_BLOCK(MASK,11,1,'MASK  ',IELBH,1,2)
C
      IF(MSK) THEN
        CALL ALLVEC(1,MASKEL,'MASKEL',IELM0,1,1)
        CALL ALLVEC(1,MASKPT,'MASKPT',IELMX,1,1)
      ELSE
        CALL ALLVEC(1,MASKEL,'MASKEL',    0,1,0)
        CALL ALLVEC(1,MASKPT,'MASKPT',    0,1,0)
      ENDIF
C
C  TABLEAU AJOUTES SI IL Y A DES TRACEURS
C
      CALL ALLBLO(T      ,'T     ')
      CALL ALLBLO(TTILD  ,'TTILD ')
      CALL ALLBLO(TN     ,'TN    ')
      CALL ALLBLO(TEXP   ,'TEXP  ')
      CALL ALLBLO(TIMP   ,'TIMP  ')
      CALL ALLBLO(TSCEXP ,'TSCEXP')
      CALL ALLBLO(VISCT  ,'VISCT ')
      CALL ALLBLO(MASKTR ,'MASKTR')
      CALL ALLBLO(TBOR   ,'TBOR  ')
      CALL ALLBLO(ATBOR  ,'ATBOR ')
      CALL ALLBLO(BTBOR  ,'BTBOR ')
      CALL ALLBLO(LITBOR ,'LITBOR')
      IF(NTRAC.GT.0) THEN
        CALL ALLVEC_IN_BLOCK(T     ,NTRAC,1,'T     ',IELMT,1,1)
        CALL ALLVEC_IN_BLOCK(TTILD ,NTRAC,1,'TTILD ',IELMT,1,1)
        CALL ALLVEC_IN_BLOCK(TN    ,NTRAC,1,'TN    ',IELMT,1,1)
        CALL ALLVEC_IN_BLOCK(TEXP  ,NTRAC,1,'TEXP  ',IELMT,1,1)
        CALL ALLVEC_IN_BLOCK(TIMP  ,NTRAC,1,'TIMP  ',IELMT,1,1)
        CALL ALLVEC_IN_BLOCK(TSCEXP,NTRAC,1,'TSCEXP',IELMT,1,1)
        IF(ITURB.EQ.2) THEN
          CALL ALLVEC_IN_BLOCK(VISCT,NTRAC,1,'VISCT ',IELMT,3,1)
        ELSE
          CALL ALLVEC_IN_BLOCK(VISCT,NTRAC,1,'VISCT ',IELMT,1,1)
        ENDIF
        CALL ALLVEC_IN_BLOCK(MASKTR,4,1,'MSKTR ',IELBH,1,2)
        IF(THOMFR) THEN
C         DEUXIEME DIMENSION UTILISEE COMME TABLEAU DE TRAVAIL
C         DANS THOMPS
          CALL ALLVEC_IN_BLOCK(TBOR,NTRAC,1,'TBOR  ',IELBT,2,1)
        ELSE
          CALL ALLVEC_IN_BLOCK(TBOR,NTRAC,1,'TBOR  ',IELBT,1,1)
        ENDIF
        CALL ALLVEC_IN_BLOCK(ATBOR  ,NTRAC,1,'ATBOR ',IELBT,1,1)
        CALL ALLVEC_IN_BLOCK(BTBOR  ,NTRAC,1,'BTBOR ',IELBT,1,1)
        CALL ALLVEC_IN_BLOCK(LITBOR ,NTRAC,2,'LITBOR',IELBT,1,1)
      ELSE
C       AT LEAST ONE ELEMENT IN BLOCKS, NOT NTRAC
        CALL ALLVEC_IN_BLOCK(T     ,1,1,'T     ',0,1,0)
        CALL ALLVEC_IN_BLOCK(TTILD ,1,1,'TTILD ',0,1,0)
        CALL ALLVEC_IN_BLOCK(TN    ,1,1,'TN    ',0,1,0)
        CALL ALLVEC_IN_BLOCK(TEXP  ,1,1,'TEXP  ',0,1,0)
        CALL ALLVEC_IN_BLOCK(TIMP  ,1,1,'TIMP  ',0,1,0)
        CALL ALLVEC_IN_BLOCK(TSCEXP,1,1,'TSCEXP',0,1,0)
        IF(ITURB.EQ.2) THEN
          CALL ALLVEC_IN_BLOCK(VISCT,1,1,'VISCT ',0,3,0)
        ELSE
          CALL ALLVEC_IN_BLOCK(VISCT,1,1,'VISCT ',0,1,0)
        ENDIF
        CALL ALLVEC_IN_BLOCK(MASKTR,4,1,'MSKTR ',0,1,0)
        CALL ALLVEC_IN_BLOCK(TBOR   ,1,1,'TBOR  ',0,1,0)
        CALL ALLVEC_IN_BLOCK(ATBOR  ,1,1,'ATBOR ',0,1,0)
        CALL ALLVEC_IN_BLOCK(BTBOR  ,1,1,'BTBOR ',0,1,0)
        CALL ALLVEC_IN_BLOCK(LITBOR ,1,2,'LITBOR',0,1,0)
      ENDIF
C
C-----------------------------------------------------------------------
C
C  COEFFICIENT DE FROTTEMENT CF
C
      CALL ALLVEC(1,CF    ,'CF    ',IELMU,1,1)
C
C  DATA FOR FRICTION SET PER ZONE
C
C  FRICTION LAW USED
C
      CALL ALLVEC(2,NKFROT,'NKFROT',IELMU,1,1)
C
C  CHESTR SUR LE BORD
C
      CALL ALLVEC(1,CHBORD,'CHBORD',IELBT,1,1)
C
      IF(FRICTB) THEN
         ALLOCATE(FRTAB%ADR(NZONMX))
         DO I=1,NZONMX
            ALLOCATE(FRTAB%ADR(I)%P)
         ENDDO
         CALL ALLVEC(2,KFROPT,'KFROPT',IELMU,1,1)
         CALL ALLVEC(1,NDEFMA,'NDEFMA',IELMU,1,1)
         IF(LINDNER) THEN
           CALL ALLVEC(1,LINDDP,'LINDDP',IELMU,1,1)
           CALL ALLVEC(1,LINDSP,'LINDSP',IELMU,1,1)
         ELSE
           CALL ALLVEC(1,LINDDP,'LINDDP',0,1,0)
           CALL ALLVEC(1,LINDSP,'LINDSP',0,1,0)
         ENDIF
         CALL ALLVEC(1,NDEF_B,'NDEF_B',IELBT,1,1)
         CALL ALLVEC(1,KFRO_B,'KFRO_B',IELBT,1,1)
      ELSE
         CALL ALLVEC(2,KFROPT,'KFROPT',0,1,0)
         CALL ALLVEC(1,NDEFMA,'NDEFMA',0,1,0)
         CALL ALLVEC(1,LINDDP,'LINDDP',0,1,0)
         CALL ALLVEC(1,LINDSP,'LINDSP',0,1,0)
         CALL ALLVEC(1,NDEF_B,'NDEF_B',0,1,0)
         CALL ALLVEC(1,KFRO_B,'KFRO_B',0,1,0)
      ENDIF
C
C  END OF DATA FOR FRICTION SET PER ZONE
C
C  TABLEAU AJOUTES EN CAS DE MODELE K-EPSILON
C
      IF(ITURB.EQ.3) THEN
        CALL ALLVEC(1,AK     ,'AK    ',IELMK,1,1)
        CALL ALLVEC(1,EP     ,'EP    ',IELME,1,1)
        CALL ALLVEC(1,AKN    ,'AKN   ',IELMK,1,1)
        CALL ALLVEC(1,EPN    ,'EPN   ',IELME,1,1)
        CALL ALLVEC(1,AKTILD ,'AKTILD',IELMK,1,1)
        CALL ALLVEC(1,EPTILD ,'EPTILD',IELME,1,1)
        CALL ALLVEC(1,KBOR   ,'KBOR  ',IELBK,1,1)
        CALL ALLVEC(1,EBOR   ,'EBOR  ',IELBE,1,1)
        CALL ALLVEC(1,CFBOR  ,'CFBOR ',IELB1,1,1)
      ELSE
        CALL ALLVEC(1,AK     ,'AK    ',0,1,0)
        CALL ALLVEC(1,EP     ,'EP    ',0,1,0)
        CALL ALLVEC(1,AKN    ,'AKN   ',0,1,0)
        CALL ALLVEC(1,EPN    ,'EPN   ',0,1,0)
        CALL ALLVEC(1,AKTILD ,'AKTILD',0,1,0)
        CALL ALLVEC(1,EPTILD ,'EPTILD',0,1,0)
        CALL ALLVEC(1,KBOR   ,'KBOR  ',0,1,0)
        CALL ALLVEC(1,EBOR   ,'EBOR  ',0,1,0)
        CALL ALLVEC(1,CFBOR  ,'CFBOR ',0,1,0)
      ENDIF
C
      CALL ALLVEC(1,UDEL   ,'UDEL  ',    IELMU,1,1)
      CALL ALLVEC(1,VDEL   ,'VDEL  ',    IELMU,1,1)
      CALL ALLVEC(1,DM1    ,'DM1   ',    IELMU,1,2)
      CALL ALLVEC(1,ZCONV  ,'ZCONV ',       10,3,1)
      CALL ALLVEC(1,FLODEL ,'FLODEL',MESH%NSEG,1,0)
      CALL ALLVEC(1,FLULIM ,'FLULIM',MESH%NSEG,1,0)
C
C-----------------------------------------------------------------------
C
C ALLOCATION DE BLOCS
C
C     FONCTIONS A CONVECTER PAR CARACTERISTIQUES
C
      CALL ALLBLO(FN    , 'FN    ')
      CALL ALLBLO(F     , 'F     ')
      CALL ALLBLO(FTILD , 'FTILD ')
      CALL ALLBLO(FNCAR , 'FNCAR ')
C
      CALL ADDBLO(FN,UN)
      CALL ADDBLO(FN,VN)
      CALL ADDBLO(FN,HN)
      CALL ADDBLO(F ,U )
      CALL ADDBLO(F ,V )
      CALL ADDBLO(F ,H )
      IF(NTRAC.GT.0) THEN
        DO ITRAC=1,NTRAC
          CALL ADDBLO(FN   ,TN%ADR(ITRAC)%P   )
          CALL ADDBLO(F    ,T%ADR(ITRAC)%P    )
        ENDDO
      ENDIF
      IF(ITURB.EQ.3) THEN
        CALL ADDBLO(FN    ,AKN)
        CALL ADDBLO(FN    ,EPN)
        CALL ADDBLO(F     ,AK )
        CALL ADDBLO(F     ,EP )
      ENDIF
C
C-----------------------------------------------------------------------
C
C     WITH FINITE VOLUMES OR KINETIC SCHEMES ADVECTION IS DONE
C     IN VOLFIN
C
      IF(EQUA(1:15).NE.'SAINT-VENANT VF') THEN
        IF(CONVV(1).AND.ICONVF(1).EQ.ADV_CAR) THEN
          CALL ADDBLO(FTILD,UTILD)
          CALL ADDBLO(FTILD,VTILD)
          CALL ADDBLO(FNCAR,UN   )
          CALL ADDBLO(FNCAR,VN   )
        ENDIF
        IF(CONVV(3).AND.NTRAC.GT.0.AND.ICONVF(3).EQ.ADV_CAR) THEN
          DO ITRAC=1,NTRAC
            CALL ADDBLO(FTILD,TTILD%ADR(ITRAC)%P)
            CALL ADDBLO(FNCAR,TN%ADR(ITRAC)%P)
          ENDDO
        ENDIF
      ENDIF
C
      IF(CONVV(4).AND.ITURB.EQ.3.AND.ICONVF(4).EQ.ADV_CAR) THEN
        CALL ADDBLO(FTILD,AKTILD)
        CALL ADDBLO(FTILD,EPTILD)
        CALL ADDBLO(FNCAR,AKN   )
        CALL ADDBLO(FNCAR,EPN   )
      ENDIF
C
C_______________________________________________________________________
C
C  TABLEAUX UTILISES POUR LE SUIVI DE FLOTTEURS
C
C_______________________________________________________________________
C
      IF(NFLOT.NE.0) THEN
        CALL ALLVEC(1,XFLOT ,'XFLOT ',NITFLO*NFLOT      ,1,0)
        CALL ALLVEC(1,YFLOT ,'YFLOT ',NITFLO*NFLOT      ,1,0)
        CALL ALLVEC(1,SHPFLO,'SHPFLO',NBPEL(IELM1)*NFLOT,1,0)
      ELSE
        CALL ALLVEC(1,XFLOT ,'XFLOT ',0,1,0)
        CALL ALLVEC(1,YFLOT ,'YFLOT ',0,1,0)
        CALL ALLVEC(1,SHPFLO,'SHPFLO',0,1,0)
      ENDIF
C
C_______________________________________________________________________
C
C  TABLEAUX UTILISES POUR LES SEUILS.
C
C-----------------------------------------------------------------------
C
      IF(NWEIRS.NE.0) THEN
C       EN FAIT TABLEAUX (NWEIRS,NPSMAX) OR NPOIN>NWEIRS*NPSMAX
        CALL ALLVEC(1,ZDIG  ,'ZDIG  ',IELM1,1,1)
        CALL ALLVEC(1,PHIDIG,'PHIDIG',IELM1,1,1)
      ELSE
        CALL ALLVEC(1,ZDIG  ,'ZDIG  ',0,1,0)
        CALL ALLVEC(1,PHIDIG,'PHIDIG',0,1,0)
      ENDIF
C
C-----------------------------------------------------------------------
C
C  TABLEAUX A LA DISPOSITION DE L'UTILISATEUR
C
      CALL ALLBLO(PRIVE ,'PRIVE ')
C
      IF(NPRIV.GT.0) THEN
C       CES TABLEAUX DOIVENT EXISTER MAIS PEUVENT ETRE VIDES
        CALL ALLVEC_IN_BLOCK(PRIVE,NPRIV,1,'PRIV  ',IELMX,1,2)
      ENDIF
C     IL FAUT AU MOINS 4 TABLEAUX MAIS ILS PEUVENT ETRE VIDES
      IF(NPRIV.LT.4) THEN
        CALL ALLVEC_IN_BLOCK(PRIVE,4-NPRIV,1,'PRIV  ',    0,1,2)
      ENDIF
C
C     ALIAS POUR LES 4 PREMIERS TABLEAUX PRIVES
C
      PRIVE1 => PRIVE%ADR(1)%P%R
      PRIVE2 => PRIVE%ADR(2)%P%R
      PRIVE3 => PRIVE%ADR(3)%P%R
      PRIVE4 => PRIVE%ADR(4)%P%R
C
C  BLOC DES VARIABLES CLANDESTINES
C
      CALL ALLBLO(VARCL,'VARCL ')
      CALL ALLVEC_IN_BLOCK(VARCL,NVARCL,1,'CL    ',IELMX,1,2)
C
C     INITIALISATION A 0 POUR LES CAS AVEC SORTIE DES CONDITIONS
C     INITIALES
C
      IF(NVARCL.GT.0) THEN
        DO I=1,NVARCL
          CALL OS('X=C     ',VARCL%ADR(I)%P,VARCL%ADR(I)%P,
     *                       VARCL%ADR(I)%P,0.D0)                
        ENDDO
      ENDIF
C
C_______________________________________________________________________
C
C                         * TABLEAUX ENTIERS *
C_______________________________________________________________________
C
      IF(MSK) THEN
        CALL ALLVEC(2,IFAMAS,'IFAMAS',IELM0,NBFEL(IELM0),1)
      ELSE
        CALL ALLVEC(2,IFAMAS,'IFAMAS',0,1,0)
      ENDIF
      CALL ALLVEC(2,LIUBOR,'LIUBOR',IELBU,1,1)
      CALL ALLVEC(2,LIVBOR,'LIVBOR',IELBU,1,1)
      CALL ALLVEC(2,LIHBOR,'LIHBOR',IELBH,1,1)
C     CLU,CLV ET CLH SONT DES TABLEAUX DE TRAVAIL DANS PROPIN
      CALL ALLVEC(2,CLU            ,'CLU   ',IELBU,1,1)
      CALL ALLVEC(2,CLV            ,'CLV   ',IELBU,1,1)
      CALL ALLVEC(2,CLH            ,'CLH   ',IELBH,1,1)
      CALL ALLVEC(2,BOUNDARY_COLOUR,'BNDCOL',IELB1,1,1)
C
      CALL ALLVEC(2,NUMLIQ,'NUMLIQ',IELB1,1,1)
      IF(ITURB.EQ.3) THEN
        CALL ALLVEC(2,LIMKEP,'LIMKEP',IELB1,2,1)
      ELSE
        CALL ALLVEC(2,LIMKEP,'LIMKEP',    0,2,0)
      ENDIF
      CALL ALLVEC(2,LIMPRO,'LIMPRO',MAX(IELBH,IELBU),6,1)
      CALL ALLVEC(2,LIMTRA,'LIMTRA',IELBT,1,1)
      CALL ALLVEC(2,SECMOU,'SECMOU',IELM0,1,1)
C
C     TABLEAU DE TRAVAIL ENTIER (DE TAILLE MINIMUM NELEM)
C
      IF(IELMX.GT.11) THEN
        CALL ALLVEC(2,IT1,'IT1   ',IELMX,1,2)
        CALL ALLVEC(2,IT2,'IT2   ',IELMX,1,2)
        CALL ALLVEC(2,IT3,'IT3   ',IELMX,1,2)
        CALL ALLVEC(2,IT4,'IT4   ',IELMX,1,2)
      ELSE
        CALL ALLVEC(2,IT1,'IT1   ',   10,1,2)
        CALL ALLVEC(2,IT2,'IT2   ',   10,1,2)
        CALL ALLVEC(2,IT3,'IT3   ',   10,1,2)
        CALL ALLVEC(2,IT4,'IT4   ',   10,1,2)
      ENDIF
C
C_______________________________________________________________________
C
C  TABLEAUX UTILISES POUR LE SUIVI DE FLOTTEURS
C
C_______________________________________________________________________
C
C     PAS DE TEST SUR NFLOT, SI IL N'Y A PAS DE FLOTTEUR
C     SI NFLOT EST NUL, LES VECTEURS SERONT A TAILLE NULLE
      CALL ALLVEC(2,DEBFLO,'DEBFLO',NFLOT         ,1,0)
      CALL ALLVEC(2,FINFLO,'FINFLO',NFLOT         ,1,0)
      CALL ALLVEC(2,ELTFLO,'ELTFLO',NFLOT         ,1,0)
      CALL ALLVEC(2,IKLFLO,'IKLFLO',NFLOT*NITFLO*3,1,0)
C
C_______________________________________________________________________
C
C  TABLEAUX UTILISES POUR LES DERIVES LAGRANGIENNES
C
C-----------------------------------------------------------------------
C
      IF(NLAG.NE.0) THEN
        CALL ALLVEC(2,DEBLAG,'DEBLAG',NLAG      ,1,0)
        CALL ALLVEC(2,FINLAG,'FINLAG',NLAG      ,1,0)
        CALL ALLVEC(2,ELTLAG,'ELTLAG',NLAG*NPOIN,1,0)
      ELSE
        CALL ALLVEC(2,DEBLAG,'DEBLAG',0         ,1,0)
        CALL ALLVEC(2,FINLAG,'FINLAG',0         ,1,0)
        CALL ALLVEC(2,ELTLAG,'ELTLAG',0         ,1,0)
      ENDIF
C
C_______________________________________________________________________
C
C  TABLEAUX UTILISES POUR LES SEUILS.
C
C-----------------------------------------------------------------------
C
C     NUMDIG A EN FAIT LA TAILLE NUMDIG(2,NWEIRS,NPSMAX)
C     NPOIN EST UN MAJORANT DE NWEIRS * NPSMAX, NOMBRES QUI SONT
C     LUS DANS LES FICHIERS DES SINGULARITES
      IF(NWEIRS.NE.0) THEN
        CALL ALLVEC(2,NUMDIG,'NUMDIG',2*NPOIN,1,0)
      ELSE
        CALL ALLVEC(2,NUMDIG,'NUMDIG',    0  ,1,0)
      ENDIF
C
C_______________________________________________________________________
C
C  TABLEAU UTILISE POUR LES NUMEROS DE ZONES.
C
C-----------------------------------------------------------------------
C
      IF(DEFZON) THEN
        CALL ALLVEC(2,ZONE,'ZONE  ',IELM1,1,1)
      ELSE
        CALL ALLVEC(2,ZONE,'ZONE  ',0    ,1,0)
      ENDIF
C
C_______________________________________________________________________
C
C  TABLEAUX NON COMMUNS A TOUS LES TYPES D'EQUATIONS RESOLUES
C_______________________________________________________________________
C
      CALL ALLBLO(SMTR     ,'SMTR  ')
      CALL ALLBLO(FLUXT    ,'FLUXT ')
      CALL ALLBLO(FLUXTEMP ,'FLUXTE')
      CALL ALLBLO(FLUHBTEMP,'FLUHBT')
      CALL ALLBLO(FLUHBOR  ,'FLUHB ')
      CALL ALLBLO(HT       ,'HT    ')
      IF(EQUA(1:15).EQ.'SAINT-VENANT VF') THEN
        CALL ALLVEC(1,QU       ,'QU    ',IELM1    ,1        ,1)
        CALL ALLVEC(1,QV       ,'QV    ',IELM1    ,1        ,1)
        CALL ALLVEC(1,HSTOK    ,'HSTOK ',IELM1    ,1        ,1)
        CALL ALLVEC(1,HCSTOK   ,'HCSTOK',2        ,MESH%NSEG,0)
        CALL ALLVEC(1,SMTR     ,'SMTR  ',IELM1    ,1        ,1)
        CALL ALLVEC(2,LOGFR    ,'LOGFR ',IELM1    ,1        ,1)
        CALL ALLVEC(1,HC       ,'HC    ',2        ,MESH%NSEG,0)
        CALL ALLVEC(1,DSZ      ,'DSZ   ',2        ,MESH%NSEG,0)
        IF(NTRAC.GT.0) THEN
          CALL ALLVEC_IN_BLOCK(FLUXT    ,NTRAC,1,'FLUXT ',MESH%NSEG,1,0)
          CALL ALLVEC_IN_BLOCK(FLUXTEMP ,NTRAC,1,'FLUXTE',MESH%NSEG,1,0)
          CALL ALLVEC_IN_BLOCK(FLUHBTEMP,NTRAC,1,'FLUHBT',IELBH    ,1,1)
          CALL ALLVEC_IN_BLOCK(FLUHBOR  ,NTRAC,1,'FLUHB ',IELBH    ,1,1)
          CALL ALLVEC_IN_BLOCK(HT       ,NTRAC,1,'HT    ',IELM1    ,1,1)
          CALL ALLVEC_IN_BLOCK(SMTR     ,NTRAC,1,'SMTR  ',IELM1    ,1,1)
        ELSE
          CALL ALLVEC_IN_BLOCK(FLUXT    ,1,1,'FLUXT ',0,1,0)
          CALL ALLVEC_IN_BLOCK(FLUXTEMP ,1,1,'FLUXTE',0,1,0)
          CALL ALLVEC_IN_BLOCK(FLUHBTEMP,1,1,'FLUHBT',0,1,1)
          CALL ALLVEC_IN_BLOCK(FLUHBOR  ,1,1,'FLUHB ',0,1,1)
          CALL ALLVEC_IN_BLOCK(HT       ,1,1,'HT    ',0,1,1) 
          CALL ALLVEC_IN_BLOCK(SMTR     ,1,1,'SMTR  ',0,1,1)        
        ENDIF
      ELSE
        CALL ALLVEC(1,QU       ,'QU    ',0 , 1,0)
        CALL ALLVEC(1,QV       ,'QV    ',0 , 1,0)
        CALL ALLVEC(1,HSTOK    ,'HSTOK ',0 , 1,0)
        CALL ALLVEC(1,HCSTOK   ,'HCSTOK',0 , 1,0)
        CALL ALLVEC(2,LOGFR    ,'LOGFR ',0 , 1,0)
        CALL ALLVEC(1,HC       ,'HC    ',0 , 1,0)
        CALL ALLVEC(1,DSZ      ,'DSZ   ',0 , 1,0)
        CALL ALLVEC_IN_BLOCK(FLUXT    ,1,1,'FLUXT ',0,1,0)
        CALL ALLVEC_IN_BLOCK(FLUXTEMP ,1,1,'FLUXTE',0,1,0)
        CALL ALLVEC_IN_BLOCK(FLUHBTEMP,1,1,'FLUHBT',0,1,1)
        CALL ALLVEC_IN_BLOCK(FLUHBOR  ,1,1,'FLUHB ',0,1,1)
        CALL ALLVEC_IN_BLOCK(HT       ,1,1,'HT    ',0,1,1) 
        CALL ALLVEC_IN_BLOCK(SMTR     ,1,1,'SMTR  ',0,1,1)  
      ENDIF
C
      IF(EQUA(1:10).EQ.'BOUSSINESQ') THEN
        CALL ALLVEC(1,H0  ,'H0    ',IELMH,1,1 )
      ELSE
        CALL ALLVEC(1,H0  ,'H0    ',0 , 1,0 )
      ENDIF
C
C-----------------------------------------------------------------------
C
C    POUR LES SURFACES LIBRES MAX, LES VITESSES MAX
C    ET LES TEMPS CORRESPONDANTS
C
      IF(SORLEO(27).OR.SORIMP(27)) THEN
        CALL ALLVEC(1,MAXZ,'MAXZ  ',IELM1,1,1 )
      ELSE
        CALL ALLVEC(1,MAXZ,'MAXZ  ',0    ,1,0 )      
      ENDIF
      IF(SORLEO(28).OR.SORIMP(28)) THEN
        CALL ALLVEC(1,TMAXZ,'TMAXZ ',IELM1,1,1 )
      ELSE
        CALL ALLVEC(1,TMAXZ,'TMAXZ ',0    ,1,0 )     
      ENDIF
      IF(SORLEO(29).OR.SORIMP(29)) THEN
        CALL ALLVEC(1,MAXV,'MAXV  ',IELM1,1,1 )
      ELSE
        CALL ALLVEC(1,MAXV,'MAXV  ',0    ,1,0 )      
      ENDIF
      IF(SORLEO(30).OR.SORIMP(30)) THEN
        CALL ALLVEC(1,TMAXV,'TMAXV ',IELM1,1,1 )
      ELSE
        CALL ALLVEC(1,TMAXV,'TMAXV ',0    ,1,0 )     
      ENDIF
C
C    POUR LES ANALYSES DE FOURIER
C
      CALL ALLBLO(AMPL,'AMPL  ')
      CALL ALLBLO(PHAS,'PHAS  ')
      IF(NPERIAF.GT.0) THEN
        CALL ALLVEC_IN_BLOCK(AMPL,NPERIAF,1,'AMPL  ',IELM1,1,2)
        CALL ALLVEC_IN_BLOCK(PHAS,NPERIAF,1,'PHAS  ',IELM1,1,2)
      ENDIF
C   
C-----------------------------------------------------------------------
C
C CONSTRUCTION DU BLOC QUI PERMET DE RELIER UN NOM DE VARIABLE
C A SON TABLEAU
C
      CALL ALLBLO(VARSOR ,'VARSOR')   
C 01
      CALL ADDBLO(VARSOR,U)
C 02
      CALL ADDBLO(VARSOR,V)
C 03
      CALL ADDBLO(VARSOR,FU)
C 04
      CALL ADDBLO(VARSOR,H)
C 05
      CALL ADDBLO(VARSOR,FV)
C 06
      CALL ADDBLO(VARSOR,ZF)
C 07
      CALL ADDBLO(VARSOR,T2)
C 08
      CALL ADDBLO(VARSOR,T3)
C 09  ANCIEN TRACEUR
C     REPETE ICI MAIS NON UTILISE CAR REMIS PLUS LOIN
      CALL ADDBLO(VARSOR,T%ADR(1)%P)
C 10
      CALL ADDBLO(VARSOR,AK)
C 11
      CALL ADDBLO(VARSOR,EP)
C 12
      CALL ADDBLO(VARSOR,VISC)
C 13 
      CALL ADDBLO(VARSOR,T4)
C 14 
      CALL ADDBLO(VARSOR,T5)
C 15 
      CALL ADDBLO(VARSOR,T6)
C 16 
      CALL ADDBLO(VARSOR,WINDX)
C 17 
      CALL ADDBLO(VARSOR,WINDY)
C 18 
      CALL ADDBLO(VARSOR,PATMOS)
C 19 
      CALL ADDBLO(VARSOR,CHESTR)
C 20 
      CALL ADDBLO(VARSOR,T7)
C 21 
      CALL ADDBLO(VARSOR,T8)
C 22 
      CALL ADDBLO(VARSOR,T9)
C 23 
      CALL ADDBLO(VARSOR,PRIVE%ADR(1)%P)
C 24 
      CALL ADDBLO(VARSOR,PRIVE%ADR(2)%P)
C 25 
      CALL ADDBLO(VARSOR,PRIVE%ADR(3)%P)
C 26 
      CALL ADDBLO(VARSOR,PRIVE%ADR(4)%P)
C 27 
      CALL ADDBLO(VARSOR,MAXZ)
C 28 
      CALL ADDBLO(VARSOR,TMAXZ)
C 29 
      CALL ADDBLO(VARSOR,MAXV)
C 30 
      CALL ADDBLO(VARSOR,TMAXV)
C 31  VITESSE DE FROTTEMENT 
      CALL ADDBLO(VARSOR,T7)
C
C     TRACERS
C
      IF(NTRAC.GT.0) THEN
        DO ITRAC=1,NTRAC
          CALL ADDBLO(VARSOR,T%ADR(ITRAC)%P)
        ENDDO
      ENDIF
C
C     FOURIER ANALYSIS
C
      IF(NPERIAF.GT.0) THEN
        DO I=1,NPERIAF
C         VARIABLES FORCEES EN SORTIE (A VOIR)
          SORLEO(32+NTRAC+2*(I-1))=.TRUE.
          SORLEO(33+NTRAC+2*(I-1))=.TRUE.
C         FIN DE VARIABLES FORCEES EN SORTIE (A VOIR)
          CALL ADDBLO(VARSOR,AMPL%ADR(I)%P)
          CALL ADDBLO(VARSOR,PHAS%ADR(I)%P)
        ENDDO
      ENDIF
C
C     OTHER POSSIBLE VARIABLES ADDED BY USER
C
      J=32+NTRAC+2*NPERIAF
900   CONTINUE
      IF(SORLEO(J).OR.SORIMP(J)) THEN
        IF(NPRIV.LT.J-27-NTRAC-2*NPERIAF) THEN
          IF(LNG.EQ.1) THEN 
            WRITE(LU,*) 'POINT : AUGMENTER LE NOMBRE'
            WRITE(LU,*) '        DE TABLEAUX PRIVES'
          ENDIF
          IF(LNG.EQ.2) THEN 
            WRITE(LU,*) 'POINT : NUMBER OF PRIVATE ARRAYS'
            WRITE(LU,*) '        TOO SMALL'
          ENDIF
          CALL PLANTE(1)
          STOP
        ENDIF
        CALL ADDBLO(VARSOR,PRIVE%ADR(J-27-NTRAC-2*NPERIAF)%P)
        J=J+1
        IF(J.LE.MAXVAR) GO TO 900
      ENDIF
C
C     CLANDESTINE VARIABLES
C
      IF(VARCL%N.NE.0) THEN
        DO I=1,VARCL%N
          CALL ADDBLO(VARSOR,VARCL%ADR(I)%P)
          SORLEO(J+I-1)=.TRUE.
          TEXTE(J+I-1)=VARCLA(I)
        ENDDO
      ENDIF
C
C=======================================================================
C
C IMPRESSIONS :
C
      IF(LISTIN) THEN
         IF(LNG.EQ.1) WRITE(LU,22)
         IF(LNG.EQ.2) WRITE(LU,23)
      ENDIF
22    FORMAT(1X,///,21X,'****************************************',/,
     *21X,              '* FIN DE L''ALLOCATION DE LA MEMOIRE  : *',/,
     *21X,              '****************************************',/)
23    FORMAT(1X,///,21X,'*************************************',/,
     *21X,              '*    END OF MEMORY ORGANIZATION:    *',/,
     *21X,              '*************************************',/)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
