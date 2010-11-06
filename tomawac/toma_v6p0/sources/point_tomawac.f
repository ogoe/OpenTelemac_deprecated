C                       ************************
                        SUBROUTINE POINT_TOMAWAC
C                       ************************
C
C**********************************************************************
C  TOMAWAC - 6.0    Michel BENOIT (EDF R&D LNHE)            06/12/2004
C**********************************************************************
C
C      FONCTION:
C      =========
C
C    ALLOCATION DE LA MEMOIRE
C
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE PAR : HOMERE_TOMAWAC
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TOMAWAC
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
C
      INTEGER MEMT5,ITAMP, DIMXK0
      INTEGER IELBT,IELM0,IELB1
      INTEGER CFG(2)
C
      INTEGER NC,NS
C
C***********************************************************************
C
      IF (LNG.EQ.1) WRITE(LU,20)
      IF (LNG.EQ.2) WRITE(LU,40)
20    FORMAT(1X,///,21X,'*******************************',/,
     *21X,              '* CONSTRUCTION DES POINTEURS: *',/,
     *21X,              '*******************************',/)
40    FORMAT(1X,///,21X,'*******************************',/,
     *21X,              '*     MEMORY ORGANISATION     *',/,
     *21X,              '*******************************',/)
C
C-----------------------------------------------------------------------
C
C                     ***********************
C                     * POINTEURS TABLEAU A *
C                     ***********************
C
C-----------------------------------------------------------------------
C
      DIMXK0 = 1000000
      IELM0  = 10
      IELM2  = 11
      IELB1 = IELBOR(IELM2,1)
      CFG(1) = 1
      CFG(2) = 1
C
C-----------------------------------------------------------------------
C
C     ALLOCATION DE LA STRUCTURE DE MAILLAGE
C
      CALL ALMESH(MESH,'MESH  ',IELM2,SPHE,CFG,WAC_FILES(WACGEO)%LU,
     *            EQUA)
C
C     ALIAS POUR CERTAINES COMPOSANTES DU MAILLAGE MESH
C
      IKLE2  => MESH%IKLE%I
      IFABOR => MESH%IFABOR%I
      NBOR   => MESH%NBOR%I
C
      X     => MESH%X%R
      Y     => MESH%Y%R
      XEL   => MESH%XEL%R
      YEL   => MESH%YEL%R
      SURDET=> MESH%SURDET%R
C
      NELEM2=>MESH%NELEM
      NPOIN2=>MESH%NPOIN
      NPTFR =>MESH%NPTFR
C
      IELBT = IELBOR(IELM2,1)
C
      NPOIN3=NPOIN2*NPLAN
C
C
C     VARIABLES 4D A CONVECTER
      CALL ALLVEC(1,SF,'SF    ',NPOIN3*NF , 1 , 0 )
C
C     COEFFICIENT B POUR CONVECTION
      CALL ALLVEC(1,SB,'SB    ',NPOIN2*NF , 1 , 0 )
C
C     TABLEAU DES FREQUENCES DE DISCRETISATION PUIS DES DELTA F
      CALL ALLVEC(1,SFR,'SFR   ' ,NF , 1 , 0 )
      CALL ALLVEC(1,SDFR,'SDFR  ',NF , 1 , 0 )
      FREQ     =>SFR%R
      DFREQ    =>SDFR%R
C
C     VARIABLES "PHYSIQUES" DE DIM NPOIN3
      CALL ALLVEC(1,SXK,'SXK   ',NPOIN2*NF, 1 , 0 )
      CALL ALLVEC(1,SCG,'SCG   ',NPOIN2*NF, 1 , 0 )
C
      IF (TSOU) THEN
        CALL ALLVEC(1,STSDER,'STSDER',NF*NPOIN3 , 1 , 0 )
        CALL ALLVEC(1,STSTOT,'STSTOT',NF*NPOIN3 , 1 , 0 )
	CALL ALLVEC(1,SDF_LIM,'SDF_LIM',NPOIN2*NF, 1 , 0 )
      ELSE
        CALL ALLVEC(1,STSDER,'STSDER', 1, 1, 0 )
        CALL ALLVEC(1,STSTOT,'STSTOT', 1, 1, 0 )
	CALL ALLVEC(1,SDF_LIM,'SDF_LIM',1, 1 , 0 )
      ENDIF
      TSDER   => STSDER%R
      TSTOT   => STSTOT%R
      DF_LIM  => SDF_LIM%R
C
C     POINTEURS DES TABLEAUX POUR LES CONDITIONS AUX LIMITES
      CALL ALLVEC(1,SFBOR,'SFBOR ',IELBT, NPLAN*NF , 2 )
C
C     TABLEAUX POUR LES INTERACTIONS NON LINEAIRES
      IF (STRIF.NE.0) THEN
        CALL ALLVEC(1,SCOEF,'SCOEF ',16   , 1 , 0 )
      ELSE
        CALL ALLVEC(1,SCOEF ,'SCOEF ', 1, 1, 0 )
      ENDIF
      COEFNL   =>SCOEF%R
C
C     CHAMP CONVECTEUR
      IF (COUSTA .OR. MAREE) THEN
        NC=NPOIN3*NF
        NS=NPOIN3*NF
      ELSE
        NS=NPOIN3*NF
        NC=NPOIN3
      ENDIF
C
      F     =>SF%R
      B     =>SB%R
      XK    =>SXK%R
      CG    =>SCG%R
      FBOR  =>SFBOR%R
C
      IF (PROP) THEN
C       PIEDS DES CARACTERISTIQUES
        CALL ALLVEC(1,SSHP1,'SSHP1 ',NS , 1 , 0 )
        CALL ALLVEC(1,SSHP2,'SSHP2 ',NS , 1 , 0 )
        CALL ALLVEC(1,SSHP3,'SSHP3 ',NS , 1 , 0 )
        CALL ALLVEC(1,SSHZ ,'SSHZ  ',NS , 1 , 0 )
C
        CALL ALLVEC(1,SCX,'SCX   ',NC , 1 , 0 )
        CALL ALLVEC(1,SCY,'SCY   ',NC , 1 , 0 )
        CALL ALLVEC(1,SCT,'SCT   ',NC , 1 , 0 )
C
        IF (COURAN) THEN
          CALL ALLVEC(1,SSHF,'SSHF  ',NS , 1 , 0 )
          CALL ALLVEC(1,SCF ,'SCF   ',NC , 1 , 0 )
        ELSE
          CALL ALLVEC(1,SSHF  ,'SSHF  ', 1, 1, 0 )
          CALL ALLVEC(1,SCF   ,'SCF   ', 1, 1, 0 )
        ENDIF
      ELSE
        CALL ALLVEC(1,SSHP1 ,'SSHP1 ', 1, 1, 0 )
        CALL ALLVEC(1,SSHP2 ,'SSHP2 ', 1, 1, 0 )
        CALL ALLVEC(1,SSHP3 ,'SSHP3 ', 1, 1, 0 )
        CALL ALLVEC(1,SSHZ  ,'SSHZ  ', 1, 1, 0 )
        CALL ALLVEC(1,SCX   ,'SCX   ', 1, 1, 0 )
        CALL ALLVEC(1,SCY   ,'SCY   ', 1, 1, 0 )
        CALL ALLVEC(1,SCT   ,'SCT   ', 1, 1, 0 )
        CALL ALLVEC(1,SSHF  ,'SSHF  ', 1, 1, 0 )
        CALL ALLVEC(1,SCF   ,'SCF   ', 1, 1, 0 )
      ENDIF
      SHP1  =>SSHP1%R
      SHP2  =>SSHP2%R
      SHP3  =>SSHP3%R
      SHZ   =>SSHZ%R
      CX    =>SCX%R
      CY    =>SCY%R
      CT    =>SCT%R
      SHF   =>SSHF%R
      CF    =>SCF%R
C
C TABLEAUX DE DIMENSION NPOIN2
C
      CALL ALLVEC(1,SZF,'SZF   ',IELM2 , 1 , 2 )
      ZF    =>SZF%R
      CALL ALLVEC(1,SDEPTH,'SDEPTH',IELM2 , 1 , 2 )
      DEPTH =>SDEPTH%R
C    
C     AJOUTE PAR JMH 16/12/2008 (PEUT-ETRE PAS TOUJOURS UTILISE, A VOIR)
C
      CALL ALLVEC(1,SBETA,'SBETA ',IELM2,1,2)
      BETA => SBETA%R      
C
C     FIN JMH
C
      IF (.NOT.PROINF) THEN
        CALL ALLVEC(1,SDZX  ,'SDZX  ',IELM2 , 1 , 2 )
        CALL ALLVEC(1,SDZY  ,'SDZY  ',IELM2 , 1 , 2 )
      ELSE
        CALL ALLVEC(1,SDZX  ,'SDZX  ', 1, 1, 0 )
        CALL ALLVEC(1,SDZY  ,'SDZY  ', 1, 1, 0 )
      ENDIF
      DZX     =>SDZX%R
      DZY     =>SDZY%R
C
C
C
C
      IF(COURAN.OR.DONTEL) THEN
         CALL ALLVEC(1,SUC ,'SUC   ',IELM2 , 1 , 2 )
         CALL ALLVEC(1,SVC ,'SVC   ',IELM2 , 1 , 2 )
         CALL ALLVEC(1,SDUX,'SDUX  ',IELM2 , 1 , 2 )
         CALL ALLVEC(1,SDUY,'SDUY  ',IELM2 , 1 , 2 )
         CALL ALLVEC(1,SDVX,'SDVX  ',IELM2 , 1 , 2 )
         CALL ALLVEC(1,SDVY,'SDVY  ',IELM2 , 1 , 2 )
         CALL ALLVEC(1,SDZHDT,'SDZHDT',IELM2 , 1 , 2 )
         IF (WAC_FILES(WACCOB)%NAME.NE.' '.OR.
     *       WAC_FILES(WACCOF)%NAME.NE.' ') THEN
!BD_INCKA dans le cas parallel on redimensionne les tableaux
! afin de bien pouvoir relire toutes les valeur des champs 
!qui peuvent etre interpole, effectivement avant on avant 
!NPOIN max lu = NPOIN3 = NPOIN2 (Npt sous domaine) * 
!NTETA mais si NCSIZE>NTETA alors le nombre de points 
!totaux est surment superieur a  npoin3
           IF (NCSIZE.GT.1) THEN
             CALL ALLVEC(1,SXRELC,'SXRELC',NPOIN3*NCSIZE*3 , 1 , 0 )
             CALL ALLVEC(1,SYRELC,'SYRELC',NPOIN3*NCSIZE*3 , 1 , 0 )
           ELSE
             CALL ALLVEC(1,SXRELC,'SXRELC',NPOIN3 , 1 , 0 )
             CALL ALLVEC(1,SYRELC,'SYRELC',NPOIN3 , 1 , 0 )
           ENDIF
!BD_INCKA fin modif
         ELSE
           CALL ALLVEC(1,SXRELC,'SXRELC', 1, 1, 0 )
           CALL ALLVEC(1,SYRELC,'SYRELC', 1, 1, 0 )
         ENDIF
      ELSE
         CALL ALLVEC(1,SUC   ,'SUC   ', 1, 1, 0 )
         CALL ALLVEC(1,SVC   ,'SVC   ', 1, 1, 0 )
         CALL ALLVEC(1,SDUX  ,'SDUX  ', 1, 1, 0 )
         CALL ALLVEC(1,SDUY  ,'SDUY  ', 1, 1, 0 )
         CALL ALLVEC(1,SDVX  ,'SDVX  ', 1, 1, 0 )
         CALL ALLVEC(1,SDVY  ,'SDVY  ', 1, 1, 0 )
         CALL ALLVEC(1,SDZHDT,'SDZHDT', 1, 1, 0 )
         CALL ALLVEC(1,SXRELC,'SXRELC', 1, 1, 0 )
         CALL ALLVEC(1,SYRELC,'SYRELC', 1, 1, 0 )
      ENDIF
      UC      =>SUC%R
      VC      =>SVC%R
      DUX     =>SDUX%R
      DUY     =>SDUY%R
      DVX     =>SDVX%R
      DVY     =>SDVY%R
      XRELC   =>SXRELC%R
      YRELC   =>SYRELC%R
C
      IF (MAREE) THEN
         IF (WAC_FILES(WACCOB)%NAME.NE.' '.OR.
     *       WAC_FILES(WACCOF)%NAME.NE.' ') THEN
           CALL ALLVEC(1,SVC1,'SVC1  ',IELM2 , 1 , 2 )
           CALL ALLVEC(1,SUC1,'SUC1  ',IELM2 , 1 , 2 )
           CALL ALLVEC(1,SVC2,'SVC2  ',IELM2 , 1 , 2 )
           CALL ALLVEC(1,SUC2,'SUC2  ',IELM2 , 1 , 2 )
         ELSE
           CALL ALLVEC(1,SVC1  ,'SVC1  ', 1, 1, 0 )
           CALL ALLVEC(1,SUC1  ,'SUC1  ', 1, 1, 0 )
           CALL ALLVEC(1,SVC2  ,'SVC2  ', 1, 1, 0 )
           CALL ALLVEC(1,SUC2  ,'SUC2  ', 1, 1, 0 )
         ENDIF
         CALL ALLVEC(1,SZM1,'SZM1  ',IELM2 , 1 , 2 )
         CALL ALLVEC(1,SZM2,'SZM2  ',IELM2 , 1 , 2 )
         IF (WAC_FILES(WACMAB)%NAME.NE.' '.OR.
     *       WAC_FILES(WACMAF)%NAME.NE.' ') THEN
!BD_INCKA dans le cas parallel on redimensionne les tableaux
! afin de bien pouvoir relire toutes les valeur des champs 
!qui peuvent etre interpole, effectivement avant on avant 
!NPOIN max lu = NPOIN3 = NPOIN2 (Npt sous domaine) * 
!NTETA mais si NCSIZE>NTETA alors le nombre de points 
!totaux est surment superieur Ã  npoin3
           IF (NCSIZE.GT.1) THEN
             CALL ALLVEC(1,SXRELM,'SXRELM',NPOIN3*NCSIZE*3 , 1 , 0 )
             CALL ALLVEC(1,SYRELM,'SYRELM',NPOIN3*NCSIZE*3 , 1 , 0 )
           ELSE
             CALL ALLVEC(1,SXRELM,'SXRELM',NPOIN3 , 1 , 0 )
             CALL ALLVEC(1,SYRELM,'SYRELM',NPOIN3 , 1 , 0 )
           ENDIF
!BD_INCKA fin modif
         ELSE
           CALL ALLVEC(1,SXRELM,'SXRELM', 1, 1, 0 )
           CALL ALLVEC(1,SYRELM,'SYRELM', 1, 1, 0 )
         ENDIF
      ELSE
         CALL ALLVEC(1,SVC1  ,'SVC1  ', 1, 1, 0 )
         CALL ALLVEC(1,SUC1  ,'SUC1  ', 1, 1, 0 )
         CALL ALLVEC(1,SVC2  ,'SVC2  ', 1, 1, 0 )
         CALL ALLVEC(1,SUC2  ,'SUC2  ', 1, 1, 0 )
         CALL ALLVEC(1,SZM1  ,'SZM1  ', 1, 1, 0 )
         CALL ALLVEC(1,SZM2  ,'SZM2  ', 1, 1, 0 )
         CALL ALLVEC(1,SXRELM,'SXRELM', 1, 1, 0 )
         CALL ALLVEC(1,SYRELM,'SYRELM', 1, 1, 0 )
      ENDIF
      DZHDT   =>SDZHDT%R
      UC1     =>SUC1%R
      VC1     =>SVC1%R
      UC2     =>SUC2%R
      VC2     =>SVC2%R
      ZM1     =>SZM1%R
      ZM2     =>SZM2%R
      XRELM   =>SXRELM%R
      YRELM   =>SYRELM%R
C
      IF (VENT) THEN
         CALL ALLVEC(1,SUV,'SUV   ',IELM2 , 1 , 2 )
         CALL ALLVEC(1,SVV,'SVV   ',IELM2 , 1 , 2 )
         IF (WAC_FILES(WACVEB)%NAME.NE.' '.OR.
     *       WAC_FILES(WACVEF)%NAME.NE.' ') THEN
!BD_INCKA dans le cas parallel on redimensionne les tableaux
! afin de bien pouvoir relire toutes les valeur des champs 
!qui peuvent etre interpole, effectivement avant on avant 
!NPOIN max lu = NPOIN3 = NPOIN2 (Npt sous domaine) * 
!NTETA mais si NCSIZE>NTETA alors le nombre de points 
!totaux est surment superieur a  npoin3
           IF (NCSIZE.GT.1) THEN
             CALL ALLVEC(1,SXRELV,'SXRELV',NPOIN3*NCSIZE*3 , 1 , 0 )
             CALL ALLVEC(1,SYRELV,'SYRELV',NPOIN3*NCSIZE*3 , 1 , 0 )
           ELSE
             CALL ALLVEC(1,SXRELV,'SXRELV',NPOIN3 , 1 , 0 )
             CALL ALLVEC(1,SYRELV,'SYRELV',NPOIN3 , 1 , 0 )
           ENDIF
!BD_INCKA fin modif
           CALL ALLVEC(1,SVV1,'SVV1  ',IELM2 , 1 , 2 )
           CALL ALLVEC(1,SUV1,'SUV1  ',IELM2 , 1 , 2 )
           CALL ALLVEC(1,SVV2,'SVV2  ',IELM2 , 1 , 2 )
           CALL ALLVEC(1,SUV2,'SUV2  ',IELM2 , 1 , 2 )
         ELSE
           CALL ALLVEC(1,SXRELV,'SXRELV', 1, 1, 0 )
           CALL ALLVEC(1,SYRELV,'SYRELV', 1, 1, 0 )
           CALL ALLVEC(1,SVV1  ,'SVV1  ', 1, 1, 0 )
           CALL ALLVEC(1,SUV1  ,'SUV1  ', 1, 1, 0 )
           CALL ALLVEC(1,SVV2  ,'SVV2  ', 1, 1, 0 )
           CALL ALLVEC(1,SUV2  ,'SUV2  ', 1, 1, 0 )
         ENDIF
      ELSE
         CALL ALLVEC(1,SUV   ,'SUV   ', 1, 1, 0 )
         CALL ALLVEC(1,SVV   ,'SVV   ', 1, 1, 0 )
         CALL ALLVEC(1,SXRELV,'SXRELV', 1, 1, 0 )
         CALL ALLVEC(1,SYRELV,'SYRELV', 1, 1, 0 )
         CALL ALLVEC(1,SVV1  ,'SVV1  ', 1, 1, 0 )
         CALL ALLVEC(1,SUV1  ,'SUV1  ', 1, 1, 0 )
         CALL ALLVEC(1,SVV2  ,'SVV2  ', 1, 1, 0 )
         CALL ALLVEC(1,SUV2  ,'SUV2  ', 1, 1, 0 )
      ENDIF
      UV      =>SUV%R
      VV      =>SVV%R
      VV1     =>SVV1%R
      UV1     =>SUV1%R
      VV2     =>SVV2%R
      UV2     =>SUV2%R
      XRELV   =>SXRELV%R
      YRELV   =>SYRELV%R
C
C
C
C
      IF (SPHE) THEN
         CALL ALLVEC(1,SCOSF,'SCOSF ',IELM2 , 1 , 2 )
         CALL ALLVEC(1,STGF,'STGF  ',IELM2 , 1 , 2 )
      ELSE
         CALL ALLVEC(1,SCOSF ,'SCOSF ', 1, 1, 0 )
         CALL ALLVEC(1,STGF  ,'STGF  ', 1, 1, 0 )
      ENDIF
      COSF    =>SCOSF%R
      TGF     =>STGF%R
C
C
C
C
C     TABLEAU DES POSITIONS RELATIVES DES PLANS DE DIRECTION PUIS
C     DES COS ET SIN TETA
      CALL ALLVEC(1,STETA,'STETA ',NPLAN+1 , 1 , 0 )
      CALL ALLVEC(1,SCOSTE,'SCOSTE',NPLAN , 1 , 0 )
      CALL ALLVEC(1,SSINTE,'SSINTE',NPLAN , 1 , 0 )
      TETA     =>STETA%R
      COSTET   =>SCOSTE%R
      SINTET   =>SSINTE%R
C
C
C
C     POINTEURS DES TABLEAUX DE TRAVAIL (PAR POINTS ET PAR ELEMENTS)
      CALL ALLVEC(1,ST0, 'ST0   ',IELM2 , 1 , 2 )
      CALL ALLVEC(1,ST1, 'ST1   ',IELM2 , 1 , 2 )
      CALL ALLVEC(1,ST2, 'ST2   ',IELM2 , 1 , 2 )
      CALL ALLVEC(1,ST3, 'ST3   ',IELM2 , 1 , 2 )
      CALL ALLVEC(1,ST4, 'ST4   ',IELM2 , 1 , 2 )
      T0      =>ST0%R
      T1      =>ST1%R
      T2      =>ST2%R
      T3      =>ST3%R
      T4      =>ST4%R
C
C
C
C     POINTEURS DE MATRICES, AM1 MATRICE SYMETRIQUE
      CALL ALLMAT(AM1,'AM1   ',IELM2,IELM2,CFG,'Q','S')
C
C
C
C     TABLEAUX DE TRAVAIL DIVERS
      CALL ALLVEC(1,SW1,   'SW1   ',NELEM2 , 4 , 0 )
      W1      =>SW1%R
      CALL ALLVEC(1,STRA31,'STRA31',NPOIN2 , 1 , 0 )
      CALL ALLVEC(1,STRA32,'STRA32',NPOIN2 , 1 , 0 )
      CALL ALLVEC(1,STRA33,'STRA33',NPOIN2 , 1 , 0 )
      CALL ALLVEC(1,STRA34,'STRA34',NPOIN2 , 1 , 0 )
      CALL ALLVEC(1,STRA35,'STRA35',NPOIN2 , 1 , 0 )
      CALL ALLVEC(1,STRA36,'STRA36',NPOIN2 , 1 , 0 )
      CALL ALLVEC(1,STRA37,'STRA37',NPOIN2 , 1 , 0 )
      CALL ALLVEC(1,STRA38,'STRA38',NPOIN2 , 1 , 0 )
      CALL ALLVEC(1,STRA39,'STRA39',NPOIN2 , 1 , 0 )
      CALL ALLVEC(1,STRA40,'STRA40',NPOIN2 , 1 , 0 )
      CALL ALLVEC(1,STRA41,'STRA41',NPOIN2 , 1 , 0 )
      CALL ALLVEC(1,STRA42,'STRA42',NPOIN2 , 1 , 0 )
      TRA31   => STRA31%R
      TRA32   => STRA32%R
      TRA33   => STRA33%R
      TRA34   => STRA34%R
      TRA35   => STRA35%R
      TRA36   => STRA36%R
      TRA37   => STRA37%R
      TRA38   => STRA38%R
      TRA39   => STRA39%R
      TRA40   => STRA40%R
      TRA41   => STRA41%R
      TRA42   => STRA42%R
C
C
C
        CALL ALLVEC(1,STRA02,'STRA02',NPOIN3*NF,2, 0 )
        TRA02 =>STRA02%R
C
C.......TABLEAUX DE TRAVAIL DIVERS
!BD_INCKA dans le cas parallel on redimensionne les tableaux
! afin de bien pouvoir relire toutes les valeur des champs 
!qui peuvent etre interpole, effectivement avant on avant 
!NPOIN max lu = NPOIN3 = NPOIN2 (Npt sous domaine) * 
!NTETA mais si NCSIZE>NTETA alors le nombre de points 
!totaux est surment superieur Ã  npoin3
        IF (NCSIZE.GT.1) THEN
          CALL ALLVEC(1,STRA01,'STRA01',NPOIN3*NCSIZE*3 , 8 , 0  )
        ELSE
          CALL ALLVEC(1,STRA01,'STRA01',NPOIN3 , 8 , 0  )
        ENDIF
!BD_INCKA fin modif
C
        CALL ALLVEC(1,STRAB1,'STRAB1', 1, 1, 0 )
        TRA01   =>STRA01%R
        TRAB1   =>STRAB1%R
C
        IF (TSOU) THEN
          CALL ALLVEC(1,STOLD,'STOLD ',NPOIN3 , 1 , 0 )
          CALL ALLVEC(1,STNEW,'STNEW ',NPOIN3 , 1 , 0 )
        ELSE
          CALL ALLVEC(1,STOLD ,'STOLD ', 1, 1, 0 )
          CALL ALLVEC(1,STNEW ,'STNEW ', 1, 1, 0 )
        ENDIF
        TOLD    => STOLD%R
        TNEW    => STNEW%R
C
        CALL ALLVEC(1,STRA43,'STRA43',NPOIN2 , 1 , 0 )
        CALL ALLVEC(1,STRA44,'STRA44',NPOIN2 , 1 , 0 )
        TRA43   => STRA43%R
        TRA44   => STRA44%R
C
        IF (SORLEO(22).OR.SORLEO(32)) THEN
          CALL ALLVEC(1,STRA56,'STRA56',NPOIN2 , 1 , 0 )
        ELSE
          CALL ALLVEC(1,STRA56,'STRA56', 1, 1, 0 )
        ENDIF
        IF (SORLEO(23).OR.SORLEO(33)) THEN
          CALL ALLVEC(1,STRA57,'STRA57',NPOIN2 , 1 , 0 )
        ELSE
          CALL ALLVEC(1,STRA57,'STRA57', 1, 1, 0 )
        ENDIF
        IF (SORLEO(25)) THEN
          CALL ALLVEC(1,STRA58,'STRA58',NPOIN2 , 1 , 0 )
        ELSE
          CALL ALLVEC(1,STRA58,'STRA58', 1, 1, 0 )
        ENDIF
        IF (SORLEO(28)) THEN
          CALL ALLVEC(1,STRA61,'STRA61',NPOIN2 , 1 , 0 )
        ELSE
          CALL ALLVEC(1,STRA61,'STRA61', 1, 1, 0 )
        ENDIF
        IF (SORLEO(29)) THEN
          CALL ALLVEC(1,STRA62,'STRA62',NPOIN2 , 1 , 0 )
        ELSE
          CALL ALLVEC(1,STRA62,'STRA62', 1, 1, 0 )
        ENDIF
        IF (SORLEO(30)) THEN
          CALL ALLVEC(1,STRA63,'STRA63',NPOIN2 , 1 , 0 )
        ELSE
          CALL ALLVEC(1,STRA63,'STRA63', 1, 1, 0 )
        ENDIF
        IF (SORLEO(31)) THEN
          CALL ALLVEC(1,STRA64,'STRA64',NPOIN2 , 1 , 0 )
        ELSE
          CALL ALLVEC(1,STRA64,'STRA64', 1, 1, 0 )
        ENDIF
        IF (SORLEO(32)) THEN
          CALL ALLVEC(1,STRA65,'STRA65',NPOIN2 , 1 , 0 )
        ELSE
          CALL ALLVEC(1,STRA65,'STRA65', 1, 1, 0 )
        ENDIF
        IF (SORLEO(33)) THEN
          CALL ALLVEC(1,STRA66,'STRA66',NPOIN2 , 1 , 0 )
        ELSE
          CALL ALLVEC(1,STRA66,'STRA66', 1, 1, 0 )
        ENDIF
        IF (SORLEO(34)) THEN
          CALL ALLVEC(1,STRA60,'STRA60',NPOIN2 , 1 , 0 )
        ELSE
          CALL ALLVEC(1,STRA60,'STRA60', 1, 1, 0 )
        ENDIF
        IF (.NOT.PROINF) THEN
          IF ( SORLEO(11).OR.SORLEO(12).OR.SORLEO(13).OR.
     *         SORLEO(14).OR.SORLEO(15) ) THEN
            CALL ALLVEC(1,STRA51,'STRA51',NPOIN2 , 1 , 0 )
            CALL ALLVEC(1,STRA52,'STRA52',NPOIN2 , 1 , 0 )
            CALL ALLVEC(1,STRA53,'STRA53',NPOIN2 , 1 , 0 )
            CALL ALLVEC(1,STRA54,'STRA54',NPOIN2 , 1 , 0 )
            CALL ALLVEC(1,STRA55,'STRA55',NPOIN2 , 1 , 0 )
          ELSE
            CALL ALLVEC(1,STRA51,'STRA51', 1, 1, 0 )
            CALL ALLVEC(1,STRA52,'STRA52', 1, 1, 0 )
            CALL ALLVEC(1,STRA53,'STRA53', 1, 1, 0 )
            CALL ALLVEC(1,STRA54,'STRA54', 1, 1, 0 )
            CALL ALLVEC(1,STRA55,'STRA55', 1, 1, 0 )
          ENDIF
          IF (SORLEO(16)) THEN
            CALL ALLVEC(1,STRA59,'STRA59',NPOIN2 , 1 , 0 )
          ELSE
            CALL ALLVEC(1,STRA59,'STRA59', 1, 1, 0 )
          ENDIF
        ELSE
          CALL ALLVEC(1,STRA51,'STRA51', 1, 1, 0 )
          CALL ALLVEC(1,STRA52,'STRA52', 1, 1, 0 )
          CALL ALLVEC(1,STRA53,'STRA53', 1, 1, 0 )
          CALL ALLVEC(1,STRA54,'STRA54', 1, 1, 0 )
          CALL ALLVEC(1,STRA55,'STRA55', 1, 1, 0 )
          CALL ALLVEC(1,STRA59,'STRA59', 1, 1, 0 )
        ENDIF
        TRA51   => STRA51%R
        TRA52   => STRA52%R
        TRA53   => STRA53%R
        TRA54   => STRA54%R
        TRA55   => STRA55%R
        TRA56   => STRA56%R
        TRA57   => STRA57%R
        TRA58   => STRA58%R
        TRA59   => STRA59%R
        TRA60   => STRA60%R
        TRA61   => STRA61%R
        TRA62   => STRA62%R
        TRA63   => STRA63%R
        TRA64   => STRA64%R
        TRA65   => STRA65%R
        TRA66   => STRA66%R
C
C...... TABLEAU DEDIE A L'UTILISATEUR
        CALL ALLVEC(1,SPRIVE,'SPRIVE',NPOIN3*NPRIV*NF, 1 , 0 )
        PRIVE   => SPRIVE%R
C
C.......TABLEAUX NON-UTILISES DANS CE MODE
        CALL ALLVEC(1,STRA15,'STRA15', 1, 1, 0 )
        CALL ALLVEC(1,STRA16,'STRA16', 1, 1, 0 )
C
C.......BLOC POUR LES SORTIES GRAPHIQUES : VARSOR
        CALL ALLBLO(VARSOR,'VARSOR')
        CALL ADDBLO(VARSOR,STRA37)
        CALL ADDBLO(VARSOR,STRA38)
        CALL ADDBLO(VARSOR,STRA32)
        CALL ADDBLO(VARSOR,STRA31)
        CALL ADDBLO(VARSOR,SZF)
        CALL ADDBLO(VARSOR,SDEPTH)
        CALL ADDBLO(VARSOR,SUC)
        CALL ADDBLO(VARSOR,SVC)
        CALL ADDBLO(VARSOR,SUV)
        CALL ADDBLO(VARSOR,SVV)
        CALL ADDBLO(VARSOR,STRA51)
        CALL ADDBLO(VARSOR,STRA52)
        CALL ADDBLO(VARSOR,STRA53)
        CALL ADDBLO(VARSOR,STRA54)
        CALL ADDBLO(VARSOR,STRA55)
        CALL ADDBLO(VARSOR,STRA59)
        CALL ADDBLO(VARSOR,SPRIVE)
        CALL ADDBLO(VARSOR,STRA33)
        CALL ADDBLO(VARSOR,STRA34)
        CALL ADDBLO(VARSOR,STRA35)
        CALL ADDBLO(VARSOR,STRA36)
        CALL ADDBLO(VARSOR,STRA56)
        CALL ADDBLO(VARSOR,STRA57)
        CALL ADDBLO(VARSOR,STRA42)
        CALL ADDBLO(VARSOR,STRA58)
        CALL ADDBLO(VARSOR,STRA44)
        CALL ADDBLO(VARSOR,STRA41)
        CALL ADDBLO(VARSOR,STRA61)
        CALL ADDBLO(VARSOR,STRA62)
        CALL ADDBLO(VARSOR,STRA63)
        CALL ADDBLO(VARSOR,STRA64)
        CALL ADDBLO(VARSOR,STRA65)
        CALL ADDBLO(VARSOR,STRA66)
        CALL ADDBLO(VARSOR,STRA60)
C       VARIABLE 35
        CALL ADDBLO(VARSOR,SBETA)
C
C
C.....BLOC POUR LA VALIDATION
C
      IF (VALID) THEN
C
C       POINTEURS DES TABLEAUX DE TRAVAIL (PAR POINTS ET PAR ELEMENTS)
        CALL ALLVEC(1,ST5, 'ST5   ',IELM2 , 1 , 2 )
        CALL ALLVEC(1,ST6, 'ST6   ',IELM2 , 1 , 2 )
        CALL ALLVEC(1,ST7, 'ST7   ',IELM2 , 1 , 2 )
        T5      =>ST5%R
        T6      =>ST6%R
        T7      =>ST7%R
C
C       ATTENTION A LA COHERENCE AVEC LE VECTEUR ALIRE DECLARE
C       DANS declarations_tomawac.f
C       7 Variables sont retenues pour la validation
C          Hauteur significative      HM0       ( 2)
C          Direction moyenne          DMOY      ( 3)
C          Etalement directionnel     SPD       ( 4)
C          Force motrice suivant x    FX        (11)
C          Force motrice suivant y    FY        (12)
C          Frequence moyenne Fm-10    FMOY      (18)
C          Frequence moyenne Fm01     FM01      (19)
        CALL ALLBLO(BST1,'BST1  ')
        CALL ADDBLO(BST1,ST0)
        CALL ADDBLO(BST1,ST1)
        CALL ADDBLO(BST1,ST2)
        CALL ADDBLO(BST1,ST3)
        CALL ADDBLO(BST1,ST0)
        CALL ADDBLO(BST1,ST0)
        CALL ADDBLO(BST1,ST0)
        CALL ADDBLO(BST1,ST0)
        CALL ADDBLO(BST1,ST0)
        CALL ADDBLO(BST1,ST0)
        CALL ADDBLO(BST1,ST4)
        CALL ADDBLO(BST1,ST5)
        CALL ADDBLO(BST1,ST0)
        CALL ADDBLO(BST1,ST0)
        CALL ADDBLO(BST1,ST0)
        CALL ADDBLO(BST1,ST0)
        CALL ADDBLO(BST1,ST0)
        CALL ADDBLO(BST1,ST6)
        CALL ADDBLO(BST1,ST7)
      ELSE
        CALL ALLVEC(1,ST5   ,'ST5   ', 1, 1, 0 )
        CALL ALLVEC(1,ST6   ,'ST6   ', 1, 1, 0 )
        CALL ALLVEC(1,ST7   ,'ST7   ', 1, 1, 0 )
      ENDIF
C
C-----------------------------------------------------------------------
C
C                     **********************
C                     * POINTEUR TABLEAU IA*
C                     **********************
C
C-----------------------------------------------------------------------
C
C
      CALL ALLVEC(2,SLIFBR,'SLIFBR',IELB1, 1 , 1 )
      CALL ALLVEC(2,SIBOR,'SIBOR ',NELEM2, 7 , 0 )
      CALL ALLVEC(2,BOUNDARY_COLOUR,'BNDCOL',IELB1,1,1)
      LIFBOR  => SLIFBR%I
      IBOR    => SIBOR%I
C
C PIEDS DES CARACTERISTIQUES
C
      IF (PROP) THEN
        CALL ALLVEC(2,SELT,'SELT  ',NS , 1 , 0 )
        CALL ALLVEC(2,SETA,'SETA  ',NS , 1 , 0 )
        IF (COURAN) THEN
          CALL ALLVEC(2,SFRE,'SFRE  ',NS , 1 , 0 )
        ELSE
          CALL ALLVEC(2,SFRE  ,'SFRE  ', 1, 1, 0 )
        ENDIF
      ELSE
        CALL ALLVEC(2,SELT  ,'SELT  ', 1, 1, 0 )
        CALL ALLVEC(2,SETA  ,'SETA  ', 1, 1, 0 )
        CALL ALLVEC(2,SFRE  ,'SFRE  ', 1, 1, 0 )
      ENDIF
      ELT   => SELT%I
      ETA   => SETA%I
      FRE   => SFRE%I
C
C TABLEAU UTILE POUR LES CARACTERISTIQUES
C
      CALL ALLVEC(2,SETAP1,'SETAP1',NPLAN, 1 , 0 )
      ETAP1  => SETAP1%I
C
C TABLEAUX DE TRAVAIL UTILISES DANS L'APPEL A INBIEF ET INIPIE
C
      CALL ALLVEC(2,SITR31,'SITR31',IELM0 , 1 , 1 )
      CALL ALLVEC(2,SITR32,'SITR32',IELM0 , 1 , 1 )
      CALL ALLVEC(2,SITR33,'SITR33',IELM0 , 1 , 1 )
      ITR31   => SITR31%I
      ITR32   => SITR32%I
      ITR33   => SITR33%I
C
C
C
C.......TABLEAUX DE TRAVAIL ENTIERS
        CALL ALLVEC(2,SITR03,'SITR03', 1, 1, 0 )
        CALL ALLVEC(2,SITRB1,'SITRB1', 1, 1, 0 )
        ITR03 => SITR03%I
        ITRB1 => SITRB1%I
        CALL ALLVEC(2,SITR01,'SITR01',NPOIN3 ,3, 0 )
        ITR01 => SITR01%I
C
C.......INTERACTIONS NON-LINEAIRES
C
        IF (STRIF.NE.0) THEN
          CALL ALLVEC(2,SIAGNL,'SIAGNL',8*NPLAN, 1 , 0 )
        ELSE
          CALL ALLVEC(2,SIAGNL,'SIAGNL', 1, 1, 0 )
        ENDIF
        IANGNL => SIAGNL%I
C
C.......PASSAGE SPECTRE RELATIF <--> ABSOLU (TRANSF)
        IF (COUSTA .OR. MAREE) THEN
          CALL ALLVEC(2,SITR11,'SITR11',NPOIN2, 1 , 0 )
          CALL ALLVEC(2,SITR12,'SITR12',NPOIN2, 1 , 0 )
          CALL ALLVEC(2,SITR13,'SITR13',NPOIN2, 1 , 0 )
        ELSE
          CALL ALLVEC(2,SITR11,'SITR11', 1, 1, 0 )
          CALL ALLVEC(2,SITR12,'SITR12', 1, 1, 0 )
          CALL ALLVEC(2,SITR13,'SITR13', 1, 1, 0 )
        ENDIF
        ITR11   => SITR11%I
        ITR12   => SITR12%I
        ITR13   => SITR13%I
C
        CALL ALLVEC(2,SKNI  ,'SKNI  ', 1, 1, 0 )
        CALL ALLVEC(2,SKNOGL,'SKNOGL', 1, 1, 0 )
        CALL ALLVEC(2,SELI  ,'SELI  ', 1, 1, 0 )
        CALL ALLVEC(2,SKELGL,'SKELGL', 1, 1, 0 )
        KNI   => SKNI%I
        KNOGL => SKNOGL%I
        ELI   => SELI%I
        KELGL => SKELGL%I
C
C***********************************************************************
C
C CONTROLES ET IMPRESSIONS
C
      IF (LNG.EQ.1) WRITE(LU,22)
      IF (LNG.EQ.2) WRITE(LU,23)
   22 FORMAT(1X,///,21X,'****************************************',/,
     *21X,              '* FIN DE L''ALLOCATION DE LA MEMOIRE  : *',/,
     *21X,              '****************************************',/)
   23 FORMAT(1X,///,21X,'*************************************',/,
     *21X,              '*    END OF MEMORY ORGANIZATION:    *',/,
     *21X,              '*************************************',/)
C
      RETURN
      END
