C                       *****************
                        SUBROUTINE CONDIW
C                       *****************
     *( AT, LT , DPI, TC1, TC2, NPC , TV1, TV2, NPV, TM1, TM2 , NPM ,
     *  NVHMA  , NVCOU )
C
C***********************************************************************
C TOMAWAC     V1.0         01/02/95          F.MARCOS (LNH) 30 87 72 66
C             V5.0         25/08/00
C***********************************************************************
C
C      FONCTION:
C      =========
C
C    INITIALISATION DES TABLEAUX DES GRANDEURS PHYSIQUES
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    F           !<-- ! DENSITE SPECTRALE D'ENERGIE                  !
C !    TETA        !<-- ! DIRECTIONS DE PROPAGATION                    !
C !    FREQ        !<-- ! FREQUENCES DISCRETISEES                      !
C !    UC,VC       !<-- ! COMPOSANTES DU CHAMP DE COURANT              !
C !    UC1,VC1     !<-- ! COMPOSANTES DU CHAMP DE COURANT INFERIEUR    !
C !    UC2,VC2     !<-- ! COMPOSANTES DU CHAMP DE COURANT SUPERIEUR    !
C !    UV,VV       !<-- ! COMPOSANTES DU CHAMP DE VENT INITIAL         !
C !    UV1,VV1     !<-- ! COMPOSANTES DU CHAMP DE VENT INFERIEUR       !
C !    UV2,VV2     !<-- ! COMPOSANTES DU CHAMP DE VENT SUPERIEUR       !
C !    ZM          !<-- ! HAUTEUR DE LA MAREE INITIALE                 !
C !    ZM1,ZM2     !<-- ! HAUTEURS DE LA MAREE INFERIEURE ET SUPERIEURE!
C !    X,Y         ! -->! COORDONNEES DES POINTS DU MAILLAGE 2D        !
C !    XRELV,YRELV ! -->! COORDONNEES DES POINTS DES VENTS RELEVES     !
C !    XRELC,YRELC ! -->! COORDONNEES DES POINTS DES COURANTS RELEVES  !
C !    XRELM,YRELM ! -->! COORDONNEES DES POINTS DES HAUTEURS RELEVEES !
C !    TRA31       !<-->! TABLEAUX DE TRAVAIL REELS                    !
C !    AT          !<-- ! TEMPS DU CALCUL                              !
C !    DDC         ! -->! DATE DE DEBUT DU CALCUL                      !
C !    TV1         !<-- ! TEMPS CORRESPONDANT AU VENT 1                !
C !    TV2         !<-- ! TEMPS CORRESPONDANT AU VENT 2                !
C !    TC1         !<-- ! TEMPS CORRESPONDANT AU COURANT 1             !
C !    TC2         !<-- ! TEMPS CORRESPONDANT AU COURANT 2             !
C !    TM1         !<-- ! TEMPS CORRESP. A LA HAUTEUR DE LA MAREE 1    !
C !    TM2         !<-- ! TEMPS CORRESP. A LA HAUTEUR DE LA MAREE 2    !
C !    COUSTA      ! -->! LOGIQUE INDIQUANT UN COURANT STATIONNAIRE    !
C !    VENT        ! -->! LOGIQUE INDIQUANT SI ON CONSIDERE UN VENT    !
C !    MAREE       ! -->! LOGIQUE INDIQUANT SI ON CONSIDERE LA MAREE   !
C !    NOMCOB-F    ! -->! NOM DU FICHIER DES COURANTS (BIN./FORMAT.)   !
C !    NOMVEB-F    ! -->! NOM DU FICHIER DES VENTS                     !
C !    NOMMAB-F    ! -->! NOM DU FICHIER DE HAUTEUR DE LA MAREE        !
C !    NCOB-F      ! -->! NUMERO LOGIQUE DU FICHIER DES COURANTS       !
C !    NVEB-F      ! -->! NUMERO LOGIQUE DU FICHIER DES VENTS          !
C !    NMAB-F      ! -->! NUMERO LOGIQUE DU FIC. DE HAUTEUR DE LA MAREE!
C !    BINCOU      ! -->! BINAIRE DU FICHIER DES COURANTS              !
C !    BINVEN      ! -->! BINAIRE DU FICHIER DES VENTS                 !
C !    BINMAR      ! -->! BINAIRE DU FICHIER DES HAUTEURS DE LA MAREE  !
C !    NBOR        ! -->! NUMERO GLOBAUX DES POINTS DE BORD            !
C !    NPTFR       ! -->! NOMBRE DE POINTS DE BORD                     !
C !    IDTEL       ! -->! RANG DE LA DONNEE TELEMAC A RECUPERER        !
C !    NPTT        ! -->! NUMERO DU PAS DE TEMPS DU FICHIER TELEMAC    !
C !    DONTEL      ! -->! LOGIQUE INDIQUANT SI ON RECUPERE UNE DON.TEL.!
C !    INDIC       ! -->! FORMAT DU FICHIER DES COURANTS               !
C !    INDIV       ! -->! FORMAT DU FICHIER DES VENTS                  !
C !    INDIV       ! -->! FORMAT DU FICHIER DES HAUTEURS               !
C !    NPOIN3      ! -->! NOMBRE DE POINTS DU MAILLAGE 3D              !
C !    NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE 2D              !
C !    NPLAN       ! -->! NOMBRE DE DIRECTIONS DE PROPAGATION          !
C !    NPV         ! -->! NOMBRE DE POINTS DU FICHIER DES VENTS        !
C !    NPC         ! -->! NOMBRE DE POINTS DU FICHIER DES COURANTS     !
C !    NPM         ! -->! NOMBRE DE POINTS DU FICHIER DES HAUTEURS     !
C !    NF          ! -->! NOMBRE DE FREQUENCES                         !
C !    F1          ! -->! FREQUENCE MINIMALE                           !
C !    RAISF       ! -->! RAISON FREQUENTIELLE                         !
C !    INISPE      ! -->! INDICATEUR D'INITIALISATION DU SPECTRE       !
C !    HM0         ! -->! HAUTEUR SIGNIFICATIVE JONSWAP                !
C !    FPIC        ! -->! FREQUENCE DE PIC JONSWAP                     !
C !    GAMMA       ! -->! FACTEUR DE FORME DE PIC JONSWAP              !
C !    SIGMAA      ! -->! VALEUR DE SIGMA JONSWAP POUR F < FP          !
C !    SIGMAB      ! -->! VALEUR DE SIGMA JONSWAP POUR F > FP          !
C !    ALPHIL      ! -->! CONSTANTE DE PHILLIPS (ALPHA)                !
C !    FETCH       ! -->! FETCH MOYEN                                  !
C !    FREMAX      ! -->! VALEUR MAXIMUM DE LA FREQUENCE DE PIC        !
C !    TETA1       ! -->! DIRECTION PRINCIPALE 1 POUR FRA              !
C !    SPRED1      ! -->! ETALEMENT DIRECTIONNEL 1 POUR FRA            !
C !    TETA2       ! -->! DIRECTION PRINCIPALE 2 POUR FRA              !
C !    SPRED2      ! -->! ETALEMENT DIRECTIONNEL 2 POUR FRA            !
C !    XLAMDA      ! -->! FACTEUR DE PONDERATION POUR LA FRA           !
C !    GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR                 !
C !    E2FMIN      ! -->! VALEUR MINIMALE D'ENERGIE                    !
C !    PRIVE       ! -->! TABLEAU POUR L'UTILISATEUR DE DIMENSION      !
C !    NPRIV       !    ! NPOIN3*NPRIV                                 !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE PAR : WAC
C SOUS-PROGRAMMES APPELES : OV, ANAMAR, LECDOI, LECDON, LECHAM, SPEINI
C
C***********************************************************************
C
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TOMAWAC
      USE BIEF
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
      DOUBLE PRECISION AT , DPI, TC1, TC2, TV1, TV2, TM1, TM2
      INTEGER          NPC, NPV, NPM
      INTEGER          LT
C
C
      INTEGER          NP0, NP1, NP2, NP3, NP4, NP5
      INTEGER          IPLAN, IFREQ , NVHMA, NVCOU, IBID
C
      CHARACTER*7   CHDON
!BD_INCKA modif pour //
      INTEGER          NPOIN_ALL
!BD_INCKA fin modif
C
C***********************************************************************
C
      AT  = 0.D0
!BD_INCKA modif pour // 
      NPOIN_ALL = MAX(NPOIN3,NPOIN3*NCSIZE*3)
!BD_INCKA fin modif pour //
!      NP0 = NPOIN3+1
!      NP1 = 2*NPOIN3
!      NP2 = NP1+1
!      NP3 = 3*NPOIN3
!      NP4 = NP3+1
!      NP5 = 4*NPOIN3
      NP0 = NPOIN_ALL+1
      NP1 = 2*NPOIN_ALL
      NP2 = NP1+1
      NP3 = 3*NPOIN_ALL
      NP4 = NP3+1
      NP5 = 4*NPOIN_ALL

C
C-----------------------------------------------------------------------
C
C   INITIALISATION DU COURANT ET DE LA HAUTEUR DE MAREE
C
C
      IF (MAREE) THEN
        IF(LNG.EQ.1) THEN
          CHDON='COURANT'
        ELSE
          CHDON='CURRENT'
        ENDIF
C       LECTURE DU COURANT DE L'ONDE DE MAREE
!BD_INCKA modif pour les bons noms de fichier
!        IF((NOMCOF(1:1).EQ.' ').AND.(NOMCOB(1:1).EQ.' ')) THEN
        IF((WAC_FILES(WACCOF)%NAME(1:1).EQ.' ').AND.
     *                (WAC_FILES(WACCOB)%NAME(1:1).EQ.' ')) THEN
!BD_INCKA fin modif
            write(*,*)'fichier analytique pour courant'
          CALL ANAMAR
     *    ( SUC%R  , SVC%R   , STRA31%R, SZM1%R   ,
     *      SZM2%R , SDZHDT%R, MESH%X%R, MESH%Y%R ,
     *      NPOIN2    , AT  , DDC  , LT         ) 
          WRITE(LU,*)' ' 
          IF (LNG.EQ.1) THEN
            WRITE(LU,*)'PRISE EN COMPTE D''UN COURANT DE MAREE'
            WRITE(LU,*)
     *      'MAIS PAS DE FICHIER DES COURANTS (OU DONNEES TELEMAC)'
            WRITE(LU,*)
     *      '==> LE COURANT DE MAREE EST INITIALISE DANS ANAMAR'
          ELSE
            WRITE(LU,*)'USE OF TIDAL CURRENT VELOCITIES'
            WRITE(LU,*)'BUT NO CURRENT FILE (NEITHER TELEMAC DATA FILE)'
            WRITE(LU,*)
     *      '==> INITIALISATION OF TIDAL CURRENT VELOCITIES IN ANAMAR'
          ENDIF
        ELSE
!BD_INCKA pour les bons noms de fichier
!          IF (NOMCOF(1:1).NE.' ') THEN
          IF (WAC_FILES(WACCOF)%NAME(1:1).NE.' ') THEN
!BD_INCKA fin modif
C           LECTURE DU COURANT DE MAREE DANS FICHIER FORMATTE
            CALL LECDOI
     *      ( SUC%R  , SVC%R   , MESH%X%R, MESH%Y%R ,
!BD_INKCA modif pour la bonne unite logique
!     *        NPOIN2, NCOF , BINCOU, NBOR , NPTFR,
     *        NPOIN2, WAC_FILES(WACCOF)%LU , BINCOU, NBOR , NPTFR,
!BD_INCKA fin modif
     *        AT , DDC , TC1, TC2, NPC   ,
     *        SXRELC%R, SYRELC%R,
!BD_INCKA modif => NPOIN3=> nombre de points totaux du domaine=> NPOIN_ALL
     *        TRA01(1:NPOIN_ALL) , TRA01(NP0:NP1) , TRA01(NP2:NP3),
     *        SUC1%R , SVC1%R, SUC2%R, SVC2%R,
     *        INDIC , NPOIN_ALL, CHDON, NVCOU)
          ELSE
C           LECTURE DU COURANT DE MAREE DANS FICHIER BINAIRE
            CALL LECDOI
     *      ( SUC%R  , SVC%R   , MESH%X%R, MESH%Y%R ,
!BD_INKCA modif pour la bonne unite logique
!     *        NPOIN2, NCOB , BINCOU, NBOR , NPTFR,
     *        NPOIN2, WAC_FILES(WACCOB)%LU , BINCOU, NBOR , NPTFR,
!BD_INCKA fin modif
     *        AT , DDC , TC1, TC2, NPC   ,
     *        SXRELC%R, SYRELC%R,
!BD_INCKA modif => NPOIN3=> nombre de points totaux du domaine=> NPOIN_ALL
     *        TRA01(1:NPOIN_ALL) , TRA01(NP0:NP1) , TRA01(NP2:NP3),
     *        SUC1%R , SVC1%R, SUC2%R, SVC2%R,
     *        INDIC , NPOIN_ALL, CHDON, NVCOU )
          ENDIF
        ENDIF
C
C       LECTURE DE LA HAUTEUR DE L'ONDE DE MAREE
!BD_INCKA modif pour lire le bon fichier de maree
!        IF((NOMMAF(1:1).EQ.' ').AND.(NOMMAB(1:1).EQ.' ')) THEN
!          IF((NOMCOF(1:1).NE.' ').OR.(NOMCOB(1:1).NE.' ')) THEN
        IF((WAC_FILES(WACMAF)%NAME(1:1).EQ.' ').AND.
     *                       (WAC_FILES(WACMAB)%NAME(1:1).EQ.' ')) THEN
          IF((WAC_FILES(WACCOF)%NAME.NE.' ').OR.
     *                            (WAC_FILES(WACCOB)%NAME.NE.' ')) THEN
!BD_INCKA fin modif
            CALL ANAMAR
     *    ( SUC%R  , SVC%R   , STRA31%R, SZM1%R   ,
     *      SZM2%R , SDZHDT%R, MESH%X%R, MESH%Y%R ,
     *      NPOIN2,   AT  , DDC , LT    )
          ENDIF
          IF(LNG.EQ.1) THEN
            WRITE(LU,*)
     *      '==> LA HAUTEUR DE LA MAREE EST INITIALISEE DANS ANAMAR'
          ELSE
            WRITE(LU,*)
     *      '==> INITIALISATION OF TIDAL WATER LEVEL IN ANAMAR'
          ENDIF
        ELSE
!BD_INCKA modif pour le cas maree lire le bon fichier wac_file etc..
!          IF (NOMMAF(1:1).NE.' ') THEN
          IF (WAC_FILES(WACMAF)%NAME(1:1).NE.' ') THEN
            write(*,*)'lecture courant dans fichier formate LECHAM MAF'
            CALL LECHAM
     *      ( STRA31%R, SDZHDT%R, MESH%X%R, MESH%Y%R ,
!BD_INCKA pour lir la bon unite logique NMAF-> WAC_FILE(...)%LU
!     *        NPOIN2, NMAF , BINMAR, NBOR  ,
     *        NPOIN2, WAC_FILES(WACMAF)%LU, BINMAR, NBOR  ,
!BD_INCKA fin modif 
     *        NPTFR, AT   , DDC , TM1  , TM2   , NPM  ,
     *        SXRELM%R , SYRELM%R ,
!BD_INCKA modif => NPOIN3=> nombre de points totaux du domaine=> NPOIN_ALL
     *        TRA01(1:NPOIN_ALL)   , SZM1%R, SZM2%R,
     *        INDIM, NPOIN_ALL, IDHMA ,
     *        NVHMA )
           ELSE
            CALL LECHAM
     *      ( STRA31%R, SDZHDT%R, MESH%X%R, MESH%Y%R ,
!BD_INCKA pour lire la bon unite logique
!     *        NPOIN2, NMAB , BINMAR, NBOR  ,
     *        NPOIN2, WAC_FILES(WACMAB)%LU , BINMAR, NBOR  ,
!BD_INCKA fin modif
     *        NPTFR, AT   , DDC , TM1  , TM2   , NPM  ,
     *        SXRELM%R , SYRELM%R ,
!BD_INCKA modif => NPOIN3=> nombre de points totaux du domaine=> NPOIN_ALL
     *        TRA01(1:NPOIN_ALL)   , SZM1%R, SZM2%R,
     *        INDIM, NPOIN_ALL, IDHMA ,
     *        NVHMA )
           ENDIF
        ENDIF
        CALL OV('X=X+Y   ', SDEPTH%R, STRA31%R, ST1%R,
     *                      0.D0, NPOIN2)
      ENDIF
C
C   INITIALISATION DU COURANT ET LECTURE EVENTUELLE D'UNE VARIABLE
C   TELEMAC
C
      IF ((COUSTA).OR.(DONTEL)) THEN
!BD_INCKA modif pour lire le bon fichier de courant
!        IF ((NOMCOF(1:1).EQ.' ').AND.(NOMCOB(1:1).EQ.' ')) THEN
        IF ((WAC_FILES(WACCOF)%NAME(1:1).EQ.' ').AND.
     *                 (WAC_FILES(WACCOB)%NAME(1:1).EQ.' ')) THEN
!BD_INCKA fin modif
          IF(COUSTA) THEN
             CALL ANACOS
     *      ( SUC%R, SVC%R, MESH%X%R, MESH%Y%R, NPOIN2)
             WRITE(LU,*)' ' 
             IF (LNG.EQ.1) THEN
               WRITE(LU,*)'PRISE EN COMPTE D''UN COURANT'
               WRITE(LU,*)
     *         'MAIS PAS DE FICHIER DES COURANTS (OU DONNEES TELEMAC)'
               WRITE(LU,*)'==> LE COURANT EST INITIALISE DANS ANACOS'
             ELSE
               WRITE(LU,*)'USE OF CURRENT VELOCITIES'
               WRITE(LU,*)
     *         'BUT NO CURRENT FILE (NEITHER TELEMAC DATA FILE)'
               WRITE(LU,*)
     *         '==> INITIALISATION OF CURRENT VELOCITIES IN ANACOS'
             ENDIF
          ELSE
             IF (LNG.EQ.1) THEN
               WRITE(LU,*)'RELECTURE D''UNE VARIABLE TELEMAC IMPOSSIBLE'
             ELSE
               WRITE(LU,*)' READING OF A TELEMAC DATA IMPOSSIBLE '
             ENDIF
             CALL PLANTE(0)
          ENDIF
        ELSE
          IF(LNG.EQ.1) THEN
            CHDON='COURANT'
          ELSE
            CHDON='CURRENT'
          ENDIF
!BD_INCKA modif pour lire le bon nom de fichier
!          IF (NOMCOF(1:1).NE.' ') THEN
          IF (WAC_FILES(WACCOF)%NAME(1:1).NE.' ') THEN
!BD_INCKA fin modif
             CALL LECDON
     *      ( SUC%R  , SVC%R   , MESH%X%R, MESH%Y%R ,
!BD_INCKA pour la bonne unite logique
!     *        NPOIN2 , NCOF   , BINCOU ,
     *        NPOIN2 , WAC_FILES(WACCOF)%LU   , BINCOU ,
!BD_INCKA fin modif
!BD_INCKA modif => NPOIN3=> nombre de points totaux du domaine=> NPOIN_ALL
     *        NBOR , NPTFR, SXRELC%R, SYRELC%R, TRA01(1:NPOIN_ALL),
     *        TRA01(NP0:NP1),TRA01(NP2:NP3),TRA01(NP4:NP5),STRA31%R,
     *        IDTEL, NPTT , DONTEL, COUSTA, INDIC  , NPOIN_ALL , CHDON)
          ELSE
             CALL LECDON
     *      ( SUC%R  , SVC%R   , MESH%X%R, MESH%Y%R ,
!BD_INCKA pour la bonne unite logique
!     *        NPOIN2 , NCOB   , BINCOU ,
     *        NPOIN2 , WAC_FILES(WACCOB)%LU   , BINCOU ,
!BD_INCKA fin modif
!BD_INCKA modif => NPOIN3=> nombre de points totaux du domaine=> NPOIN_ALL
     *        NBOR , NPTFR, SXRELC%R, SYRELC%R, TRA01(1:NPOIN_ALL),
     *        TRA01(NP0:NP1),TRA01(NP2:NP3),TRA01(NP4:NP5),STRA31%R,
     *        IDTEL, NPTT , DONTEL, COUSTA, INDIC  , NPOIN_ALL , CHDON)
          ENDIF
        ENDIF
        DZHDT = 0.D0
      ENDIF
C
C-----------------------------------------------------------------------
C
C   INITIALISATION DU VENT
C
      IF (VENT) THEN
        IF(LNG.EQ.1) THEN
          CHDON='VENT   '
        ELSE
          CHDON='WIND   '
        ENDIF
C
!BD_INCKA modif pour bon noms
!        IF ((NOMVEF(1:1).EQ.' ').AND.(NOMVEB(1:1).EQ.' ')) THEN
        IF ((WAC_FILES(WACVEF)%NAME(1:1).EQ.' ').AND.
     *                 (WAC_FILES(WACVEB)%NAME(1:1).EQ.' ')) THEN
!BD_INCKA fin modif
          CALL ANAVEN
     *   ( SUV%R, SVV%R, MESH%X%R, MESH%Y%R ,
     *     NPOIN2,AT,DDC,VX_CTE,VY_CTE)
          WRITE(LU,*)' '
          IF (LNG.EQ.1) THEN
            WRITE(LU,*)'PRISE EN COMPTE D''UN VENT'
            WRITE(LU,*)'MAIS PAS DE FICHIER DE VENT'
            WRITE(LU,*)'==> LE VENT EST INITIALISE DANS ANAVEN'
          ELSE
            WRITE(LU,*)'USE OF WIND VELOCITIES'
            WRITE(LU,*)'BUT NO WIND FILE '
            WRITE(LU,*)'==> INITIALISATION OF WIND VELOCITIES IN ANAVEN'
          ENDIF 
        ELSE
!BD_INCA pour les bon noms
!          IF (NOMVEF(1:1).NE.' ') THEN
          IF (WAC_FILES(WACVEF)%NAME(1:1).NE.' ') THEN
!BD_INCKA fin modif
           IF (VENSTA) THEN
             CALL LECDON
     *      ( SUV%R  , SVV%R   , MESH%X%R, MESH%Y%R ,
!BD_INCKA pour la bonne untie logique
!     *        NPOIN2 , NVEF   , BINVEN ,
     *        NPOIN2 , WAC_FILES(WACVEF)%LU   , BINVEN ,
!BD_INCA fin modif
!BD_INCKA modif => NPOIN3=> nombre de points totaux du domaine=> NPOIN_ALL
     *        NBOR , NPTFR, SXRELV%R, SYRELV%R, TRA01(1:NPOIN_ALL),
     *        TRA01(NP0:NP1),TRA01(NP2:NP3),TRA01(NP4:NP5),STRA31%R,
     *        IDTEL, NPTT , .FALSE., VENSTA, INDIV  , NPOIN_ALL ,CHDON)
           ELSE  
            CALL LECDOI
     *      ( SUV%R , SVV%R , MESH%X%R , MESH%Y%R ,
!BD_INCKA pour la bonne unite logique
!     *        NPOIN2, NVEF , BINVEN, NBOR , NPTFR,
     *        NPOIN2, WAC_FILES(WACVEF)%LU , BINVEN, NBOR , NPTFR,
!BD_INCKA fin modif
!BD_INCKA modif => NPOIN3=> nombre de points totaux du domaine=> NPOIN_ALL
     *        AT , DDC , TV1, TV2, NPV   , SXRELV%R, SYRELV%R ,
     *        TRA01(1:NPOIN_ALL) , TRA01(NP0:NP1) , TRA01(NP2:NP3),
     *        SUV1%R , SVV1%R, SUV2%R, SVV2%R,
     *        INDIV , NPOIN_ALL, CHDON, IBID )
           ENDIF
          ELSE
           IF (VENSTA) THEN
             CALL LECDON
     *      ( SUV%R  , SVV%R   , MESH%X%R, MESH%Y%R ,
!BD_INCA modif pour bonne unite logique
!     *        NPOIN2 , NVEB   , BINVEN ,
     *        NPOIN2 , WAC_FILES(WACVEB)%LU   , BINVEN ,
!BD_INCKA fin modif
!BD_INCKA modif => NPOIN3=> nombre de points totaux du domaine=> NPOIN_ALL
     *        NBOR , NPTFR, SXRELV%R, SYRELV%R, TRA01(1:NPOIN_ALL),
     *        TRA01(NP0:NP1),TRA01(NP2:NP3),TRA01(NP4:NP5),STRA31%R,
     *        IDTEL, NPTT , .FALSE., VENSTA, INDIV  , NPOIN_ALL , CHDON)
           ELSE  
            CALL LECDOI
     *      ( SUV%R , SVV%R , MESH%X%R , MESH%Y%R ,
!BD_INCKA modif pour bonne untie logique
!     *        NPOIN2, NVEB , BINVEN, NBOR , NPTFR,
     *        NPOIN2, WAC_FILES(WACVEB)%LU, BINVEN, NBOR , NPTFR,
!BD_INCKA fin modif
     *        AT , DDC , TV1, TV2, NPV   , SXRELV%R, SYRELV%R ,
!BD_INCKA modif => NPOIN3=> nombre de points totaux du domaine=> NPOIN_ALL
     *        TRA01(1:NPOIN_ALL) , TRA01(NP0:NP1) , TRA01(NP2:NP3),
     *        SUV1%R , SVV1%R, SUV2%R, SVV2%R,
     *        INDIV , NPOIN_ALL, CHDON, IBID )
           ENDIF
          ENDIF
         ENDIF
      ENDIF      
C
C-----------------------------------------------------------------------
C
C   INITIALISATION DE TETA
C     PAR DEFAUT ON REPARTIE UNIFORMEMENT LES DIRECTIONS DE PROPAGATION
C
      DO IPLAN = 1,NPLAN+1
         TETA(IPLAN) = (IPLAN-1)*DPI/NPLAN
      ENDDO
C
C-----------------------------------------------------------------------
C
C
C     INITIALISATION DE FREQ ET DFREQ, ON REPARTIE SELON UNE LOI EXP
C     LES FREQUENCES DE PROPAGATION
C
      DO IFREQ = 1,NF
        FREQ(IFREQ) = F1*RAISF**(IFREQ-1)
      ENDDO
C
C-----------------------------------------------------------------------
C
C     INITIALISATION DE F
C
      CALL SPEINI
     *  ( SF%R  , TRA01(1:NF)   , TRA01(NP0:NP1),
     *    SUV%R , SVV%R      , SFR%R  , STETA%R  , GRAVIT,
     *    FREMAX   , FETCH , SIGMAA, SIGMAB , GAMMA  , FPIC  , HM0   ,
     *    ALPHIL   , TETA1 , SPRED1, TETA2  , SPRED2 , XLAMDA, NPOIN2,
     *    NPLAN    , NF    , INISPE, E2FMIN , DEPTH  , FRABI  )
C
C-----------------------------------------------------------------------
C
      RETURN
      END
