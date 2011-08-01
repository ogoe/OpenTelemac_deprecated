C                       *****************
                        SUBROUTINE CORMAR
C                       *****************
     *( AT    , LT    , TC1   , TC2   , TV1   , TV2   , TM1   , TM2   ,
     *  NPC   , NPM   , NVHMA , NVCOU )
C
C***********************************************************************
C TOMAWAC     V5.0         25/08/00
C***********************************************************************
C
C      FONCTION:
C      =========
C
C    ACTUALISATION DES TABLEAUX DES GRANDEURS PHYSIQUES
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    F           !<-- ! DENSITE SPECTRALE D'ENERGIE                  !
C !    PRIVE       ! -->! TABLEAU POUR L'UTILISATEUR DE DIMENSION      !
C !    NPRIV       !    ! NPOIN3*NPRIV                                 !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE PAR : WAC
C SOUS-PROGRAMMES APPELES :   NOUDON , NOUMAR , ANAMAR , OV
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
      INTEGER          NPC , NPM, NVHMA, NVCOU
      INTEGER          LT
      DOUBLE PRECISION AT, TC1, TC2 , TV1, TV2, TM1 , TM2
C
C     VARIABLES LOCALES
      INTEGER N1,N2,N3,N4
!BD_INCKA // on modifier la taille des tableaux, on agrandit car
! certain tableau sont sousestimé du à ncsize
      INTEGER NPOIN4
C
C-----------------------------------------------------------------------
C         MISE A JOUR DES TABLEAUX DE COURANT ET DE HAUTEUR DE MAREE
C       ==============================================================
C
C            MISE A JOUR DU COURANT POUR LA DATE AT
C         ---------------------------------------------
C
!BD_INCKA modif //
      NPOIN4 = MAX(NPOIN3,NPOIN3*NCSIZE*3)
      N1=NPOIN4+1
      N2=2*NPOIN4
      N3=N2+1
      N4=3*NPOIN4
!      N1=NPOIN3+1
!      N2=2*NPOIN3
!      N3=N2+1
!      N4=3*NPOIN3
!BD_INCKA dans toute cette subroutine on replace NPOIN3 par NPOIN4
C
!BD_INCKA modif pour bon format
!      IF (NOMCOB(1:1).NE.' ') THEN
      IF (WAC_FILES(WACCOB)%NAME(1:1).NE.' ') THEN
!BD_INCKA fin modif
!        write(*,*)'on lit des courant dans le fichier U'
        CALL NOUDON
     * ( SUC%R , SVC%R  , MESH%X%R, MESH%Y%R, NPOIN2     ,
!BD_INCKA bon format
!     *   NCOB     , BINCOU    , NBOR       , NPTFR , AT , DDC        ,
     *   WAC_FILES(WACCOB)%LU , BINCOU , NBOR      , NPTFR , AT , DDC ,
!BD_INCA fin midf
     *   TC1      , TC2       , NPC        , SXRELC%R, SYRELC%R,
     *   TRA01(1:NPOIN4),TRA01(N1:N2),TRA01(N3:N4) ,
     *   SUC1%R, SVC1%R , SUC2%R  , SVC2%R  , INDIC  ,
     *   'COURANT', NVCOU     )

!BD_INCKA modif bon nom de fichier
!      ELSEIF (NOMCOF(1:1).NE.' ') THEN
      ELSEIF (WAC_FILES(WACCOF)%NAME(1:1).NE.' ') THEN
!BD_INCKA fin modif
!        write(*,*)'on lit des courant dans le fichier F'
        CALL NOUDON
     * ( SUC%R , SVC%R  , MESH%X%R, MESH%Y%R, NPOIN2     ,
!BD_INCKA bon format
!     *   NCOF     , BINCOU    , NBOR       , NPTFR , AT , DDC        ,
     *   WAC_FILES(WACCOF)%LU  , BINCOU , NBOR   , NPTFR , AT , DDC   ,
!BD_INCKA fin modif
     *   TC1      , TC2       , NPC        , SXRELC%R, SYRELC%R,
     *   TRA01(1:NPOIN4),TRA01(N1:N2),TRA01(N3:N4) , 
     *   SUC1%R, SVC1%R , SUC2%R  , SVC2%R  , INDIC  ,
     *   'COURANT', NVCOU )
      ELSE
        CALL ANAMAR
     * ( SUC%R   , SVC%R   , TRA01(1:NPOIN2), ZM1 , ZM2 ,
     *   SDZHDT%R, MESH%X%R, MESH%Y%R    ,
     *   NPOIN2     , AT   , DDC , LT )
      ENDIF
C
C            MISE A JOUR DE LA HAUTEUR D'EAU POUR LA DATE AT
C         ------------------------------------------------------
C
!BD_INCKA modif pour bon format
!      IF (NOMMAB(1:1).NE.' ') THEN
      IF (WAC_FILES(WACMAB)%NAME(1:1).NE.' ') THEN
!BD_INCKA fin modif 
        CALL NOUMAR
     * (TRA01(1:NPOIN2) , SDZHDT%R, MESH%X%R , MESH%Y%R ,
!BD_INCKA bon format
!     *  NPOIN2  , NMAB  , BINMAR     , NBOR , NPTFR, AT    , DDC ,
     *  NPOIN2,WAC_FILES(WACMAB)%LU,BINMAR,NBOR,NPTFR,AT    , DDC ,
!BD_INCKA fin modif
     *  TM1     , TM2   , NPM        , SXRELM%R , SYRELM%R ,
     *  TRA01(N1:N2), TRA01(N3:N4), ZM1, ZM2 , INDIM , IDHMA , NVHMA)
!BD_INCKA modif bon nom
!      ELSEIF (NOMMAF(1:1).NE.' ') THEN
      ELSEIF (WAC_FILES(WACMAF)%NAME(1:1).NE.' ') THEN
!BD_INCKA fin modif 
        CALL NOUMAR
     * (TRA01(1:NPOIN2) , SDZHDT%R, MESH%X%R , MESH%Y%R ,
!BD_INCKA modif bon format
!     *  NPOIN2  , NMAF  , BINMAR     , NBOR , NPTFR, AT    , DDC ,
     *  NPOIN2 ,WAC_FILES(WACMAF)%LU,BINMAR,NBOR,NPTFR,AT   , DDC ,
!BD_INCKA fin modif
     *  TM1     , TM2   , NPM        , SXRELM%R , SYRELM%R ,
     *  TRA01(N1:N2), TRA01(N3:N4), ZM1, ZM2 , INDIM , IDHMA , NVHMA)
      ELSE
!BD_INCKA modif bon nom
!        IF((NOMCOF(1:1).NE.' ').OR.(NOMCOB(1:1).NE.' ')) THEN
        IF((WAC_FILES(WACCOF)%NAME(1:1).NE.' ').OR.
     *                    (WAC_FILES(WACCOB)%NAME(1:1).NE.' ')) THEN
!BD_INCKA fin modif
!         write(*,*)'on lit de la maree dans le fichier analytique'
          CALL ANAMAR 
     *   ( SUC%R   , SVC%R   , TRA01(1:NPOIN2), ZM1 , ZM2 ,
     *     SDZHDT%R, MESH%X%R, MESH%Y%R    , NPOIN2    ,
     *     AT   , DDC , LT )
       ENDIF
      ENDIF
C
      CALL OV('X=X+Y   ', SDEPTH%R , TRA01(1:NPOIN2) , ST0%R ,
     *         0.D0 , NPOIN2)
C
C
C            MISE A JOUR DES GRADIENTS DE COURANT ET
C                DE PROFONDEUR D'EAU A LA DATE AT
C         ------------------------------------------------------
C
C W1 ( ex MASKEL) EST MIS A 1 POUR GRADF
C
      CALL OV ( 'X=C     ' , SW1%R , ST0%R , ST1%R ,
     *          1.D0 , NELEM2 )
C
      CALL VECTOR(ST1,'=','GRADF          X',IELM2,1.D0,SDEPTH,
     * ST0,ST0,ST0,ST0,ST0,MESH,.FALSE.,SW1)
C
      CALL VECTOR(ST2,'=','GRADF          X',IELM2,1.D0,SUC,
     * ST0,ST0,ST0,ST0,ST0,MESH,.FALSE.,SW1)
C
      CALL VECTOR(ST3,'=','GRADF          X',IELM2,1.D0,SVC,
     * ST0,ST0,ST0,ST0,ST0,MESH,.FALSE.,SW1)
C
      CALL VECTOR(ST4,'=','GRADF          X',IELM2,1.D0,MESH%X,
     * ST0,ST0,ST0,ST0,ST0,MESH,.FALSE.,SW1)
!BD_INCKA modif //
       IF(NCSIZE.GT.1) THEN
          CALL PARCOM(ST1,2,MESH)
          CALL PARCOM(ST2,2,MESH) 
          CALL PARCOM(ST3,2,MESH)
          CALL PARCOM(ST4,2,MESH)      
       ENDIF
!BD_INCKA fin modif //
C
      CALL OV('X=Y/Z   ',SDZX%R,ST1%R,ST4%R,0.D0,NPOIN2)
      CALL OV('X=Y/Z   ',SDUX%R,ST2%R,ST4%R,0.D0,NPOIN2)
      CALL OV('X=Y/Z   ',SDVX%R,ST3%R,ST4%R,0.D0,NPOIN2)
C
      CALL VECTOR(ST1,'=','GRADF          Y',IELM2,1.D0,SDEPTH,
     * ST0,ST0,ST0,ST0,ST0,MESH,.FALSE.,SW1)
C
      CALL VECTOR(ST2,'=','GRADF          Y',IELM2,1.D0,SUC,
     * ST0,ST0,ST0,ST0,ST0,MESH,.FALSE.,SW1)
C
      CALL VECTOR(ST3,'=','GRADF          Y',IELM2,1.D0,SVC,
     * ST0,ST0,ST0,ST0,ST0,MESH,.FALSE.,SW1)
C
      CALL VECTOR(ST4,'=','GRADF          Y',IELM2,1.D0,MESH%Y,
     * ST0,ST0,ST0,ST0,ST0,MESH,.FALSE.,SW1)
!BD_INCKA modif //
       IF(NCSIZE.GT.1) THEN
          CALL PARCOM(ST1,2,MESH)
          CALL PARCOM(ST2,2,MESH) 
          CALL PARCOM(ST3,2,MESH)
          CALL PARCOM(ST4,2,MESH)      
       ENDIF
!BD_INCKA fin modif //
C
      CALL OV('X=Y/Z   ',SDZY%R,ST1%R,ST4%R,0.D0,NPOIN2)
      CALL OV('X=Y/Z   ',SDUY%R,ST2%R,ST4%R,0.D0,NPOIN2)
      CALL OV('X=Y/Z   ',SDVY%R,ST3%R,ST4%R,0.D0,NPOIN2)
C     
C
C-----------------------------------------------------------------------
C
      RETURN
      END
