C                       *****************
                        SUBROUTINE LECSUI
C                       *****************
C
     *(F,NPLAN,NF,TETA,FREQ,NELEM2,NPOIN2,AT,UC,VC,UC1,VC1,UC2,VC2,
     * UV,VV,UV1,VV1,UV2,VV2,VENT,TV1,TV2,
     * COURAN,NPRE,BINPRE,TRA02,
     * DEPTH,TC1,TC2,ZM1,ZM2,DZHDT,TM1,TM2,MAREE)
C
C***********************************************************************
C TOMAWAC   V1.0            01/02/95       F MARCOS   (LNH) 30 87 72 66
C***********************************************************************
C
C      FONCTION:
C      =========
C
C    LIT LES DONNEES POUR UNE SUITE DE CALCUL
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    F           !<-- !  DENSITE SPECTRALE D'ENERGIE                 !
C !    NPLAN       ! -->!  NOMBRE DE PLANS OU DE DIRECTIONS            !
C !    NF          ! -->!  NOMBRE DE FREQUENCES                        !
C !    TETA        !<-- !  DISTRIBUTION DES DIRECTIONS                 !
C !    FREQ        !<-- !  DISTRIBUTION DES FREQUENCES                 !
C !    NELEM2      ! -->!  NOMBRE D'ELEMENTS 2D                        !
C !    NPOIN2      ! -->!  NOMBRE DE POINTS DU MAILLAGE 2D             !
C !    AT          ! -->!  TEMPS                                       !
C !    U,V         !<-- !  COMPOSANTES DU COURANT                      !
C !    U1,V1       !<-- !  COMPOSANTES DU VENT                         !
C !    U2,V2       !<-- !  COMPOSANTES DU VENT                         !
C !    UV,VV       !<-- !  COMPOSANTES DU VENT                         !
C !    VENT        ! -->!  LOGIQUE INDIQUANT SI IL YA UN VENT          !
C !    TV1,2       ! -->!  TEMPS DES ENREGISTREMENTS DE VENT 1 ET 2    !
C !    COURAN      ! -->!  LOGIQUE INDIQUANT SI IL YA UN COURANT       !
C !    NPRE        ! -->!  NUMERO D'UNITE LOGIQUE DU FICHIER DU CALCUL !
C !                !    !  PRECEDENT                                   !
C !    BINPRE      ! -->!  BINAIRE DU FICHIER DU CALCUL PRECEDENT      !
C !    TRA02       !<-->!  TABLEAU DE TRAVAIL                          !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE PAR : WAC
C SOUS-PROGRAMMES APPELES : LIT , OV
C
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
      INTEGER NPRE,NF,NPLAN,NELEM2,NPOIN2
      INTEGER I,ISTAT,IB(2),NTOT
C
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF),AT,ATT(1)
      DOUBLE PRECISION TETA(NPLAN+1),FREQ(NF),PI,TV1,TV2,Z(1)
      DOUBLE PRECISION UC(NPOIN2),VC(NPOIN2),UV(NPOIN2),VV(NPOIN2)
      DOUBLE PRECISION UV1(NPOIN2),VV1(NPOIN2),UV2(NPOIN2),VV2(NPOIN2)
      DOUBLE PRECISION UC1(NPOIN2),VC1(NPOIN2),UC2(NPOIN2),VC2(NPOIN2)
      DOUBLE PRECISION DEPTH(NPOIN2),ZM1(NPOIN2),ZM2(NPOIN2)
      DOUBLE PRECISION DZHDT(NPOIN2)
      DOUBLE PRECISION TC1,TC2,TM1,TM2
C
      DOUBLE PRECISION TRA02(*)
C
      LOGICAL COURAN,VENT,MAREE
C
      CHARACTER*3 BINPRE
      CHARACTER*72 CAR
C
      REAL, ALLOCATABLE :: W(:)
      ALLOCATE(W(NPOIN2*NPLAN*NF))
C
C***********************************************************************
C
      PI=3.141592D0
      REWIND NPRE
C
C LECTURE DU TITRE
C
      CALL LIT(F,W,IB,CAR,72,'CH',NPRE,BINPRE,ISTAT)
      WRITE(LU,*) ' '
      IF(LNG.EQ.1) THEN
         WRITE(LU,*) '**** SUITE DE CALCUL ****'
         WRITE(LU,*) ' '
         WRITE(LU,*) 'TITRE DU CALCUL PRECEDENT :'
         WRITE(LU,*) '     ',CAR
      ELSE
         WRITE(LU,*) '**** FOLLOWING COMPUTATION ****'
         WRITE(LU,*) ' '
         WRITE(LU,*) 'TITLE OF THE PREVIOUS COMPUTATION :'
      ENDIF
      WRITE(LU,*) '     ',CAR
C
C LECTURE DE NPLAN, NF ET VERIFICATION
C
      CALL LIT(F,W,IB,CAR,2,'I ',NPRE,BINPRE,ISTAT)
      IF ((IB(1).NE.NPLAN).OR.(IB(2).NE.NF)) THEN
       IF(LNG.EQ.1) THEN
         WRITE(LU,*) '**** ERREUR DANS LECSUI : ****'
         WRITE(LU,*) 'LE NOMBRE DE DIRECTIONS ET/OU CELUI DE FREQUENCES'
         WRITE(LU,*) ' NE CORRESPOND PAS'
         WRITE(LU,*) 'VALEURS LUES : NDIR=',IB(1),' NF=',IB(2)
         WRITE(LU,*) 'VALEURS ATTENDUES : NDIR=',NPLAN,' NF=',NF
       ELSE
         WRITE(LU,*) '**** ERROR IN LECSUI : ****'
         WRITE(LU,*) 'THE NUMBER OF DIRECTIONS AND/OR FREQUENCIES'
         WRITE(LU,*) '   IS NOT CORRESPONDING '
         WRITE(LU,*) 'READ VALUES : NDIR=',IB(1),' NF=',IB(2)
         WRITE(LU,*) 'EXPECTED VALUES : NDIR=',NPLAN,' NF=',NF
       ENDIF
       CALL PLANTE(0)
      ENDIF
C
C LECTURE DE NELEM2, NPOIN2 ET VERIFICATION
C
      CALL LIT(F,W,IB,CAR,2,'I ',NPRE,BINPRE,ISTAT)
      IF ((IB(1).NE.NELEM2).OR.(IB(2).NE.NPOIN2)) THEN
       IF(LNG.EQ.1) THEN
         WRITE(LU,*) '**** ERREUR DANS LECSUI ****'
         WRITE(LU,*) 'LE NOMBRE DE POINTS ET/OU CELUI D''ELEMENTS 2D NE'
         WRITE(LU,*) 'CORRESPOND PAS'
         WRITE(LU,*) 'VALEURS LUES : NELEM2=',IB(1),' NPOIN2=',IB(2)
         WRITE(LU,*) 'VALEURS ATTENDUES : NELEM2=',NELEM2,
     *               ' NPOIN2=',NPOIN2
       ELSE
         WRITE(LU,*) '**** ERROR IN LECSUI : ****'
         WRITE(LU,*) 'THE NUMBER OF POINTS AND/OR 2D ELEMENTS '
         WRITE(LU,*) '   IS NOT CORRESPONDING '
         WRITE(LU,*) 'READ VALUES     : NELEM2=',IB(1),' NPOIN2=',IB(2)
         WRITE(LU,*) 'EXPECTED VALUES : NELEM2=',NELEM2,
     *               ' NPOIN2=',NPOIN2
       ENDIF
       CALL PLANTE(0)
      ENDIF
C
C LECTURE DU TEMPS
C
      CALL LIT(ATT,W,IB,CAR,1,'R4',NPRE,BINPRE,ISTAT)
      AT = ATT(1)
      IF(LNG.EQ.1) THEN
         WRITE(LU,*) '- REPRISE DE CALCUL AU TEMPS  ',AT
      ELSE
         WRITE(LU,*) '- COMPUTATIONAL RESUMPTION AT TIME ',AT
      ENDIF
C
C LECTURE DE TETA
C
      CALL LIT(TETA,W,IB,CAR,NPLAN,'R4',NPRE,BINPRE,ISTAT)
      TETA(NPLAN+1)=2.D0*PI+TETA(1)
      IF(LNG.EQ.1) THEN
         WRITE(LU,*) 'DISTRIBUTION DES DIRECTIONS :'
      ELSE
         WRITE(LU,*) 'DISTRIBUTION OF THE DIRECTIONS:'
      ENDIF
      DO I=1,NPLAN
        WRITE(LU,*) '       ',TETA(I)*180/PI,' Degres'
      ENDDO
C
C LECTURE DE FREQ
C
      CALL LIT(FREQ,W,IB,CAR,NF,'R4',NPRE,BINPRE,ISTAT)
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'DISTRIBUTION DES FREQUENCES :'
      ELSE
        WRITE(LU,*) 'DISTRIBUTION OF THE FREQUENCIES:'
      ENDIF
      DO I=1,NF
        WRITE(LU,*) '       ',FREQ(I),' Hertz'
      ENDDO
C
C LECTURE DE F
C
      NTOT=NPOIN2*NPLAN*NF
      CALL LIT(F,W,IB,CAR,NTOT,'R4',NPRE,BINPRE,ISTAT)
C
C LECTURE CONDITIONNELLE DE U,V,UV,VV
C
      IF (COURAN) THEN
      CALL LIT(UC ,W,IB,CAR,NPOIN2,'R4',NPRE,BINPRE,ISTAT)
      CALL LIT(VC ,W,IB,CAR,NPOIN2,'R4',NPRE,BINPRE,ISTAT)
C     
C       ON MET LES TRIPLETS U,V,TV1 ET 2 a UV,VV AT
C
      TC1=AT
      TC2=AT
        CALL OV( 'X=Y     ' , UC1 , UC , Z , 0.D0 , NPOIN2)
        CALL OV( 'X=Y     ' , UC2 , UC , Z , 0.D0 , NPOIN2)
        CALL OV( 'X=Y     ' , VC1 , VC , Z , 0.D0 , NPOIN2)
        CALL OV( 'X=Y     ' , VC2 , VC , Z , 0.D0 , NPOIN2)
      ENDIF
C
      IF (VENT) THEN
      CALL LIT(UV,W,IB,CAR,NPOIN2,'R4',NPRE,BINPRE,ISTAT)
      CALL LIT(VV,W,IB,CAR,NPOIN2,'R4',NPRE,BINPRE,ISTAT)
C     
C       ON MET LES TRIPLETS U,V,TV1 ET 2 a UV,VV AT
C
      TV1=AT
      TV2=AT
        CALL OV( 'X=Y     ' , UV1 , UV , Z , 0.D0 , NPOIN2)
        CALL OV( 'X=Y     ' , UV2 , UV , Z , 0.D0 , NPOIN2)
        CALL OV( 'X=Y     ' , VV1 , VV , Z , 0.D0 , NPOIN2)
        CALL OV( 'X=Y     ' , VV2 , VV , Z , 0.D0 , NPOIN2)
      ENDIF
C
      IF (MAREE) THEN
      CALL LIT(DEPTH ,W,IB,CAR,NPOIN2,'R4',NPRE,BINPRE,ISTAT)
C     
C       ON MET LES TRIPLETS U,V,TV1 ET 2 a UV,VV AT
C
      TM1=AT
      TM2=AT
        CALL OV( 'X=Y     ' , ZM1 , DEPTH , Z , 0.D0 , NPOIN2)
        CALL OV( 'X=Y     ' , ZM2 , DEPTH , Z , 0.D0 , NPOIN2)
        CALL OV( 'X=C     ' , DZHDT , DEPTH , Z , 0.D0 , NPOIN2)
      ENDIF
C
C-----------------------------------------------------------------------
C
      DEALLOCATE(W)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
