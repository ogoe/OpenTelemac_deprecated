C                       ****************
                        SUBROUTINE SOR3D
C                       ****************
C
     *(F,NPLAN,NF,TETA,FREQ,NELEM2,NPOIN2,AT,U,V,UV,VV,DEPTH,VENT,
     * COURAN,MAREE,TITRE,NR3D,BINR3D)
C
C***********************************************************************
C TOMAWAC   V1.0            01/02/95       F MARCOS   (LNH) 30 87 72 66
C***********************************************************************
C
C      FONCTION:
C      =========
C
C    ECRIT LES DONNEES POUR UNE SUITE DE CALCUL ULTERIEURE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    F           !<-- !  DENSITE SPECTRALE D'ENERGIE                 !
C !    NPLAN       ! -->!  NOMBRE DE PLANS OU DE DIRECTIONS            !
C !    NF          ! -->!  NOMBRE DE FREQUENCES                        !
C !    TETA        ! -->!  DISTRIBUTION DES DIRECTIONS                 !
C !    FREQ        ! -->!  DISTRIBUTION DES FREQUENCES                 !
C !    NELEM2      ! -->!  NOMBRE D'ELEMENTS 2D                        !
C !    NPOIN2      ! -->!  NOMBRE DE POINTS DU MAILLAGE 2D             !
C !    AT          ! -->!  TEMPS                                       !
C !    U,V         ! -->!  COMPOSANTES DU COURANT                      !
C !    UV,VV       ! -->!  COMPOSANTES DU VENT                         !
C !    VENT        ! -->!  LOGIQUE INDIQUANT SI IL YA UN VENT          !
C !    COURAN      ! -->!  LOGIQUE INDIQUANT SI IL YA UN COURANT       !
C !    TITRE       ! -->!  TITRE DU CAS                                !
C !    NR3D        ! -->!  NUMERO D'UNITE LOGIQUE DU FICHIER DES       !
C !                !    !  RESULTATS GLOBAUX                           !
C !    BINR3D      ! -->!  BINAIRE DU FICHIER DES RESULTATS GLOBAUX    !
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
      IMPLICIT NONE
C
      INTEGER NR3D,NF,NPLAN,NELEM2,NPOIN2
C
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF),AT, ATT(1)
      DOUBLE PRECISION FREQ(NF),TETA(NPLAN)
      DOUBLE PRECISION U(NPOIN2),V(NPOIN2),UV(NPOIN2),VV(NPOIN2)
      DOUBLE PRECISION DEPTH(NPOIN2)
C
      INTEGER ISTAT,IB(2),NTOT
C
      LOGICAL COURAN,VENT,MAREE
C
      CHARACTER*3 BINR3D,CAR
      CHARACTER*72 TITRE
C
C***********************************************************************
C
C
C ECRITURE DU TITRE
C
      CALL ECRI2(F,IB,TITRE,80,'CH',NR3D,BINR3D,ISTAT)
C
C ECRITURE DE NPLAN, NF
C
      IB(1)=NPLAN
      IB(2)=NF
      CALL ECRI2(F,IB,CAR,2,'I ',NR3D,BINR3D,ISTAT)
C
C ECRITURE DE NELEM2, NPOIN2
C
      IB(1)=NELEM2
      IB(2)=NPOIN2
      CALL ECRI2(F,IB,CAR,2,'I ',NR3D,BINR3D,ISTAT)
C
C ECRITURE DU TEMPS
C
      ATT(1)=AT
      CALL ECRI2(ATT,IB,CAR,1,'R4',NR3D,BINR3D,ISTAT)
C
C ECRITURE DE TETA
C
      CALL ECRI2(TETA,IB,CAR,NPLAN,'R4',NR3D,BINR3D,ISTAT)
C
C ECRITURE DE FREQ
C
      CALL ECRI2(FREQ,IB,CAR,NF,'R4',NR3D,BINR3D,ISTAT)
C
C ECRITURE DE F
C
      NTOT=NPOIN2*NPLAN*NF
      CALL ECRI2(F,IB,CAR,NTOT,'R4',NR3D,BINR3D,ISTAT)
C
C ECRITURE CONDITIONNELLE DE U,V,UV,VV
C
      IF (COURAN) THEN
      CALL ECRI2(U ,IB,CAR,NPOIN2,'R4',NR3D,BINR3D,ISTAT)
      CALL ECRI2(V ,IB,CAR,NPOIN2,'R4',NR3D,BINR3D,ISTAT)
      ENDIF
      IF (VENT) THEN
      CALL ECRI2(UV,IB,CAR,NPOIN2,'R4',NR3D,BINR3D,ISTAT)
      CALL ECRI2(VV,IB,CAR,NPOIN2,'R4',NR3D,BINR3D,ISTAT)
      ENDIF
      IF (MAREE) THEN
      CALL ECRI2(DEPTH,IB,CAR,NPOIN2,'R4',NR3D,BINR3D,ISTAT)
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
