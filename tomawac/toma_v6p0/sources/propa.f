C                       *****************
                        SUBROUTINE PROPA
C                       *****************
C
     *   (F,B,SHP1,SHP2,SHP3,SHZ,SHF,ELT,ETA,FRE,IKLE2,ETAP1,
     *    NPOIN3,NPOIN2,NELEM2,NPLAN,NF,COURAN,TRA01,TRA02)
C
C***********************************************************************
C  TOMAWAC VERSION 1.0    5/12/95           F. MARCOS (LNH) 30 87 72 66
C***********************************************************************
C
C     FONCTION  : ETAPE DE CONVECTION 
C         INTERPOLE AUX PIEDS DES CARACTERISTIQUES
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    F           !<-->!  FONCTION A CONVECTER                        !
C !    B           ! -->!  FACTEUR B                                   !
C !    SHP1-2-3    !<-->! COORDONNEES BARYCENTRIQUES 2D AU PIED DES    !
C !                !    ! COURBES CARACTERISTIQUES.                    !
C !    SHT         !<-->! COORDONNEES BARYCENTRIQUES SUIVANT TETA DES  !
C !                !    ! NOEUDS DANS LEURS ETAGES "ETA" ASSOCIES.     !
C !    SHF         !<-->! COORDONNEES BARYCENTRIQUES SUIVANT F DES     !
C !                !    ! NOEUDS DANS LEURS FREQUENCES "FRE" ASSOCIEES.!
C !    ELT         !<-->! NUMEROS DES ELEMENTS 2D CHOISIS POUR CHAQUE  !
C !                !    ! NOEUD.                                       !
C !    ETA         !<-->! NUMEROS DES ETAGES CHOISIS POUR CHAQUE NOEUD.!
C !    FRE         !<-->! NUMEROS DES FREQ. CHOISIES POUR CHAQUE NOEUD.!
C !    IKLE2       ! -->! TRANSITION ENTRE LES NUMEROTATIONS LOCALE    !
C !                !    ! ET GLOBALE DU MAILLAGE 2D.                   !
C !    ETAP1       !<-->! TABLEAU DE TRAVAIL DONNANT LE NUMERO DE      !
C !                !    ! L'ETAGE SUPERIEUR                            !
C !    NPOIN3      ! -->!  NOMBRE DE POINTS DU MAILLAGE 3D             !
C !    NPOIN2      ! -->!  NOMBRE DE POINTS DU MAILLAGE 2D             !
C !    NELEM2      ! -->!  NOMBRE D ELEMENTS DU MAILLAGE 2D            !
C !    NPLAN       ! -->!  NOMBRE DE PLANS OU DE DIRECTIONS            !
C !    NF          ! -->!  NOMBRE DE FREQUENCES                        !
C !    COURAN      ! -->!  LOGIQUE INDIQUANT SI IL Y A DU COURANT      !
C !    TRA01       !<-->!  TABLEAU DE TRAVAIL                          !
C !    TRA02       !<-->!  TABLEAU DE TRAVAIL                          !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : WAC
C
C  SOUS-PROGRAMME APPELE : INTERP,INTER4D,LIT
C
C***********************************************************************
C
      USE BIEF
      USE TOMAWAC_MPI
C
      IMPLICIT NONE
C
      INTEGER NPOIN3,NPOIN2,NELEM2,NPLAN,NF
C
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF) 
      DOUBLE PRECISION SHP1(NPOIN3,NF) , SHP2(NPOIN3,NF)
      DOUBLE PRECISION SHP3(NPOIN3,NF) , SHZ(NPOIN3,NF)
      DOUBLE PRECISION SHF(NPOIN3,NF)
      DOUBLE PRECISION B(NPOIN2,NF)
      DOUBLE PRECISION TRA01(NPOIN3,8),TRA02(NPOIN2,NPLAN,NF)
      INTEGER ELT(NPOIN3,NF),ETA(NPOIN3,NF),FRE(NPOIN3,NF)
      INTEGER IKLE2(NELEM2,3),ETAP1(NPLAN)
      LOGICAL COURAN
      REAL WW(1)
C
      DOUBLE PRECISION X(1)
      INTEGER IFF,I,ISTAT,LU, IB(1)
      CHARACTER*3 CAR
C
C----------------------------------------------------------------------
C
      LU=6
C
        IF (.NOT.COURAN) THEN
C
         DO 300 IFF=1,NF
C
            IFREQ = IFF
            CALL INTERP_TOMAWAC
     *        (F(1,1,IFF),B(1,IFF),SHP1(1,IFF),SHP2(1,IFF),
     *             SHP3(1,IFF),SHZ(1,IFF),ELT(1,IFF),ETA(1,IFF),IKLE2,
     *         ETAP1,NPOIN2,NELEM2,NPLAN,TRA01)
C
300      CONTINUE
C
        ELSE
C
            CALL INTER4D
     *       (F,B,SHP1,SHP2,SHP3,SHZ,SHF,ELT,ETA,
     *        FRE,IKLE2,ETAP1,NPOIN2,NELEM2,NPLAN,NF,TRA02)
C
        ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
