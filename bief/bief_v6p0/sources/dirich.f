C                       *****************
                        SUBROUTINE DIRICH
C                       *****************
C
     *(F, S, SM , FBOR,LIMDIR,WORK,MESH,KDIR,MSK,MASKPT)
C
C***********************************************************************
C BIEF VERSION 5.9      11/07/2008    J-M HERVOUET (LNHE) 01 30 87 80 18
C
C
C 11/07/2008 : DIMENSION DE LIMPRO MODIFIEE SI ELEMENTS QUADRATIQUES
C
C***********************************************************************
C
C FONCTION :
C
C     PRISE EN COMPTE DES POINTS DE TYPE DIRICHLET DANS UN SYSTEME
C     D'EQUATIONS LINEAIRES AVEC MATRICE SYMETRIQUE.
C
C     DANS LES EQUATIONS DES POINTS QUI NE SONT PAS DE TYPE DIRICHLET :
C     LES VALEURS DIRICHLETS, DONC CONNUES SONT ELIMINEES
C
C     DANS LES EQUATIONS DES POINTS QUI SONT DE TYPE DIRICHLET :
C     ON MET UNE EQUATION FIXANT LA VALEUR IMPOSEE.
C
C     ATTENTION : POUR LES SYSTEMES DE BLOCS DE MATRICES
C                 LES MATRICES EXTRA-DIAGONALES DOIVENT ETRE
C                 NON SYMETRIQUES CAR LA PRISE EN COMPTE DES POINTS
C                 DE TYPE DIRICHLET LES DISSYMETRISE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   F            |<-->|  VALEURS A L'ETAPE N+1 ET INITIALISATION
C |   S            |<-->|  MATRICE DU SYSTEME
C |   SM           | -->|  SECOND MEMBRE DU SYSTEME.
C |   FBOR         | -->|  CONDITIONS AUX LIMITES DES POINTS DIRICHLET.
C |   LIMDIR       | -->|  TYPES DE CONDITIONS AUX LIMITES .
C |                |    |  SI LIMDIR(K) = KDIR LE KIEME POINT DE BORD
C |                |    |  EST DE TYPE DIRICHLET.
C |                |    |  DIMENSION LIMDIR(LIMDIM,6)
C |                |    |  LIMDIM EST NPTFR OU 2*NPTFR (EN QUADRATIQUE)
C |   WORK         | -->|  BLOC DE TABLEAUX DE TRAVAIL.
C |   MESH         | -->|  BLOC DES TABLEAUX D'ENTIERS DU MAILLAGE.
C |   KDIR         | -->|  CONVENTION POUR LES CONDITIONS DE DIRICHLET
C |   MSK          | -->| SI OUI, PRESENCE D'ELEMENTS MASQUES.         |
C |   MASKPT       | -->| TABLEAU DE MASQUAGE DES POINTS               |
C |                |    |  =1. : NORMAL   =0. : POINT MASQUE.          |
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C----------------------------------------------------------------------
C
C  PROGRAMMES APPELES : OV , MATMAT
C
C**********************************************************************
C
      USE BIEF, EX_DIRICH => DIRICH
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C                                                   DIMLIM,6
      INTEGER        , INTENT(IN)    :: KDIR,LIMDIR(*)
      LOGICAL        , INTENT(IN)    :: MSK
      TYPE(BIEF_OBJ) , INTENT(INOUT) :: WORK
      TYPE(BIEF_OBJ) , INTENT(INOUT) :: F,SM,S
      TYPE(BIEF_OBJ) , INTENT(IN)    :: FBOR,MASKPT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER DIMLIM
C
C----------------------------------------------------------------------
C
C  SI S EST UNE MATRICE
C
      IF(S%TYPE.EQ.3) THEN
C
      CALL DIRI01(F, S, SM , FBOR,LIMDIR,WORK%ADR(1)%P,WORK%ADR(2)%P,
     *            MESH,KDIR,MSK,MASKPT)
C
C  SI S EST UN BLOC DE 4 MATRICES
C
      ELSEIF(S%TYPE.EQ.4.AND.S%N.EQ.4) THEN
C
      DIMLIM=MAX(FBOR%ADR(1)%P%DIM1,
     *           FBOR%ADR(2)%P%DIM1)
C
      CALL DIRI04(F%ADR(1)%P,F%ADR(2)%P,
     *     S%ADR(1)%P,S%ADR(2)%P,S%ADR(3)%P,S%ADR(4)%P,
     *     SM%ADR(1)%P,SM%ADR(2)%P,
     *     WORK%ADR(1)%P,WORK%ADR(2)%P,WORK%ADR(3)%P,WORK%ADR(4)%P,
     *     FBOR%ADR(1)%P,FBOR%ADR(2)%P,
     *     LIMDIR(1:DIMLIM),LIMDIR(DIMLIM+1:2*DIMLIM),
     *     MESH,KDIR,MSK,MASKPT)
C
C  SI S EST UN BLOC DE 9 MATRICES
C
      ELSEIF(S%TYPE.EQ.4.AND.S%N.EQ.9) THEN
C
      DIMLIM=MAX(FBOR%ADR(1)%P%DIM1,
     *           FBOR%ADR(2)%P%DIM1,
     *           FBOR%ADR(3)%P%DIM1)
C
      CALL DIRI09(F%ADR(1)%P,F%ADR(2)%P,F%ADR(3)%P,
     *            S%ADR(1)%P,S%ADR(2)%P,S%ADR(3)%P,
     *            S%ADR(4)%P,S%ADR(5)%P,S%ADR(6)%P,
     *            S%ADR(7)%P,S%ADR(8)%P,S%ADR(9)%P,
     *            SM%ADR(1)%P,SM%ADR(2)%P,SM%ADR(3)%P,
     *            WORK%ADR(1)%P,WORK%ADR(2)%P,WORK%ADR(3)%P,
     *            WORK%ADR(4)%P,WORK%ADR(5)%P,WORK%ADR(6)%P,
     *            FBOR%ADR(1)%P,FBOR%ADR(2)%P,FBOR%ADR(3)%P,
     *            LIMDIR(         1:  DIMLIM),
     *            LIMDIR(  DIMLIM+1:2*DIMLIM),
     *            LIMDIR(2*DIMLIM+1:3*DIMLIM),
     *            MESH,KDIR,MSK,MASKPT)
C
C  ERREUR
C
      ELSE
C
         IF (LNG.EQ.1) WRITE(LU,1000) S%TYPE
         IF (LNG.EQ.2) WRITE(LU,1001) S%TYPE
1000     FORMAT(1X,'DIRICH (BIEF) : MAUVAIS TYPE POUR S :',1I6)
1001     FORMAT(1X,'DIRICH (BIEF): WRONG TYPE FOR S:',1I6)
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C----------------------------------------------------------------------
C
      RETURN
      END
