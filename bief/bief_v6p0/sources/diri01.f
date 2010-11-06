C                       *****************
                        SUBROUTINE DIRI01
C                       *****************
C
     *(F, S, SM ,FBOR,LIMDIR,WORK1,WORK2,MESH,KDIR,MSK,MASKPT)
C
C***********************************************************************
C BIEF VERSION 6.0      07/08/2009    J-M HERVOUET (LNHE) 01 30 87 80 18
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
C
C     BEWARE: THIS SUBROUTINE IS NOT PROTECTED AGAINST DIAGONAL EQUAL
C             TO ZERO ON DIRICHLET POINTS, IT WILL SET THEN A EQUATION
C             0 X = 0 ON SUCH POINTS.
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
C |   WORK1,2      | -->|  TABLEAUX DE TRAVAIL.
C |   KDIR         | -->|  CONVENTION POUR LES CONDITIONS DE DIRICHLET
C |   MSK          | -->| SI OUI, PRESENCE D'ELEMENTS MASQUES.         |
C |   MASKPT       | -->| TABLEAU DE MASQUAGE DES POINTS               |
C |                |    |  =1. : NORMAL   =0. : POINT MASQUE.          |
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C----------------------------------------------------------------------
C
C  PROGRAMMES APPELES : OV , MATVEC , MATMAT
C
C**********************************************************************
C
      USE BIEF, EX_DIRI01 => DIRI01
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: F,S,SM,WORK1,WORK2
      TYPE(BIEF_OBJ), INTENT(IN)    :: FBOR,MASKPT
      INTEGER, INTENT(IN) :: LIMDIR(*), KDIR
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      LOGICAL, INTENT(IN) :: MSK
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION C,Z(1)
C
      INTEGER IELMSM,IELMFB
C
      CHARACTER*1 OLDDIA
C
C----------------------------------------------------------------------
C
C  DEPLOIEMENT DE LA STRUCTURE DE MAILLAGE
C
C----------------------------------------------------------------------
C
C  CONSTRUCTION D'UN TABLEAU QUI VAUT ZERO PARTOUT SAUF AUX POINTS DE
C  TYPE DIRICHLET OU ON MET LA VALEUR SOUHAITEE PRISE DANS FBOR
C  FBOR DOIT ETRE NUL QUAND LE POINT N'EST PAS DIRICHLET.
C
      CALL CPSTVC(SM,WORK1)
C
      IELMSM=SM%ELM
      IELMFB=FBOR%ELM
      IF(IELMSM.EQ.IELMFB) THEN
        CALL MATVEC( 'X=AY    ' ,WORK2,S,FBOR,C, MESH )
      ELSE
        CALL OS( 'X=0     ' , X=WORK1 )
        CALL OSDBIF( 'X=Y     ' ,WORK1,FBOR,LIMDIR,KDIR,MESH)
        CALL MATVEC( 'X=AY    ' ,WORK2,S,WORK1,C, MESH )
      ENDIF
C
C----------------------------------------------------------------------
C
C  LE PRODUIT S WORK1 EST RETRANCHE AU SECOND MEMBRE
C  CECI SIGNIFIE QUE LES VALEURS DES POINTS DIRICHLET NE SONT PLUS
C  CONSIDEREES COMME DES INCONNUES DANS LES EQUATIONS DES AUTRES POINTS.
C
      CALL OS( 'X=X-Y   ' , X=SM , Y=WORK2 )
C
C----------------------------------------------------------------------
C
C  CONSTRUCTION D'UN TABLEAU QUI VAUT 1. PARTOUT SAUF POUR LES POINTS DE
C  TYPE DIRICHLET OU IL VAUT 0.
C
C  DE PLUS POUR LES POINTS DIRICHLET ON MET DANS LA MATRICE UNE EQUATION
C  DE LA FORME DS(N) * X = DS(N) * FBOR QUI DONNERA X=FBOR
C  ET ON INITIALISE F A SA VALEUR CONNUE.
C  THIS ASSUMES THAT DS(N) IS NOT ZERO
C
      CALL DIRAUX(SM,S%D,FBOR,WORK2,F,LIMDIR,KDIR,MESH )
C
C  MASQUAGE : POUR LES POINTS DES ELEMENTS MASQUES ON MET
C             UNE EQUATION X=0. AU COEFFICIENT DIAGONAL PRES
C
      IF(MSK) THEN
        CALL OV( 'X=XY    ', SM%R   ,MASKPT%R ,Z,C,   SM%DIM1)
        CALL OV( 'X=XY    ', F%R    ,MASKPT%R ,Z,C,    F%DIM1)
        CALL OV( 'X=XY    ', WORK2%R,MASKPT%R ,Z,C,WORK2%DIM1)
      ENDIF
C
C----------------------------------------------------------------------
C
C  MULTIPLICATION DE S A DROITE ET A GAUCHE PAR WORK2 : L'EFFET EST
C  D'EFFACER LES LIGNES ET LES COLONNES DE S QUI CORRESPONDENT A DES
C  POINTS DIRICHLET. COMME ON NE VEUT PAS TOUCHER A LA DIAGONALE ON LA
C  DECLARE ICI COMME NULLE.
C
      OLDDIA=S%TYPDIA
      S%TYPDIA='0'
      CALL OM( 'M=DMD   ' , S , S , WORK2 , C , MESH )
      S%TYPDIA=OLDDIA
C
C----------------------------------------------------------------------
C
      RETURN
      END 
