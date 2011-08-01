C                       *****************
                        SUBROUTINE DIRI04
C                       *****************
C
     *(X1,X2,
     * A11,A12,A21,A22,
     * SM1,SM2,T1,T2,T3,T4,
     * XBOR1,XBOR2,LIDIR1,LIDIR2,
     * MESH,KDIR,MSK,MASKPT)
C
C***********************************************************************
C BIEF VERSION 5.1           30/01/95    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C   FONCTION: TRAITEMENT DES POINTS DE TYPE DIRICHLET POUR UN
C             SYSTEME COMPOSE COMME SUIT (BLOC DE 4 MATRICES)
C
C         (     A11          A12              )  ( X1 )   ( SM1 )
C         (                                   )  (    ) = (     )
C         (     A21          A22              )  ( X2 )   ( SM2 )
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   X1 ,X2       |<-- |  VALEURS A L'ETAPE N+1 (INITIALISATION)
C |   A11,A12 ETC  | -->|  MATRICES
C |   SM1,SM2      | -->|  SECONDS MEMBRES DU SYSTEME.
C |   T1,..4       | -->|  TABLEAUX DE TRAVAIL DU SYSTEME
C |   XBOR1,2      | -->|  CONDITIONS AUX LIMITES SUR X1,2,3.
C |   LIDIR1,2     | -->|  TYPES DE CONDITIONS AUX LIMITES POUR X1,2,3
C |   MESH         | -->|  BLOC DES TABLEAUX D'ENTIERS DU MAILLAGE.
C |   KDIR         | -->|  CONDITION A LA LIMITE DE TYPE DIRICHLET
C |   MSK          | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |   MASKPT       | -->|  TABLEAU DE MASQUAGE DES POINTS
C |                |    |  =1. : NORMAL   =0. : POINT MASQUE.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C**********************************************************************
C
C SOUS PROGRAMMES APPELES :
C
C**********************************************************************
C
      USE BIEF, EX_DIRI04 => DIRI04
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: X1,X2,SM1,SM2,T1,T2,T3,T4
      TYPE(BIEF_OBJ), INTENT(INOUT) :: A11,A12,A21,A22
      TYPE(BIEF_OBJ), INTENT(IN)    :: XBOR1,XBOR2,MASKPT
      INTEGER, INTENT(IN)           :: KDIR,LIDIR1(*),LIDIR2(*)
      TYPE(BIEF_MESH), INTENT(INOUT):: MESH
      LOGICAL, INTENT(IN)           :: MSK
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION C,Z(1)
C
      CHARACTER*1 STODIA
C
C-----------------------------------------------------------------------
C
C 1) ON CONSTRUIT DES TABLEAUX T1,T2    QUI CONTIENNENT
C    LES VALEURS DE X1 ET X2 IMPOSEES SI LE POINT EST DE TYPE DIRICHLET
C    ET 0 SINON.
C
C    PAR AILLEURS X1,X2  3 RECOIVENT LEUR VALEUR DIRICHLET.
C
C=======================================================================
C
C   CONDITION AUX LIMITES SUR X1 : "XBOR1" IMPOSE
C
      CALL CPSTVC(X1,T1)
      CALL OS ( 'X=C     ' , T1 , T1 , T1 , 0.D0 )
      CALL OSDBIF ( 'X=Y     ',T1,XBOR1,LIDIR1,KDIR,MESH)
C
C-----------------------------------------------------------------------
C
C   CONDITION AUX LIMITES SUR X2 : "XBOR2" IMPOSE
C
      CALL CPSTVC(X2,T2)
      CALL OS  ( 'X=C     ' , T2 , T2 , T2 , 0.D0 )
      CALL OSDBIF ( 'X=Y     ',T2,XBOR2,LIDIR2,KDIR,MESH)
C
C=======================================================================
C
C   2) ON CALCULE LE PRODUIT DE LA MATRICE DU SYSTEME A RESOUDRE PAR
C      T1,T2 ET LE RESULTAT EST RETRANCHE DES SECONDS
C      MEMBRES.
C
      CALL MATVEC('X=AY    ',T3,A11,T1,C,MESH,LEGO=.FALSE.)
      CALL MATVEC('X=X+AY  ',T3,A12,T2,C,MESH,LEGO=.TRUE. )
      CALL MATVEC('X=AY    ',T4,A21,T1,C,MESH,LEGO=.FALSE.)
      CALL MATVEC('X=X+AY  ',T4,A22,T2,C,MESH,LEGO=.TRUE. )
C
      CALL CPSTVC(X1,SM1)
      CALL CPSTVC(X2,SM2)
      CALL OS( 'X=X-Y   ' , SM1 , T3 , T3 , C )
      CALL OS( 'X=X-Y   ' , SM2 , T4 , T4 , C )
C
C=======================================================================
C
C  SECONDS MEMBRES DES EQUATIONS PORTANT SUR DES POINTS DIRICHLETS
C  ET PREPARATION DU SYSTEME LINEAIRE.
C
      CALL DIRAUX(SM1,A11%D,XBOR1,T1,X1,LIDIR1,KDIR,MESH )
      CALL DIRAUX(SM2,A22%D,XBOR2,T2,X2,LIDIR2,KDIR,MESH )
C
      IF(MSK) THEN
        CALL OV( 'X=XY    ',SM1%R,MASKPT%R,Z,C,SM1%DIM1)
        CALL OV( 'X=XY    ', X1%R,MASKPT%R,Z,C,X1%DIM1)
        CALL OV( 'X=XY    ', T1%R,MASKPT%R,Z,C,T1%DIM1)
        CALL OV( 'X=XY    ',SM2%R,MASKPT%R,Z,C,SM2%DIM1)
        CALL OV( 'X=XY    ', X2%R,MASKPT%R,Z,C,X2%DIM1)
        CALL OV( 'X=XY    ', T2%R,MASKPT%R,Z,C,T2%DIM1)
      ENDIF
C
C=======================================================================
C
C   EFFACEMENT DES LIGNES ET DES COLONNES DES POINTS DIRICHLETS
C
C   C'EST L'EQUIVALENT D'UN PRECONDITIONNEMENT DIAGONAL
C   AVEC LES TABLEAUX T1,T2,T3
C
C   MAIS ON NE VEUT PAS TOUCHER AUX DIAGONALES DE A11,A22,A33
C   ON LEUR MET ALORS UN FAUX TYPE : '0' CAR DANS CE CAS
C   LE SOUS-PROGRAMME OM N'Y TOUCHE PAS.
C
C=======================================================================
C PRECONDITIONNEMENT DE A11 :
C=======================================================================
C
      STODIA = A11%TYPDIA
      A11%TYPDIA='0'
      CALL OM( 'M=DMD   ' , A11,A11 ,T1,C,MESH)
      A11%TYPDIA=STODIA
C
C=======================================================================
C PRECONDITIONNEMENT DE A12 :
C=======================================================================
C
      CALL OM( 'M=DM    ' , A12,A12 ,T1,C,MESH)
      CALL OM( 'M=MD    ' , A12,A12 ,T2,C,MESH)
C
C=======================================================================
C PRECONDITIONNEMENT DE A21 :
C=======================================================================
C
      CALL OM( 'M=DM    ' , A21,A21 ,T2,C,MESH)
      CALL OM( 'M=MD    ' , A21,A21 ,T1,C,MESH)
C
C=======================================================================
C PRECONDITIONNEMENT DE A22 :
C=======================================================================
C
      STODIA = A22%TYPDIA
      A22%TYPDIA='0'
      CALL OM( 'M=DMD   ' , A22,A22 ,T2,C,MESH)
      A22%TYPDIA=STODIA
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
 
