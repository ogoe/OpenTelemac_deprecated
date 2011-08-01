C                       *****************
                        SUBROUTINE DIRI09
C                       *****************
C
     *(X1,X2,X3,
     * A11,A12,A13,A21,A22,A23,A31,A32,A33,
     * SM1,SM2,SM3,T1,T2,T3,T4,T5,T6,
     * XBOR1,XBOR2,XBOR3,LIDIR1,LIDIR2,LIDIR3,
     * MESH,KDIR,MSK,MASKPT)
C
C***********************************************************************
C BIEF VERSION 5.1           30/01/95    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C   FONCTION: TRAITEMENT DES POINTS DE TYPE DIRICHLET POUR UN
C             SYSTEME COMPOSE COMME SUIT (BLOC DE 9 MATRICES)
C
C         (     A11          A12         A13  )  ( X1 )   ( SM1 )
C         (                                   )  (    )   (     )
C         (    T                              )  (    )   (     )
C         (     A21          A22         A23  )  ( X2 ) = ( SM2 )
C         (                                   )  (    )   (     )
C         (    T            T                 )  (    )   (     )
C         (     A31          A32         A33  )  ( X3 )   ( SM3 )
C
C
C   ATTENTION A LA TRANSPOSITION DE A21 A31 ET A32 QUI PERMET DE
C   CONFONDRE A L'APPEL A12 ET A21 , A31 ET A13 , A32 ET A23 LORSQUE
C   LE BLOC EST SYMETRIQUE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   X1 ,X2 ,X3   |<-- |  VALEURS A L'ETAPE N+1 (INITIALISATION)
C |   A11,A12 ETC  |<-->|  MATRICES
C |   SM1,SM2,SM3  | -->|  SECONDS MEMBRES DU SYSTEME.
C |   S1,2,3       | -->|  SIGNES DEVANT LES SECONDS MEMBRES.
C |   T1,..6       | -->|  TABLEAUX DE TRAVAIL DU SYSTEME
C |   XBOR1,2,3.   | -->|  CONDITIONS AUX LIMITES SUR X1,2,3.
C |   LIDIR1,2,3   | -->|  TYPES DE CONDITIONS AUX LIMITES POUR X1,2,3
C |   KDIR         | -->|  CONDITION A LA LIMITE DE TYPE DIRICHLET
C |   MSK          | -->| SI OUI, PRESENCE D'ELEMENTS MASQUES.         |
C |   MASKPT       | -->| TABLEAU DE MASQUAGE DES POINTS               |
C |                |    |  =1. : NORMAL   =0. : POINT MASQUE.          |
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C**********************************************************************
C
C SOUS PROGRAMMES APPELES :
C
C**********************************************************************
C
      USE BIEF, EX_DIRI09 => DIRI09
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: X1,X2,X3,SM1,SM2,SM3
      TYPE(BIEF_OBJ), INTENT(INOUT) :: T1,T2,T3,T4,T5,T6
      TYPE(BIEF_OBJ), INTENT(INOUT) :: A11,A12,A13,A21,A22
      TYPE(BIEF_OBJ), INTENT(INOUT) :: A23,A31,A32,A33
      TYPE(BIEF_OBJ), INTENT(IN)    :: XBOR1,XBOR2,XBOR3,MASKPT
      INTEGER, INTENT(IN)           :: LIDIR1(*),LIDIR2(*),LIDIR3(*)
      INTEGER, INTENT(IN)           :: KDIR
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
C 1) ON CONSTRUIT DES TABLEAUX T1,T2,T3 QUI CONTIENNENT
C    LES VALEURS DE X1,2,3 IMPOSEES SI LE POINT EST DE TYPE DIRICHLET
C    ET 0 SINON.
C
C    PAR AILLEURS X1,X2,X3 RECOIVENT LEUR VALEUR DIRICHLET.
C
C=======================================================================
C
C   CONDITION AUX LIMITES SUR X1 : "XBOR1" IMPOSE
C
      CALL CPSTVC(X1,T1)
      CALL OS  ( 'X=C     ' , T1 , T1 , T1 , 0.D0 )
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
C-----------------------------------------------------------------------
C
C   CONDITIONS AUX LIMITES SUR X3 : "XBOR3" IMPOSE
C
      CALL CPSTVC(X3,T3)
      CALL OS  ( 'X=C     ' , T3 , T3 , T3 , 0.D0 )
      CALL OSDBIF ( 'X=Y     ',T3,XBOR3,LIDIR3,KDIR,MESH)
C
C=======================================================================
C
C   2) ON CALCULE LE PRODUIT DE LA MATRICE DU SYSTEME A RESOUDRE PAR
C      T1,T2,T3 ET LE RESULTAT EST RETRANCHE DES SECONDS
C      MEMBRES.
C
      CALL MATVEC('X=AY    ',T4,A11,T1,C,MESH,LEGO=.FALSE.)
      CALL MATVEC('X=X+AY  ',T4,A12,T2,C,MESH,LEGO=.FALSE.)
      CALL MATVEC('X=X+AY  ',T4,A13,T3,C,MESH,LEGO=.TRUE. )
      CALL MATVEC('X=AY    ',T5,A21,T1,C,MESH,LEGO=.FALSE.)
      CALL MATVEC('X=X+AY  ',T5,A22,T2,C,MESH,LEGO=.FALSE.)
      CALL MATVEC('X=X+AY  ',T5,A23,T3,C,MESH,LEGO=.TRUE. )
      CALL MATVEC('X=AY    ',T6,A31,T1,C,MESH,LEGO=.FALSE.)
      CALL MATVEC('X=X+AY  ',T6,A32,T2,C,MESH,LEGO=.FALSE.)
      CALL MATVEC('X=X+AY  ',T6,A33,T3,C,MESH,LEGO=.TRUE. )
C
      CALL CPSTVC(X1,SM1)
      CALL CPSTVC(X2,SM2)
      CALL CPSTVC(X3,SM3)
      CALL OS( 'X=X-Y   ' , SM1 , T4 , T4 , C )
      CALL OS( 'X=X-Y   ' , SM2 , T5 , T5 , C )
      CALL OS( 'X=X-Y   ' , SM3 , T6 , T6 , C )
C
C=======================================================================
C
C  SECONDS MEMBRES DES EQUATIONS PORTANT SUR DES POINTS DIRICHLETS
C  ET PREPARATION DU SYSTEME LINEAIRE.
C
      CALL DIRAUX(SM1,A11%D,XBOR1,T1,X1,LIDIR1,KDIR,MESH)
      CALL DIRAUX(SM2,A22%D,XBOR2,T2,X2,LIDIR2,KDIR,MESH)
      CALL DIRAUX(SM3,A33%D,XBOR3,T3,X3,LIDIR3,KDIR,MESH)
C
C ON APPELLE OV PLUTOT QUE OS CAR SM1 ET MASKPT N'ONT PAS TOUJOURS
C LA MEME LONGUEUR.
C
      IF(MSK) THEN
        CALL OV( 'X=XY    ',SM1%R,MASKPT%R,Z,C,SM1%DIM1)
        CALL OV( 'X=XY    ', X1%R,MASKPT%R,Z,C,X1%DIM1)
        CALL OV( 'X=XY    ', T1%R,MASKPT%R,Z,C,T1%DIM1)
        CALL OV( 'X=XY    ',SM2%R,MASKPT%R,Z,C,SM2%DIM1)
        CALL OV( 'X=XY    ', X2%R,MASKPT%R,Z,C,X2%DIM1)
        CALL OV( 'X=XY    ', T2%R,MASKPT%R,Z,C,T2%DIM1)
        CALL OV( 'X=XY    ',SM3%R,MASKPT%R,Z,C,SM3%DIM1)
        CALL OV( 'X=XY    ', X3%R,MASKPT%R,Z,C,X3%DIM1)
        CALL OV( 'X=XY    ', T3%R,MASKPT%R,Z,C,T3%DIM1)
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
C PRECONDITIONNEMENT DE A13 :
C=======================================================================
C
      CALL OM( 'M=DM    ' , A13,A13 ,T1,C,MESH)
      CALL OM( 'M=MD    ' , A13,A13 ,T3,C,MESH)
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
C=======================================================================
C PRECONDITIONNEMENT DE A23 :
C=======================================================================
C
      CALL OM( 'M=DM    ' , A23,A23 ,T2,C,MESH)
      CALL OM( 'M=MD    ' , A23,A23 ,T3,C,MESH)
C
C=======================================================================
C PRECONDITIONNEMENT DE A31 :
C=======================================================================
C
      CALL OM( 'M=DM    ' , A31,A31 ,T3,C,MESH)
      CALL OM( 'M=MD    ' , A31,A31 ,T1,C,MESH)
C
C=======================================================================
C PRECONDITIONNEMENT DE A32 :
C=======================================================================
C
      CALL OM( 'M=DM    ' , A32,A32 ,T3,C,MESH)
      CALL OM( 'M=MD    ' , A32,A32 ,T2,C,MESH)
C
C=======================================================================
C PRECONDITIONNEMENT DE A33 :
C=======================================================================
C
      STODIA = A33%TYPDIA
      A33%TYPDIA='0'
      CALL OM( 'M=DMD   ' , A33,A33 ,T3,C,MESH)
      A33%TYPDIA=STODIA
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
 
