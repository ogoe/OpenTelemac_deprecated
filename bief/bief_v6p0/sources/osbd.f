C                       ***************
                        SUBROUTINE OSBD
C                       ***************
C
     * ( OP , X , Y , Z , C , MESH )
C
C***********************************************************************
C BIEF VERSION 5.5           06/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C  FONCTION : OPERATIONS SUR LES VECTEURS
C
C             X,Y ET Z DOIVENT ETRE DES STRUCTURES
C
C             ICI X EST UN VECTEUR DEFINI SUR LE BORD
C                 Y ET Z SONT DES VECTEURS DEFINIS SUR LE DOMAINE
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES VECTEURS X,Y ET Z ET LA CONSTANTE C. LE RESULTAT
C   EST LE VECTEUR X.
C
C   OP = 'X=Y     '     :  VALEURS DE BORD DE Y MISES DANS X
C   OP = 'X=+Y    '     :  IDEM
C   OP = 'X=X+Y   '     :  VALEURS DE BORD DE Y AJOUTEES A X
C   OP = 'X=X-Y   '     :  VALEURS DE BORD DE Y RETRANCHEES A X
C   OP = 'X=CY    '     :  VALEURS DE BORD DE CY MISES DANS X
C   OP = 'X=X+CY  '     :  VALEURS DE BORD DE CY AJOUTEES A X
C   OP = 'X=CXY   '     :
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |      OP        | -->| CHAINE DE CARACTERES INDIQUANT L'OPERATION
C |                |   >| A EFFECTUER.
C |      X         |<-- | VECTEUR RESULTAT
C |      Y         | -->| VECTEUR OPERANDE
C |      Z         | -->| VECTEUR OPERANDE
C |      C         | -->| CONSTANTE DONNEE
C |      NPOIN     | -->| DIMENSION DES VECTEURS
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELANTS :
C PROGRAMMES APPELES   : NEANT
C
C**********************************************************************
C
      USE BIEF, EX_OSBD => OSBD
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DOUBLE PRECISION, INTENT(IN) :: C
C
      CHARACTER(LEN=8), INTENT(IN) :: OP
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: X
      TYPE(BIEF_OBJ), INTENT(IN)    :: Y,Z
      TYPE(BIEF_MESH), INTENT(IN)   :: MESH
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER NPTFR,IELMX,IELMY
C
C-----------------------------------------------------------------------
C
      IF(X%TYPE.NE.2.OR.Y%TYPE.NE.2) THEN
        IF (LNG.EQ.1) WRITE(LU,100)
        IF (LNG.EQ.2) WRITE(LU,101)
100     FORMAT(1X,'OSBD (BIEF) : X ET Y NE SONT PAS DES VECTEURS')
101     FORMAT(1X,'OSBD (BIEF) : X AND Y ARE NOT VECTORS')
        CALL PLANTE(1)
        STOP
      ENDIF
C
      IELMX = X%ELM
      IELMY = Y%ELM
C
      IF((DIMENS(IELMX).NE.1.OR.DIMENS(IELMY).NE.2).AND.
     *   (DIMENS(IELMX).NE.2.OR.DIMENS(IELMY).NE.3)) THEN
        IF (LNG.EQ.1) WRITE(LU,102)
        IF (LNG.EQ.2) WRITE(LU,103)
102     FORMAT(1X,'OSBD (BIEF) : X ET Y ONT DE MAUVAISES DIMENSIONS')
103     FORMAT(1X,'OSBD (BIEF) : X AND Y HAVE WRONG DIMENSIONS')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      NPTFR = X%DIM1
C
      IF(  IELMY.EQ.11.OR.
     *     IELMY.EQ.61.OR.
     *     IELMY.EQ.21.OR.
     *     IELMY.EQ.71.OR.
     *     IELMY.EQ.12.OR.
     *     IELMY.EQ.31.OR.
     *     IELMY.EQ.41.OR. 
     *     IELMY.EQ.51      ) THEN
C       TABLEAU NBOR
        CALL OVBD( OP , X%R , Y%R , Z%R , C , MESH%NBOR%I , NPTFR )
!
!  JMH 23/06/2008 : MIS PROVISOIREMENT POUR VOIR QUI EN A BESOIN
!                   C'EST UNE ERREUR EN PARALLELISME CAR NELBOR
!                   PEUT ETRE NUL
!
!     ELSEIF(IELMY.EQ.10.OR.
!    *       IELMY.EQ.20.OR.
!    *       IELMY.EQ.30.OR.
!    *       IELMY.EQ.70.OR.
!    *       IELMY.EQ.40.OR.
!    *       IELMY.EQ.50     ) THEN
C       TABLEAU NELBOR
!       CALL OVBD( OP , X%R , Y%R , Z%R , C , MESH%NELBOR%I , NPTFR )
      ELSE
        IF (LNG.EQ.1) WRITE(LU,104)
        IF (LNG.EQ.2) WRITE(LU,105)
104     FORMAT(1X,'OSBD (BIEF) : DISCRETISATIONS NON PREVUES')
105     FORMAT(1X,'OSBD (BIEF) : DISCRETIZATIONS NOT IMPLEMENTED')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
