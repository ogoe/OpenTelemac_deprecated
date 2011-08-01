C                       ***************
                        SUBROUTINE OSDB
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
C             ICI X EST UN VECTEUR DEFINI SUR LE DOMAINE
C                 Y ET Z SONT DES VECTEURS DE BORD
C
C   Y NE DOIT PAS ETRE UNE STRUCTURE FACTICE.
C   Z NON PROGRAMME POUR L'INSTANT.
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES VECTEURS X,Y ET Z ET LA CONSTANTE C. LE RESULTAT
C   EST LE VECTEUR X.
C
C   OP = 'X=Y     '     :  Y COPIE DANS X
C   OP = 'X=+Y    '     :  IDEM
C   OP = 'X=X+Y   '     :  Y AJOUTE A X
C   OP = 'X=X-Y   '     :  Y RETRANCHE DE X
C   OP = 'X=CY    '     :  CY MIS DANS X
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
      USE BIEF, EX_OSDB => OSDB
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DOUBLE PRECISION, INTENT(IN)  :: C
      CHARACTER(LEN=8), INTENT(IN)  :: OP
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
100     FORMAT(1X,'OSDB (BIEF) : X ET Y NE SONT PAS DES VECTEURS')
101     FORMAT(1X,'OSDB (BIEF) : X AND Y ARE NOT VECTORS')
        CALL PLANTE(1)
        STOP
      ENDIF
C
      IELMX = X%ELM
      IELMY = Y%ELM
C
      IF((DIMENS(IELMX).NE.2.OR.DIMENS(IELMY).NE.1).AND.
     *   (DIMENS(IELMX).NE.3.OR.DIMENS(IELMY).NE.2)) THEN
        IF (LNG.EQ.1) WRITE(LU,102)
        IF (LNG.EQ.2) WRITE(LU,103)
102     FORMAT(1X,'OSDB (BIEF) : X ET Y ONT DE MAUVAISES DIMENSIONS')
103     FORMAT(1X,'OSDB (BIEF) : X AND Y HAVE WRONG DIMENSIONS')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      NPTFR = Y%DIM1
C
CC MAILLAGE-3D
      IF(IELMX.EQ.11.OR.IELMX.EQ.21.OR.IELMX.EQ.31.OR.IELMX.EQ.61.OR.
     *   IELMX.EQ.12.OR.IELMX.EQ.41.OR.IELMX.EQ.51.OR.IELMX.EQ.81) THEN
C       TABLEAU NBOR
        CALL OVDB( OP , X%R , Y%R , Z%R , C , MESH%NBOR%I , NPTFR )
      ELSEIF(IELMX.EQ.10.OR.IELMX.EQ.20.OR.
     *       IELMX.EQ.40.OR.IELMX.EQ.50.OR.IELMX.EQ.80) THEN
C       TABLEAU NELBOR
        CALL OVDB( OP , X%R , Y%R , Z%R , C , MESH%NELBOR%I , NPTFR )
      ELSE
        IF (LNG.EQ.1) WRITE(LU,104)
        IF (LNG.EQ.2) WRITE(LU,105)
104     FORMAT(1X,'OSDB (BIEF) : DISCRETISATIONS NON PREVUES')
105     FORMAT(1X,'OSDB (BIEF) : DISCRETIZATIONS NOT IMPLEMENTED')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
