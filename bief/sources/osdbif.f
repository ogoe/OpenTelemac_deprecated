C                       *****************
                        SUBROUTINE OSDBIF
C                       *****************
C
     * ( OP , X , Y , INDIC , CRITER , MESH )
C
C***********************************************************************
C BIEF VERSION 5.6        22/08/05    J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C  FONCTION : OPERATIONS CONDITIONNELLES SUR LES VECTEURS
C
C             X,Y ET Z DOIVENT ETRE DES STRUCTURES
C
C             ICI X EST UN VECTEUR DEFINI SUR LE DOMAINE
C                 Y EST UN VECTEUR DE BORD
C
C             INDIC EST UN TABLEAU : PAS UNE STRUCTURE ||||||||
C
C  |||||||| : ICI L'OPERATION N'EST EFFECTUEE QUE SI LA CONDITION
C             INDIC(K)=CRITERE POUR UN NUMERO K DE BORD.
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES VECTEURS X,Y ET Z ET LA CONSTANTE C. LE RESULTAT
C   EST LE VECTEUR X.
C
C   OP = 'X=Y     '     :  Y COPIE DANS X
C   OP = 'X=+Y    '     :  IDEM
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
C |      INDIC     | -->| TABLEAU D'INDICATEURS.
C |      CRITER    | -->| CRITERE POUR FAIRE L'OPERATION.
C |      MESH      | -->| STRUCTURE DU MAILLAGE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELANTS :
C PROGRAMMES APPELES   : NEANT
C
C**********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      TYPE(BIEF_OBJ) :: X,Y
      TYPE(BIEF_MESH) :: MESH
C
      INTEGER K,NPTFR,IELMX,IELMY
      INTEGER INDIC(*),CRITER
C
      CHARACTER*8 OP
C
C-----------------------------------------------------------------------
C
      IF(X%TYPE.NE.2.OR.Y%TYPE.NE.2) THEN
        IF (LNG.EQ.1) WRITE(LU,100)
        IF (LNG.EQ.2) WRITE(LU,101)
100     FORMAT(1X,'OSDBIF (BIEF) : X ET Y NE SONT PAS DES VECTEURS')
101     FORMAT(1X,'OSDBIF (BIEF) : X AND Y ARE NOT VECTORS')
        CALL PLANTE(1)
        STOP
      ENDIF
C
      IELMX = X%ELM
      IELMY = Y%ELM
!
! JP Renaud 18/08/2005
! Modification for 3D meshes: the domain vector dimension is 3 
! and the boundary vector dimension is 2. So the possible
! combinations are: 
!     -2D: dimesn(ielmx)==2 _and_ dimesn(ielmy)==1
!     -3D: dimesn(ielmx)==3 _and_ dimesn(ielmy)==2
!
      IF( .NOT. (DIMENS(IELMX).EQ.3 .AND.DIMENS(IELMY).EQ.2 )
     &    .AND.
     &    .NOT. (DIMENS(IELMX).EQ.2 .AND. DIMENS(IELMY).EQ.1 ) ) THEN
!
        IF (LNG.EQ.1) WRITE(LU,102)
        IF (LNG.EQ.2) WRITE(LU,103)
102     FORMAT(1X,'OSDBIF (BIEF) : X ET Y MAUVAISES DIMENSIONS')
103     FORMAT(1X,'OSDBIF (BIEF) : X AND Y WRONG DIMENSIONS')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      NPTFR = Y%DIM1
C
C-----------------------------------------------------------------------
C
      IF(OP(1:8).EQ.'X=Y     '.OR.
     *   OP(1:8).EQ.'X=+Y    ') THEN
C
        DO K=1,NPTFR
          IF(INDIC(K).EQ.CRITER) X%R(MESH%NBOR%I(K)) = Y%R(K)
        ENDDO
C
C-----------------------------------------------------------------------
C
      ELSE
C
         IF (LNG.EQ.1) WRITE(LU,1000) OP
         IF (LNG.EQ.2) WRITE(LU,1001) OP
1000     FORMAT(1X,'OSDBIF (BIEF) : OPERATION INCONNUE: ',A8)
1001     FORMAT(1X,'OSDBIF (BIEF) : UNKNOWN OPERATION: ',A8)
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
