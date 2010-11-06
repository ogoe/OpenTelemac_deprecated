C                       ***************
                        SUBROUTINE OVBD
C                       ***************
C
     * ( OP , X , Y , Z , C , NBOR , NPTFR )
C
C***********************************************************************
C BIEF VERSION 5.1           06/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C  FONCTION : OPERATIONS SUR LES VECTEURS
C
C             ICI X EST UN VECTEUR DEFINI SUR LE BORD
C                 Y ET Z SONT DES VECTEURS DEFINIS SUR TOUT LE DOMAINE.
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
C |      NBOR      | -->| NUMEROTATION GLOBALE DES POINTS DE BORD.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELANTS :
C PROGRAMMES APPELES   : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NPTFR,NBOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: X(*)
      DOUBLE PRECISION, INTENT(IN)    :: Y(*),Z(*),C      
      CHARACTER(LEN=8), INTENT(IN)    :: OP
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K
C
C-----------------------------------------------------------------------
C
      IF(OP(1:8).EQ.'X=Y     '.OR.
     *   OP(1:8).EQ.'X=+Y    ') THEN
C
        DO K=1,NPTFR
          X(K) = Y(NBOR(K))
        ENDDO
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=X+Y   ') THEN
C
        DO K=1,NPTFR
          X(K) = X(K) + Y(NBOR(K))
        ENDDO
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=Y+Z   ') THEN
C
        DO K=1,NPTFR
          X(K) = Y(NBOR(K)) + Z(NBOR(K))
        ENDDO
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=X-Y   ') THEN
C
        DO K=1,NPTFR
          X(K) = X(K) - Y(NBOR(K))
        ENDDO
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=CY    ') THEN
C
        DO K=1,NPTFR
          X(K) = C * Y(NBOR(K))
        ENDDO
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=X+CY  ') THEN
C
        DO K=1,NPTFR
          X(K) = X(K) + C * Y(NBOR(K))
        ENDDO
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=CXY   ') THEN
C
        DO K=1,NPTFR
          X(K) = C * X(K) * Y(NBOR(K))
        ENDDO
C
C-----------------------------------------------------------------------
C
      ELSE
C
         IF (LNG.EQ.1) WRITE(LU,1000) OP
         IF (LNG.EQ.2) WRITE(LU,1001) OP
1000     FORMAT(1X,'OVBD (BIEF) : OPERATION INCONNUE: ',A8)
1001     FORMAT(1X,'OVBD (BIEF) : UNKNOWN OPERATION: ',A8)
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
