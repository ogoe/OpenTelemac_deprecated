C                       *****************
                        SUBROUTINE DIRAUX
C                       *****************
C
     * ( X , Y , Z , W , F , INDIC , CRITER , MESH )
C
C***********************************************************************
C BIEF VERSION 5.1           06/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C  FONCTION : AIDE A LA PREPARATION D'UN SYSTEME LINEAIRE AVEC DES
C             CONDITIONS DE DIRICHLET.
C
C             X,Y ET Z DOIVENT ETRE DES STRUCTURES
C
C             ICI X EST UN VECTEUR DEFINI SUR LE DOMAINE
C                 Y EST UN VECTEUR DEFINI SUR LE DOMAINE
C                 Z EST UN VECTEUR DEFINI SUR LE DOMAINE OU DE BORD.
C
C             INDIC EST UN TABLEAU, PAS UNE STRUCTURE ||||||||||
C
C  |||||||| : ICI L'OPERATION N'EST EFFECTUEE QUE SI LA CONDITION
C             INDIC(K)=CRITER POUR UN NUMERO K GLOBAL OU DE BORD.
C
C  OPERATIONS EFFECTUEES :
C
C             W EST MIS A 0.D0 POUR LES POINTS OU INDIC(K) = CRITER
C                    ET A 1.D0 POUR LES AUTRES.
C
C             X EST MIS EGAL A Y MULTIPLIE PAR Z SI INDIC(K) = CRITER
C
C             F EST MIS EGAL Z SI INDIC(K) = CRITER
C
C  CES OPERATIONS SERVENT A TRAITER LES POINTS DE TYPE DIRICHLET
C
C             X EST LE SECOND MEMBRE QUI SERA EGAL A LA DIAGONALE
C             MULTIPLIEE PAR LA VALEUR DIRICHLET Z.
C
C             F EST L'INCONNUE MISE A SA VALEUR DIRICHLET
C
C             W EST UN TABLEAU DE TRAVAIL QUI VA SERVIR A ANNULER
C             LES TERMES DES MATRICES QUI TOUCHENT AUX POINTS DIRICHLET.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |      X         |<-- | VECTEUR RESULTAT
C |      Y         | -->| VECTEUR OPERANDE
C |      Z         | -->| VECTEUR OPERANDE
C |      W         | -->| TABLEAU DE TRAVAIL
C |      F         | -->| INCONNUE MISE A SA VALEUR DIRICHLET
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
      USE BIEF, EX_DIRAUX => DIRAUX
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ) , INTENT(INOUT) :: X,W,F
      TYPE(BIEF_OBJ) , INTENT(IN)    :: Y,Z
      INTEGER        , INTENT(IN)    :: INDIC(*),CRITER
      TYPE(BIEF_MESH), INTENT(IN)    :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K,NPOIN,IELMX,IELMZ,N
C
C-----------------------------------------------------------------------
C
      NPOIN = Z%DIM1
C
C-----------------------------------------------------------------------
C
C  W MIS A 1.
C
      CALL OS( 'X=C     ' , X=W , C=1.D0 )
C
C-----------------------------------------------------------------------
C
      IELMX=X%ELM
      IELMZ=Z%ELM
C
      IF(IELMX.NE.IELMZ) THEN
C
        DO K=1,NPOIN
          IF(INDIC(K).EQ.CRITER) THEN
            N = MESH%NBOR%I(K)
            X%R(N) = Y%R(N) * Z%R(K)
            W%R(N) = 0.D0
            F%R(N) = Z%R(K)
          ENDIF
        ENDDO
C
      ELSE
C
        DO K=1,NPOIN
          IF(INDIC(K).EQ.CRITER) THEN
            X%R(K) = Y%R(K) * Z%R(K)
            W%R(K) = 0.D0
            F%R(K) = Z%R(K)
          ENDIF
        ENDDO
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
