C                       ***********************
                        LOGICAL FUNCTION CMPOBJ
C                       ***********************
C
     *( OBJ1 , OBJ2 )
C
C***********************************************************************
C BIEF VERSION 5.1            01/03/90    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C  FONCTION  : COMPARAISON DE DEUX OBJETS
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   OBJ1,2       |<-->| LES DEUX STRUCTURES A COMPARER
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : AUCUN
C
C***********************************************************************
C
      USE BIEF, EX_CMPOBJ => CMPOBJ
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(IN) ::  OBJ1,OBJ2
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELM1,IELM2,TYP1,TYP2
C
C-----------------------------------------------------------------------
C
      CMPOBJ = .FALSE.
C
      TYP1 = OBJ1%TYPE
      TYP2 = OBJ2%TYPE
C
      IF(TYP1.EQ.TYP2) THEN
C
        IF(TYP1.EQ.2) THEN
C
C         VECTEURS : ON VERIFIE LA DISCRETISATION
C
          IELM1 = OBJ1%ELM
          IELM2 = OBJ2%ELM
          IF(IELM1.EQ.IELM2) CMPOBJ = .TRUE.
C
        ELSEIF(TYP1.EQ.4) THEN
C
C         BLOCS : ON VERIFIE LE NOMBRE D'OBJETS
C
          IF(OBJ1%N.EQ.OBJ2%N) CMPOBJ=.TRUE.
C
        ELSE
C
          IF(LNG.EQ.1) WRITE(LU,100)
          IF(LNG.EQ.2) WRITE(LU,101)
100       FORMAT(1X,'CMPOBJ (BIEF) : POUR BLOCS ET VECTEURS SEULEMENT')
101       FORMAT(1X,'CMPOBJ (BIEF) : FOR BLOCS AND VECTORS ONLY')
          CALL PLANTE(0)
          STOP
        ENDIF
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
