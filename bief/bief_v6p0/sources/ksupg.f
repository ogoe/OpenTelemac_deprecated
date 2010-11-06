C                       ****************
                        SUBROUTINE KSUPG
C                       ****************
C
     *(KX,KY,XMUL,U,V,MESH)
C
C***********************************************************************
C BIEF VERSION 5.6           08/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C  FONCTION : CALCUL D'UN VECTEUR QUI SERA UTILISE
C             PAR LA METHODE :
C
C             STREAMLINE UPWIND PETROV GALERKIN
C
C             AVEC UN DECENTREMENT EGAL A UN
C
C
C                    DX   U
C             KX = -----------
C                  2 NORME(U)
C
C                    DY   V
C             KY = -----------
C                  2 NORME(U)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      KX,KY     | -->|  COORDONNEES DU VECTEUR UNITAIRE.
C |      XMUL      | -->|  COEFFICIENT MULTIPLICATEUR.
C |      U,V       | -->|  COMPOSANTES DE LA VITESSE.
C |      MESH      | -->|  MAILLAGE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C  PROGRAMMES APPELES :
C
C**********************************************************************
C
      USE BIEF, EX_KSUPG => KSUPG
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      TYPE(BIEF_MESH), INTENT(IN)   :: MESH
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: KX,KY
      TYPE(BIEF_OBJ), INTENT(IN)    :: U,V
C
      DOUBLE PRECISION, INTENT(IN)  :: XMUL
C
C-----------------------------------------------------------------------
C
      IF(U%ELM.EQ.11.OR.U%ELM.EQ.12.OR.U%ELM.EQ.13) THEN
C
        CALL KSPG11(KX%R,KY%R,MESH%XEL%R,MESH%YEL%R,U%R,V%R,
     *              MESH%IKLE%I,MESH%NELEM,MESH%NELMAX,XMUL)
C
C  ELEMENT NON PREVU : ERREUR
C
      ELSE
        IF (LNG.EQ.1) WRITE(LU,100) U%ELM
        IF (LNG.EQ.2) WRITE(LU,101) U%ELM
100     FORMAT(1X,'KSUPG (BIEF) : U%ELM = ',1I6,' ELEMENT NON PREVU')
101     FORMAT(1X,'KSUPG (BIEF): U%ELM = ',1I6,' ELEMENT NOT AVAILABLE')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
