C                       *****************
                        SUBROUTINE PARMOY
C                       *****************
C
     *( X , MESH )
C
C***********************************************************************
C BIEF VERSION 5.1           24/04/97    J-M HERVOUET (LNH) 30 87 80 18
C COPYRIGHT 1997              AFTER REINHARD HINKELMANN (HANNOVER UNI.)
C***********************************************************************
C
C   FONCTION : MOYENNE D'UN VECTEUR AUX INTERFACES ENTRE
C              SOUS-DOMAINES.
C
C              X PEUT ETRE UN BLOC DE VECTEURS, EN CE CAS, TOUS LES
C              VECTEURS DU BLOC SONT TRAITES.
C
C              ATTENTION ||||
C
C              SI LES VECTEURS ONT UNE DEUXIEME DIMENSION
C              ELLE EST POUR L'INSTANT IGNOREE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C |      X         |<-->| VECTEUR OU BLOC DE VECTEURS.
C |      MESH      | -->| MAILLAGE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   :  PLANTE
C
C***********************************************************************
C
      USE BIEF, EX_PARMOY => PARMOY
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     STRUCTURES : MAILLAGE, VECTEURS OU BLOCS
C
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: X
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,TYPX
C
C  COMPLEMENT AUX INTERFACES :
C
      CALL PARCOM( X , 2 , MESH )
C
C-----------------------------------------------------------------------
C
      TYPX = X%TYPE
C
C-----------------------------------------------------------------------
C
C  CAS OU LES STRUCTURES SONT DES BLOCS
C
      IF(TYPX.EQ.4) THEN
C
        DO 10 I=1,X%N
          CALL OS('X=XY    ',X=X%ADR(I)%P,Y=MESH%FAC)
10      CONTINUE
C
C-----------------------------------------------------------------------
C
C  CAS OU LA STRUCTURE EST UN VECTEUR
C
      ELSEIF(TYPX.EQ.2) THEN
C
        CALL OS('X=XY    ',X=X,Y=MESH%FAC)
C
C-----------------------------------------------------------------------
C
C  ERREUR SUR LA STRUCTURE
C
      ELSE
C
         IF (LNG.EQ.1) WRITE(LU,50) X%NAME,X%TYPE
         IF (LNG.EQ.1) WRITE(LU,53)
50       FORMAT(1X,'PARMOY (BIEF) : NOM DE X : ',A6,'  TYPE : ',1I6)
53       FORMAT(1X,'                CAS NON PREVU')
         IF (LNG.EQ.2) WRITE(LU,51) X%NAME,X%TYPE
         IF (LNG.EQ.2) WRITE(LU,54)
51       FORMAT(1X,'PARMOY (BIEF) : NAME OF X: ',A6,'  TYPE : ',1I6)
54       FORMAT(1X,'                UNEXPECTED CASE')
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
 
