C                       ******************************
                        DOUBLE PRECISION FUNCTION DOTS
C                       ******************************
C
     *( X , Y )
C
C***********************************************************************
C BIEF VERSION 5.1           08/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C   FONCTION : PRODUIT SCALAIRE DE DEUX OBJETS QUI PEUVENT ETRE :
C
C              DEUX STRUCTURES DE VECTEURS
C
C              OU
C
C              DEUX STRUCTURES DE BLOCS DE VECTEURS EN NOMBRE ET
C              CARACTERISTIQUES IDENTIQUES
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
C |      X ET Y    | -->| LES DEUX STRUCTURES DONT ON VEUT LE PRODUIT
C |                |    | SCALAIRE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   :  PLANTE
C
C***********************************************************************
C
      USE BIEF, EX_DOTS => DOTS
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(IN) :: X,Y
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IBL
C
C-----------------------------------------------------------------------
C
C  CAS OU LES STRUCTURES SONT DES BLOCS
C
      IF(X%TYPE.EQ.4) THEN
C
       DOTS = 0.D0
       DO IBL = 1 , X%N
         DOTS=DOTS+DOT(X%ADR(IBL)%P%DIM1,X%ADR(IBL)%P%R,Y%ADR(IBL)%P%R)
       ENDDO
C
C-----------------------------------------------------------------------
C
C  CAS OU LES STRUCTURES NE SONT PAS DES BLOCS
C  (ON SUPPOSE ICI QUE Y EST DU MEME TYPE QUE X)
C
      ELSEIF(X%TYPE.EQ.2) THEN
C
        DOTS = DOT(X%DIM1,X%R,Y%R)
C
C-----------------------------------------------------------------------
C
C  ERREUR
C
      ELSE
C
         IF (LNG.EQ.1) WRITE(LU,50) X%NAME,X%TYPE
         IF (LNG.EQ.1) WRITE(LU,51) Y%NAME,Y%TYPE
         IF (LNG.EQ.1) WRITE(LU,53)
         IF (LNG.EQ.2) WRITE(LU,60) X%NAME,X%TYPE
         IF (LNG.EQ.2) WRITE(LU,61) Y%NAME,Y%TYPE
         IF (LNG.EQ.2) WRITE(LU,63)
50       FORMAT(1X,'DOTS (BIEF) : NOM DE X : ',A6,'  TYPE : ',1I6)
51       FORMAT(1X,'              NOM DE Y : ',A6,'  TYPE : ',1I6)
53       FORMAT(1X,'              CAS NON PREVU')
60       FORMAT(1X,'DOTS (BIEF) : NAME OF X : ',A6,'  TYPE : ',1I6)
61       FORMAT(1X,'              NAME OF Y : ',A6,'  TYPE : ',1I6)
63       FORMAT(1X,'              NOT IMPLEMENTED')
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
