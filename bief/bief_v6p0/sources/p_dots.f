C                       ********************************
                        DOUBLE PRECISION FUNCTION P_DOTS
C                       ********************************
C
     *( X , Y , MESH )
C
C***********************************************************************
C BIEF VERSION 5.1           24/04/97    J-M HERVOUET (LNH) 30 87 80 18
C COPYRIGHT 1997              AFTER REINHARD HINKELMANN (HANNOVER UNI.)
C***********************************************************************
C
C   FONCTION : COMME DOTS MAIS EN TENANT COMPTE DU PARALLELISME.
C
C              PRODUIT SCALAIRE DE DEUX OBJETS QUI PEUVENT ETRE :
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
C |      MESH      | -->| MAILLAGE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   :  PLANTE
C
C***********************************************************************
C
      USE BIEF, EX_P_DOTS => P_DOTS
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      TYPE(BIEF_MESH), INTENT(IN) :: MESH
      TYPE(BIEF_OBJ), INTENT(IN)  :: X,Y
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER NPX,IBL,TYPX
C
      DOUBLE PRECISION P_DSUM
      EXTERNAL         P_DSUM
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
       P_DOTS = 0.D0
C
       IF(NCSIZE.LE.1.OR.NPTIR.EQ.0) THEN
         DO 99 IBL = 1 , X%N
           P_DOTS=P_DOTS+DOT(X%ADR(IBL)%P%DIM1,X%ADR(IBL)%P%R,
     *                                         Y%ADR(IBL)%P%R)
99       CONTINUE
       ELSE
         DO 100 IBL = 1 , X%N
           P_DOTS=P_DOTS+P_DOT(X%ADR(IBL)%P%DIM1,X%ADR(IBL)%P%R,
     *                                           Y%ADR(IBL)%P%R,
     *                                           MESH%FAC%R)
100      CONTINUE
       ENDIF
C
C-----------------------------------------------------------------------
C
C  CAS OU LES STRUCTURES NE SONT PAS DES BLOCS
C  (ON SUPPOSE ICI QUE Y EST DU MEME TYPE QUE X)
C
      ELSEIF(TYPX.EQ.2) THEN
C
        NPX = X%DIM1
C
        IF(Y%DIM1.NE.NPX) THEN
          IF (LNG.EQ.1) WRITE(LU,50) X%NAME,X%TYPE
          IF (LNG.EQ.1) WRITE(LU,51) Y%NAME,Y%TYPE
          IF (LNG.EQ.1) WRITE(LU,52) X%DIM1,Y%DIM1
          IF (LNG.EQ.2) WRITE(LU,60) X%NAME,X%TYPE
          IF (LNG.EQ.2) WRITE(LU,61) Y%NAME,Y%TYPE
          IF (LNG.EQ.2) WRITE(LU,62) X%DIM1,Y%DIM1
52        FORMAT(1X,'TAILLES DIFFERENTES : ',1I6,' ET ',1I6)
62        FORMAT(1X,'DIFFERENT SIZES: ',1I6,' AND ',1I6)
          CALL PLANTE(1)
          STOP
        ENDIF
C
        IF(NCSIZE.LE.1.OR.NPTIR.EQ.0) THEN
          P_DOTS=DOT(NPX,X%R,Y%R)
        ELSE
          P_DOTS=P_DOT(NPX,X%R,Y%R,MESH%FAC%R)
        ENDIF
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
50       FORMAT(1X,'P_DOTS (BIEF) : NOM DE X : ',A6,'  TYPE : ',1I6)
51       FORMAT(1X,'                NOM DE Y : ',A6,'  TYPE : ',1I6)
53       FORMAT(1X,'                CAS NON PREVU')
60       FORMAT(1X,'P_DOTS (BIEF) : NAME OF X : ',A6,'  TYPE : ',1I6)
61       FORMAT(1X,'                NAME OF Y : ',A6,'  TYPE : ',1I6)
63       FORMAT(1X,'                NOT IMPLEMENTED')
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C SOMME FINALE SUR TOUS LES SOUS-DOMAINES
C
      IF(NCSIZE.GT.1) P_DOTS = P_DSUM(P_DOTS)
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
