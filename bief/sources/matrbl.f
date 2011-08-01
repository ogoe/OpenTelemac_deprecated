C                       *****************
                        SUBROUTINE MATRBL
C                       *****************
C
     *( OP , X , A , Y , C , MESH )
C
C***********************************************************************
C BIEF VERSION 5.1           24/04/97    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C   FONCTION : OPERATIONS MATRICE VECTEUR
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES VECTEURS X,Y ET LA MATRICE M. LE RESULTAT
C   EST LE VECTEUR X (DONT UNE PARTIE NON ASSEMBLEE PEUT ETRE DANS
C   LE TABLEAU W SI LEGO = .FALSE.)
C
C   CES OPERATIONS SONT DIFFERENTES SUIVANT LE TYPE DE DIAGONALE
C   ET LE TYPE DES TERMES EXTRADIAGONAUX.
C
C   OPERATIONS PROGRAMMEES :
C
C      OP = 'X=AY    '  : X = AY
C      OP = 'X=X+AY  '  : X = X + AY
C      OP = 'X=X-AY  '  : X = X - AY
C      OP = 'X=X+CAY '  : X = X + C AY
C      OP = 'X=TAY   '  : X = TA Y (TRANSPOSEE DE A)
C      OP = 'X=X+TAY '  : X = X + TA Y
C      OP = 'X=X-TAY '  : X = X - TA Y
C      OP = 'X=X+CTAY'  : X = X + C TA Y
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C |      OP        | -->|OPERATION A EFFECTUER
C |      X         |<-- |VECTEUR IMAGE
C |      DA        | -->| DIAGONALE DE LA MATRICE.
C |      TYPDIA    | -->| TYPE DE LA DIAGONALE (CHAINE DE CARACTERES)
C |                |    | TYPDIA = 'Q' : DIAGONALE QUELCONQUE
C |                |    | TYPDIA = 'I' : DIAGONALE IDENTITE.
C |                |    | TYPDIA = '0' : DIAGONALE NULLE.
C |      XA        | -->| TERMES EXTRA-DIAGONAUX DE LA MATRICE
C |      TYPEXT    | -->| TYPE DES TERMES EXTRADIAGONAUX
C |                |    | TYPEXT = 'Q' : QUELCONQUES.
C |                |    | TYPEXT = 'S' : SYMETRIQUES.
C |                |    | TYPEXT = '0' : NULS.
C |      Y         | -->| VECTEUR OPERANDE
C |      C         | -->| CONSTANTE DONNEE
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   : MATVEC
C
C***********************************************************************
C
      USE BIEF, EX_MATRBL => MATRBL
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER(LEN=8), INTENT(IN)   :: OP
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: X
      TYPE(BIEF_OBJ), INTENT(IN)     :: A,Y
      DOUBLE PRECISION, INTENT(IN)   :: C
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER S   
C
C-----------------------------------------------------------------------
C
C     CAS OU LES STRUCTURES SONT DES BLOCS
C
      IF(A%TYPE.EQ.4) THEN
C
        S = X%N
C
        IF(S.EQ.1) THEN
C
         CALL MATVEC( OP,X%ADR(1)%P,A%ADR(1)%P,Y%ADR(1)%P,C,MESH)
C
        ELSEIF(S.EQ.2) THEN
C
          IF(OP(1:8).EQ.'X=AY    ') THEN
            CALL MATVEC('X=AY    ',
     *      X%ADR(1)%P,A%ADR(1)%P,Y%ADR(1)%P,C,MESH,LEGO=.FALSE.)
            CALL MATVEC('X=X+AY  ',
     *      X%ADR(1)%P,A%ADR(2)%P,Y%ADR(2)%P,C,MESH,LEGO=.TRUE.)
            CALL MATVEC('X=AY    ',
     *      X%ADR(2)%P,A%ADR(3)%P,Y%ADR(1)%P,C,MESH,LEGO=.FALSE.)
            CALL MATVEC('X=X+AY  ',
     *      X%ADR(2)%P,A%ADR(4)%P,Y%ADR(2)%P,C,MESH,LEGO=.TRUE.)
          ELSEIF(OP(1:8).EQ.'X=TAY   ') THEN
            CALL MATVEC('X=TAY   ',
     *      X%ADR(1)%P,A%ADR(1)%P,Y%ADR(1)%P,C,MESH,LEGO=.FALSE.)
            CALL MATVEC('X=X+TAY ',
     *      X%ADR(1)%P,A%ADR(3)%P,Y%ADR(2)%P,C,MESH,LEGO=.TRUE.)
            CALL MATVEC('X=TAY   ',
     *      X%ADR(2)%P,A%ADR(2)%P,Y%ADR(1)%P,C,MESH,LEGO=.FALSE.)
            CALL MATVEC('X=X+TAY ',
     *      X%ADR(2)%P,A%ADR(4)%P,Y%ADR(2)%P,C,MESH,LEGO=.TRUE.)
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,10) OP
            IF (LNG.EQ.2) WRITE(LU,11) OP
            CALL PLANTE(1)
            STOP
          ENDIF
C
        ELSEIF(S.EQ.3) THEN
C
          IF(OP(1:8).EQ.'X=AY    ') THEN
            CALL MATVEC('X=AY    ',
     *      X%ADR(1)%P,A%ADR(1)%P,Y%ADR(1)%P,C,MESH,LEGO=.FALSE.)
            CALL MATVEC('X=X+AY  ',
     *      X%ADR(1)%P,A%ADR(2)%P,Y%ADR(2)%P,C,MESH,LEGO=.FALSE.)
            CALL MATVEC('X=X+AY  ',
     *      X%ADR(1)%P,A%ADR(3)%P,Y%ADR(3)%P,C,MESH,LEGO=.TRUE.)
            CALL MATVEC('X=AY    ',
     *      X%ADR(2)%P,A%ADR(4)%P,Y%ADR(1)%P,C,MESH,LEGO=.FALSE.)
            CALL MATVEC('X=X+AY  ',
     *      X%ADR(2)%P,A%ADR(5)%P,Y%ADR(2)%P,C,MESH,LEGO=.FALSE.)
            CALL MATVEC('X=X+AY  ',
     *      X%ADR(2)%P,A%ADR(6)%P,Y%ADR(3)%P,C,MESH,LEGO=.TRUE. )
            CALL MATVEC('X=AY    ',
     *      X%ADR(3)%P,A%ADR(7)%P,Y%ADR(1)%P,C,MESH,LEGO=.FALSE.)
            CALL MATVEC('X=X+AY  ',
     *      X%ADR(3)%P,A%ADR(8)%P,Y%ADR(2)%P,C,MESH,LEGO=.FALSE.)
            CALL MATVEC('X=X+AY  ',
     *      X%ADR(3)%P,A%ADR(9)%P,Y%ADR(3)%P,C,MESH,LEGO=.TRUE.)                                      
          ELSEIF(OP(1:8).EQ.'X=TAY   ') THEN
            CALL MATVEC('X=TAY   ',
     *      X%ADR(1)%P,A%ADR(1)%P,Y%ADR(1)%P,C,MESH,LEGO=.FALSE.)
            CALL MATVEC('X=X+TAY ',
     *      X%ADR(1)%P,A%ADR(4)%P,Y%ADR(2)%P,C,MESH,LEGO=.FALSE.)
            CALL MATVEC('X=X+TAY ',
     *      X%ADR(1)%P,A%ADR(7)%P,Y%ADR(3)%P,C,MESH,LEGO=.TRUE.)
            CALL MATVEC('X=TAY   ',
     *      X%ADR(2)%P,A%ADR(2)%P,Y%ADR(1)%P,C,MESH,LEGO=.FALSE.)
            CALL MATVEC('X=X+TAY ',
     *      X%ADR(2)%P,A%ADR(5)%P,Y%ADR(2)%P,C,MESH,LEGO=.FALSE.)
            CALL MATVEC('X=X+TAY ',
     *      X%ADR(2)%P,A%ADR(8)%P,Y%ADR(3)%P,C,MESH,LEGO=.TRUE.)
            CALL MATVEC('X=TAY   ',
     *      X%ADR(3)%P,A%ADR(3)%P,Y%ADR(1)%P,C,MESH,LEGO=.FALSE.)
            CALL MATVEC('X=X+TAY ',
     *      X%ADR(3)%P,A%ADR(6)%P,Y%ADR(2)%P,C,MESH,LEGO=.FALSE.)
            CALL MATVEC('X=X+TAY ',
     *      X%ADR(3)%P,A%ADR(9)%P,Y%ADR(3)%P,C,MESH,LEGO=.TRUE.)
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,10) OP
            IF (LNG.EQ.2) WRITE(LU,11) OP
10          FORMAT(1X,'MATRBL (BIEF) : OPERATION INCONNUE : ',A8)
11          FORMAT(1X,'MATRBL (BIEF) : UNKNOWN OPERATION  : ',A8)
            CALL PLANTE(0)
            STOP
          ENDIF
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,150) S
          IF (LNG.EQ.2) WRITE(LU,151) S
          IF (LNG.EQ.1) WRITE(LU,50) X%NAME,X%TYPE
          IF (LNG.EQ.1) WRITE(LU,51) Y%NAME,Y%TYPE
          IF (LNG.EQ.1) WRITE(LU,52) A%NAME,A%TYPE
          IF (LNG.EQ.1) WRITE(LU,53)
          IF (LNG.EQ.2) WRITE(LU,60) X%NAME,X%TYPE
          IF (LNG.EQ.2) WRITE(LU,61) Y%NAME,Y%TYPE
          IF (LNG.EQ.2) WRITE(LU,62) A%NAME,A%TYPE
          IF (LNG.EQ.2) WRITE(LU,63)
150       FORMAT(1X,'MATRBL (BIEF) : TROP DE VECTEURS INCONNUS :',1I6)
151       FORMAT(1X,'MATRBL (BIEF) : TOO MANY VECTORS          :',1I6)
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
C-----------------------------------------------------------------------
C
C  CAS OU LES STRUCTURES NE SONT PAS DES BLOCS
C
      ELSEIF(A%TYPE.EQ.3.AND.X%TYPE.EQ.4.AND.Y%TYPE.EQ.4) THEN
C
        CALL MATVEC( OP , X%ADR(1)%P , A , Y%ADR(1)%P , C , MESH )
C
C-----------------------------------------------------------------------
C
      ELSEIF(A%TYPE.EQ.3.AND.X%TYPE.EQ.2.AND.Y%TYPE.EQ.2) THEN
C
        CALL MATVEC( OP , X          , A , Y          , C , MESH )
C
C-----------------------------------------------------------------------
C
C  ERREUR
C
      ELSE
C
         IF (LNG.EQ.1) WRITE(LU,50) X%NAME,X%TYPE
         IF (LNG.EQ.1) WRITE(LU,51) Y%NAME,Y%TYPE
         IF (LNG.EQ.1) WRITE(LU,52) A%NAME,A%TYPE
         IF (LNG.EQ.1) WRITE(LU,53)
         IF (LNG.EQ.2) WRITE(LU,60) X%NAME,X%TYPE
         IF (LNG.EQ.2) WRITE(LU,61) Y%NAME,Y%TYPE
         IF (LNG.EQ.2) WRITE(LU,62) A%NAME,A%TYPE
         IF (LNG.EQ.2) WRITE(LU,63)
50       FORMAT(1X,'MATRBL (BIEF) : NOM DE X : ',A6,'  TYPE : ',1I6)
51       FORMAT(1X,'                NOM DE Y : ',A6,'  TYPE : ',1I6)
52       FORMAT(1X,'                NOM DE A : ',A6,'  TYPE : ',1I6)
53       FORMAT(1X,'                CAS NON PREVU')
60       FORMAT(1X,'MATRBL (BIEF) : NAME OF X : ',A6,'  TYPE : ',1I6)
61       FORMAT(1X,'                NAME OF Y : ',A6,'  TYPE : ',1I6)
62       FORMAT(1X,'                NAME OF A : ',A6,'  TYPE : ',1I6)
63       FORMAT(1X,'                NOT IMPLEMENTED')
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  COMPLEMENT DU VECTEUR EN CAS DE PARALLELISME
C
      IF(NCSIZE.GT.1) THEN
        CALL PARCOM(X,2,MESH)
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
