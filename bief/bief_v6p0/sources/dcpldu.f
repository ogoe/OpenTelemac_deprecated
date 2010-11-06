C                       *****************
                        SUBROUTINE DCPLDU
C                       *****************
C
     *(B,A,MESH,COPY,LV)
C
C***********************************************************************
C BIEF VERSION 5.1           23/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FONCTION : DECOMPOSITION L D U DES MATRICES ELEMENTAIRES CONTENUES
C            DANS LA MATRICE A.
C
C            (A PEUT AUSSI ETRE UN BLOC DE MATRICES, DANS CE CAS ON
C             TRAITE LES MATRICES DIAGONALES DU BLOC).
C
C            ON EXIGE ICI QUE LA DIAGONALE DE A SOIT L'IDENTITE
C
C            CHAQUE MATRICE ELEMENTAIRE EST DECOMPOSEE SOUS LA FORME :
C
C            LE X DE X UE
C
C            LE : TRIANGULAIRE INFERIEURE AVEC DES 1 SUR LA DIAGONALE.
C            DE : DIAGONALE
C            UE : TRIANGULAIRE SUPERIEURE AVEC DES 1 SUR LA DIAGONALE.
C
C                                                T
C            SI LA MATRICE EST SYMETRIQUE : LE =  UE
C
C            LES MATRICES "DE" SONT CONSIDEREES COMME DES DIAGONALES
C            DE TAILLE NPOIN X NPOIN QU'IL FAUT IMAGINER COMPLETEES
C            AVEC DES 1 POUR LES POINTS QUI N'APPARTIENNENT PAS A
C            L'ELEMENT CONSIDERE.
C
C            ON EFFECTUE ENSUITE LE PRODUIT DE TOUTES CES DIAGONALES
C            QUI DONNE LA DIAGONALE DB.
C
C   ||||||   ENFIN : DB EST INVERSEE CAR C'EST SOUS CETTE FORME QU'ELLE
C                    SERA UTILISEE DANS DOWNUP.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      B         |<-- |  MATRICE RESULTAT.
C |      A         |<-- |  MATRICE A.
C |      MESH      | -->|  BLOC DES TABLEAUX D'ENTIERS DU MAILLAGE.
C |      COPY      | -->|  SI .TRUE. A EST COPIEE DANS B.
C |                |    |  SINON ON CONSIDERE QUE B EST DEJA REMPLIE.
C |      LV        | -->|  LONGUEUR DU VECTEUR POUR LA VECTORISATION.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : DECLDU
C
C**********************************************************************
C
      USE BIEF, EX_DCPLDU => DCPLDU
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ) , INTENT(INOUT) :: B
      TYPE(BIEF_OBJ) , INTENT(IN)    :: A
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      LOGICAL        , INTENT(IN)    :: COPY
      INTEGER        , INTENT(IN)    :: LV
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER SA,SB,I
C
C-----------------------------------------------------------------------
C
      IF(A%TYPE.EQ.3) THEN
        SA = 0
      ELSEIF(A%TYPE.EQ.4) THEN
        SA = A%N
      ELSE
        IF (LNG.EQ.1) WRITE(LU,300) A%TYPE
        IF (LNG.EQ.2) WRITE(LU,400) A%TYPE
300     FORMAT(1X,'DCPLDU (BIEF) :',1I6,' TYPE DE A NON PREVU.')
400     FORMAT(1X,'DCPLDU (BIEF) :',1I6,' UNEXPECTED TYPE FOR A.')
        CALL PLANTE(0)
        STOP
      ENDIF
C
      IF(B%TYPE.EQ.3) THEN
        SB = 0
      ELSEIF(B%TYPE.EQ.4) THEN
        SB = B%N
      ELSE
        IF (LNG.EQ.1) WRITE(LU,301) B%TYPE
        IF (LNG.EQ.2) WRITE(LU,401) B%TYPE
301     FORMAT(1X,'DCPLDU (BIEF) :',1I6,' TYPE DE B NON PREVU.')
401     FORMAT(1X,'DCPLDU (BIEF) :',1I6,' UNEXPECTED TYPE FOR B.')
        CALL PLANTE(0)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(SA.EQ.0.AND.SB.EQ.0) THEN
C
         CALL DECLDU(B,A,MESH,COPY,LV)
C
      ELSEIF(SB.GT.0.AND.SA.GT.0) THEN
C
C       ON PREND LES DIAGONALES DU BLOC A.
C
        DO I=1,SB
          CALL DECLDU(B%ADR(I)%P,A%ADR(1+(SB+1)*(I-1))%P,
     *                MESH,COPY,LV)
        ENDDO
C
      ELSEIF(SA.NE.0.AND.SB.EQ.0) THEN
C
C       ON PREND LA PREMIERE DIAGONALE DU BLOC A.
C
        CALL DECLDU(B,A%ADR(1)%P,MESH,COPY,LV)
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,302)
        IF (LNG.EQ.2) WRITE(LU,402)
302     FORMAT(1X,'DCPLDU (BIEF) : CAS NON PREVU')
402     FORMAT(1X,'DCPLDU (BIEF) : UNEXPECTED CASE')
        CALL PLANTE(0)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
 
