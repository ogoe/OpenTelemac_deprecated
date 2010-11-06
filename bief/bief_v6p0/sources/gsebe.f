C                       ****************
                        SUBROUTINE GSEBE
C                       ****************
C
     *(B,A,MESH)
C
C***********************************************************************
C BIEF VERSION 5.1           23/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FONCTION : DECOMPOSITION DES MATRICES ELEMENTAIRES CONTENUES
C            DANS LA MATRICE A SUIVANT LA METHODE GAUSS-SEIDEL EBE.
C
C            (A PEUT AUSSI ETRE UN BLOC DE MATRICES, DANS CE CAS ON
C             TRAITE TOUTES LES MATRICES DU BLOC).
C
C            ON EXIGE ICI QUE LA DIAGONALE DE A SOIT L'IDENTITE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      B         |<-- |  MATRICE RESULTAT.
C |      A         |<-- |  MATRICE A.
C |      MESH      | -->|  BLOC DES TABLEAUX D'ENTIERS DU MAILLAGE.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : DECLDU
C
C**********************************************************************
C
      USE BIEF, EX_GSEBE => GSEBE
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(IN) :: A
      TYPE(BIEF_OBJ), INTENT(INOUT) :: B
      TYPE(BIEF_MESH), INTENT(IN) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER SA,SB,I
C
      DOUBLE PRECISION C
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
300     FORMAT(1X,'GSEBE (BIEF) :',1I6,' TYPE DE A NON PREVU.')
400     FORMAT(1X,'GSEBE (BIEF) :',1I6,' UNEXPECTED TYPE FOR A.')
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
301     FORMAT(1X,'GSEBE (BIEF) :',1I6,' TYPE DE B NON PREVU.')
401     FORMAT(1X,'GSEBE (BIEF) :',1I6,' UNEXPECTED TYPE FOR B.')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(SA.EQ.0.AND.SB.EQ.0) THEN
C
C         B%D EST ICI UNE STRUCTURE DE VECTEUR
C         QU'ON MET COMME DIAGONALE BIDON.
          CALL OM( 'M=N     ' , B , A , B%D , C , MESH )
          B%TYPDIA='I'
C
      ELSEIF(SA.GT.0.AND.SB.GT.0) THEN
C
C       ON PREND LES DIAGONALES DU BLOC A.
C
        DO 10 I=1,SB
          CALL OM( 'M=N     ' ,  B%ADR(I)%P ,
     *              A%ADR(1+(SB+1)*(I-1))%P , B%ADR(I)%P%D ,
     *              C , MESH )
          B%ADR(I)%P%TYPDIA='I'
10      CONTINUE
C
      ELSEIF(SA.NE.0.AND.SB.EQ.0) THEN
C
C       ON PREND LA PREMIERE DIAGONALE DU BLOC A.
C
        CALL OM( 'M=N     ' ,B,A%ADR(1)%P,B%D,C,MESH)
        B%TYPDIA='I'
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,302)
        IF (LNG.EQ.2) WRITE(LU,402)
302     FORMAT(1X,'GSEBE (BIEF) : CAS NON PREVU')
402     FORMAT(1X,'GSEBE (BIEF) : UNEXPECTED CASE')
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
 
