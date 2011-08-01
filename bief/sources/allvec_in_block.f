C                       **************************
                        SUBROUTINE ALLVEC_IN_BLOCK
C                       **************************
C
     *( BLO , N , NAT , NOMGEN , IELM , NDIM , STATUT )
C
C***********************************************************************
C BIEF VERSION 5.1            11/07/95    J-M HERVOUET (LNH) 30 87 80 18
C                             march 1999  Jacek A. Jankowski pinxit
C***********************************************************************
C
C  FONCTION  : ALLOCATION EN MEMOIRE DE N VECTEURS QUI SERONT PLACES
C              DANS UN BLOC DONNE.
C
c this modification of allvec_in_block allows adding a number of 
c identic numbered vectors to an already existing block without 
c destroying the previous structure
c
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   BLO          |<-->| BLOC OU ON VA ALLOUER LES VECTEURS
C |   N            | -->| NOMBRE DE VECTEURS A ALLOUER
C |   NAT          |<-- | 1: VECTEUR REEL   2:VECTEUR ENTIER
C |   NOMGEN       | -->| NOM GENERIQUE FORTRAN DES VECTEURS
C |   IELM         | -->| TYPE D'ELEMENT DU VECTEUR, OU DIMENSION
C |                |    | (SUIVANT LE STATUT, VOIR PLUS BAS)
C |   NDIM         | -->| DEUXIEME DIMENSION DU VECTEUR
C |   STATUT       | -->| STATUT DU VECTEUR :
C |                |    | 0 : VECTEUR LIBRE, IELM EST ALORS SA DIMENSION
C |                |    | 1 : VECTEUR DEFINI SUR LE MAILLAGE
C |                |    |     IELM EST ALORS LE TYPE D'ELEMENT
C |                |    |     CHANGEMENT DE DISCRETISATION INTERDIT
C |                |    | 2 : COMME 1 MAIS CHANGEMENTS AUTORISES
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_ALLVEC_IN_BLOCK => ALLVEC_IN_BLOCK
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: BLO
      INTEGER         , INTENT(IN)    :: IELM,NDIM,STATUT,NAT,N
      CHARACTER(LEN=6), INTENT(IN)    :: NOMGEN
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IDEB,I,II
C
      CHARACTER(LEN=6) :: NOM
      CHARACTER*1 CHIFFRE(0:9)
      DATA CHIFFRE/'0','1','2','3','4','5','6','7','8','9'/
      SAVE CHIFFRE
C
C-----------------------------------------------------------------------
C
      IDEB = 6
      DO 5 I=5,2,-1
        IF(NOMGEN(I:I).EQ.' ') IDEB = I
5     CONTINUE
C
C-----------------------------------------------------------------------
C
      IF(BLO%N+N.LE.BLO%MAXBLOCK) THEN
C
      IF(N.GT.0) THEN
C
      DO 10 I = BLO%N+1 , BLO%N+N
C
C  NAME OF THE VECTOR
C
        NOM=NOMGEN
        IF(I.LT.10) THEN
          IDEB = MIN(6,IDEB)
          NOM(IDEB:IDEB) = CHIFFRE(I)
        ELSEIF(I.LT.100) THEN
          IDEB = MIN(5,IDEB)
          NOM(IDEB  :IDEB  ) = CHIFFRE(I/10)
          NOM(IDEB+1:IDEB+1) = CHIFFRE(I-10*(I/10))
        ELSEIF(I.LT.1000) THEN
          IDEB = MIN(4,IDEB)
          NOM(IDEB  :IDEB  ) = CHIFFRE(I/100)
          II=I-100*(I/100)
          NOM(IDEB+1:IDEB+1) = CHIFFRE(II/10)
          NOM(IDEB+2:IDEB+2) = CHIFFRE(II-10*(II/10))
        ELSE
          STOP 'MORE THAN 999 VECTORS ASKED IN ALLVEC_IN_BLOCK'
        ENDIF
C
C  ALLOCATION OF THE VECTOR
C
        ALLOCATE(BLO%ADR(I)%P)
        CALL ALLVEC(NAT,BLO%ADR(I)%P,NOM,IELM,NDIM,STATUT)
C
10    CONTINUE
C
      BLO%N=BLO%N+N
C
      ENDIF
C
      ELSE
C
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'ALLVEC_IN_BLOCK : PLUS DE ',BLO%MAXBLOCK,' (',N,')'
        WRITE(LU,*) '                  VECTEURS DEMANDES'
        WRITE(LU,*) '                  CHANGER MAXBLOCK DANS ALLBLO'
      ENDIF
      IF(LNG.EQ.2) THEN
        WRITE(LU,*) 'ALLVEC_IN_BLOCK : MORE THAN '
        WRITE(LU,*) '                  ',BLO%MAXBLOCK,'(',N,')'
        WRITE(LU,*) '                  VECTORS TO BE ALLOCATED'
        WRITE(LU,*) '                  CHANGE MAXBLOCK IN ALLBLO'
      ENDIF
      STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
