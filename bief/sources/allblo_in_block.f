C                       **************************
                        SUBROUTINE ALLBLO_IN_BLOCK
C                       **************************
C
     *( BLO , N , NOMGEN )
C
C***********************************************************************
C BIEF VERSION 5.1            11/07/95    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C  FONCTION  : ALLOCATION EN MEMOIRE DE N BLOCS QUI SERONT PLACES
C              DANS UN BLOC DONNE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   BLO          |<-->| BLOC OU ON VA ALLOUER LES VECTEURS
C |   N            | -->| NOMBRE DE VECTEURS A ALLOUER
C |   NOMGEN       | -->| NOM GENERIQUE FORTRAN DES VECTEURS
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_ALLBLO_IN_BLOCK => ALLBLO_IN_BLOCK
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: BLO
      INTEGER         , INTENT(IN)    :: N
      CHARACTER(LEN=6), INTENT(IN)    :: NOMGEN
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IDEB,I,II
C
      CHARACTER(LEN=6) :: NOM
      CHARACTER(LEN=1) :: CHIFFRE(0:9)
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
      IF(N.LE.BLO%MAXBLOCK) THEN
C
      DO I = 1 , N
C
C  NAME OF THE BLOC
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
          STOP 'TOO MANY BLOCKS IN ALLBLO_IN_BLOCK'
        ENDIF
C
C  ALLOCATION OF THE BLOC
C
        ALLOCATE(BLO%ADR(I)%P)
        CALL ALLBLO(BLO%ADR(I)%P,NOM)
        BLO%N=BLO%N+1
C
      ENDDO
C
      ELSE
C
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'ALLBLO_IN_BLOCK : PLUS DE ',BLO%MAXBLOCK,' (',N,')'
        WRITE(LU,*) '                  BLOCS DEMANDES'
        WRITE(LU,*) '                  CHANGER MAXBLOCK DANS ALLBLO'
      ENDIF
      IF(LNG.EQ.2) THEN
        WRITE(LU,*) 'ALLBLO_IN_BLOCK : MORE THAN '
        WRITE(LU,*) '                 ',BLO%MAXBLOCK,' (',N,')'
        WRITE(LU,*) '                  BLOCKS TO BE ALLOCATED'
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
 
