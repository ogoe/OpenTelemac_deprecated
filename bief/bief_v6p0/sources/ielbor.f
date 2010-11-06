C                       ***********************
                        INTEGER FUNCTION IELBOR
C                       ***********************
C
     *( IELM , I )
C
C***********************************************************************
C BIEF VERSION 5.9             06/02/08      J-M HERVOUET 01 30 87 80 18
C***********************************************************************
C
C  FONCTION  : DONNE LE TYPE D'ELEMENT DE BORD CORRESPONDANT A UN
C              TYPE D'ELEMENT DONNE DANS LE DOMAINE.
C
C              LORSQU'IL Y A PLUSIEURS TYPES (CAS DES PRISMES PAR
C              EXEMPLE) ON UTILISE L'INDICE I POUR LES DISTINGUER.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   IELM         | -->| TYPE D'ELEMENT SUR LE DOMAINE.
C |   I            | -->| CAS DE PLUSIEURS TYPES D'ELEMENTS DE BORD
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER IELM,I
C
C-----------------------------------------------------------------------
C
      IF(IELM.EQ.10.OR.IELM.EQ.20) THEN
        IELBOR = 0
      ELSEIF(IELM.EQ.11.OR.IELM.EQ.12.OR.IELM.EQ.21) THEN
        IELBOR = 1
      ELSEIF(IELM.EQ.13) THEN
        IELBOR = 2
      ELSEIF(IELM.EQ.30) THEN
        IELBOR = 10
      ELSEIF(IELM.EQ.31) THEN
        IELBOR = 81
      ELSEIF(IELM.EQ.51.AND.I.EQ.1) THEN
        IELBOR = 11
      ELSEIF(IELM.EQ.51.AND.I.EQ.2) THEN
        IELBOR = 61
      ELSEIF(IELM.EQ.50.AND.I.EQ.1) THEN
        IELBOR = 10
      ELSEIF(IELM.EQ.50.AND.I.EQ.2) THEN
        IELBOR = 60
      ELSEIF(IELM.EQ.40.AND.I.EQ.1) THEN
        IELBOR = 10
      ELSEIF(IELM.EQ.40.AND.I.EQ.2) THEN
        IELBOR = 70
      ELSEIF(IELM.EQ.41.AND.I.EQ.1) THEN
        IELBOR = 11
      ELSEIF(IELM.EQ.41.AND.I.EQ.2) THEN
        IELBOR = 71
      ELSE
        IF(LNG.EQ.1) WRITE(LU,100) IELM
        IF(LNG.EQ.2) WRITE(LU,101) IELM
100     FORMAT(1X,'IELBOR (BIEF) : ',1I6,' ELEMENT NON PREVU')
101     FORMAT(1X,'IELBOR (BIEF) : ',1I6,' ELEMENT NOT IMPLEMENTED')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
