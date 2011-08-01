C                       **********************
                        INTEGER FUNCTION NBFEL
C                       **********************
C
     *(IELM)
C
C***********************************************************************
C BIEF VERSION 5.5         08/04/04    J-M HERVOUET (LNH) 01 30 87 80 18
C***********************************************************************
C
C FONCTION : DONNE LE NOMBRE DE FACES PAR ELEMENT, 3 POUR UN TRIANGLE
C            PAR EXEMPLE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  IELM          | -->| TYPE D'ELEMENT
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER IELM,NDS(0:81,7)
C
      COMMON/NODES/NDS
C
C-----------------------------------------------------------------------
C
      IF(IELM.LT.0.OR.IELM.GT.81) THEN
        IF(LNG.EQ.1) WRITE(LU,200) IELM
        IF(LNG.EQ.2) WRITE(LU,201) IELM
 200    FORMAT(1X,'NBFEL (BIEF) : MAUVAIS ARGUMENT : ',1I6)
 201    FORMAT(1X,'NBFEL (BIEF) : WRONG ARGUMENT: ',1I6)
        CALL PLANTE(1)
        STOP
      ENDIF
C
      NBFEL = NDS(IELM,4)
C
C-----------------------------------------------------------------------
C
      RETURN
      END

 