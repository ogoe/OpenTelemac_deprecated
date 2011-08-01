C                       ********************
                        SUBROUTINE COMP_IKLE
C                       ********************
C
     *(IKLE,IKLBOR,ELTSEG,NBOR,IELM,NELEM,NELMAX,NPOIN,NPTFR)
C
C***********************************************************************
C BIEF VERSION 5.9      20/03/2008    J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C FONCTION : EXTENSION DES TABLEAUX DES CONNECTIVITES ET DU TABLEAU NBOR
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      IKLE      |<-->|  CONNECTIVITY TABLE FOR ALL POINTS 
C |      IKLBOR    |<-->|  CONNECTIVITY TABLE FOR BOUNDARY POINTS
C |      ELTSEG    | -->|  SEGMENT NUMBERS OF AN ELEMENT 
C |      NBOR      |<-->|  GLOBAL NUMBERS OF BOUNDARY POINTS                       
C |      IELM      | -->|  TYPE OF ELEMENT                    
C |      NELEM     | -->|  NOMBRE D'ELEMENTS
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS
C |      NPOIN     | -->|  NOMBRE DE SOMMETS DU MAILLAGE
C |      NPTFR     | -->|  NUMBER OF (LINEAR) BOUNDARY POINTS
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : INBIEF
C
C  SOUS-PROGRAMME APPELE : CPIK12
C
C**********************************************************************
C
      USE BIEF   !, EX_COMP_IKLE => COMP_IKLE
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NELEM,NELMAX,IELM,NPOIN,NPTFR
      INTEGER, INTENT(IN)    :: ELTSEG(NELMAX,3)
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,*),IKLBOR(NPTFR,*),NBOR(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      IF(IELM.EQ.12) THEN
C
        CALL CPIK12(IKLE,NELEM,NELMAX,NPOIN)
C
      ELSEIF(IELM.EQ.13.OR.IELM.EQ.14) THEN
C
        CALL CPIK13(IKLE,IKLBOR,ELTSEG,NBOR,NELEM,NELMAX,NPOIN,NPTFR)
C
      ELSE
C
        IF(LNG.EQ.1) WRITE(LU,10) IELM
        IF(LNG.EQ.2) WRITE(LU,11) IELM
10      FORMAT(1X,'CPIKLE : DISCRETISATION NON PREVUE :'    ,I6)
11      FORMAT(1X,'CPIKLE: DISCRETIZATION NOT IMPLEMENTED:',I6)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
