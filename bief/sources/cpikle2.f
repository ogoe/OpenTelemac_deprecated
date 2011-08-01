C                       ******************
                        SUBROUTINE CPIKLE2
C                       ******************
C
     *(IKLE3,KLEI3,IKLES,NELEM2,NELMAX2,NPOIN2,NPLAN)
C
C***********************************************************************
C BIEF VERSION 5.9           23/08/99    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FONCTION : EXTENSION DU TABLEAU DES CONNECTIVITES
C            CAS D'UN MAILLAGE DE PRISMES, IKLE CONSTRUIT
C            A PARTIR DE LA CONNECTIVITE DU MAILLAGE DE TRIANGLES.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      IKLE      |<-->|  TABLEAU DES CONNECTIVITES                   |
C |      NELEM     | -->|  NOMBRE D'ELEMENTS
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS
C |      NPOIN     | -->|  NOMBRE DE SOMMETS DU MAILLAGE
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : INBIEF
C
C  SOUS-PROGRAMME APPELE : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NELEM2,NELMAX2,NPOIN2,NPLAN
      INTEGER, INTENT(INOUT) :: IKLES(3,NELEM2)
      INTEGER, INTENT(INOUT) :: IKLE3(NELMAX2,NPLAN-1,6)
      INTEGER, INTENT(INOUT) :: KLEI3(6,NELMAX2,NPLAN-1)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,I
C
C-----------------------------------------------------------------------
C
C     BOTTOM AND TOP OF ALL LAYERS
C
      IF(NPLAN.GE.2) THEN
        DO I = 1,NPLAN-1
          DO IELEM = 1,NELEM2          
            IKLE3(IELEM,I,1) = IKLES(1,IELEM) + (I-1)*NPOIN2
            IKLE3(IELEM,I,2) = IKLES(2,IELEM) + (I-1)*NPOIN2
            IKLE3(IELEM,I,3) = IKLES(3,IELEM) + (I-1)*NPOIN2
            IKLE3(IELEM,I,4) = IKLES(1,IELEM) +  I   *NPOIN2
            IKLE3(IELEM,I,5) = IKLES(2,IELEM) +  I   *NPOIN2
            IKLE3(IELEM,I,6) = IKLES(3,IELEM) +  I   *NPOIN2
            KLEI3(1,IELEM,I) = IKLES(1,IELEM) + (I-1)*NPOIN2
            KLEI3(2,IELEM,I) = IKLES(2,IELEM) + (I-1)*NPOIN2
            KLEI3(3,IELEM,I) = IKLES(3,IELEM) + (I-1)*NPOIN2
            KLEI3(4,IELEM,I) = IKLES(1,IELEM) +  I   *NPOIN2
            KLEI3(5,IELEM,I) = IKLES(2,IELEM) +  I   *NPOIN2
            KLEI3(6,IELEM,I) = IKLES(3,IELEM) +  I   *NPOIN2    
          ENDDO
        ENDDO
      ELSE
        IF(LNG.EQ.1) WRITE(LU,*) 'CPIKLE2 : IL FAUT AU MOINS 2 PLANS'
        IF(LNG.EQ.2) WRITE(LU,*) 'CPIKLE2 : MINIMUM OF 2 PLANES NEEDED'
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
