C                       *****************
                        SUBROUTINE STEADY
C                       *****************
C
     *(H1,H2,NPH,U1,U2,NPU,V1,V2,NPV,NTRAC,T1,T2,NPT,CRIPER,ARRET)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.8   05/09/07   J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C     FONCTION  : CONTROLE SI UN ETAT PERMANENT EST ATTEINT
C
C     NOTE : ARRET N'EST PAS INITIALISE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  H1,H2,NPH     | -->| HAUTEURS A COMPARER ET NOMBRE DE POINTS
C |  U1,U2,NPU     | -->| VITESSES A COMPARER ET NOMBRE DE POINTS
C |  V1,V2,NPU     | -->| VITESSES A COMPARER ET NOMBRE DE POINTS
C |    TRAC        | -->| LOGIQUE INDIQUANT S'IL Y A UN TRACEUR.
C |  T1,T2,NPT     | -->| TRACEURS A COMPARER ET NOMBRE DE POINTS
C |  CRIPER        | -->| CRITERES D'ARRET
C |                |    | DANS L'ORDRE SUIVANT : H , U , V , T
C |    ARRET       |<-- | LOGIQUE MIS A TRUE SI PERMANENT ATTEINT
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)          :: NPH,NPU,NPV,NPT,NTRAC
      LOGICAL, INTENT(INOUT)       :: ARRET
      DOUBLE PRECISION, INTENT(IN) :: H1(NPH),H2(NPH),U1(NPU),U2(NPU)
      DOUBLE PRECISION, INTENT(IN) :: V1(NPV),V2(NPV)
      DOUBLE PRECISION, INTENT(IN) :: CRIPER(3)
      TYPE(BIEF_OBJ)  , INTENT(IN) :: T1,T2
C   
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,ITRAC
C
C-----------------------------------------------------------------------
C
C  CONTROLE DE LA HAUTEUR
C
      DO 10 I = 1 , NPH
        IF(ABS(H1(I)-H2(I)).GT.CRIPER(1)) GO TO 1000
10    CONTINUE
C
C-----------------------------------------------------------------------
C
C  CONTROLE DE U
C
      DO 20 I = 1 , NPU
        IF(ABS(U1(I)-U2(I)).GT.CRIPER(2)) GO TO 1000
20    CONTINUE
C
C-----------------------------------------------------------------------
C
C  CONTROLE DE V
C
      DO 30 I = 1 , NPV
        IF(ABS(V1(I)-V2(I)).GT.CRIPER(2)) GO TO 1000
30    CONTINUE
C
C-----------------------------------------------------------------------
C
C  CONTROLES DU TRACEUR
C
      IF(NTRAC.GT.0) THEN
C
      DO ITRAC=1,NTRAC
        DO I = 1 , NPT
          IF(ABS(T1%ADR(ITRAC)%P%R(I)
     *          -T2%ADR(ITRAC)%P%R(I)).GT.CRIPER(3)) GO TO 1000
        ENDDO
      ENDDO
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      ARRET=.TRUE.
      IF(LNG.EQ.1) WRITE(LU,100)
      IF(LNG.EQ.2) WRITE(LU,200)
100   FORMAT(/,1X,'ETAT PERMANENT ATTEINT')
200   FORMAT(/,1X,'THE STEADY STATE HAS BEEN REACHED')
C
C-----------------------------------------------------------------------
C
1000  CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
