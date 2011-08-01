C                       *****************
                        SUBROUTINE HPROPA
C                       *****************
C
     *(HPROP ,HN,H,PROLIN,HAULIN,TETA,NSOUSI)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.8   16/07/07    J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C      FONCTION:  CALCUL DE LA HAUTEUR DE PROPAGATION
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   HPROP        |<-- |  HAUTEUR DE PROPAGATION
C |   HN           |<-- |  HAUTEUR AU PAS DE TEMPS PRECEDENT
C |   H            |<-- |  HAUTEUR
C |   PROLIN       | -->| CORRESPOND AU MOT CLE:"PROPAGATON LINEARISEE"
C |   HAULIN       | -->| PROFONDEUR MOYENNE POUR LA LINEARISATION
C |   TETA         | -->| SEMI-IMPLICITATION SUR H.
C |   NSOUSI       | -->| NOMBRE DE SOUS ITERATIONS
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMME APPELE : OV
C
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)           :: NSOUSI
      LOGICAL, INTENT(IN)           :: PROLIN
      DOUBLE PRECISION, INTENT(IN)  :: TETA,HAULIN
      TYPE(BIEF_OBJ), INTENT(IN)    :: HN,H
      TYPE(BIEF_OBJ), INTENT(INOUT) :: HPROP
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      IF(PROLIN) THEN
        CALL OS( 'X=C     ' , X=HPROP , C=HAULIN    )
      ELSEIF(NSOUSI.EQ.1) THEN
        CALL OS( 'X=Y     ' , X=HPROP , Y=HN )
      ELSE
        CALL OS( 'X=CY    ' , X=HPROP , Y=HN , C=1.D0-TETA )
        CALL OS( 'X=X+CY  ' , X=HPROP , Y=H  , C= TETA )
      ENDIF
C
C-----------------------------------------------------------------------
C
C     CLIPPING HPROP
C
      IF(.NOT.PROLIN) THEN
        CALL OS('X=+(Y,C)',X=HPROP,Y=HPROP,C=0.D0)
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
