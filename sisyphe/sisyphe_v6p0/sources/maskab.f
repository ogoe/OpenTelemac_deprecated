C                       *****************
                        SUBROUTINE MASKAB
C                       *****************
C
     *(HN , Q , QU , QV , NPOIN)
C
C***********************************************************************
C SISYPHE VERSION 5.1                             E. PELTIER    11/09/95
C                                                 C. LENORMANT
C                                                 J.-M. HERVOUET
C***********************************************************************
C
C     FONCTION : SOUS-PROGRAMME UTILISATEUR
C                SERVANT A ELIMINER LES HAUTEURS D'EAU NEGATIVES
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   HN           |<-->| HAUTEUR D' EAU
C |   Q            | -->| DEBIT LIQUIDE
C |   QU , QV      | -->| COMPOSANTES DU DEBIT VECTORIEL
C |   NPOIN        | -->| NOMBRE DE POINTS
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C APPELE PAR: SISYPHE
C***********************************************************************
C
      IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN):: NPOIN
C
      DOUBLE PRECISION, INTENT(IN)    :: HN(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: Q(NPOIN),QU(NPOIN),QV(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C
C-----------------------------------------------------------------------
C
C
C ECRETEMENT DES HAUTEURS D'EAU
C
C
      DO I=1,NPOIN
C
C  TRAITEMENT DES VALEURS NEGATIVES DANS LE DOMAINE
C
         IF(HN(I).LE.0.D0) THEN
C
            Q(I)  = 0.D0
            QU(I) = 0.D0
            QV(I) = 0.D0
C
         ENDIF
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE MASKAB
