C                       *****************
                        SUBROUTINE ENTART
C                       *****************
C
     *(ITITRE,X,LT,NBR,NBRTOT,ALEMON,ALEMUL,BALAYE)
C
C***********************************************************************
C
C  ARTEMIS VERSION 5.1   02/06/99    D. AELBRECHT (LNH) 01 30 87 74 12
C
C  LINKED TO BIEF VERS. 5.0          J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C      FONCTION: ECRIT LES ENTETES DES DIFFERENTS CALCULS D'AGITATION
C                SUR LE LISTING.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   ITITRE       | -->| TYPE DE TITRE A IMPRIMER                     |
C |   X            | -->| REEL A IMPRIMER                              |
C |   LT           | -->| NUMERO DU CALCUL                             |
C |   NBR          | -->| NUMERO DE LA PERIODE OU DIRECTION EN COURS   |
C |   NBRTOT       | -->| NOMBRES TOTAL DES PERIODES OU DIRECTIONS     |
C |   ALEMON       | -->| VRAI SI HOULE ALEATOIRE MONODIRECTIONNELLE   |
C |   ALEMUL       | -->| VRAI SI HOULE ALEATOIRE MULTIDIRECTIONNELLE  |
C |   BALAYE       | -->| VRAI SI BALAYAGE EN PERIODES                 |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : ARTEMI
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      USE INTERFACE_ARTEMIS, EX_ENTART => ENTART
C                  
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER LT,ITITRE,NBR,NBRTOT
C
      DOUBLE PRECISION X
C
      LOGICAL ALEMON,ALEMUL,BALAYE
C
      CHARACTER*32 TEXTFR(5),TEXTGB(5)
C
C-----------------------------------------------------------------------
C
      DATA TEXTFR / 'PERIODE     ' ,
     *              ' SECONDES   ' ,
     *              'DIRECTION   ' ,
     *              ' DEGRES     ' ,
     *              '            ' /
      DATA TEXTGB / 'PERIOD      ' ,
     *              ' SECONDS    ' ,
     *              'DIRECTION   ' ,
     *              ' DEGREES    ' ,
     *              '            ' /
C
C CAS DE LA HOULE REGULIERE
C
      IF (.NOT.ALEMON .AND. .NOT.ALEMUL .AND. .NOT.BALAYE) THEN
         NBRTOT = 1
      ENDIF
C
C-----------------------------------------------------------------------
C
C   IMPRESSION : PERIODE DE LA HOULE CALCULEE
C
C
      IF (ITITRE.EQ.1) THEN
         IF(LNG.EQ.1) WRITE(LU,100) TEXTFR(1),NBR,NBRTOT,X,TEXTFR(2)
         IF(LNG.EQ.2) WRITE(LU,100) TEXTGB(1),NBR,NBRTOT,X,TEXTGB(2)
      ENDIF
C
100   FORMAT(/,80('='),/,7X,A8,I2,'/',I2,' : ',F12.4,A10,/)
C
C
C-----------------------------------------------------------------------
C
C   IMPRESSION : DIRECTION DE LA HOULE CALCULEE
C
C
      IF (ITITRE.EQ.2) THEN
         IF(LNG.EQ.1) WRITE(LU,110) TEXTFR(3),NBR,NBRTOT,X,TEXTFR(4)
         IF(LNG.EQ.2) WRITE(LU,110) TEXTGB(3),NBR,NBRTOT,X,TEXTGB(4)
      ENDIF
C
110   FORMAT(/,7X,A10,I2,'/',I2,' : ',F12.4,A10,/)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
