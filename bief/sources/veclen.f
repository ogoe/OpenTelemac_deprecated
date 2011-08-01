C                       *****************
                        SUBROUTINE VECLEN
C                       *****************
C
     * (LV,NDP,IKLE,NELEM,NELMAX,NPOIN,V)
C
C***********************************************************************
C BIEF VERSION 5.1           11/03/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        IDEE DE BASE DE J.-P. GREGOIRE
C***********************************************************************
C
C     FONCTION : RECHERCHE D'UNE LONGUEUR DE VECTEUR SANS DEPENDANCE
C                ARRIERE POUR DES BOUCLES PORTANT SUR LES ELEMENTS
C
C                ON NE RECHERCHE QUE LES VALEURS 1, 64, 128,
C                                                256, 512, OU 1024
C
C                LE PRINCIPE EST D'EXECUTER EN MODE SCALAIRE ET
C                VECTORIEL UN ALGORITHME QUI CALCULE LE NOMBRE
C                D'ELEMENTS ADJACENTS A CHAQUE POINT.
C
C                EN MODE VECTORIEL AVEC DEPENDANCES, LE RESULTAT EST
C                FAUX.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   LV           |<-- | LONGUEUR DE VECTEUR ADMISSIBLE
C |   NDP          | -->| NOMBRE DE SOMMETS PAR ELEMENT
C |   IKLE         | -->| TABLE DE CONNECTIVITE
C |   NELEM        | -->| NOMBRE D'ELEMENTS
C |   NELMAX       | -->| NOMBRE TOTAL D'ELEMENTS
C |   NPOIN        | -->| NOMBRE DE POINTS DU MAILLAGE
C |   V            | -->| TABLEAU DE TRAVAIL REEL, DE DIMENSION NPOIN
C |________________|____|______________________________________________
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PREDON
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      USE BIEF !, EX_VECLEN => VECLEN
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(INOUT) :: LV
      INTEGER, INTENT(IN)    :: NELEM,NELMAX,NDP,NPOIN
      INTEGER, INTENT(IN)    :: IKLE(NELMAX,NDP)
C
      DOUBLE PRECISION, INTENT(INOUT) :: V(NPOIN)
C
C-----------------------------------------------------------------------
C
      IF(NDP.EQ.3) THEN
        CALL VECLE3(LV,IKLE,NELEM,NELMAX,NPOIN,V)
      ELSEIF(NDP.EQ.4) THEN
        CALL VECLE4(LV,IKLE,NELEM,NELMAX,NPOIN,V)
      ELSEIF(NDP.EQ.6) THEN
        CALL VECLE6(LV,IKLE,NELEM,NELMAX,NPOIN,V)
      ELSE
        IF(LNG.EQ.1) WRITE(LU,50) NDP
        IF(LNG.EQ.2) WRITE(LU,60) NDP
50      FORMAT(1X,'VECLEN : VALEUR DE NDP NON PREVUE : ',1I6)
60      FORMAT(1X,'VECLEN : UNEXPECTED VALUE OF NDP: ',1I6)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
