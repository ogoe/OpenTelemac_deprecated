C                       *****************
                        SUBROUTINE VECLE3
C                       *****************
C
     * (LV,IKLE,NELEM,NELMAX,NPOIN,V)
C
C***********************************************************************
C BIEF VERSION 5.1           19/08/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        IDEE DE BASE DE J.-P. GREGOIRE
C***********************************************************************
C
C     FONCTION : RECHERCHE D'UNE LONGUEUR DE VECTEUR SANS DEPENDANCE
C                ARRIERE POUR DES BOUCLES PORTANT SUR LES ELEMENTS
C
C                ON NE RECHERCHE QUE LES VALEURS :
C                1 , 64 , 128 , 256 , 512 , 1024
C
C                LE PRINCIPE EST D'EXECUTER EN MODE SCALAIRE ET
C                VECTORIEL UN ALGORITHME QUI CALCULE LE NOMBRE
C                D'ELEMENTS ADJACENTS A CHAQUE POINT.
C
C                EN MODE VECTORIEL AVEC DEPENDANCES, LE RESULTAT EST
C                FAUX.
C
C     ELEMENT TRAITE : TRIANGLE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   LV           |<-- | LONGUEUR DE VECTEUR ADMISSIBLE
C |   IELM         | -->| TYPE D'ELEMENT
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
      USE BIEF !, EX_VECLE3 => VECLE3
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(INOUT)          :: LV
      INTEGER, INTENT(IN)             :: NELEM,NELMAX,NPOIN
      INTEGER, INTENT(IN)             :: IKLE(NELMAX,3)
      DOUBLE PRECISION, INTENT(INOUT) :: V(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IB
      DOUBLE PRECISION Y(1),Z(1)
C
C-----------------------------------------------------------------------
C
      LV = 1024
C
5     CONTINUE
C
C  INITIALISATION DE V A ZERO
C
      CALL OV( 'X=C     ' , V , Y , Z , 0.D0 , NPOIN )
C
C  MODE SCALAIRE
C
      DO 10 IELEM = 1 , NELEM
        V(IKLE(IELEM,1)) = V(IKLE(IELEM,1)) + 1.D0
        V(IKLE(IELEM,2)) = V(IKLE(IELEM,2)) + 1.D0
        V(IKLE(IELEM,3)) = V(IKLE(IELEM,3)) + 1.D0
10    CONTINUE
C
C  MODE VECTORIEL AVEC VECTORISATION FORCEE
C (COMMANDES FUJITSU, PUIS COMMANDES CRAY)
C
      DO 20 IB = 1,(NELEM+LV-1)/LV
*VOCL LOOP,NOVREC
CDIR$ IVDEP
      DO 30 IELEM = 1+(IB-1)*LV , MIN(NELEM,IB*LV)
        V(IKLE(IELEM,1)) = V(IKLE(IELEM,1)) - 1.D0
        V(IKLE(IELEM,2)) = V(IKLE(IELEM,2)) - 1.D0
        V(IKLE(IELEM,3)) = V(IKLE(IELEM,3)) - 1.D0
30    CONTINUE
20    CONTINUE
C
C-----------------------------------------------------------------------
C
      IF(DOT(NPOIN,V,V).GT.0.5D0.AND.LV.NE.1) THEN
        LV = LV/2
        IF(LV.LT.64) THEN
          LV = 1
          GO TO 1000
        ENDIF
        GO TO 5
      ENDIF
C
1000  CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
 
