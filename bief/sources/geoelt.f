C                       *****************
                        SUBROUTINE GEOELT
C                       *****************
C
     *(SURDET,SURFAC,XEL,YEL,NELEM,NELMAX,IELM)
C
C***********************************************************************
C BIEF VERSION 5.1           10/01/95    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C     FONCTION  : CALCUL DE DETERMINANTS ET DE CERTAINES VALEURS POUR
C                 LES COORDONNEES ISO-PARAMETRIQUES
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    SURDET      |<-- | INVERSE DU DETERMINANT DE LA TRANSFORMEE     |
C |    SURFAC      |<-- | SURFACES DES ELEMENTS                        |
C |    TAILLE      |<-- | TAILLE DES ELEMENTS (LONGUEUR CARACTERISTIQUE)
C |    XEL,YEL     | -->| COORDONNEES DES NOEUDS PAR ELEMENT           |
C |    NELEM       | -->| NOMBRE D'ELEMENTS                            |
C |    IELM        | -->| TYPE D'ELEMENT                               |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PREDON
C
C SOUS-PROGRAMME APPELE : SURVOL
C
C***********************************************************************
C
      USE BIEF, EX_GEOELT => GEOELT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: IELM,NELEM,NELMAX
      DOUBLE PRECISION, INTENT(OUT) :: SURDET(NELEM),SURFAC(NELEM)
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,*),YEL(NELMAX,*)      
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM
C
      DOUBLE PRECISION XSOM(4,2)
C
      DOUBLE PRECISION T12,T13,T22,T23,DET,Z(1)
C
C-----------------------------------------------------------------------
C
      CALL SURVOL(SURFAC, XEL,YEL,Z,NELEM,NELMAX,IELM)
C
      IF(IELM.EQ.11) THEN
C
         DO 400 IELEM = 1 , NELEM
C
         XSOM(1,1) = XEL(IELEM,1)
         XSOM(2,1) = XEL(IELEM,2)
         XSOM(3,1) = XEL(IELEM,3)
         XSOM(1,2) = YEL(IELEM,1)
         XSOM(2,2) = YEL(IELEM,2)
         XSOM(3,2) = YEL(IELEM,3)
C
         T12 = - XSOM(1,1) + XSOM(2,1)
         T13 = - XSOM(1,1) + XSOM(3,1)
         T22 = - XSOM(1,2) + XSOM(2,2)
         T23 = - XSOM(1,2) + XSOM(3,2)
C
         DET = T12*T23 - T22*T13
C
         IF(DET.LT.1.D-20) THEN
           IF(LNG.EQ.1) WRITE(LU,98) IELEM
           IF(LNG.EQ.2) WRITE(LU,99) IELEM
98         FORMAT(1X,'GEOELT: ELEMENT ',1I6,' : DETERMINANT NEGATIF')
99         FORMAT(1X,'GEOELT: ELEMENT ',1I6,' : NEGATIVE DETERMINANT')
           STOP
         ENDIF
C
         SURDET(IELEM) = 1.D0/DET
C
400      CONTINUE
C
      ELSE
            IF(LNG.EQ.1) WRITE(LU,10) IELM
            IF(LNG.EQ.2) WRITE(LU,11) IELM
10          FORMAT(1X,'GEOELT: TYPE D''ELEMENT INCONNU :',1I6)
11          FORMAT(1X,'GEOELT: UNKNOWN TYPE OF ELEMENT :',1I6)
            CALL PLANTE(1)
            STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
