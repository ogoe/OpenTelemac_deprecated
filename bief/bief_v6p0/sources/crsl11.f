C                       *****************
                        SUBROUTINE CRSL11
C                       *****************
C
     *(NEWSL,OLDSL,ZF,IKLE,NELEM,NELMAX)
C
C***********************************************************************
C  BIEF VERSION 5.5          27/11/92    J-M JANIN    (LNH) 30 87 72 84
C
C***********************************************************************
C
C FONCTION : CALCUL D'UNE SURFACE LIBRE CORRIGEE PAR ELEMENTS
C            POUR TENIR COMPTE DES BANCS DECOUVRANTS
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      NEWSL     | -- |  SURFACE LIBRE MODIFIEE, PAR ELEMENTS
C |      OLDSL     | -->|  SURFACE LIBRE REELLE.                       |
C |      ZF        | -->|  COTE DU FOND                                |
C |      IKLE      | -->|  TABLES DE CONNECTIVITE
C |      NELEM     | -->|  NOMBRE D'ELEMENTS
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : PROPAG
C
C  SOUS-PROGRAMME APPELE :
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NELEM,NELMAX
      DOUBLE PRECISION, INTENT(INOUT) :: NEWSL(NELMAX,3)
      DOUBLE PRECISION, INTENT(IN)    :: OLDSL(*),ZF(*)
      INTEGER, INTENT(IN)             :: IKLE(NELMAX,3)      
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IK(3),J(3)
C
C-----------------------------------------------------------------------
C
      INTRINSIC MAX
C
C-----------------------------------------------------------------------
C
C  1) CLASSEMENT PAR ORDRE CROISSANT DES COTES DE FONDS ET CORRECTION
C     EVENTUELLE DES COTES DE SURFACE LIBRE DES ELEMENTS DECOUVERTS
C
C-----------------------------------------------------------------------
C
         DO 4 IELEM = 1 , NELEM
C
           IK(1) = IKLE(IELEM,1)
           IK(2) = IKLE(IELEM,2)
           IK(3) = IKLE(IELEM,3)
           J(1) = 1           
           J(2) = 2           
           J(3) = 3
C
C          TRI DES 3 POINTS
C
           IF(ZF(IK(2)).LT.ZF(IK(1)))  THEN
              J(2)=1
              J(1)=2
           ENDIF
           IF(ZF(IK(3)).LT.ZF(IK(J(2)))) THEN
              J(3)=J(2)
              J(2)=3
              IF(ZF(IK(3)).LT.ZF(IK(J(1)))) THEN
                 J(2)=J(1)
                 J(1)=3
              ENDIF
           ENDIF
C
C          CORRECTIONS
C
           NEWSL(IELEM,J(1))=MAX(ZF(IK(J(2))),OLDSL(IK(J(1))))
           NEWSL(IELEM,J(2))=OLDSL(IK(J(2)))
           NEWSL(IELEM,J(3))=
     *      OLDSL(IK(J(3)))-MAX(0.D0,ZF(IK(J(3)))-OLDSL(IK(J(2))))
C
4        CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
