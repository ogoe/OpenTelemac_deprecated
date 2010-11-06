C                       *****************
                        SUBROUTINE DECV21
C                       *****************
C
     *(TETA,SL,ZF,IKLE,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.1           09/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FONCTION : MARQUAGE DES BANCS DECOUVRANTS
C
C            ELEMENT DECOUVRANT : TETA = 0.
C            ELEMENT NORMAL     : TETA = 1.
C
C            LE CRITERE DE DECOUVREMENT EST CELUI DE J.-M. JANIN :
C            LE FOND D'UN POINT DE L'ELEMENT EST PLUS HAUT QUE LA
C            SURFACE LIBRE D'UN AUTRE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      TETA      |<-- |  INDICATEUR (PAR ELEMENT)                    |
C |      SL,ZF     | -->|  SURFACE LIBRE ET FOND                       |
C |      MESH      | -->|  STRUCTURE DE MAILLAGE
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : DECVRT
C
C  SOUS-PROGRAMME APPELE : OS
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(IN)  :: NELEM,NELMAX
      INTEGER         , INTENT(IN)  :: IKLE(NELMAX,*)
      DOUBLE PRECISION, INTENT(OUT) :: TETA(NELEM)
      DOUBLE PRECISION, INTENT(IN)  :: SL(*),ZF(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM
C
      DOUBLE PRECISION SL1,SL2,SL3,SL4,ZF1,ZF2,ZF3,ZF4
C
      INTRINSIC MAX,MIN
C
C-----------------------------------------------------------------------
C
      CALL OV( 'X=C     ' , TETA , TETA , TETA , 1.D0 , NELEM )
C
C-----------------------------------------------------------------------
C
         DO 1 IELEM = 1 , NELEM
C
           SL1 = SL(IKLE(IELEM,1))
           SL2 = SL(IKLE(IELEM,2))
           SL3 = SL(IKLE(IELEM,3))
           SL4 = SL(IKLE(IELEM,4))
C
           ZF1 = ZF(IKLE(IELEM,1))
           ZF2 = ZF(IKLE(IELEM,2))
           ZF3 = ZF(IKLE(IELEM,3))
           ZF4 = ZF(IKLE(IELEM,4))
C
           IF(MAX(ZF1,ZF2,ZF3,ZF4).GT.MIN(SL1,SL2,SL3,SL4)) THEN
             TETA(IELEM) = 0.D0
           ENDIF
C
1        CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
