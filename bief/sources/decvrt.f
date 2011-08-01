C                       *****************
                        SUBROUTINE DECVRT
C                       *****************
C
     *(TETA,SL,ZF,MESH)
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
C  APPELE PAR : PROPAG ET CVDFTR
C
C  SOUS-PROGRAMME APPELE : DECV11
C
C**********************************************************************
C
      USE BIEF, EX_DECVRT => DECVRT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ) , INTENT(INOUT) :: TETA
      TYPE(BIEF_OBJ) , INTENT(IN)    :: SL,ZF
      TYPE(BIEF_MESH), INTENT(IN)    :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NELEM,NELMAX,IELMS,IELMZ
C
C-----------------------------------------------------------------------
C
      IELMS=SL%ELM
      IELMZ=ZF%ELM
C
C-----------------------------------------------------------------------
C
C  DEPLOIEMENT DE LA STRUCTURE DE MAILLAGE
C
      NELEM = MESH%NELEM
      NELMAX= MESH%NELMAX
C
C-----------------------------------------------------------------------
C
      IF(IELMS.EQ.11.AND.IELMZ.EQ.11) THEN
C
        CALL DECV11(TETA%R,SL%R,ZF%R,MESH%IKLE%I,NELEM,NELMAX)
C
      ELSEIF((IELMS.EQ.21.AND.IELMZ.EQ.21).OR.
     *       (IELMS.EQ.12.AND.IELMZ.EQ.12)      ) THEN
C
        CALL DECV21(TETA%R,SL%R,ZF%R,MESH%IKLE%I,NELEM,NELMAX)
C
      ELSE
C
        IF(LNG.EQ.1) WRITE(LU,10) IELMS,IELMZ
        IF(LNG.EQ.2) WRITE(LU,11) IELMS,IELMZ
10      FORMAT(1X,'DECVRT : DISCRETISATION NON PREVUE :'    ,I6,' ',I6)
11      FORMAT(1X,'DECVRT : DISCRETIZATION NOT IMPLEMENTED:',I6,' ',I6)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
