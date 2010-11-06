C                       *****************
                        SUBROUTINE CORRSL
C                       *****************
C
     *(NEWSL,OLDSL,ZF,MESH)
C
C***********************************************************************
C BIEF VERSION 5.1           27/11/92    J-M JANIN    (LNH) 30 87 72 84
C                                        J-M HERVOUET (LNH) 30 87 80 18
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
C |      MESH      | -->|  STRUCTURE DE MAILLAGE
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
      USE BIEF, EX_CORRSL => CORRSL
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ) , INTENT(INOUT) :: NEWSL
      TYPE(BIEF_OBJ) , INTENT(IN)    :: OLDSL,ZF
      TYPE(BIEF_MESH), INTENT(IN)    :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NELEM,NELMAX,IELM
C
C-----------------------------------------------------------------------
C
      NELEM = MESH%NELEM
      NELMAX= MESH%NELMAX
C
      IELM=OLDSL%ELM
C
C-----------------------------------------------------------------------
C
C     LE VECTEUR NEWSL EST MARQUE COMME ETANT DISCONTINU
      NEWSL%DIMDISC=IELM
C
C-----------------------------------------------------------------------
C
      IF(IELM.EQ.11) THEN
C
        CALL CRSL11(NEWSL%R,OLDSL%R,ZF%R,MESH%IKLE%I,NELEM,NELMAX)
        NEWSL%ELM=10
C
C-----------------------------------------------------------------------
C
      ELSEIF(IELM.EQ.12) THEN
C
        CALL CRSL12(NEWSL%R,OLDSL%R,ZF%R,MESH%IKLE%I,NELEM,NELMAX)
        NEWSL%ELM=10
C
C-----------------------------------------------------------------------
C
      ELSE
C
         IF(LNG.EQ.1) WRITE(LU,10) IELM
         IF(LNG.EQ.2) WRITE(LU,11) IELM
10       FORMAT(1X,'CORRSL : DISCRETISATION INCONNUE :',I6)
11       FORMAT(1X,'CORRSL : UNKNOWN DISCRETIZATION :',I6)
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
