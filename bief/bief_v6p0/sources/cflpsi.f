C                       *****************
                        SUBROUTINE CFLPSI
C                       *****************
C
     *(SYGMA,U,V,DT,IELM,MESH,MSK,MASKEL)
C
C***********************************************************************
C BIEF VERSION 5.3           17/08/94      C MOULIN   (LNH) 30 87 83 81
C                                          + MODIFS JMH LE 17/08/94
C***********************************************************************
C
C  FONCTION  : CALCULE LE NOMBRE DE COURANT EN CHAQUE POINT DU MAILLAGE
C              POUR CHAQUE PAS DE TEMPS.
C
C  ATTENTION : LES COORDONNEES SONT DONNEES ICI PAR ELEMENTS
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      SYGMA     |<-- | NOMBRE DE COURANT.
C |      U         | -->| VITESSE SUIVANT X.                           |
C |      V         | -->| VITESSE SUIVANT Y.                           |
C |      DT        | -->| PAS DE TEMPS DU CALCUL.                      |
C |      IELM      | -->| TYPE D'ELEMENT.                              |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : TELMAC
C
C SOUS-PROGRAMME APPELE : LISSAG
C
C***********************************************************************
C
      USE BIEF, EX_CFLPSI => CFLPSI
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: SYGMA
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: U,V,MASKEL
      DOUBLE PRECISION, INTENT(IN)    :: DT
      INTEGER         , INTENT(IN)    :: IELM
      TYPE(BIEF_MESH) , INTENT(INOUT) :: MESH
      LOGICAL         , INTENT(IN)    :: MSK
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C TRIANGLES P1
C
      IF(IELM.EQ.11) THEN
C
        CALL CFLP11(U%R,V%R,MESH%XEL%R,MESH%YEL%R,
     *              MESH%IKLE%I,MESH%NELEM,MESH%NELMAX,MESH%W%R)
C
C-----------------------------------------------------------------------
C
C TRIANGLES QUASI-BULLE
C
      ELSEIF(IELM.EQ.12) THEN
C
        CALL CFLP12(U%R,V%R,MESH%XEL%R,MESH%YEL%R,
     *              MESH%IKLE%I,MESH%NELEM,MESH%NELMAX,MESH%W%R)
C
C-----------------------------------------------------------------------
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,100) IELM
        IF (LNG.EQ.2) WRITE(LU,101) IELM
100     FORMAT(1X,'CFLPSI : IELM = ',1I6,'  CAS NON PREVU |')
101     FORMAT(1X,'CFLPSI: IELM = ',1I6,' COMBINATION NOT AVAILABLE |')
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C ASSEMBLAGE DES LIJ
C
      CALL ASSVEC(SYGMA%R,MESH%IKLE%I,NBPTS(IELM),
     *            MESH%NELEM,MESH%NELMAX,IELM,
     *            MESH%W%R,.TRUE.,MESH%LV,MSK,MASKEL%R)
C
C-----------------------------------------------------------------------
C
C RESULTAT FINAL
C
C     MASSE DES BASES MISE DANS LE TABLEAU DE TRAVAIL DE BIEF
C
      CALL VECTOR(MESH%T,'=','MASBAS          ',
     *            IELM,1.D0,U,U,U,U,U,U,MESH,MSK,MASKEL)
C     CORRECTION JMH 27/01/2003
      IF(NCSIZE.GT.1) CALL PARCOM(MESH%T,2,MESH)
C
C DIVISION PAR LA MASSE DES BASES
C
      CALL CPSTVC(MESH%T,SYGMA)
      CALL OS( 'X=CY/Z  ' , SYGMA , SYGMA , MESH%T , DT ,2,0.D0,1.D-6)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
