C                       *****************
                        SUBROUTINE ASSEX3
C                       *****************
C
     *(XM,STO,NAME,IELM1,IELM2,TYPEXT,XMT,DIM1XMT,DIM2XMT,STOXMT,
     * MESH,NELMAX,ELTSEG,ORISEG)
C
C***********************************************************************
C BIEF VERSION 6.0      05/02/2010    J-M HERVOUET (LNHE) 01 30 87 80 18
C                                       
C***********************************************************************
C
C FONCTION : ASSEMBLING EXTRA-DIAGONAL TERMS OF MATRICES
C            IN THE CASE OF EDGE-BASED STORAGE
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C |      XM        | -->| ASSEMBLED OFF-DIAGONAL TERMS  
C |     STO        | -->| STORAGE REQUIRED IN XM 1: EBE  3: EDGE-BASED  
C |      XMT       | -->| OFF-DIAGONAl TERMS OF THE WORK MATRIX  
C |    DIM1XMT     | -->| FIRST DIMENSION OF XMT 
C |    DIM2XMT     | -->| SECOND DIMENSION OF XMT 
C |    STOXMT      | -->| STORAGE OF OFF-DIAGONAL TERMS
C |                |    | 1: XMT(NELMAX,*)  2: XMT(*,NELMAX) 
C |                | -->|  
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   : 
C
C***********************************************************************
C
      USE BIEF, EX_ASSEX3 => ASSEX3
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(INOUT) :: STO
      CHARACTER(LEN=6), INTENT(IN)    :: NAME
      INTEGER         , INTENT(IN)    :: IELM1,IELM2,NELMAX
      INTEGER         , INTENT(IN)    :: DIM1XMT,DIM2XMT,STOXMT
      INTEGER         , INTENT(IN)    :: ELTSEG(NELMAX,*)
      INTEGER         , INTENT(IN)    :: ORISEG(NELMAX,3)
      CHARACTER(LEN=1), INTENT(IN)    :: TYPEXT
      DOUBLE PRECISION, INTENT(INOUT) :: XMT(DIM1XMT,DIM2XMT)
      DOUBLE PRECISION, INTENT(INOUT) :: XM(*)
      TYPE(BIEF_MESH) , INTENT(IN)    :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NELEM,STOM
C
C-----------------------------------------------------------------------
C
C  EXTRACTION DES CARACTERISTIQUES DE LA MATRICE M
C
      STOM = STO
      IF(STOM.NE.1) THEN
        IF (LNG.EQ.1) WRITE(LU,500) NAME,STOM
        IF (LNG.EQ.2) WRITE(LU,501) NAME,STOM
500     FORMAT(1X,'ASSEX3 (BIEF) : MATRICE M (NOM REEL : ',A6,')',/,1X,
     *            '                STOCKAGE NON PREVU : ',1I6)
501     FORMAT(1X,'ASSEX3 (BIEF) : MATRIX  M (REAL NAME:',A6,')',/,1X,
     *            '                UNEXPECTED STORAGE: ',1I6)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(DIMENS(IELM1).NE.MESH%DIM) THEN
C        MATRICE DE BORD : NON TRAITE ICI
         IF (LNG.EQ.1) WRITE(LU,100) NAME
         IF (LNG.EQ.2) WRITE(LU,101) NAME
         IF (LNG.EQ.1) WRITE(LU,200) IELM1,IELM2
         IF (LNG.EQ.2) WRITE(LU,201) IELM1,IELM2
         IF (LNG.EQ.1) WRITE(LU,300)
         IF (LNG.EQ.2) WRITE(LU,301)
         CALL PLANTE(0)
         STOP
      ENDIF
C
      NELEM  = MESH%NELEM
C
C-----------------------------------------------------------------------
C
      IF(IELM1.EQ.11.AND.IELM2.EQ.11) THEN
C
C       MATRICE DE TRIANGLES P1-P1
C
        IF(TYPEXT.EQ.'S') THEN
          CALL AS3_1111_S(XM,NBSEG(11),
     *                    XMT,NELMAX,NELEM,
     *                    ELTSEG(1,1),ELTSEG(1,2),ELTSEG(1,3))
        ELSEIF(TYPEXT.EQ.'Q') THEN
          CALL AS3_1111_Q(XM,NBSEG(11),
     *                    XMT,NELMAX,NELEM,
     *                    ELTSEG(1,1),ELTSEG(1,2),ELTSEG(1,3),
     *                    ORISEG(1,1),ORISEG(1,2),ORISEG(1,3))     
        ENDIF
C
      ELSEIF(IELM1.EQ.11.AND.IELM2.EQ.12) THEN
C
C       MATRICE DE TRIANGLES P1-QB
C
          CALL AS3_1112(XM,NBSEG(IELM1),NBSEG(IELM2),
     *                  XMT,NELMAX,NELEM,
     *                  ELTSEG(1,1),ELTSEG(1,2),ELTSEG(1,3),
     *                  ELTSEG(1,4),ELTSEG(1,5),ELTSEG(1,6),
     *                  ORISEG(1,1),ORISEG(1,2),ORISEG(1,3))
C 
      ELSEIF(IELM1.EQ.11.AND.IELM2.EQ.13) THEN
C
C       MATRICE DE TRIANGLES P1-QUADRATIQUE
C
          CALL AS3_1113(XM,NBSEG(IELM1),NBSEG(IELM2),
     *                  XMT,NELMAX,NELEM,ELTSEG,ORISEG)
C 
      ELSEIF(IELM1.EQ.13.AND.IELM2.EQ.11) THEN
C
C       MATRICE DE TRIANGLES QUADRATIQUE-P1
C
          CALL AS3_1311(XM,NBSEG(IELM2),NBSEG(IELM1),
     *                  XMT,NELMAX,NELEM,ELTSEG,ORISEG)     
C
      ELSEIF(IELM1.EQ.12.AND.IELM2.EQ.11) THEN
C
C       MATRICE DE TRIANGLES P1-QB
C
          CALL AS3_1211(XM,NBSEG(11),NBSEG(12),
     *                  XMT,NELMAX,NELEM,
     *                  ELTSEG(1,1),ELTSEG(1,2),ELTSEG(1,3),
     *                  ELTSEG(1,4),ELTSEG(1,5),ELTSEG(1,6),
     *                  ORISEG(1,1),ORISEG(1,2),ORISEG(1,3))
C 
      ELSEIF(IELM1.EQ.12.AND.IELM2.EQ.12) THEN
C
C       MATRICE DE TRIANGLES QB-QB
C
        IF(TYPEXT.EQ.'S') THEN
          CALL AS3_1212_S(XM,NBSEG(11),NBSEG(12),XMT,NELMAX,NELEM,
     *                    ELTSEG(1,1),ELTSEG(1,2),ELTSEG(1,3),
     *                    ELTSEG(1,4),ELTSEG(1,5),ELTSEG(1,6))
        ELSEIF(TYPEXT.EQ.'Q') THEN
          CALL AS3_1212_Q(XM,NBSEG(11),NBSEG(12),
     *                    XMT,NELMAX,NELEM,
     *                    ELTSEG(1,1),ELTSEG(1,2),ELTSEG(1,3),
     *                    ELTSEG(1,4),ELTSEG(1,5),ELTSEG(1,6),
     *                    ORISEG(1,1),ORISEG(1,2),ORISEG(1,3))     
        ENDIF
C
      ELSEIF(IELM1.EQ.13.AND.IELM2.EQ.13) THEN
C
C       MATRICE DE TRIANGLE QUADRATIQUE
C
        IF(TYPEXT.EQ.'S') THEN
          CALL AS3_1313_S(XM,NBSEG(IELM1),
     *                    XMT,DIM1XMT,DIM2XMT,STOXMT,
     *                    NELMAX,NELEM,ELTSEG)
        ELSEIF(TYPEXT.EQ.'Q') THEN
          CALL AS3_1313_Q(XM,NBSEG(IELM1),
     *                    XMT,DIM1XMT,DIM2XMT,STOXMT,
     *                    NELMAX,NELEM,ELTSEG,ORISEG)     
        ENDIF
C
      ELSEIF(IELM1.EQ.41.AND.IELM2.EQ.41) THEN
C
C       MATRICE DE PRISMES 
C
        IF(TYPEXT.EQ.'S') THEN
          CALL AS3_4141_S(XM,NBSEG(IELM1),
     *                    XMT,DIM1XMT,DIM2XMT,STOXMT,
     *                    NELMAX,NELEM,ELTSEG)
        ELSEIF(TYPEXT.EQ.'Q') THEN
          CALL AS3_4141_Q(XM,NBSEG(IELM1),
     *                    XMT,DIM1XMT,DIM2XMT,STOXMT,
     *                    NELMAX,NELEM,ELTSEG,ORISEG)     
        ENDIF
C
      ELSE
C
C        COMBINAISON DE IELM1 ET IELM2 NON PREVUE : ERREUR
C
         IF (LNG.EQ.1) WRITE(LU,100) NAME
         IF (LNG.EQ.2) WRITE(LU,101) NAME
         IF (LNG.EQ.1) WRITE(LU,200) IELM1,IELM2
         IF (LNG.EQ.2) WRITE(LU,201) IELM1,IELM2
         IF (LNG.EQ.1) WRITE(LU,300)
         IF (LNG.EQ.2) WRITE(LU,301)
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  CODAGE DU NOUVEAU TYPE DE STOCKAGE
C
      STO=3
C
C-----------------------------------------------------------------------
C
100   FORMAT(1X,'ASSEX3 (BIEF) : MATRICE M (NOM REEL : ',A6,')')
200   FORMAT(1X,'                IELM1 = ',1I6,' IELM2 = ',1I6)
300   FORMAT(1X,'                CAS NON PREVU')
C
101   FORMAT(1X,'ASSEX3 (BIEF) : MATRIX  M (REAL NAME:',A6,')')
201   FORMAT(1X,'                IELM1 = ',1I6,' IELM2 = ',1I6)
301   FORMAT(1X,'                THIS CASE IS NOT IMPLEMENTED')
C
C-----------------------------------------------------------------------
C
      RETURN
      END
