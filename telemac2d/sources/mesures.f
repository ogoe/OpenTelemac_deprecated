C                       ******************
                        SUBROUTINE MESURES
C                       ******************
C
     *(ITER,TT)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.2          17/08/01    J-M HERVOUET
C***********************************************************************
C
C  USER SUBROUTINE
C
C  FUNCTION  : READS MEASUREMENTS H, U AND V AT TIME AT
C              GIVES THE CORRESPONDING WEIGHTS ALPHA1, ALPHA2 AND ALPHA3
C 
C-----------------------------------------------------------------------
C  ARGUMENTS 
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|_______________________________________________
C |      AT        | -->| TEMPS DE LA MESURE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT :
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC2D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: ITER
      DOUBLE PRECISION, INTENT(IN) :: TT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION TPS,C
      LOGICAL OKH,OKU,OKV
      INTEGER I
C
C-----------------------------------------------------------------------
C
      IF(T2D_FILES(T2DREF)%NAME(1:1).NE.' ') THEN
C
C-----------------------------------------------------------------------
C
C       WHEN MEASUREMENTS ARE IN A SELAFIN FILE
C
        CALL FIND_IN_SEL(HD,TEXTE(4)(1:16),T2D_FILES(T2DREF)%LU,
     *                   W,OKH,RECORD=ITER,TIME=TPS)
        CALL FIND_IN_SEL(UD,TEXTE(1)(1:16),T2D_FILES(T2DREF)%LU,
     *                   W,OKU,RECORD=ITER,TIME=TPS)
        CALL FIND_IN_SEL(VD,TEXTE(2)(1:16),T2D_FILES(T2DREF)%LU,
     *                   W,OKV,RECORD=ITER,TIME=TPS)
C
        IF(.NOT.OKH.OR..NOT.OKU.OR..NOT.OKV) THEN
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'MESURES : PROBLEME DE LECTURE DE HD, UD OU VD'
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'MESURES : PROBLEM WHEN READIND HD, UD, OR VD'
          ENDIF
          CALL PLANTE(1)
          STOP
        ENDIF
        IF(ABS(TT-TPS).GT.1.D-3) THEN
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'MESURES : PROBLEME DE LECTURE DU TEMPS'
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'MESURES : PROBLEM WHEN READIND TIME'
          ENDIF
          CALL PLANTE(1)
          STOP      
        ENDIF
C       UD AND VD MAY BE QUASI-BUBBLE
C       (BUT ALPHA2 AND ALPHA3 WILL BE SET TO 0)
        IF(UD%ELM.EQ.12) THEN
          CALL CHGDIS(UD,11,12,MESH)
          CALL CHGDIS(VD,11,12,MESH)
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
C
C-----------------------------------------------------------------------
C
C      CASE TO BE IMPLEMENTED HERE (OTHER FILE FORMAT, ETC.)
C
       IF(LNG.EQ.1) WRITE(LU,*) 'MESURES A PROGRAMMER DANS MESURES'
       IF(LNG.EQ.2) WRITE(LU,*) 'MEASUREMENTS TO IMPLEMENT IN MESURES'
       CALL PLANTE(1)
       STOP
C
C-----------------------------------------------------------------------
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C     WEIGHT FUNCTIONS FOR ALL THE TIME-STEPS
C
      CALL VECTOR(T1,'=','MASBAS          ',
     *            HD%ELM,1.D0,T3,T3,T3,T3,T3,
     *            T3,MESH,MSK,MASKEL)
      CALL OS( 'X=Y     ' , ALPHA1 , T1 , T1 , C )
C
C     CASE OF QUASI-BUBBLE ELEMENT FOR UD
      IF(HD%ELM.NE.UD%ELM) THEN
        CALL VECTOR(T1,'=','MASBAS          ',
     *              UD%ELM,1.D0,T3,T3,T3,T3,T3,
     *              T3,MESH,MSK,MASKEL)
      ENDIF
C
      CALL OS( 'X=Y     ' , ALPHA2 , T1 , T1 , C )
      CALL OS( 'X=Y     ' , ALPHA3 , T1 , T1 , C )
C
C     CANCELLING WEIGHTS FOR QUASI-BUBBLE POINTS
C
      IF(UD%ELM.EQ.12) THEN
        DO I=NPOIN+1,NPOIN+NELEM
          ALPHA2%R(I)=0.D0
          ALPHA3%R(I)=0.D0
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
