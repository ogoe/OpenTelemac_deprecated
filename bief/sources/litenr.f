C                       *****************
                        SUBROUTINE LITENR
C                       *****************
C
     *(VARSOR,CLAND,
     * NPRE,STD,HIST,NHIST,NPOIN,AT,TEXTPR,TEXTLU,
     * NVAR,VARCLA,NVARCL,TROUVE,ALIRE,W,LISTIN,MAXVAR,
     * NPOIN_PREV,NPLAN_PREV,WD)
C
C***********************************************************************
C BIEF VERSION 5.9        09/12/08    J-M HERVOUET (LNHE) 01 30 71 80 18
C
C***********************************************************************
C
C   FUNCTION: READING A TIME-STEP IN A SELAFIN FILE
C
C             THE FILE IS ASSUMED OPENED AND THE TIME-STEP READY TO READ
C
C             IN 3D : INTERPOLATION IF THE NUMBER OF PLANES IS DIFFERENT
C
C
C
C  09/12/2008 : STD IS NOW A STRING OF ANY SIZE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   VARSOR       |<-- | BLOC DES TABLEAUX CONTENANT LES VARIABLES
C |   CLAND        |<-- | VARIABLES CLANDESTI-NES
C |   NPRE         | -->| NUMERO DE CANAL DU FICHIER
C |   STD          | -->| FORMAT DU FICHIER
C |   HIST         | -->| TABLEAU DE VALEURS MISES DANS L'ENREGISTREMENT
C |                |    | DU TEMPS.
C |   NHIST        | -->| NOMBRE DE VALEURS DANS LE TABLEAU HIST.
C |   NPOIN        | -->| NOMBRE DE POINTS DANS LE MAILLAGE
C |   AT           | -->| TEMPS
C |   TEXTPR       | -->| NOMS ET UNITES DES VARIABLES DANS LE LOGICIEL
C |   TEXTLU       | -->| NOMS ET UNITES DES VARIABLES LUES DANS LE
C |                |    | FICHIER ET VARIABLES CLANDEST.
C |   NVAR         | -->| NOMBRE DE VARIABLES NORMALES
C |   VARCLA       | -->| TABLEAU OU L'ON RANGE LES VARIABLES
C |                |    | CLANDESTINES.
C |   NVARCL       | -->| NOMBRE DE VARIABLES CLANDESTINES.
C |                |    | NVAR + NVARCL FERA LE NOMBRE TOTAL.
C |   TROUVE       |<-- | INDIQUE (TROUVE(K)=1) LES VARIABLES TROUVEES
C |                |    | DANS LE FICHIER.
C |                |    | DE K =  1 A MAXVAR VARIABLES NORMALES
C |                |    | DE K = MAXVAR+1 A MAXVAR+10 VARIABLES CLANDESTINES.
C |   ALIRE        | -->| VARIABLES QU'IL FAUT LIRE (POUR LES AUTRES ON
C |                |    | SAUTE L'ENREGISTREMENT CORRESPONDANT)
C |                |    | LES VARIABLES CLANDESTINES SONT LUES
C |                |    | SYSTEMATIQUEMENT.
C |   W            | -->| TABLEAU DE TRAVAIL REEL, DE TAILLE NPOIN.
C |   LISTIN       | -->| SI OUI, IMPRESSION D'INFORMATIONS SUR LISTING
C |   MAXVAR       | -->| NOMBRE MAXIMUM DE VARIABLES
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : LIT
C
C***********************************************************************
C
      USE BIEF_DEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: VARSOR,CLAND 
      INTEGER, INTENT(IN)             :: NPRE,NHIST,NPOIN,MAXVAR,NVARCL
      INTEGER, INTENT(IN)             :: NVAR,ALIRE(MAXVAR)
      INTEGER, INTENT(OUT)            :: TROUVE(MAXVAR)
      INTEGER, INTENT(IN), OPTIONAL   :: NPOIN_PREV,NPLAN_PREV
      CHARACTER(LEN=*)                :: STD
      CHARACTER(LEN=32)               :: TEXTPR(MAXVAR),TEXTLU(MAXVAR)
      CHARACTER(LEN=32)               :: VARCLA(NVARCL)
      DOUBLE PRECISION, INTENT(INOUT) :: HIST(*) 
      DOUBLE PRECISION, INTENT(OUT)   :: AT
      REAL                            :: W(NPOIN)
      LOGICAL, INTENT(IN)             :: LISTIN 
C
      DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: WD(*)    
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER ISTAT,L,IBID(1),NVAL,K,IHIST,NPLAN,NPOIN2,IPLAN,IP1,IP2,I     
      DOUBLE PRECISION AAT(20),TETA,ARG   
      CHARACTER(LEN=1) CBID      
      LOGICAL OK,INTERP
C
C-----------------------------------------------------------------------
C
C     INTERPOLATION ?
C
      INTERP=.FALSE.
      IF(PRESENT(NPOIN_PREV).AND.PRESENT(NPLAN_PREV)) THEN
        IF(NPOIN_PREV.NE.NPOIN) THEN
          INTERP=.TRUE.
          NPOIN2=NPOIN_PREV/NPLAN_PREV
          NPLAN=NPOIN/NPOIN2
          IF(.NOT.PRESENT(WD)) THEN
            CALL PLANTE(1)
            STOP 'WD NOT PRESENT IN LITENR'
          ENDIF
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C
      DO L=1,MAXVAR
        TROUVE(L) = 0
      ENDDO
C
C-----------------------------------------------------------------------
C
C  LECTURE DE L'ENREGISTREMENT DU TEMPS
C
      NVAL = 1 + NHIST
      CALL LIT(AAT ,W,IBID,CBID,NVAL,'R4',NPRE,STD,ISTAT)
      AT = AAT(1)
C
C  RECUPERATION DE VALEURS ECRITES DANS L'ENREGISTREMENT DU TEMPS
C
      IF(NHIST.NE.0) THEN
        DO 139 IHIST = 1, NHIST
          HIST(IHIST) = AAT(1+IHIST)
139     CONTINUE
      ENDIF
C
      IF(LISTIN.AND.LNG.EQ.1) WRITE(LU,140) AT
      IF(LISTIN.AND.LNG.EQ.2) WRITE(LU,141) AT
140   FORMAT(//,1X,'TEMPS DE L''ENREGISTREMENT : ',G16.7,' S')
141   FORMAT(//,1X,'TIME OF RECORD: ',G16.7,' S')
C
C-----------------------------------------------------------------------
C
C  LECTURE DES VARIABLES DEMANDEES SI ON LES TROUVE
C
C     CORRECTION JMH 16/12/03: NVARCL IS TAKEN INTO ACCOUNT HERE IN NVAR
C                              AS AN OUPUT OF SKIPGEO
C     DO 80 K=1,NVAR+NVARCL
      DO 80 K=1,NVAR
C
        OK = .FALSE.
C
      DO 81 L=1,MAXVAR
C
      IF (TEXTLU(K)(1:32).EQ.TEXTPR(L)(1:32) ) THEN
C
        OK = .TRUE.
C
        IF(ALIRE(L).EQ.1) THEN
C
          IF(.NOT.INTERP) THEN
C
            CALL LIT(VARSOR%ADR(L)%P%R,
     *               W,IBID,CBID,NPOIN,'R4',NPRE,STD,ISTAT)
C
          ELSE
            CALL LIT(WD,W,IBID,CBID,NPOIN_PREV,'R4',NPRE,STD,ISTAT)
C           COPYING BOTTOM AND FREE SURFACE
            CALL OV('X=Y     ',VARSOR%ADR(L)%P%R,WD,WD,0.D0,NPOIN2)
            CALL OV('X=Y     ',VARSOR%ADR(L)%P%R(NPOIN-NPOIN2+1:NPOIN),
     *                         WD(NPOIN_PREV-NPOIN2+1:NPOIN_PREV),
     *                         WD,0.D0,NPOIN2)
C           INTERPOLATING OTHER PLANES
            IF(NPLAN.GT.2) THEN
              DO IPLAN=2,NPLAN-1
                ARG=(NPLAN_PREV-1)*FLOAT(IPLAN-1)/FLOAT(NPLAN-1)
                TETA=ARG-INT(ARG)
C               IP1 : LOWER PLANE NUMBER - 1
                IP1=INT(ARG)
C               IP2 : UPPER PLANE NUMBER - 1
                IP2=IP1+1
                DO I=1,NPOIN2
                  VARSOR%ADR(L)%P%R(I+NPOIN2*(IPLAN-1))=
     *            TETA *WD(I+NPOIN2*IP2)+(1.D0-TETA)*WD(I+NPOIN2*IP1)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
C
        ELSE
C
          IF(LISTIN.AND.LNG.EQ.1) WRITE(LU,75) TEXTLU(K)
          IF(LISTIN.AND.LNG.EQ.2) WRITE(LU,76) TEXTLU(K)
75        FORMAT(/,1X,'LA VARIABLE : ',A32,/,1X,
     *                'EST DANS LE FICHIER MAIS ELLE N''EST PAS LUE')
76        FORMAT(/,1X,'VARIABLE : ',A32,/,1X,
     *                'IS IN THE FILE BUT WILL NOT BE READ')
          CALL LIT(AAT,W,IBID,CBID ,2,'R4',NPRE,STD,ISTAT)
C
        ENDIF
C
        TROUVE(L)=1
C
      ENDIF
C
81    CONTINUE
C
C     THIS SHOULD NEVER HAPPEN NOW ?????
C
      IF(NVARCL.NE.0) THEN
C
      DO 82 L=1,NVARCL
C
      IF(TEXTLU(K)(1:32).EQ.VARCLA(L)(1:32) ) THEN
        OK = .TRUE.
        CALL LIT(CLAND%ADR(L)%P%R,
     *           W,IBID,CBID,NPOIN,'R4',NPRE,STD,ISTAT)
        TROUVE(MAXVAR+L)=1
      ENDIF
C
82    CONTINUE
      ENDIF
C
        IF(.NOT.OK) THEN
          IF(LISTIN.AND.LNG.EQ.1) WRITE(LU,77) TEXTLU(K)
          IF(LISTIN.AND.LNG.EQ.2) WRITE(LU,78) TEXTLU(K)
77        FORMAT(/,1X,'LA VARIABLE : ',A32,/,1X,
     *                'EST INCONNUE, ELLE NE SERA PAS CONSERVEE')
78        FORMAT(/,1X,'VARIABLE : ',A32,/,1X,
     *                'UNKNOWN, IT WILL NOT BE KEPT')
          CALL LIT(AAT,W,IBID,CBID ,2,'R4',NPRE,STD,ISTAT)
        ENDIF
C
80    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
