C                       *****************
                        SUBROUTINE INITAB
C                       *****************
C
C 
C$DC$ : Ajout Arg NPOIN2 pour dimensionnement tableaux
C
c    * (IBOR1,IFABOR1)
     * (IBOR1,IFABOR1,NELEM2_DIM)
C
C***********************************************************************
C  TOMAWAC VERSION 1.2    23/05/96        F.MARCOS     (LNH) 30 87 72 66
C***********************************************************************
C
C   FONCTION : CE SOUS-PROGRAMME INITIALISE DES TABLEAUX UTILES
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : COW
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TOMAWAC
C
      IMPLICIT NONE
C
      DOUBLE PRECISION DEGRAD,Z,C
C
C$DC$ : NPOIN2 -> NPOIN2_DIM pour dimensionnement tableaux
C
c     INTEGER IBOR1(NELEM2,7),IFABOR1(NELEM2,3)
      INTEGER NELEM2_DIM
      INTEGER IBOR1(NELEM2_DIM,7),IFABOR1(NELEM2_DIM,3)

      INTEGER          IPLAN, IPOIN, IELEM2, IFREQ
      DOUBLE PRECISION AUXI
C
C-----------------------------------------------------------------------
C
      DO IPLAN = 1,NPLAN
         COSTET(IPLAN) = COS(TETA(IPLAN))
         SINTET(IPLAN) = SIN(TETA(IPLAN))
         IF (ABS(COSTET(IPLAN)).LT.1.D-10) COSTET(IPLAN)=0.D0
         IF (ABS(SINTET(IPLAN)).LT.1.D-10) SINTET(IPLAN)=0.D0
         IF (IPLAN.LT.NPLAN) THEN
            ETAP1(IPLAN)=IPLAN+1
         ELSE
            ETAP1(IPLAN)=1
         ENDIF
      ENDDO
C
       AUXI=(RAISF-1.D0)/2.D0
       DFREQ(1)=AUXI*FREQ(1)
       DFREQ(NF)=AUXI*FREQ(NF-1)
       DO IFREQ = 2,NF-1
         DFREQ(IFREQ) = AUXI*(FREQ(IFREQ)+FREQ(IFREQ-1))
         DO IPOIN=1,NPOIN2
           B(IPOIN+(IFREQ-1)*NPOIN2)=0.D0
         ENDDO
       ENDDO
C
      IF (SPHE) THEN
         DEGRAD=1.745329252D-2
         DO 30 IPOIN=1,NPOIN2
           COSF(IPOIN)=COS(Y(IPOIN)*DEGRAD)
           TGF(IPOIN)=TAN(Y(IPOIN)*DEGRAD)
30       CONTINUE
      ENDIF
C
      DO 40 IELEM2=1,NELEM2
         IBOR1(IELEM2,1)=IFABOR1(IELEM2,1)
         IBOR1(IELEM2,2)=IFABOR1(IELEM2,2)
         IBOR1(IELEM2,3)=IFABOR1(IELEM2,3)
         IBOR1(IELEM2,4)=1
         IBOR1(IELEM2,5)=1
         IBOR1(IELEM2,6)=1
         IBOR1(IELEM2,7)=1
40    CONTINUE
C
C INITIALISATION DES GRADIENTS DE DEPTH, U ET V
C
C W1 ( ex MASKEL) EST MIS A 1 POUR GRADF
C
      CALL OV ( 'X=C     ' , SW1%R, ST1%R, ST2%R,
     *                       1.D0 , NELEM2 )
C
      IF (.NOT.PROINF)
     *CALL VECTOR(ST1,'=','GRADF          X',IELM2,1.D0,SDEPTH,
     * ST0,ST0,ST0,ST0,ST0,MESH,.FALSE.,SW1)
C
      IF (COURAN) THEN
      CALL VECTOR(ST2,'=','GRADF          X',IELM2,1.D0,SUC,
     * ST0,ST0,ST0,ST0,ST0,MESH,.FALSE.,SW1)
C
      CALL VECTOR(ST3,'=','GRADF          X',IELM2,1.D0,SVC,
     * ST0,ST0,ST0,ST0,ST0,MESH,.FALSE.,SW1)
      ENDIF
C
      CALL VECTOR(ST4,'=','GRADF          X',IELM2,1.D0,MESH%X,
     * ST0,ST0,ST0,ST0,ST0,MESH,.FALSE.,SW1)
C
!BD_INCKA modif //
       IF(NCSIZE.GT.1) THEN
          IF (.NOT.PROINF) CALL PARCOM(ST1,2,MESH)
          CALL PARCOM(ST4,2,MESH) 
          IF (COURAN) THEN
            CALL PARCOM(ST2,2,MESH) 
            CALL PARCOM(ST3,2,MESH) 
          ENDIF       
       ENDIF
!BD_INCKA fin modif //
      IF (.NOT.PROINF)
     * CALL OV('X=Y/Z   ',SDZX%R,ST1%R,ST4%R,C,NPOIN2)
      IF (COURAN) THEN
       CALL OV('X=Y/Z   ',SDUX%R,ST2%R,ST4%R,C,NPOIN2)
       CALL OV('X=Y/Z   ',SDVX%R,ST3%R,ST4%R,C,NPOIN2)
      ENDIF
C
      IF (.NOT.PROINF)
     * CALL VECTOR(ST1,'=','GRADF          Y',IELM2,1.D0,SDEPTH,
     *  ST0,ST0,ST0,ST0,ST0,MESH,.FALSE.,SW1)
C
      IF (COURAN) THEN
       CALL VECTOR(ST2,'=','GRADF          Y',IELM2,1.D0,SUC,
     *  ST0,ST0,ST0,ST0,ST0,MESH,.FALSE.,SW1)
C
      CALL VECTOR(ST3,'=','GRADF          Y',IELM2,1.D0,SVC,
     * ST0,ST0,ST0,ST0,ST0,MESH,.FALSE.,SW1)
      ENDIF
C
      CALL VECTOR(ST4,'=','GRADF          Y',IELM2,1.D0,MESH%Y,
     * ST0,ST0,ST0,ST0,ST0,MESH,.FALSE.,SW1)
C
!BD_INCKA modif //
       IF(NCSIZE.GT.1) THEN
          IF (.NOT.PROINF) CALL PARCOM(ST1,2,MESH)
          CALL PARCOM(ST4,2,MESH)   
          IF (COURAN) THEN
            CALL PARCOM(ST2,2,MESH) 
            CALL PARCOM(ST3,2,MESH) 
          ENDIF        
       ENDIF
!BD_INCKA fin modif // 
      IF (.NOT.PROINF)
     * CALL OV('X=Y/Z   ',SDZY%R,ST1%R,ST4%R,C,NPOIN2)
      IF (COURAN) THEN
       CALL OV('X=Y/Z   ',SDUY%R,ST2%R,ST4%R,C,NPOIN2)
       CALL OV('X=Y/Z   ',SDVY%R,ST3%R,ST4%R,C,NPOIN2)
      ENDIF
C
C-----------------------------------------------------------------------
      RETURN
      END
