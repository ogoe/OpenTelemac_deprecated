C                       ***************************
                        SUBROUTINE PRERES_TELEMAC2D
C                       ***************************
C
C***********************************************************************
C TELEMAC 2D VERSION 6.0  24/11/2009  J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C     FONCTION  : PREPARATION DE VARIABLES QUI SERONT ECRITES SUR
C                 LE FICHIER DE RESULTATS OU SUR LE LISTING.
C
C-----------------------------------------------------------------------
C
C     FUNCTION  : PREPARING VARIABLES WHICH WILL BE WRITTEN IN THE
C                 RESULTS FILE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      LT        | -->| ITERATION NUMBER
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : TELMAC
C
C  SOUS-PROGRAMME APPELE : OV
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC2D
      USE INTERFACE_TELEMAC2D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C     
      LOGICAL IMP,LEO,DEJA1,DEJA2,DEJA3
C
      INTEGER LTT,N,IMAX,I
C
      DOUBLE PRECISION HHH,XMAX,NF,PI,AMP,PHA
C
      INTRINSIC MAX,SQRT
C
      DOUBLE PRECISION P_DMAX,P_DMIN
      EXTERNAL         P_DMAX,P_DMIN
C
C-----------------------------------------------------------------------
C
      DATA DEJA1/.FALSE./
      DATA DEJA2/.FALSE./
      DATA DEJA3/.FALSE./
      SAVE DEJA1,DEJA2,DEJA3,NF
C
C-----------------------------------------------------------------------
C
C     THE OUTPUT VARIABLES ARE BUILD ONLY IF NECESSARY, HENCE THE TESTS
C     BELOW, WHICH MUST BE THE SAME THAN IN BIEF_DESIMP (BIEF LIBRARY)
C
C     THIS WILL TRIGGER AN OUTPUT OF LAST TIME-STEP
C     BUT NOT WITH PARAMETER ESTIMATION (LISPRD WOULD STAY AT 1
C     FOR FURTHER COMPUTATIONS)
      IF(LT.EQ.NIT.AND.ESTIME(1:1).EQ.' ') THEN
        LISPRD=1
        LEOPRD=1
      ENDIF
C
      IMP=.FALSE.
      LEO=.FALSE.
      LTT=(LT/LISPRD)*LISPRD
      IF(LT.EQ.LTT.AND.LT.GE.PTINIL) IMP=.TRUE.
      LTT=(LT/LEOPRD)*LEOPRD
      IF(LT.EQ.LTT.AND.LT.GE.PTINIG) LEO=.TRUE.
C
      IF(LT.EQ.0) THEN
        IMP=OUTINI
        LEO=OUTINI
      ENDIF
C
C-----------------------------------------------------------------------
C
C 1)  PART WHICH MUST BE DONE EVEN IF THERE IS NO OUTPUT FOR THIS STEP
C     BUT ONLY AFTER FIRST TIME STEP FOR GRAPHIC PRINTOUTS      
C
C----------------------------------------------------------------------- 
C
      IF(LT.GE.PTINIG) THEN
C
C=======================================================================
C CALCUL DE LA COTE MAXIMUM ET TEMPS ASSOCIE
C=======================================================================
C
      IF(SORLEO(27).OR.SORIMP(27)) THEN
        IF(.NOT.DEJA1) THEN
          CALL OS('X=Y     ',X=MAXZ ,Y=ZF)
          CALL OS('X=C     ',X=TMAXZ,C=AT)
          DEJA1=.TRUE.
        ELSE
          DO N=1,NPOIN
            XMAX=H%R(N)+ZF%R(N)
C           DRY LAND EXCLUDED (TO AVOID RANDOM TIMES)
            IF(XMAX.GT.MAXZ%R(N).AND.H%R(N).GT.0.01D0) THEN
              MAXZ%R(N)=XMAX
              IF(SORLEO(28).OR.SORIMP(28)) TMAXZ%R(N)=AT
            ENDIF
          ENDDO
        ENDIF
      ENDIF
C
C=======================================================================
C CALCUL DE LA VITESSE MAXIMUM ET TEMPS ASSOCIE
C=======================================================================
C
      IF(SORLEO(29).OR.SORIMP(29)) THEN
        IF(.NOT.DEJA2) THEN
          CALL OS('X=C     ',X=MAXV ,C=0.D0)
          CALL OS('X=C     ',X=TMAXV,C=AT)
          DEJA2=.TRUE.
        ELSE
          DO N=1,NPOIN
            XMAX=SQRT(U%R(N)**2+V%R(N)**2)
C           DRY LAND EXCLUDED (TO AVOID RANDOM TIMES)
            IF(XMAX.GT.MAXV%R(N).AND.H%R(N).GT.0.01D0) THEN
              MAXV%R(N)=XMAX
              IF(SORLEO(30).OR.SORIMP(30)) TMAXV%R(N)=AT
            ENDIF
          ENDDO
        ENDIF
      ENDIF
C
C=======================================================================
C IMPRESSIONS POUR LES POINTS REMARQUABLES
C=======================================================================
C
      IF(LT.EQ.NIT.AND.NPTS.GT.0) THEN
        DO I=27,30
C         CAUTION : HERE SORLEO IS USED INSTEAD OF SORIMP
          IF(SORLEO(I)) THEN
            WRITE(LU,*) ' '
            WRITE(LU,*) ' '
            WRITE(LU,*) ' '
            WRITE(LU,*) TEXTE(I)(1:16)
            WRITE(LU,*) ' '
            DO N=1,NPTS
!             IN PARALLEL POINT NOT ALWAYS EXISTING, MAYBE ELSEWHERE
              IF(NCSIZE.GT.0) THEN
                WRITE(LU,*) NAME_PTS(N),' : ',
     *                    P_DMIN(VARSOR%ADR(I)%P%R(LIST_PTS(N)))+
     *                    P_DMAX(VARSOR%ADR(I)%P%R(LIST_PTS(N)))
              ELSE
                WRITE(LU,*) NAME_PTS(N),' : ',
     *                                    VARSOR%ADR(I)%P%R(LIST_PTS(N))
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
C     
C-----------------------------------------------------------------------
C
      ELSE
C
C     CAS OU OUTINI=.TRUE. : PRIORITE SUR PTINIG, ON MET LES VALEURS
C     DE LT=0 CAR SINON ELLES RESTENT NON INITIALISEES
       IF(SORLEO(27).OR.SORIMP(27)) CALL OS('X=Y     ',X=MAXZ ,Y=ZF)
       IF(SORLEO(28).OR.SORIMP(28)) CALL OS('X=C     ',X=TMAXZ,C=AT)
       IF(SORLEO(29).OR.SORIMP(29)) CALL OS('X=C     ',X=MAXV ,C=0.D0)
       IF(SORLEO(30).OR.SORIMP(30)) CALL OS('X=C     ',X=TMAXV,C=AT)
C
C     ENDIF DE : IF(LT.GE.PTINIG) THEN
      ENDIF
C     
C-----------------------------------------------------------------------
C
C 2)  PART WHICH MUST BE DONE ONLY IF THERE IS AN OUTPUT FOR THIS STEP
C    
C-----------------------------------------------------------------------
C
C     PAS D'IMPRESSION, PAS DE SORTIE SUR FICHIER, ON RESSORT
      IF(.NOT.(LEO.OR.IMP)) GO TO 1000
C
C
C=======================================================================
C CALCUL DE LA CELERITE (MISE DANS FU, VOIR LE BLOC VARSOR)
C=======================================================================
C
      IF((LEO.AND.SORLEO(3)).OR.(IMP.AND.SORIMP(3))) THEN
        CALL CPSTVC(ZF,FU)
        DO N=1,NPOIN
          FU%R(N) = SQRT ( GRAV * MAX(H%R(N),0.D0) )
        ENDDO
      ENDIF
C
C=======================================================================
C CALCUL DE LA SURFACE LIBRE (= H + ZF, MISE DANS FV)
C=======================================================================
C
      IF((LEO.AND.SORLEO(5)).OR.(IMP.AND.SORIMP(5))) THEN
        CALL CPSTVC(ZF,FV)
        DO N=1,NPOIN
          FV%R(N) = H%R(N)+ZF%R(N)
        ENDDO
      ENDIF
C
C=======================================================================
C CALCUL DU NOMBRE DE FROUDE
C=======================================================================
C
      IF((LEO.AND.SORLEO(7)).OR.(IMP.AND.SORIMP(7))) THEN
        CALL CPSTVC(ZF,T2)
        DO N=1,NPOIN
          HHH = MAX( H%R(N) , 1.D-8 )
          T2%R(N) = SQRT (( U%R(N)**2 + V%R(N)**2 ) / ( HHH*GRAV ))
        ENDDO
      ENDIF
C
C=======================================================================
C CALCUL DU DEBIT SCALAIRE
C=======================================================================
C
      IF((LEO.AND.SORLEO(8)).OR.(IMP.AND.SORIMP(8))) THEN
        CALL CPSTVC(ZF,T3)
        DO N=1,NPOIN
         T3%R(N) = SQRT (U%R(N)**2 + V%R(N)**2) * H%R(N)
        ENDDO
      ENDIF
C
C=======================================================================
C CALCUL DU DEBIT VECTORIEL , COMPOSANTE SUIVANT X
C=======================================================================
C
      IF((LEO.AND.SORLEO(13)).OR.(IMP.AND.SORIMP(13))) THEN
        CALL CPSTVC(ZF,T4)
        DO N=1,NPOIN
          T4%R(N)=H%R(N)*U%R(N)
        ENDDO
      ENDIF
C
C=======================================================================
C CALCUL DU DEBIT VECTORIEL , COMPOSANTE SUIVANT Y
C=======================================================================
C
      IF((LEO.AND.SORLEO(14)).OR.(IMP.AND.SORIMP(14))) THEN
        CALL CPSTVC(ZF,T5)
        DO N=1,NPOIN
          T5%R(N)=H%R(N)*V%R(N)
        ENDDO
      ENDIF
C
C=======================================================================
C CALCUL DE LA VITESSE SCALAIRE
C=======================================================================
C
      IF((LEO.AND.SORLEO(15)).OR.(IMP.AND.SORIMP(15))) THEN
        CALL OS( 'X=N(Y,Z)' , X=T6 , Y=U , Z=V )
      ENDIF
C
C=======================================================================
C CALCUL DU NOMBRE DE COURANT
C=======================================================================
C
      IF((LEO.AND.SORLEO(22)).OR.(IMP.AND.SORIMP(22))) THEN
C                             IELM
        CALL CFLPSI(T9,U,V,DT,11,MESH,MSK,MASKEL)
        CALL MAXI(XMAX,IMAX,T9%R,NPOIN)
        IF(NCSIZE.GT.1) THEN
          IF(LNG.EQ.1) WRITE(LU,78) P_DMAX(XMAX)
          IF(LNG.EQ.2) WRITE(LU,79) P_DMAX(XMAX)
        ELSE
          IF(LNG.EQ.1) WRITE(LU,78) XMAX
          IF(LNG.EQ.2) WRITE(LU,79) XMAX          
        ENDIF
78      FORMAT(1X,'PRERES : NOMBRE DE COURANT MAXIMUM :',G16.7)
79      FORMAT(1X,'PRERES: MAXIMUM COURANT NUMBER: ',G16.7)
      ENDIF
C
C=======================================================================
C CALCUL DE LA VITESSE DE FROTTEMENT
C=======================================================================
C
      IF((LEO.AND.SORLEO(31)).OR.(IMP.AND.SORIMP(31))) THEN
        CALL CPSTVC(CF,T7)
        DO N=1,NPOIN
          T7%R(N) = SQRT(0.5D0*CF%R(N)*(U%R(N)**2+V%R(N)**2))
        ENDDO
      ENDIF
C
C=======================================================================
C
1000  CONTINUE
C
C=======================================================================
C HARMONIC ANALYSIS USING LEAST MEAN ERROR SQUARE METHOD
C=======================================================================
C
      IF(NPERIAF.GT.0) CALL SPECTRE
C
C=======================================================================
C
      RETURN
      END
