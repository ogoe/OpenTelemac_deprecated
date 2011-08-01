C                        *************************
                         SUBROUTINE HOMERE_ADJ_T2D
C                        *************************
C
C***********************************************************************
C TELEMAC2D VERSION 6.0                 19/09/00      A LEOPARDI (UNINA)
C***********************************************************************
C
C    FONCTIONS:
C    ==========
C
C 1) ACQUISITION DE TOUTES LES DONNEES NECESSAIRES
C    AU CALCUL DES POINTEURS: FICHIER CAS + PARTIELLEMENT LA GEOMETRIE
C
C 2) LOOP OF CALIBRATION
C
C
C-----------------------------------------------------------------------
C
C APPELE PAR :              HOMERE_TELEMAC2D
C
C SOUS-PROGRAMMES APPELES : LECDON , POINT , TELEMAC2D , COST_FUNCTION
C                           MAXIM ,
C                           GRCOUT , METGRA , INTERPOL , NEWSTR ,
C                           CLASS , SPLITINZ , ASSIGNSTR 
C
C**********************************************************************
C
C-----------------------------------------------------------------------
C 1:  DECLARATIONS
C-----------------------------------------------------------------------
C
      USE BIEF
      USE DECLARATIONS_TELEMAC2D
      USE INTERFACE_TELEMAC2D
C
      IMPLICIT NONE
      INTEGER     LNG,LU
      COMMON/INFO/LNG,LU
C
C     TYPE INTEGER:
C
      INTEGER I,NPARAM,ITER
      INTEGER IH,IU,IV,NVAR
      INTEGER TROUVE(MAXVAR+10)
      INTEGER NLAGR,ILAGR,NPOINRES
C     
C     TYPE REAL:
C     
      DOUBLE PRECISION ROX,JCOUT,JR,JCOUTN
      DOUBLE PRECISION R02,R03
      DOUBLE PRECISION C
      DOUBLE PRECISION HIST(1)
      DOUBLE PRECISION JSTEP0,JCOUT1,JCOUT2,JCOUT3
      DOUBLE PRECISION ERRH,ERRU,ERRV,AT1
C
C     TYPE LOGICAL
C
      LOGICAL RSTART
C
      CHARACTER(LEN=72) :: TITFIC
C      
C-----------------------------------------------------------------------
C  VARIABLES A LIRE :
C  0 : ON NE LIT PAS    1 : ON LIT  (VOIR SS-PG NOMVAR)
C
      INTEGER ALIRRES(MAXVAR)
C
C     ALIRRES IS TO READ U,V,H IN TELEMAC RESULTS
C
      DATA ALIRRES /1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
C
C-----------------------------------------------------------------------
C 2:  INIT PIT
C-----------------------------------------------------------------------
C
C     ENTETE SUR LISTING
C
      WRITE(LU,92)     
92    FORMAT(/////,
     *14X,'   AAAAA    DDDD       JJ        TTTTTTT  22222  DDDD ',/,
     *14X,'   A   A    D   D      JJ           T         2  D   D',/,
     *14X,'   AAAAA    D   D      JJ   IIII    T     22222  D   D',/,
     *14X,'   A   A    D   D      JJ           T     2      D   D',/,
     *14X,'   A   A    DDDD    JJJJJ           T     22222  DDDD ',/,
     *14X,'                                                      ',/,
     *14X,'               VERSION 6.0   FORTRAN 90               ',/,
     *14X,/////)
C
C ALLOCATION DES VECTEURS, MATRICES ET BLOCS
C   
      CALL POINT_ADJ_T2D
C
C     SERA REFAIT DANS TELEMAC2D, MAIS SEMBLE NECESSAIRE ICI
C     POUR AVOIR NZONE
C
      IF(DEFZON) CALL DEF_ZONES
      IF(LNG.EQ.1) WRITE(LU,*) 'NOMBRE DE ZONES : ',NZONE
      IF(LNG.EQ.2) WRITE(LU,*) 'NUMBER OF ZONES: ',NZONE
      IF(NZONE.GT.NPOIN) THEN
        IF(LNG.EQ.2) WRITE(LU,*) 'ERREUR: PLUS DE ZONES QUE DE POINTS'
        IF(LNG.EQ.2) WRITE(LU,*) 'ERROR: MORE ZONES THAN POINTS'
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C 
C     NUMBER OF PARAMETERS TO BE ESTIMATED
C
C     A MODIFIER SI ON ESTIME FROTTEMENT + AUTRE CHOSE
      NPARAM = NPOIN
      IF(NZONE.GT.0) NPARAM = NZONE
C
C-----------------------------------------------------------------------
C    
C     INITIALISATIONS ET OPTIONS POUR L'IDENTIFICATION DE PARAMETRES
C
C     INIT SORLEOA AND SORIMPA
      DO I=1,MAXVAR
         SORLEOA(I)=.FALSE.
         SORIMPA(I)=.FALSE.
      ENDDO
C
      IF(INCLU2(ESTIME,'DEBUG')) THEN
C        CV1, CV2, CV3 (seconds membres des systemes adjoints)
         SORLEOA(20)=.TRUE.
         SORLEOA(21)=.TRUE.
         SORLEOA(22)=.TRUE.
C        SORTIE DES VARIABLES ADJOINTES PP, QQ, RR
         SORLEOA(23)=.TRUE.
         SORLEOA(24)=.TRUE.
         SORLEOA(25)=.TRUE.
      ELSE
C        OUTPUT OF BOTTOM TOPOGRAPHY (6) AND FRICTION (19)
         SORLEOA(6)=.TRUE.
         SORLEOA(19)=.TRUE.
      ENDIF
C     
      CALL OS('X=C     ',PRIVE,PRIVE,PRIVE,0.D0)
C
C     OPTION DE DESCENTE : OPTID = 1 METHODE DE GRADIENT SIMPLE
C                          OPTID = 2 METHODE DE GRADIENT CONJUGUE
C                          OPTID = 3 INTERPOLATION DE LAGRANGE POUR RHO
C
      WRITE(LU,*) 'OPTID = ',OPTID
      IF(OPTID.EQ.0) THEN
         WRITE(LU,*) 'PLAN D''EXPERIENCE'
      ELSEIF(OPTID.EQ.1) THEN 
         WRITE(LU,*) 'GRADIENT METHOD'
      ELSEIF(OPTID.EQ.2) THEN
         WRITE(LU,*) 'CONJUGATE GRADIENT METHOD'
      ELSEIF(OPTID.EQ.3) THEN
         WRITE(LU,*) 'LAGRANGE INTERPOLATION'
      ELSE
         WRITE(LU,*) 'WRONG OPTION FOR COMPUTATION OF RHO'
         CALL PLANTE(1)
         STOP  
      ENDIF
C         
      IF(OPTID.EQ.3) THEN
        NLAGR=3
      ELSE
        NLAGR=1
      ENDIF
C
C======================================================================
C     INIT ADJOINT
C======================================================================
C     
      JSTEP0=1.D0
      JR=0.D0
      JCOUTN = 8000000.D0
      NITERA = 0
C
      JCOUT=0.D0     
      RSTART=.TRUE.
C
      IF(OPTID.NE.0) THEN 
        CALL OV('X=C     ',GRADJ%R  , GRADJ%R  , GRADJ%R  ,0.D0,NPARAM)
        CALL OV('X=C     ',GRADJN%R , GRADJN%R , GRADJN%R ,0.D0,NPARAM)
      ENDIF
C      
C======================================================================
C      MANAGEMENT OF FILES
C======================================================================
C
C HEADER FOR ASCII OUTPUT (FORMAT SCOPT)
C
       IF(T2D_FILES(T2DRFO)%NAME.NE.' ') THEN
         WRITE(T2D_FILES(T2DRFO)%LU,300) TITCAS
300      FORMAT('''',A,'''',1I2)
         WRITE(T2D_FILES(T2DRFO)%LU,300) 'PARAMETER ESTIMATION'
         WRITE(T2D_FILES(T2DRFO)%LU,300) 'IDENTIFICATION OF FRICTION'
         WRITE(T2D_FILES(T2DRFO)%LU,300) 'ITERATION'
         WRITE(T2D_FILES(T2DRFO)%LU,300) 'COST'
         WRITE(T2D_FILES(T2DRFO)%LU,300) 'ERROR ON H (M)'
         WRITE(T2D_FILES(T2DRFO)%LU,300) 'ERROR ON U (M/S)'
         WRITE(T2D_FILES(T2DRFO)%LU,300) 'ERROR ON V (M/S)'
         IF(DEFZON) THEN
           DO I=1,NZONE
             WRITE(T2D_FILES(T2DRFO)%LU,300) 'FRICTION ZONE ',I
           ENDDO
         ELSE
           WRITE(T2D_FILES(T2DRFO)%LU,300) 'FRICTION POINT 1'
         ENDIF
       ENDIF
C
C ****************************************************************************
C LOOP OF CALIBRATION
C ****************************************************************************
C
C     INITIAL STRICKLERS ARE SET IN STRCHE
C     (IF NOT MODIFIED ALL VALUES AT FFON)
C
C     LINES BELOW WILL BE REPEATED WHEN NITERA = 1 IN TELEMAC2D (SEE CALL TO FONSTR)
C     THEY ARE KEPT HERE TO INITIALISE SETSTR2
C
      CALL OS('X=C     ',X=CHESTR,C=FFON)
      CALL STRCHE
      CALL INITSTR(CHESTR,SETSTR,ZONE%I,NZONE,NPOIN,T1)
      CALL ASSIGNSTR(CHESTR,SETSTR,ZONE%I,NZONE,NPOIN)
C
95    CONTINUE           
C
C
C------------------------------------------------------------------------------
C
C     PLAN D'EXPERIENCE : LECTURE DES COEFFICIENTS
C
C     LIGNE DE COMMENTAIRES / SKIPPING A LINE OF COMMENTS
      IF(OPTID.EQ.0) READ(T2D_FILES(T2DFO1)%LU,*)
500   CONTINUE
C    
C------------------------------------------------------------------------------    
C
C     NITERA : NUMBER OF ITERATIONS
C
      NITERA = NITERA + 1
C
      WRITE(LU,*) ' '
      WRITE(LU,*) '------------------'
      WRITE(LU,*) 'ITERATION : ',NITERA
      WRITE(LU,*) '------------------'
      WRITE(LU,*) ' '
C    
C TO PRESERVE THE VALUES OF STRICKLERS'.
C
      CALL OS( 'X=Y     ' ,X=SETSTR2 , Y=SETSTR )
C
C READING THE NEW STRICKLERS IN A FILE IF OPTID=0
C
      IF(OPTID.EQ.0) THEN
        READ(T2D_FILES(T2DFO1)%LU,*,END=999) ITER,
     *       (SETSTR%R(I),I=1,NPARAM)
        IF(ITER.NE.NITERA) THEN
          IF(LNG.EQ.1) WRITE(LU,*)'PB. DANS LE PLAN D''EXPERIENCE ITER=',ITER
          IF(LNG.EQ.2) WRITE(LU,*)'PB. IN LIST OF TESTS AT ITER=',ITER
          STOP
        ENDIF
        CALL ASSIGNSTR(CHESTR,SETSTR,ZONE%I,NZONE,NPOIN)
      ENDIF
C
C *** LOOP FOR LAGRANGE INTERPOLATION ***
C
      DO ILAGR=1,NLAGR
C      
         IF(OPTID.EQ.3) THEN      
           WRITE(LU,*) ' '
           WRITE(LU,*) '------------------'
           IF(LNG.EQ.1) WRITE(LU,*) 'SOUS-ITERATION : ',ILAGR
           IF(LNG.EQ.2) WRITE(LU,*) 'SUB-ITERATION : ',ILAGR
           WRITE(LU,*) '------------------'
           WRITE(LU,*) ' '
         ENDIF
C
C-----------------------------------------------------------------------
C
         IF(LNG.EQ.1) WRITE(LU,100)
         IF(LNG.EQ.2) WRITE(LU,101)
         WRITE(LU,102)    
100      FORMAT(/////,1X,'LISTING DE TELEMAC-2D ',78('-'))
101      FORMAT(/////,1X,'LISTING OF TELEMAC-2D ',78('-'))
102      FORMAT(/////,
     *14X,'TTTTT  EEEEE  L      EEEEE  M   M  AAAAA  CCCCC',/,
     *14X,'  T    E      L      E      MM MM  A   A  C    ',/,
     *14X,'  T    EEE    L      EEE    M M M  AAAAA  C    ',/,
     *14X,'  T    E      L      E      M   M  A   A  C    ',/,
     *14X,'  T    EEEEE  LLLLL  EEEEE  M   M  A   A  CCCCC',/,
     *14X,'                                               ',/,
     *14X,'        2D    VERSION 6.0     FORTRAN 90       ',/,
     *14X,'                                               ',/,
     *14X,'DIRECT MODE DIRECT MODE DIRECT MODE DIRECT MODE',/,
     *14X,/////)
C
      ADJO=.FALSE.
      CALL TELEMAC2D(PASS= -1,ATDEP=0.D0,NITER=0,CODE='       ')
C
C  /* TEMPORAL LOOP (COMPUTATION OF COST FUNCTION) */
C
      JCOUT=0.D0
C       
C SKIP GEOMETRY
C
      REWIND T2D_FILES(T2DRES)%LU
C
      CALL SKIPGEO(T2D_FILES(T2DRES)%LU,TITFIC,NPOINRES,NVARRES,TEXRES)
      IF(NPOINRES.NE.NPOIN) THEN
        WRITE(LU,*) 'ERROR: NPOINRES DIFFERENT FROM NPOIN'
        WRITE(LU,*) 'NPOINRES = ',NPOINRES
        WRITE(LU,*) 'NPOIN    = ',NPOIN
        CALL PLANTE(1)
        STOP
      ENDIF
C       
C SKIP INITIAL CONDITION
C
      IF(OUTINI) THEN
       CALL LITENR(VARSOR,VARCL,T2D_FILES(T2DRES)%LU,'STD',HIST,0,
     *             NPOIN,AT1,TEXTE,TEXRES,NVARRES,VARCLA,0,TROUVE,
     *             ALIRRES,W,.FALSE.,MAXVAR)  
      ENDIF
C
      ERRH=0.D0
      ERRU=0.D0
      ERRV=0.D0
      IH=0
      IU=0
      IV=0
C      
      DO LT=1,NIT
C
C     IN STEADY STATE ONLY THE LAST TIME-STEP IS CONSIDERED
C
      IF(     INCLU2(ESTIME,'PERMANENT')
     *    .OR.INCLU2(ESTIME,'STEADY'   )  ) THEN
C
        IF(LT.EQ.1) THEN
          REWIND T2D_FILES(T2DRES)%LU
          CALL BIEF_SUITE(VARSOR,VARCL,ITER,T2D_FILES(T2DRES)%LU,
     *                    T2D_FILES(T2DRES)%FMT,
     *                    HIST,0,NPOIN,AT1,TEXTE,VARCLA,
     *                    NVARCL,TROUVE,ALIRRES,LISTIN,.TRUE.,MAXVAR)
C         GETTING MEASUREMENTS AND WEIGHTS AT THE SAME TIME
C         HERE ALPHA1, ALPHA2 AND ALPHA3 ARE ALSO SET.
C         ITER OF LAST RECORD GIVEN BY THE CALL TO SUITE 
          CALL MESURES(ITER,AT1)
        ENDIF
C
      ELSE
C
C  READING TELEMAC2D RESULTS (RESULTS FILE - UNIT NRES)
C
        ITER=LT
        IF(OUTINI) ITER=ITER+1
C
        CALL LITENR(VARSOR,VARCL,T2D_FILES(T2DRES)%LU,'STD',HIST,0,
     *              NPOIN,AT1,TEXTE,TEXRES,NVARRES,VARCLA,0,TROUVE,
     *             ALIRRES,W,.FALSE.,MAXVAR)
C
C       GETTING MEASUREMENTS AND WEIGHTS AT THE SAME TIME
C
C       HERE ALPHA1, ALPHA2 AND ALPHA3 ARE ALSO SET. 
        CALL MESURES(ITER,AT1)         
C
      ENDIF
C   
C       COMPUTATION OF COST FUNCTION :   
C
        CALL COST_FUNCTION(JCOUT,OPTCOST,'FCT')
C     
C       COMPUTATION OF DIFFERENCES BETWEEN MEASUREMENTS AND COMPUTED VALUES 
C
        CALL ERRMAX(H,HD,C,I)
        IF(ERRH.LT.C) THEN
          ERRH=C
          IH=I
        ENDIF
        CALL ERRMAX(U,UD,C,I)
        IF(ERRU.LT.C) THEN
          ERRU=C
          IU=I
        ENDIF
        CALL ERRMAX(V,VD,C,I)
        IF(ERRV.LT.C) THEN
          ERRV=C
          IU=I
        ENDIF
C
C END OF TEMPORAL LOOP (COST FUNCTION COMPUTED)
C
      ENDDO
C
      IF(NITERA.EQ.1.AND.ILAGR.EQ.1) THEN
        IF(JCOUT.GT.0.D0) THEN
          JSTEP0=JCOUT
        ELSE
          JSTEP0=1.D0
        ENDIF
      ENDIF
C 
      IF(LNG.EQ.1) THEN
C
      WRITE(LU,*)'FONCTION COUT =',JCOUT,' VALEUR INITIALE:',JSTEP0
      WRITE(LU,*)'ERREUR MAXIMUM SUR H =',ERRH
      WRITE(LU,*)'ERREUR MAXIMUM SUR U =',ERRU
      WRITE(LU,*)'ERREUR MAXIMUM SUR V =',ERRV 
C 
      ELSEIF(LNG.EQ.2) THEN
C              
      WRITE(LU,*) 'COST FUNCTION =',JCOUT,' INITIAL VALUE :',JSTEP0
      WRITE(LU,*) 'MAX ERROR ON H =',ERRH
      WRITE(LU,*) 'MAX ERROR ON U =',ERRU
      WRITE(LU,*) 'MAX ERROR ON V =',ERRV
C
      ENDIF
C  
      IF(ILAGR.EQ.1) THEN
        JR = JCOUT/JSTEP0
        IF(LNG.EQ.1) WRITE(LU,*) 'FONCTION COUT RELATIVE : ',JR
        IF(LNG.EQ.2) WRITE(LU,*) 'RELATIVE COST FUNCTION: ',JR
        IF(T2D_FILES(T2DRFO)%NAME(1:1).NE.' ') THEN
          IF(DEFZON) THEN
            WRITE(T2D_FILES(T2DRFO)%LU,*) NITERA,JR,ERRH,ERRU,ERRV,
     *                    (SETSTR%R(I),I=1,NZONE)
          ELSE
            WRITE(T2D_FILES(T2DRFO)%LU,*) NITERA,JR,ERRH,ERRU,ERRV,
     *                     SETSTR%R(1)
          ENDIF
        ENDIF
      ENDIF
C
C
C
      IF(OPTID.EQ.0) GO TO 500
C
C
C
C  TEST, DECISIONAL STEP & ADJOINT SYSTEM (ONLY FOR ILAGR=1)
C
C  TEST: TWO CRITERIA
C
      IF (ILAGR.EQ.1) THEN
C       
C     DECISIONAL STEP :
C
      IF(      JR.LE.TOLEST(4).OR.
     *      (ERRH.LE.TOLEST(1).AND.
     *       ERRU.LE.TOLEST(2).AND.
     *       ERRV.LE.TOLEST(3))       ) THEN
C     
         IF(LISTIN) THEN
           IF(LNG.EQ.1) WRITE(LU,395) NITERA
           IF(LNG.EQ.2) WRITE(LU,396) NITERA
         ENDIF
C
395      FORMAT(/,1X,'------------------------------------------',/
     *           ,1X,'    SOLUTION TROUVEE EN ',1I3,' ITERATIONS',/
     *           ,1X,'------------------------------------------')     
396      FORMAT(/,1X,'-----------------------------------------',/
     *           ,1X,'    SOLUTION FOUND IN ',1I3,' ITERATIONS',/
     *           ,1X,'-----------------------------------------')
         WRITE(LU,*) 'GRADIENT OF ZONE 1 : ',GRADJ%R(1)
         WRITE(LU,*) 'STRICKLER OF POINT 10 : ',CHESTR%R(10)
         GO TO 999
C     
      ELSEIF (NITERA.GT.MAXEST) THEN
C
         IF(LNG.EQ.1) THEN        
         WRITE(LU,*) 'PAS DE CONVERGENCE EN ',NITERA,' ITERATIONS'
         WRITE(LU,*) 'STRICKLER DU POINT 10 : ',CHESTR%R(10)     
         WRITE(LU,398) MAXEST,JCOUT
398      FORMAT(1X,'SOLUTION NON TROUVEE APRES ',1I6,1X,
     *           'ITERATIONS',/,1X,
     *           'PRECISION  :',G16.7,1X,'JCOUTN/JCOUT1 :',G16.7)
         ELSEIF(LNG.EQ.2) THEN
         WRITE(LU,*) 'NO CONVERGENCE AFTER ',NITERA,' ITERATIONS'
         WRITE(LU,*) 'STRICKLER OF POINT 10 : ',CHESTR%R(10)     
         WRITE(LU,399) MAXEST,JCOUT
399      FORMAT(1X,'SOLUTION NOT FOUND AFTER ',1I6,1X,
     *           'ITERATIONS',/,1X,
     *           'PRECISION  :',G16.7,1X,'JCOUTN/JCOUT1 :',G16.7)
         ENDIF
         GO TO 999
C
      ELSEIF (JCOUT.GT.JCOUTN.AND..NOT.RSTART) THEN
C
         IF(LNG.EQ.1) THEN
           WRITE(LU,*) 'LA FONCTION COUT AUGMENTE : STOP' 
           WRITE(LU,*) 'STRICKLER DU POINT 10 : ',CHESTR%R(10)
         ELSEIF(LNG.EQ.2) THEN
           WRITE(LU,*) 'COST FUNCTION INCREASES : STOP' 
           WRITE(LU,*) 'STRICKLER OF POINT 10 : ',CHESTR%R(10)
         ENDIF
C        GO TO 999
C
      ENDIF
C
C ADJOINT SYSTEM
C
         IF(LNG.EQ.1) WRITE(LU,403)
         IF(LNG.EQ.2) WRITE(LU,404)
         WRITE(LU,405)    
403      FORMAT(/////,1X,'LISTING D" ESTIMATION',82(1H-))
404      FORMAT(/////,1X,'LISTING OF ESTIMATION',82(1H-))
405      FORMAT(/////,
     *14X,'TTTTT  EEEEE  L      EEEEE  M   M  AAAAA  CCCCC',/,
     *14X,'  T    E      L      E      MM MM  A   A  C    ',/,
     *14X,'  T    EEE    L      EEE    M M M  AAAAA  C    ',/,
     *14X,'  T    E      L      E      M   M  A   A  C    ',/,
     *14X,'  T    EEEEE  LLLLL  EEEEE  M   M  A   A  CCCCC',/,
     *14X,'                                               ',/,
     *14X,'        2D    VERSION 6.0     FORTRAN 90       ',/,
     *14X,'                                               ',/,
     *14X,' ADJOINT MODE ADJOINT MODE ADJOINT MODE ADJOINT',/,
     *14X,/////)
C
C        INITIALISING THE GRADIENT WHICH WILL BE COMPUTED
C        BY PROPAG_ADJ
C
         CALL OV('X=C     ',GRADJ%R,GRADJ%R,GRADJ%R,0.D0,NPARAM)
C
C        SERIES OF ADJOINT SYSTEMS
C
         ADJO=.TRUE.
         CALL TELEMAC2D(PASS= -1,ATDEP=0.D0,
     *                  NITER=0,CODE='       ')
C
         IF(NZONE.GT.0) THEN 
           DO I=1,NZONE
             WRITE(LU,*) 'GRADJ(',I,')= ',GRADJ%R(I)
C            STOP 'PROVISOIRE'
           ENDDO
         ENDIF
C
C END OF: IF (ILAGR.EQ.1)
      ENDIF
C
C
      IF(ILAGR.EQ.1) THEN
C
C       GRADIENT METHOD: COMPUTATION OF RHO AND DIRECTION
C       JCOUT1 IS JCOUT FOR RHO=0
        JCOUT1=JCOUT
        CALL METGRA(ROX,ESTIME,GRADJ,GRADJN,JCOUT1,DESC,NPARAM,OPTID,
     *              RSTART,R02,R03)
        IF(OPTID.EQ.3) THEN
          CALL NEWSTR(SETSTR,SETSTR2,DESC,ROX,RSTART,NPARAM,
     *                ESTIME,KFROT)
          CALL ASSIGNSTR(CHESTR ,SETSTR,ZONE%I,NZONE,NPOIN)
        ENDIF
      ELSEIF (ILAGR.EQ.2) THEN
C JCOUT2 IS JCOUT FOR RHO=ROX
         JCOUT2=JCOUT
         CALL NEWSTR(SETSTR,SETSTR2,DESC,0.5D0*ROX,RSTART,NPARAM,
     *               ESTIME,KFROT)
         CALL ASSIGNSTR(CHESTR ,SETSTR,ZONE%I,NZONE,NPOIN)
      ELSEIF (ILAGR.EQ.3) THEN         
C JCOUT3 IS JCOUT FOR RHO=1/2 ROX
         JCOUT3=JCOUT
C
      ENDIF
C
C *** END OF LAGRANGIAN LOOP ***
      ENDDO
C
C  CASE OF A NEW ITERATION:
C
C  COMPUTATION OF NEW VALUE OF ROX (IF LAGRANGE)
      IF(OPTID.EQ.3) CALL INTERPOL(ROX,R02,R03,JCOUT1,JCOUT2,JCOUT3)
C
      WRITE(LU,*) 'ITERATION = ',NITERA
      WRITE(LU,*) 'STRICKLERS = ',SETSTR%R(1)
      WRITE(LU,*) 'J         = ',JCOUT1
      WRITE(LU,*) 'JR        = ',JR
C      
C  COMPUTE THE NEW SET OF COEFFICIENTS
C
      CALL NEWSTR(SETSTR,SETSTR2,DESC,ROX,RSTART,NPARAM,ESTIME,KFROT)
      CALL ASSIGNSTR(CHESTR ,SETSTR,ZONE%I,NZONE,NPOIN)
C
C  COST FUNCTION AND GRADIENT AT PREVIOUS ITERATION
C
      JCOUTN=JCOUT1
      CALL OV( 'X=Y     ' , GRADJN%R , GRADJ%R, GRADJ%R , C , NPARAM )
C
      GOTO 95
C
C ****************************************************************************
C  END OF LOOP OF CALIBRATION
C ****************************************************************************
C
999   CONTINUE
C
C PRINT INFORMATIONS ABOUT THE LAST ITERATION
C
      WRITE(LU,*) 'ITERATION = ',NITERA
      WRITE(LU,*) 'STRICKLERS = ',SETSTR%R(1)
      WRITE(LU,*) 'J         = ',JCOUT1
      WRITE(LU,*) 'JR        = ',JR
C
C-----------------------------------------------------------------------
C
C     ECRITURE FINALE DU FICHIER DE GEOMETRIE RESULTAT QUI CONTIENDRA
C     AUSSI LE FROTTEMENT (EN MODE DEBUG, ON SORT A LA PLACE LES VARIABLES
C     ADJOINTES, VOIR DANS TELEMAC-2D)
C
      IF(.NOT.INCLU2(ESTIME,'DEBUG')) THEN
C
        CALL ECRGEO(MESH%X%R,MESH%Y%R,MESH%NPOIN,MESH%NBOR%I,
     *              T2D_FILES(T2DRBI)%LU,NVAR,TEXTE,VARCLA,NVARCL,
     *              TITCAS,SORLEOA,MAXVAR,MESH%IKLE%I,
     *              MESH%NELEM,MESH%NPTFR,3,MARDAT,MARTIM,
     *              NCSIZE,NPTIR,MESH%KNOLG%I,I3=I_ORIG,I4=J_ORIG)
        CALL BIEF_DESIMP('SERAFIN ',VARSOR,
     *                   HIST,0,NPOIN,T2D_FILES(T2DRBI)%LU,
     *                   'STD',0.D0,0,LISPRD,LEOPRD,
     *                   SORLEOA,SORIMPA,MAXVAR,TEXTE,0,0)
C
      ENDIF
C
C-----------------------------------------------------------------------
C    
      RETURN      
      END
