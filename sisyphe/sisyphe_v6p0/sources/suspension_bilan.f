      ! *************************** !
        SUBROUTINE SUSPENSION_BILAN 
      ! *************************** !

     &(MESH,CST,HN,ZFCL_S,MASKEL, 
     & IELMT,ITRA,LT,NIT,DT,CSF,
     & MASSOU,MASED0,MSK,ENTET,MASTEN,MASTOU,MASINI,T2,
     & T3,MASFIN,MASDEPT,MASDEP,AGGLOT,
     & VOLU2D,NUMLIQ,NFRLIQ,NPTFR,FLBORTRA)

C***********************************************************************
C SISYPHE VERSION 5.8  29/10/07  J-M HERVOUET  CORRECTIONS EN PARALLELE
C SISYPHE VERSION 5.6  22/12/04  F. HUVELIN                            
C SISYPHE VERSION 5.4  mai 2003  M. GONZALES DE LINARES                
C SUBIEF  VERSION 5.1  13/12/00  C. MOULIN (LNH)        01 30 87 83 81
C
C 05/05/2008 : CALCUL DE LA MASSE EN TENANT COMPTE DU MASS-LUMPING
C 28/05/2008 : FLUX DONNES PAR FRONTIERE
C 10/06/2008 : FLUX DE TRACEUR DONNE DANS FLBORTRA (SORTI DE CVDFTR)
C              13 ARGUMENTS SUPPRIMES
C 
C**********************************************************************
C
C
C               ! =============================== !
C               ! Mass-balance for the suspension !
C               ! =============================== !
C
C
C
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT
C**********************************************************************C
C                                                                      C
C                 SSSS I   SSSS Y   Y PPPP  H   H EEEEE                C
C                S     I  S      Y Y  P   P H   H E                    C
C                 SSS  I   SSS    Y   PPPP  HHHHH EEEE                 C
C                    S I      S   Y   P     H   H E                    C
C                SSSS  I  SSSS    Y   P     H   H EEEEE                C
C                                                                      C
C----------------------------------------------------------------------C
C                             ARGUMENTS                                C
C .________________.____.______________________________________________C
C |      NOM       |MODE|                   ROLE                       C
C |________________|____|______________________________________________C
C |   HPROP        | <= | Work array for the water depth               C
C |   DISP_C       | <= | Work array for the dispersion                C
C |________________|____|______________________________________________C
C                    <=  Can't be change by the user                   C
C                    =>  Can be changed by the user                    C 
C ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY SUSPENSION_COMPUTATION                                     !
!                                                                      !
! CALL      ------                                                     !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!

      ! 1/ MODULES
      ! ----------
      USE INTERFACE_SISYPHE,EX_SUSPENSION_BILAN => SUSPENSION_BILAN
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU

      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE(BIEF_MESH),  INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ),   INTENT(IN)    :: CST,HN,VOLU2D
      TYPE(BIEF_OBJ),   INTENT(IN)    :: ZFCL_S,MASKEL,FLBORTRA
      INTEGER,          INTENT(IN)    :: IELMT,ITRA,LT,NIT,NFRLIQ,NPTFR
      INTEGER,          INTENT(IN)    :: NUMLIQ(NFRLIQ)
      DOUBLE PRECISION, INTENT(IN)    :: DT,CSF
      DOUBLE PRECISION, INTENT(IN)    :: MASSOU,MASED0,AGGLOT
      LOGICAL,          INTENT(IN)    :: MSK,ENTET
      DOUBLE PRECISION, INTENT(INOUT) :: MASTEN,MASTOU,MASINI
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: T2,T3
      DOUBLE PRECISION, INTENT(INOUT) :: MASFIN,MASDEPT,MASDEP

      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER IFRLIQ,I
      DOUBLE PRECISION            :: ERREUR, PERDUE, RELATI, FLUXT
C     300 EST ICI MAXFRO, LE NOMBRE MAXIMUM DE FRONTIERES LIQUIDES
      DOUBLE PRECISION FLT_BOUND(300)      

      ! 4/ EXTERNAL FUNCTION
      ! --------------------
      DOUBLE PRECISION, EXTERNAL :: P_DSUM

!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!

      ! ************************************** !
      ! I - QUANTITY OF SEDIMENT IN SUSPENSION !
      ! ************************************** !
!
      IF(AGGLOT.GT.0.999999D0) THEN
!       ICI ON SUPPOSE QUE AGGLOT=1.D0
        CALL OS('X=YZ    ',X=T2,Y=VOLU2D,Z=CST)
      ELSE
        CALL VECTOR(T2,'=','MASVEC          ',IELMT,
     *              1.D0-AGGLOT,CST,T3,T3,T3,T3,T3,MESH,MSK,MASKEL)
        CALL VECTOR(T3,'=','MASBAS          ',IELMT,
     *                   AGGLOT,T2,T2,T2,T2,T2,T2,MESH,MSK,MASKEL)
        CALL OS('X=X+YZ  ',X=T2,Y=T3,Z=CST)
      ENDIF
C
      MASFIN = DOTS(T2,HN)
      IF(NCSIZE.GT.1) MASFIN=P_DSUM(MASFIN) 
      
      ! ************************** !
      ! II - TOTAL MASS OF DEPOSIT ! 
      ! ************************** !
      
      CALL VECTOR(T2, '=', 'MASVEC          ', IELMT, CSF, ZFCL_S, HN,
     &            HN, HN, HN, HN, MESH, MSK, MASKEL)
      MASDEPT = BIEF_SUM(T2)
      IF(NCSIZE.GT.1) MASDEPT = P_DSUM(MASDEPT)

      ! *************************************************** !
      ! III - MASSE DE SEDIMENTS DEPOSEE (OU ERODEE) TOTALE !
      ! *************************************************** !

      MASDEP = MASDEP + MASDEPT
C
C=======================================================================
C
C   CALCUL DES FLUX (IL MANQUE LE FLUX DIFFUSIF,...A VOIR)
C
C=======================================================================
C
      FLUXT=0.D0
C
      IF(NFRLIQ.GT.0) THEN
        DO IFRLIQ=1,NFRLIQ
          FLT_BOUND(IFRLIQ)=0.D0
        ENDDO
        IF(NPTFR.GT.0) THEN
          DO I=1,NPTFR
!           NOTE: ON POURRAIT DEFINIR FLUX_BOUNDARIES ENTRE 0 ET NFRLIQ
            IFRLIQ=NUMLIQ(I)
            IF(IFRLIQ.GT.0) THEN
              FLT_BOUND(IFRLIQ)=FLT_BOUND(IFRLIQ)+FLBORTRA%R(I)
            ENDIF
          ENDDO
        ENDIF
        IF(NCSIZE.GT.1) THEN
          DO IFRLIQ=1,NFRLIQ
            FLT_BOUND(IFRLIQ)=P_DSUM(FLT_BOUND(IFRLIQ))
          ENDDO
        ENDIF
        DO IFRLIQ=1,NFRLIQ
          FLUXT=FLUXT+FLT_BOUND(IFRLIQ)
        ENDDO
      ENDIF

      ! ********************************************** !
      ! VII - QUANTITY ENTERED THROUGH LIQUID BOUNDARY ! 
      ! ********************************************** !
      
      MASTEN = MASTEN - FLUXT * DT

      ! ************************************** !
      ! VIII - QUANTITY CREATED BY SOURCE TERM ! 
      ! ************************************** !
      
      MASTOU = MASTOU + MASSOU

      ! ***************************** !
      ! IX - RELATIVE ERROR ON VOLUME ! 
      ! ***************************** !
!
! CORRECTION JMH 17/03/05 : MISSING TERM                
!                                         - MASDEPT
      ERREUR = MASINI + MASSOU - DT*FLUXT - MASDEPT - MASFIN
!
      IF (MASFIN > 1.D-8) ERREUR = ERREUR / MASFIN

      ! *********** !
      ! X - LISTING ! 
      ! *********** !
      
      IF(ENTET) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,1005) ITRA,MASINI
          WRITE(LU,1100) ITRA,MASFIN
          IF(NFRLIQ.GT.0) THEN
            DO IFRLIQ=1,NFRLIQ
              WRITE(LU,1110) IFRLIQ,ITRA,-FLT_BOUND(IFRLIQ)
            ENDDO
          ENDIF
          IF (ABS(MASDEPT) > 1.D-8) WRITE(LU,1115) MASDEPT
          IF (ABS(MASSOU ) > 1.D-8) WRITE(LU,1116) MASSOU
          WRITE(LU,1120) ERREUR
        ELSEIF(LNG.EQ.2) THEN
          WRITE(LU,2005) ITRA,MASINI
          WRITE(LU,2100) ITRA,MASFIN
          IF(NFRLIQ.GT.0) THEN
            DO IFRLIQ=1,NFRLIQ
              WRITE(LU,2110) IFRLIQ,ITRA,-FLT_BOUND(IFRLIQ)
            ENDDO
          ENDIF
          IF(ABS(MASDEPT) > 1.D-8) WRITE(LU,2115) MASDEPT
          IF(ABS(MASSOU ) > 1.D-8) WRITE(LU,2116) MASSOU
          WRITE(LU,2120) ERREUR
        ENDIF
      ENDIF

      ! ************************************** !
      ! XI - LISTING OF THE FINAL MASS-BALANCE ! 
      ! ************************************** !
      
      IF(LT.EQ.NIT.AND.ENTET) THEN

         PERDUE = MASED0 + MASTEN + MASTOU - MASFIN - MASDEP
         RELATI = PERDUE
         IF(MAX(MASED0,MASFIN) > 1.D-10) THEN
           RELATI = RELATI / MAX(MASED0,MASFIN)
         ENDIF

         IF(LNG.EQ.1) THEN
            WRITE(LU,3000) ITRA
            WRITE(LU,1140) RELATI
            WRITE(LU,1160) ITRA, MASED0, MASFIN
            IF(ABS(MASTEN) > 1.D-8) WRITE(LU,1161) MASTEN
            IF(ABS(MASTOU) > 1.D-8) WRITE(LU,1164) MASTOU
            IF(ABS(MASDEP) > 1.D-8) WRITE(LU,1167) MASDEP
            WRITE(LU,1166) PERDUE
         ELSEIF(LNG.EQ.2) THEN
            WRITE(LU,3100) ITRA
            WRITE(LU,2140) RELATI
            WRITE(LU,2160) ITRA,MASED0, MASFIN
            IF(ABS(MASTEN) > 1.D-8) WRITE(LU,2161) MASTEN
            IF(ABS(MASTOU) > 1.D-8) WRITE(LU,2164) MASTOU
            IF(ABS(MASDEP) > 1.D-8) WRITE(LU,2167) MASDEP
            WRITE(LU,2166) PERDUE
         ENDIF
      ENDIF

      ! *************************** !
      ! XII - UPDATING INITIAL MASS ! 
      ! *************************** !
      
      MASINI = MASFIN
      
      !----------------------------------------------------------------!
1005  FORMAT(1X,'QUANTITE DE LA CLASSE    ',I2
     &         ,' EN SUSPENSION AU TEMPS T    : ',G16.7,' M3')
1100  FORMAT(1X,'QUANTITE DE LA CLASSE    ',I2
     &         ,' EN SUSPENSION AU TEMPS T+DT : ',G16.7,' M3')
1110  FORMAT(1X,'FRONTIERE ',1I3,' FLUX TRACEUR ',1I2,' = ',G16.7,
     *          ' ( >0 : ENTRANT  <0 : SORTANT )')
1112  FORMAT(1X,'FLUX IMPOSE DE LA CLASSE ',I2
     &         ,'                             : ',G16.7,' M3/S')
1113  FORMAT(1X,'FLUX LIBRE  DE LA CLASSE ',I2
     &         ,'                             : ',G16.7,' M3/S')
1114  FORMAT(1X,'FLUX DE LA CLASSE        ',I2
     &         ,' PAR ONDE INCIDENTE          : ',G16.7,' M3/S') 
1115  FORMAT(1X,'VOLUME DEPOSE SUR LE FOND  '
     &         ,'                             : ',G16.7,' M3')
1116  FORMAT(1X,'VOLUME CREE PAR TERME SOURCE  '
     &         ,   '                          : ',G16.7,' M3')
1120  FORMAT(1X,'ERREUR RELATIVE SUR LE VOLUME  '
     &         ,    '                         : ', G16.7)
1140  FORMAT(1X,'ERREUR RELATIVE CUMULEE SUR LE VOLUME : ', G16.7)
1160  FORMAT(1X,'QUANTITE INITIALE DE ',I2,'               : ',G16.7
     &         ,' M3',
     &     /,1X,'QUANTITE FINALE                       : ', G16.7,' M3')
1161  FORMAT(1X,'QUANTITE ENTREE AUX FRONT. LIQUID.    : ', G16.7,' M3')
1164  FORMAT(1X,'QUANTITE CREEE PAR TERME SOURCE       : ', G16.7,' M3')
1166  FORMAT(1X,'QUANTITE TOTALE PERDUE                : ', G16.7,' M3')
1167  FORMAT(1X,'VOLUME TOTAL DEPOSE SUR LE FOND       : ', G16.7,' M3')
3000  FORMAT(/,1X,'        *** ','BILAN FINAL DE LA CLASSE ',I2,' ***')
      !----------------------------------------------------------------!
2005  FORMAT(1X,'QUANTITY OF CLASS                 ',I2
     &         ,' IN SUSPENSION AT TIME T    : ',G16.7,' M3')
2100  FORMAT(1X,'QUANTITY OF CLASS                 ',I2
     &         ,' IN SUSPENSION AT TIME T+DT : ',G16.7,' M3')
2110  FORMAT(1X,'BOUNDARY ',1I3,' FLUX TRACER ',1I2,' = ',G16.7,
     *          ' ( >0 : ENTERING  <0 : EXITING )')
2112  FORMAT(1X,'PRESCRIBED SEDIMENT FLUX OF CLASS ',I2
     &         ,'                            : ',G16.7,' M3/S')
2113  FORMAT(1X,'FREE FLUX OF CLASS                ',I2
     &         ,'                            : ',G16.7,' M3/S')              
2114  FORMAT(1X,'FLUX OF SEDIMENT CLASS            ',I2
     &         ,' ADDED BY INCIDENT WAVE     : ',G16.7,' M3/S')
2115  FORMAT(1X,'VOLUME OF DEPOSIT                   '
     &         ,'                            : ',G16.7, ' M3')
2116  FORMAT(1X,'VOLUME CREATED BY SOURCE TERM       '
     &         ,'                            : ',G16.7, ' M3')
2120  FORMAT(1X,'RELATIVE ERROR ON VOLUME            '
     &         ,'                            : ',G16.7)
2140  FORMAT(1X,'RELATIVE ERROR CUMULATED ON VOLUME : ', G16.7       )
2160  FORMAT(1X,'INITIAL QUANTITY OF ',I2,'           : '  , G16.7
     &         ,' M3',
     &     /,1X,'FINAL QUANTITY                     : ', G16.7, ' M3')
2161  FORMAT(1X,'QUANTITY ENTERED THROUGH LIQ. BND. : ', G16.7, ' M3')
2164  FORMAT(1X,'QUANTITY CREATED BY SOURCE TERM    : ', G16.7, ' M3')
2166  FORMAT(1X,'TOTAL QUANTITY LOST                : ', G16.7, ' M3')
2167  FORMAT(1X,'TOTAL MASS OF DEPOSIT              : ', G16.7, ' M3')
3100  FORMAT(/,1X,'      *** ','FINAL BALANCE FOR TRACER',I2,' ***')
      !----------------------------------------------------------------!

!======================================================================!
!======================================================================!

      RETURN
      END
