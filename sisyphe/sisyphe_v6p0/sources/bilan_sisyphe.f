C                       ************************
                        SUBROUTINE BILAN_SISYPHE
C                       ************************
C
     *( E      , ESOMT  , QSX    , QSY    , MESH   , MSK    , MASKEL ,
     *  T1     , T2     , S      , IELMU  , VCUMU  , DT     , NPTFR  ,
     *  INFO   , ZFCL_C , QSCLXC , QSCLYC , NSICLA ,
     *  VOLTOT , DZF_GF , MASS_GF, LGRAFED, NUMLIQ , NFRLIQ)
C
C  T2 NON UTILISE
C
C***********************************************************************
C SISYPHE VERSION 5.9               CHANGED BY CMGDL FOR GRADED SEDIMENT       
C***********************************************************************
C
C     FONCTION  : EFFECTUE LE BILAN DE MASSE
C                 
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C PROGRAMME APPELANT : RESOLU
C PROGRAMMES APPELES :
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)          :: NPTFR,NFRLIQ,IELMU,NSICLA
      INTEGER, INTENT(IN)          :: NUMLIQ(NPTFR)
      DOUBLE PRECISION, INTENT(IN) :: DT   
      LOGICAL, INTENT(IN)          :: MSK, INFO
C
      LOGICAL,          INTENT(IN)    :: LGRAFED
      DOUBLE PRECISION, INTENT(INOUT) :: MASS_GF,VCUMU
      DOUBLE PRECISION, INTENT(IN)    :: VOLTOT(10) 
C
C-----------------------------------------------------------------------
C
C  STRUCTURES DE VECTEURS
C 
      TYPE(BIEF_OBJ), INTENT(IN)    :: MASKEL,S,ZFCL_C,QSCLXC,QSCLYC
      TYPE(BIEF_OBJ), INTENT(IN)    :: E,ESOMT,QSX,QSY,DZF_GF
      TYPE(BIEF_OBJ), INTENT(INOUT) :: T1,T2 
C
C-----------------------------------------------------------------------
C
C  STRUCTURES DE MAILLAGE
C
      TYPE(BIEF_MESH) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,IFRLIQ,IPTFR
      DOUBLE PRECISION RMASSE,RCUMU,RMASCLA(10)
      DOUBLE PRECISION VCUMUCLA(10),MASST,FLUXT
C     300 EST ICI MAXFRO, LE NOMBRE MAXIMUM DE FRONTIERES LIQUIDES
      DOUBLE PRECISION FLT_BOUND(300)  
C  
      DOUBLE PRECISION P_DSUM
      EXTERNAL         P_DSUM  
C
C-----------------------------------------------------------------------
C
C     CALCUL DE L'EVOLUTION (E)
C
      CALL VECTOR(T1,'=','MASVEC          ',IELMU,
     *            1.D0,E,S,S,S,S,S,MESH,MSK,MASKEL)
      RMASSE = BIEF_SUM(T1)
      IF(NCSIZE.GT.1) RMASSE = P_DSUM(RMASSE)
C                                                                       
C=======================================================================
C
C   CALCUL DE L'INTEGRALE DE L'EVOLUTION (ESOMT)
C
      CALL VECTOR(T1,'=','MASVEC          ',IELMU,
     *            1.D0,ESOMT,S,S,S,S,S,MESH,MSK,MASKEL)
      RCUMU = BIEF_SUM(T1)
      IF(NCSIZE.GT.1) RCUMU = P_DSUM(RCUMU)
C                                                                       
C=======================================================================
C                                                                       
C   CALCUL DES FLUX AUX FRONTIERES                               
C  
      CALL VECTOR(T1,'=','FLUBOR          ',IELBOR(IELMU,1),
     *            1.D0,S,S,S,QSX,QSY,S,MESH,MSK,MASKEL)
C
      FLUXT=0.D0
C
      IF(NFRLIQ.GT.0) THEN
        DO IFRLIQ=1,NFRLIQ
          FLT_BOUND(IFRLIQ)=0.D0
        ENDDO
        IF(NPTFR.GT.0) THEN
          DO IPTFR=1,NPTFR
            IFRLIQ=NUMLIQ(IPTFR)
            IF(IFRLIQ.GT.0) THEN
              FLT_BOUND(IFRLIQ)=FLT_BOUND(IFRLIQ)+T1%R(IPTFR)
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
C
      VCUMU = VCUMU - FLUXT*DT
C
C     BILAN EN GRANULOMETRIE ETENDUE
C
      IF(NSICLA.GT.1) THEN
C
        DO I=1,NSICLA
C
C       CALCUL DE L'EVOLUTION PAR CLASSE
C       
        CALL VECTOR(T1,'=','MASVEC          ',IELMU,
     &              1.D0,ZFCL_C%ADR(I)%P,S,S,S,S,S,
     &              MESH,MSK,MASKEL)
        RMASCLA(I) = BIEF_SUM(T1)
        IF(NCSIZE.GT.1) RMASCLA(I) = P_DSUM(RMASCLA(I))
C
C       CALCUL DES FLUX LIBRES PAR CLASSE
C
        CALL VECTOR(T1,'=','FLUBOR          ',IELBOR(IELMU,1),
     *              1.D0,S,S,S,QSCLXC%ADR(I)%P,QSCLYC%ADR(I)%P,
     *              S,MESH,MSK,MASKEL)
C
        FLUXT=0.D0
        IF(NFRLIQ.GT.0) THEN
          IF(NPTFR.GT.0) THEN
            DO IPTFR=1,NPTFR
              IFRLIQ=NUMLIQ(IPTFR)
              IF(IFRLIQ.GT.0) THEN
                FLUXT=FLUXT+T1%R(IPTFR)
              ENDIF
            ENDDO
          ENDIF
          IF(NCSIZE.GT.1) FLUXT=P_DSUM(FLUXT)
        ENDIF                               
C
        VCUMUCLA(I) = - FLUXT*DT
C
        ENDDO
C
      ENDIF
C
C=======================================================================
C
C     GRAIN-FEEDING
!     IF (LGRAFED) THEN
!        CALL VECTOR(T1,'=','MASVEC          ',IELMU,
!    &               1.D0,DZF_GF,S,S,S,S,S,MESH,MSK,MASKEL)
!        MASST = BIEF_SUM(T1)
!        IF(NCSIZE.GT.1) MASST = P_DSUM(MASST)
!        MASS_GF = MASS_GF + MASST
!     ENDIF
!     IF(DREDGESIM) ...   ?????
C
C  IMPRESSION DU BILAN
C
      IF (INFO) THEN  
C
          WRITE(LU,*)
          IF(LNG.EQ.1) THEN
            WRITE(LU,1000)                                           
            WRITE(LU,1010) RMASSE 
            IF(NFRLIQ.GT.0) THEN                                      
              DO IFRLIQ=1,NFRLIQ
                WRITE(LU,1110) IFRLIQ,-FLT_BOUND(IFRLIQ)
              ENDDO
            ENDIF
            WRITE(LU,1030) RCUMU
            WRITE(LU,1031) VCUMU                            
          ELSEIF(LNG.EQ.2) THEN
            WRITE(LU,2000)                                          
            WRITE(LU,2010) RMASSE
            IF(NFRLIQ.GT.0) THEN                                       
              DO IFRLIQ=1,NFRLIQ
                WRITE(LU,2110) IFRLIQ,-FLT_BOUND(IFRLIQ)
              ENDDO
            ENDIF
            WRITE(LU,2030) RCUMU
            WRITE(LU,2031) VCUMU                           
          ENDIF
         IF (NSICLA>1) THEN
          DO I=1,NSICLA
           IF(LNG.EQ.1) THEN
             WRITE(LU,*) '     BILAN POUR LA CLASSE DE SEDIMENT :',I
             WRITE(LU,*) '     VOLUME TOTAL DE LA CLASSE :',VOLTOT(I)
             WRITE(LU,3011) RMASCLA(I)
             WRITE(LU,3032) VCUMUCLA(I)
           ELSEIF(LNG.EQ.2) THEN
             WRITE(LU,*) '     MASS BALANCE FOR SEDIMENT CLASS :',I
             WRITE(LU,*) '     TOTAL VOLUME:',VOLTOT(I)
             WRITE(LU,3010) RMASCLA(I)
             WRITE(LU,3031) VCUMUCLA(I)
           ENDIF           
          ENDDO
         ENDIF
         IF (LGRAFED) THEN
            IF (LNG.EQ.1) THEN
               WRITE(LU, 4000) MASST
               WRITE(LU, 4010) MASS_GF
            ENDIF
            IF (LNG.EQ.2) THEN
               WRITE(LU, 4001) MASST
               WRITE(LU, 4011) MASS_GF
            ENDIF
         ENDIF
      ENDIF
C
1000  FORMAT(1X,'BILAN DE MASSE : ')
1010  FORMAT(1X,'SOMME DES EVOLUTIONS : ',G16.7)
1020  FORMAT(1X,'FLUX IMPOSE          : ', G16.7,' M3/S'
     *         ,'  ( >0 : ENTRANT  <0 : SORTANT )')                     
1021  FORMAT(1X,'FLUX LIBRE           : ', G16.7,' M3/S'
     *         ,'  ( >0 : ENTRANT  <0 : SORTANT )')           
1030  FORMAT(1X,'SOMME DES EVOLUTIONS CUMULEES : ',G16.7)  
1031  FORMAT(1X,'VOLUME ENTRE AUX FRONTIERES   : ',G16.7,' M3'
     *         ,'  ( SI <0 VOLUME SORTI )')
1110  FORMAT(1X,'FRONTIERE ',1I3,' FLUX EN CHARRIAGE = ',G16.7,
     *          ' ( >0 : ENTRANT  <0 : SORTANT )')
2000  FORMAT(1X,'MASS-BALANCE          : ')
2010  FORMAT(1X,'SUM OF THE EVOLUTIONS : ',G16.7)
2020  FORMAT(1X,'PRESCRIBED FLOW       : ',G16.7,' M3/S'
     *         ,'  ( >0 : ENTERING  <0 : EXITING )')
2021  FORMAT(1X,'FREE FLOW             : ',G16.7,' M3/S'
     *         ,'  ( >0 : ENTERING  <0 : EXITING )')
2030  FORMAT(1X,'SUM OF THE CUMULATED EVOLUTIONS : ',G16.7)
2031  FORMAT(1X,'VOLUME THAT ENTERED THE DOMAIN  : ',G16.7,' M3'
     *         ,'  ( IF <0 EXIT )')
2110  FORMAT(1X,'BOUNDARY ',1I3,' BEDLOAD FLUX = ',G16.7,
     *          ' ( >0 : ENTERING  <0 : EXITING )')
3010  FORMAT(1X,'SUM OF THE EVOLUTIONS FOR THIS CLASS: ',G16.7)
3031  FORMAT(1X,'VOLUME THAT ENTERED THE DOMAIN FOR THIS CLASS: '
     *       ,G16.7,' M3')
3011  FORMAT(1X,'SOMME DES EVOLUTIONS POUR CETTE CLASSE : ',G16.7)
3032  FORMAT(1X,'VOLUME ENTRE DANS LE DOMAINE POUR CETTE CLASSE : '
     *       ,G16.7,' M3')

4000  FORMAT(1X,'GRAIN-FEEDING A CET INSTANT       : ',G16.7)
4010  FORMAT(1X,'GRAIN-FEEDING JUSQU''A MAINTENANT : ',G16.7)
4001  FORMAT(1X,'GRAIN-FEEDING THIS MOMENT : ',G16.7)
4011  FORMAT(1X,'GRAIN-FEEDING UNTIL NOW   : ',G16.7)

C
C-----------------------------------------------------------------------
C
      RETURN
      END
