C                        ********************
                         SUBROUTINE TASSEMENT 
C                        ********************
C
     *(ZF,NPOIN,DTS,ELAY,DZF_TASS,T2,LT,AVAIL,NSICLA,ES,XMVS,
     * XKV,TRANS_MASS,CONC_VASE,NCOUCH_TASS,MS_SABLE,MS_VASE)
C
C***********************************************************************
C PROGRAMME APPELANT : SISYPHE
C PROGRAMMES APPELES : 
C***********************************************************************
C
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER,          INTENT(IN)    :: NPOIN,NSICLA
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: DZF_TASS,ZF,ELAY,T2
      DOUBLE PRECISION, INTENT(INOUT) :: MS_SABLE(NPOIN,10)
      DOUBLE PRECISION, INTENT(INOUT) :: MS_VASE(NPOIN,10)
      DOUBLE PRECISION, INTENT(IN)    :: DTS
      INTEGER, INTENT(IN)             :: LT,NCOUCH_TASS
      DOUBLE PRECISION, INTENT(INOUT) :: AVAIL(NPOIN,10,NSICLA)
      DOUBLE PRECISION, INTENT(INOUT) :: ES(NPOIN,10)
      DOUBLE PRECISION, INTENT(IN)    :: TRANS_MASS(10),CONC_VASE(10) 
      DOUBLE PRECISION, INTENT(IN)    :: XMVS,XKV
C   
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,J
      DOUBLE PRECISION CONC_SABLE
C
      DOUBLE PRECISION TAUX(10),TRANSFERT_MASSE_VASE(10)
      DOUBLE PRECISION TRANSFERT_MASSE_SABLE(10)
      DOUBLE PRECISION EPAI_SABLE(10),EPAI_VASE(10)    
C
C CALCUL EPAISSEUR TOTALE SEDIMENT (SABLE + VASE) AVANT TASSEMENT
C 
      CONC_SABLE=XMVS/XKV
C
C T2: epaisseur de vase avant tassement
C
      DO I=1,NPOIN
C
        T2%R(I)=0.D0
        DO J=1,NCOUCH_TASS
          EPAI_VASE(J)=MS_VASE(I,J)/CONC_VASE(J)
          ES(I,J)=EPAI_VASE(J)
          IF(NSICLA.GT.1) THEN
            EPAI_SABLE(J)=MS_SABLE(I,J)/XMVS
            ES(I,J)=EPAI_VASE(J)+EPAI_SABLE(J)
          ENDIF
          T2%R(I)=T2%R(I)+ES(I,J)                       
        ENDDO
C        
        DO J=1,NCOUCH_TASS
          IF(MS_VASE(I,J).GE.1.D-6) THEN
            TRANSFERT_MASSE_VASE(J)=MIN(MS_VASE(I,J),
     &              MS_VASE(I,J)*DTS*TRANS_MASS(J))
            IF(NSICLA.GT.1) THEN
              TAUX(J)=TRANSFERT_MASSE_VASE(J)/MS_VASE(I,J)
              TRANSFERT_MASSE_SABLE(J)=TAUX(J)*MS_SABLE(I,J)
            ENDIF
          ELSE
            TRANSFERT_MASSE_VASE(J)=0.D0
            IF(NSICLA.GT.1) TRANSFERT_MASSE_SABLE(J)=0.D0
          ENDIF
C**************ARRET DE TASSEMENT SI LA VASE A REMPLI LES INTERSTICES
C**************   ENTRE LES GRAINS DE SABLE                
          IF(NSICLA.GT.1.AND.EPAI_SABLE(J).GE.ES(I,J)) THEN
            TRANSFERT_MASSE_VASE(J) =0.D0
            TRANSFERT_MASSE_SABLE(J)=0.D0
          ENDIF 
        ENDDO
C
        DO J=1,NCOUCH_TASS
          IF(J.EQ.NCOUCH_TASS) THEN         
             MS_VASE(I,J)=MAX(0.D0,MS_VASE(I,J)
     &            +TRANSFERT_MASSE_VASE(J-1))
             IF(NSICLA.GT.1) THEN
                MS_SABLE(I,J)=MAX(0.D0,MS_SABLE(I,J)
     &                                  +TRANSFERT_MASSE_SABLE(J-1))
             ENDIF
          ELSEIF(J.EQ.1) THEN
             MS_VASE(I,J)=MAX(0.D0,MS_VASE(I,J)
     &            -TRANSFERT_MASSE_VASE(J))
            IF(NSICLA.GT.1) THEN
              MS_SABLE(I,J)=MAX(0.D0,MS_SABLE(I,J)
     &            -TRANSFERT_MASSE_SABLE(J))
            ENDIF              
          ELSE
             MS_VASE(I,J)=MAX(0.D0,MS_VASE(I,J)
     &            +TRANSFERT_MASSE_VASE(J-1)-TRANSFERT_MASSE_VASE(J))
             IF(NSICLA.GT.1) THEN 
               MS_SABLE(I,J)=MAX(0.D0,MS_SABLE(I,J)
     &         +TRANSFERT_MASSE_SABLE(J-1)-TRANSFERT_MASSE_SABLE(J))
             ENDIF
          ENDIF
        ENDDO
C
        ELAY%R(I)=0.D0     
C           
        DO J=1,NCOUCH_TASS
          EPAI_VASE(J)=MS_VASE(I,J)/CONC_VASE(J)
          ES(I,J) = EPAI_VASE (J)
          IF(NSICLA.GT.1) THEN
            EPAI_SABLE(J)=MS_SABLE(I,J)/XMVS
            ES(I,J)=EPAI_VASE(J)+EPAI_SABLE(J)
          ENDIF
          ELAY%R(I)=ELAY%R(I) + ES(I,J)
        ENDDO
C          
C       EVOLUTION DU FOND DUE AU TASSEMENT
C
        DZF_TASS%R(I)=ELAY%R(I)-T2%R(I)        
C          
C NOTE JMH : J'AVAIS COMPRIS QUE CLASSE 1 = VASE
C            ET A PARTIR DE 2 : SABLE, ET LA QU'EST-CE QU'ON FAIT ??
C
        IF(NSICLA.GT.1) THEN
          DO J=1,NCOUCH_TASS
           IF(ES(I,J).GE.1.D-6) THEN
             AVAIL(I,J,1)=MS_SABLE(I,J)/XMVS/ES(I,J)
             AVAIL(I,J,2)=MS_VASE(I,J)/CONC_VASE(J)/ES(I,J)             
           ELSE
             AVAIL(I,J,1)=0.D0
             AVAIL(I,J,2)=0.D0
           ENDIF
          ENDDO    
        ENDIF
C
      ENDDO 
C
C-----------------------------------------------------------------------
C
      RETURN
      END
