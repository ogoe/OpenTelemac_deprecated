C                       *********************            
                        SUBROUTINE INIT_MIXTE
C                       *********************
C
     *(XMVS,NPOIN,AVAIL,NSICLA,ES,ELAY,NCOUCH_TASS,CONC_VASE,
     * MS_SABLE,MS_VASE,ZF,ZR,AVA0)
C
C***********************************************************************
C SISYPHE VERSION 6.0                            
C                                                 
C                                                
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT      
C***********************************************************************
C
C 
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |                |<-- | 
C |________________|____|______________________________________________
C MODE : -->(INPUT ), <--(RESULT), <-->(MODIFIED VARIABLE)
C-----------------------------------------------------------------------
C CALLED BY SUBROUTINE  : SISYPHE
C 
C***********************************************************************
C
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NPOIN,NSICLA,NCOUCH_TASS 
      DOUBLE PRECISION, INTENT(IN)    :: XMVS
      DOUBLE PRECISION, INTENT(INOUT) :: AVAIL(NPOIN,10,NSICLA)
      DOUBLE PRECISION, INTENT(INOUT) :: ES(NPOIN,10),ELAY(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: ZR(NPOIN),ZF(NPOIN)
      DOUBLE PRECISION,  INTENT(INOUT) :: MS_SABLE(NPOIN,10)
      DOUBLE PRECISION,  INTENT(INOUT) :: MS_VASE(NPOIN,10)
      DOUBLE PRECISION, INTENT(IN)    :: CONC_VASE(10)  
      DOUBLE PRECISION, INTENT(IN)   :: AVA0(NSICLA)
C
C-----------------------------------------------------------------------
C     LOCAL VARIABLES
C
      INTEGER I,J
C
      DOUBLE PRECISION EPAI_VASE(10),EPAI_SABLE(10)  
      DOUBLE PRECISION DIFF
C
C ------------------------------------------------------------------------
C
C*******COMPOSITION INITIALE DU SEDIMENT identique en chaque point
C répartion uniforme des épaisseurs de sédiments 
C épaisseurs des couches de vase  (volume de chaque phase/
          EPAI_VASE(1)=0.0525D0
          EPAI_VASE(2)=0.0385D0
          EPAI_VASE(3)=0.03995D0
          EPAI_VASE(4)=0.0437D0
          EPAI_VASE(5)=0.0517D0
          EPAI_VASE(6)=0.1259D0
          EPAI_VASE(7)=0.4889D0
          EPAI_VASE(8)=1.5071D0
          EPAI_VASE(9)=0.86410D0
C          EPAI_VASE(10)=100.D0
C
C épaisseurs des couches de sable
C  (indice des vides inclus)   
C Pour sediments mixtes seulement
          IF (NSICLA.GT.1) THEN
C verification des donnees
            DO J= 1, NCOUCH_TASS-1
C              
              EPAI_SABLE(J) = AVA0(1)/AVA0(2)* EPAI_VASE(J)
C
            ENDDO        
           ENDIF
C
C          
CC EPAISSEUR DE VASES = ZF-ZR 
C Calcul de la dernière épaisseur : 
C 
        DO I=1,NPOIN
C
           ELAY(I)=0.D0
           DO J=1,NCOUCH_TASS-1
             ES(I,J)= EPAI_VASE(J)
C            si mixte    
             IF(NSICLA.GT.1) THEN
               ES(I,J)= ES(I,J)  + EPAI_SABLE(J)
             ENDIF 
             ELAY(I)=ELAY(I)+ES(I,J)
           ENDDO
C
           DIFF= (ZF(I)-ZR(I)) - ELAY(I)
C
           IF(DIFF.GT.0.D0) THEN
             ES(I,NCOUCH_TASS) = DIFF
             ELAY(I) = ZF(I)-ZR(I)
           ELSE          
C            définition des fonds rigides à zr=zf
             DO J=NCOUCH_TASS,1,-1
               ES(I,J) =0.D0
             ENDDO
             ELAY(I)=0.D0 
           ENDIF            
         ENDDO

C     
C calcul des masses de vase et de sable initiales
C 
      DO I=1,NPOIN
C
        DO J=1,NCOUCH_TASS
C
C *************** ON COMBLE LES VIDES ENTRE LES GRAINS DE SABLE
C vase pure
             IF(NSICLA.EQ.1) THEN
               MS_VASE(I,J) = ES(I,J)*CONC_VASE(J)
               AVAIL(I,J,1) = 1.D0                             
C
C SI mixte
C
            ELSE
              MS_VASE(I,J) = ES(I,J)*CONC_VASE(J)*AVA0(2)           
              MS_SABLE(I,J)= ES(I,J)*XMVS*AVA0(1)
              IF(ES(I,J).GE.1.D-6) THEN
                AVAIL(I,J,1)= AVA0(1)
                AVAIL(I,J,2)= AVA0(2)        
              ELSE
                AVAIL(I,J,1)= 0.D0
                AVAIL(I,J,2)= 0.D0        
              ENDIF 
            ENDIF                    
        ENDDO
C
      ENDDO
C 
C-----------------------------------------------------------------------
C    
      RETURN
      END
