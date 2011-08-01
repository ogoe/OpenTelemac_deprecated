C                       *****************                               
                        SUBROUTINE INTEMP                               
C                       *****************                               
C                                                                       
     *(W,FLUX,AIRS,DT,NPOIN,ZF,CF,EPS,DDMIN,KFROT,SMH)                   
C                                                                       
C***********************************************************************
C                                         N. GOUTAL 24/11/97
C-----------------------------------------------------------------------
C                                                                       
C  FONCTION  : . INTEGRATION EN TEMPS                                   
C  
C
C  NOTE JMH : DDMIN NON UTILISE
C                                                                     
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C ! . FLUX         ! -->! FLUX DE ROE                                  !
C ! . AIRS         ! -->! TABLEAU DES AIRES DES CELLULES               !
C ! . H            !<-->! ENTHALPIE                                    !
C !                !    !                                              !
C ! . DT           ! -->! PAS DE TEMPS.                                !
C ! . 3            ! -->! DIMENSION DU SYSTEME                         !
C !________________!____!______________________________________________!
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)   
C        -- (TABLEAU DE TRAVAIL)                                        
C-----------------------------------------------------------------------
C     - SOUS PROGRAMME(S) APPELANT : RESOLU                             
C     - SOUS PROGRAMME(S) APPELE   : AUCUN                              
C     - PORTABILITE: CRAY                                               
C***********************************************************************
C                                                                       
      IMPLICIT NONE
C                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NPOIN,KFROT
      DOUBLE PRECISION, INTENT(INOUT) :: W(3,NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: FLUX(NPOIN,3),DT,EPS,DDMIN
      DOUBLE PRECISION, INTENT(IN)    :: AIRS(NPOIN),ZF(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: CF(NPOIN),SMH(NPOIN) 
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                       
      INTEGER I                                                   
C                                                                                             
      DOUBLE PRECISION DELTA,KF,ST2D,ALPHAF            
C
C-----------------------------------------------------------------------
C                                                                       
C------                                                                 
C 1. CALCUL DE W A L'INSTANT N+1                                        
C------                                                                 
C                                                                       
      DO 10 I = 1 , NPOIN                                               
         W(1,I) = W(1,I) - DT *( FLUX (I,1)-SMH(I)) / AIRS(I)                     
         W(2,I) = W(2,I) - DT * FLUX (I,2) / AIRS(I)                     
         W(3,I) = W(3,I) - DT * FLUX (I,3) / AIRS(I)                     
10    CONTINUE                                                          
C                                                                       
C     PRISE EN COMPTE DU FROTTEMENT (PAS DE FROTTEMENT RITTER)    
C     *****************************
      IF (KFROT.NE.0) THEN
C                                                                       
         DO 20 I = 1,NPOIN                           
C                                                                       
C FH-FRDATA
!            IF (W(1,I).GT.EPS/10.D0) THEN
            IF ((W(1,I).GT.EPS/10.D0).AND.(CF(I).GT.1.D-12)) THEN
C FH-FRDATA
               ST2D = CF(I)
               KF = 9.81D0*DT*DSQRT(W(2,I)**2+W(3,I)**2)/
     *              (ST2D*ST2D*W(1,I)**(7.D0/3.D0))
               IF (KF.GT.1.D-6) THEN                                          
                  DELTA = (1.D0+4.D0*KF)                                      
                  ALPHAF = (-1.D0+SQRT(DELTA))/(2*KF)                         
               ELSE                                                          
                  ALPHAF = 1.D0 - KF                                          
               ENDIF
                  W(2,I) = ALPHAF * W(2,I)                                      
                  W(3,I) = ALPHAF * W(3,I)                                      
            ENDIF                                                           
C                                                                       
20       CONTINUE
C 
      ENDIF                                                           
C                                                          
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END
