C                         ******************
                          SUBROUTINE SPECTRE
C                         ******************
C   
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.9    16/10/2008                        Chun WANG
C
C  16/10/2008 JMH : VERIFICATION OF TIME RANGE FOR ANALYSIS
C
C***********************************************************************
C
C      FONCTION:    Harmonic analysis of tide waves using Leat-Squares 
C                   Fitting Method (refer to "Simulation Des
C                   Courants de Maree en Manche et Proche Atlantique", 
C                   EDF report, J. M. Janin et. al. Page 40-41.)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !   H            !--> !  water depth
C !   ZF           !--> !  bathymetry
C !   NPOIN        !--> !  number of points in the mesh.
C !   LT           !--> !  iteration time step
C !   PERIAF       !--> !  periods of waves
C !   NPERIAF      !--> !  number of waves
C !   M            !--> !  sampling points  
C !   DT           !--> !  time interval
C !   AT           !--> !  current time
C !   AFBGN        !--> !  time step at which harmonic analysis begins
C !   AFEND        !--> !  time step at which harmonic analysis ends
C !   AM           !<-- !  coefficient matrix in least square method
C !   BM           !<-- !  inverse matrix of AM
C !________________!____!_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C_______________________________________________________________________
C
C SOUS-PROGRAMME APPELE : DEBIMP
C
C-----------------------------------------------------------------------     
C
      USE BIEF
      USE DECLARATIONS_TELEMAC2D
      USE INTERFACE_TELEMAC2D
C
      IMPLICIT NONE
      INTEGER LNG,LU             
      COMMON/INFO/LNG,LU
C
      DOUBLE PRECISION PI,CFX,CFY,A      
C
      DOUBLE PRECISION, ALLOCATABLE :: AM(:,:),BM(:,:),HA(:)
C
      INTEGER M,I,J,K,N,AFBGN,AFEND
      SAVE AFBGN,AFEND,AM,BM,HA
C
      LOGICAL DEJA
      DATA DEJA/.FALSE./
C
      INTRINSIC COS,SIN,ACOS,DMOD,DATAN2 
C
      DOUBLE PRECISION P_DMIN,P_DMAX
      EXTERNAL         P_DMIN,P_DMAX   
C
C-----------------------------------------------------------------------     
C
      IF(.NOT.DEJA) THEN
C
        AFBGN=INT((TAFBGN-AT)/DT +1.D-6)
        AFEND=INT((TAFEND-AT)/DT +1.D-6)
        IF(AFEND.LE.AFBGN) THEN
          IF(LNG.EQ.1) THEN
            WRITE(LU,*)
            WRITE(LU,*) 'BORNES EN TEMPS POUR L''ANALYSE DE FOURIER'
            WRITE(LU,*) 'VERIFIER OU RENSEIGNER CE MOT-CLE'
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*)
            WRITE(LU,*) 'TIME RANGE FOR FOURIER ANALYSIS'
            WRITE(LU,*) 'CHECK OR USE THIS KEY-WORD'
          ENDIF
          CALL PLANTE(1)
          STOP
        ENDIF
        ALLOCATE(AM(2*NPERIAF,2*NPERIAF))
        ALLOCATE(BM(2*NPERIAF,2*NPERIAF))
        ALLOCATE(HA(2*NPERIAF))
        DEJA=.TRUE.
C       INITIALISATION A ZERO DES LE PREMIER APPEL (PENSER A DESIMP)
        DO I = 1,NPERIAF
          DO J = 1,NPOIN
            AMPL%ADR(I)%P%R(J) = 0.D0
            PHAS%ADR(I)%P%R(J) = 0.D0
          ENDDO	
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------     
C   
      PI = ACOS(-1.D0)
      M = AFEND-AFBGN+1      
C    
      IF(LT.GE.AFBGN.AND.LT.LE.AFEND) THEN      
         DO I=1,NPERIAF
            A=2.D0*PI*DMOD(AT/PERIAF(I),1.D0)
            CFX=COS(A)/M
            CFY=SIN(A)/M
            DO J=1,NPOIN	    
         AMPL%ADR(I)%P%R(J)=AMPL%ADR(I)%P%R(J)+(H%R(J)+ZF%R(J))*CFX
         PHAS%ADR(I)%P%R(J)=PHAS%ADR(I)%P%R(J)+(H%R(J)+ZF%R(J))*CFY
            ENDDO
         ENDDO
      ENDIF
C
C Establish the coefficient matrices and inverse it. Multiplied by 
C right hand sides, we get the unkowns of linear eqns. Now, AMPL 
C contains the amplitude of the SIN harmonic components of surface 
C elevation; PHAS contains the amplitude of the COS harmonic 
C components of surface elevation. 
C
      IF(LT.EQ.AFEND) THEN      
         CALL COEFMAT(PERIAF,DT,M,AM,NPERIAF)  
         N=2*NPERIAF         
         CALL INVMTX(AM,BM,N)                              
         DO J=1,NPOIN      
            DO I=1,NPERIAF
               HA(I) = AMPL%ADR(I)%P%R(J)
               HA(I+NPERIAF) = PHAS%ADR(I)%P%R(J)
            ENDDO   
            DO I = 1,NPERIAF
               AMPL%ADR(I)%P%R(J) = 0
               PHAS%ADR(I)%P%R(J) = 0 
               DO K=1,2*NPERIAF	       
         AMPL%ADR(I)%P%R(J) = AMPL%ADR(I)%P%R(J)+BM(I,K)*HA(K)
         PHAS%ADR(I)%P%R(J) = PHAS%ADR(I)%P%R(J)+BM(I+NPERIAF,K)*HA(K)
               ENDDO
            ENDDO
         ENDDO
C	 	 
         DO J = 1, NPOIN	 
            DO I = 1,NPERIAF
               CFX = AMPL%ADR(I)%P%R(J)
               CFY = PHAS%ADR(I)%P%R(J)	       
               AMPL%ADR(I)%P%R(J) = SQRT(CFX**2+CFY**2) 
               PHAS%ADR(I)%P%R(J) = 180.D0*DATAN2(CFY,CFX)/PI  
               IF (PHAS%ADR(I)%P%R(J).LT.0.D0) 
     *	           PHAS%ADR(I)%P%R(J)=PHAS%ADR(I)%P%R(J)+360.D0
            ENDDO
         ENDDO
      ENDIF
C
C Output the amplitudes and phases of each interested points.
C
      IF(LT.EQ.NIT) THEN
         DO I = 1, NPERIAF
            IF(NPTS.GT.0) THEN
               WRITE(LU,*) ' '
               WRITE(LU,*) ' '
               WRITE(LU,*) ' '	      
               IF(LNG.EQ.1) WRITE(LU,*) 'ANALYSE DE LA PERIODE ', 
     *	                                 PERIAF(I), ' S :'
               IF(LNG.EQ.2) WRITE(LU,*) 'ANALYSE OF PERIOD ', 
     *	                                 PERIAF(I), ' S :'  
               WRITE(LU,*) ' '      
               WRITE(LU,90) 'NOM DE POINT', 'AMPLITUDE', 'PHASE'
               WRITE(LU,*) ' '                
               DO J = 1, NPTS 
!                 IN PARALLEL POINT NOT ALWAYS EXISTING, MAYBE ELSEWHERE
                  IF(NCSIZE.GT.0) THEN
                    WRITE(LU,100) NAME_PTS(J), 
     *                            P_DMIN(AMPL%ADR(I)%P%R(LIST_PTS(J)))
     *                           +P_DMAX(AMPL%ADR(I)%P%R(LIST_PTS(J))),
     *	                          P_DMIN(PHAS%ADR(I)%P%R(LIST_PTS(J)))
     *                           +P_DMAX(PHAS%ADR(I)%P%R(LIST_PTS(J)))
                  ELSE         
                    WRITE(LU,100) NAME_PTS(J), 
     *                            AMPL%ADR(I)%P%R(LIST_PTS(J)),
     *	                          PHAS%ADR(I)%P%R(LIST_PTS(J))
                  ENDIF
               ENDDO
            ENDIF   
         ENDDO       
      ENDIF
C
 90   FORMAT(1X, A15, A16  , A11  )  
 100  FORMAT(1X, A15, F16.3, F11.2)        
C
C-----------------------------------------------------------------------     
C
      RETURN
      END
