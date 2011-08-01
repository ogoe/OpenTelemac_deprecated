C                       *******************                                                    
                        SUBROUTINE INTERPOL                                   
C                       *******************                                     
C     
     *(RO,R02,R03,JCOUT1,JCOUT2,JCOUT3)
C     
C***********************************************************************
C     TELEMAC 2D VERSION 5.2    27/04/93    E. BARROS
C             UPGRADE TO 5.1    05/10/00    A. LEOPARDI (UNINA)
C***********************************************************************
C     
C     FUNCTION  : COMPUTATION OF RO
C     MINIMUM OF THE FUNCTION :   A * (RO**2) + B * RO +C            
C     IT IS :  RO = -B / 2 * A                                       
C     
C-----------------------------------------------------------------------
C     ARGUMENTS                                         
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    RO          |<-- | COEFFICIENT OF DESC                          |
C |    JCOUT1      | -->| J(RO)                                        |
C |    JCOUT2      | -->| J(R02)                                       |
C |    JCOUT3      | -->| J(R03)                                       |
C |    R02         | -->| COEFFICIENT                                  |
C |    R03         | -->| COEFFICIENT                                  |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C     APPELE PAR :             HOMERE_PIT
C     
C     SOUS-PROGRAMME APPELE :  RIEN
C     
C***********************************************************************
C     
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION , INTENT(IN)    :: R02,R03,JCOUT1,JCOUT2,JCOUT3
      DOUBLE PRECISION , INTENT(INOUT) :: RO
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C     
      DOUBLE PRECISION COEFA,COEFB,COEFC
      DOUBLE PRECISION ROMAX
C     
      INTRINSIC ABS
C
C-----------------------------------------------------------------------
C                                                                      
      COEFA = ((JCOUT1*(R02-R03))+(R03*JCOUT2)-(JCOUT3*R02))/(R02*R03   
     *     *(R02-R03))                                              
C     
      COEFB = ((-JCOUT1*((R02*R02)-(R03*R03)))-(JCOUT2*R03*R03)
     *     + (JCOUT3*R02*R02))/(R02*R03*(R02-R03))
C     
      COEFC = JCOUT1                                                   
C         
      IF(COEFA.LE.0.D0) THEN                                           
        WRITE(LU,*) 'INTERPOL : COEFFICIENT A LESS THAN ZERO:',COEFA
        CALL PLANTE(1)
        STOP 
      ENDIF                                                            
C                                
      RO = - COEFB / (2.D0 * COEFA)
C     
C     LIMITATION OF THE VALUE RHO :
C     
      IF(ABS(R02).GE.ABS(R03)) THEN
        ROMAX = 2.D0*R02
      ELSEIF(ABS(R03).GE.ABS(R02)) THEN
        ROMAX = 2.D0*R03
      ENDIF
C     
      IF(ABS(RO).GT.ABS(ROMAX)) THEN
        WRITE(LU,*) 'INTERPOL : LIMIT VALUE OF RHO'
        RO = ROMAX
      ENDIF
C         
C-----------------------------------------------------------------------
C     
      RETURN
      END
