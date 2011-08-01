C                       *****************                                                                                         
                        SUBROUTINE METGRA                             
C                       *****************                                        
C     
     *(RO,ESTIME,GRADJ,GRADJN,JCOUT1,DESC,NPARAM,OPTID,RSTART,R02,R03)
C     
C***********************************************************************
C     PROGICIEL : TELEMAC 2D        02/08/93  E. BARROS
C             UPGRADE TO 5.2        04/10/00  A. LEOPARDI (UNINA)
C***********************************************************************
C     
C     FUNCTION: ONE STEP OF GRADIENT METHOD                                                   
C     
C-----------------------------------------------------------------------
C     ARGUMENTS                                         
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |    GRADJ       | -->| GRADIENT OF COST FUNCTION (ITERATION K)
C |    GRADJN      | -->| GRADIENT OF COST FUNCTION (ITERATION K-1)
C |    RO          |<-->| COEFFICIENT OF THE GRADIENT
C |    NPARAM      | -->| TOTAL NUMBER OF PARAMETERS TO ESTIMATE
C |    DESC        |<-- | VECTOR USED TO CHANGE THE SET OF STRICKLERS'
C |    OPTID       | -->| METHOD 1=GRADIENT, 2=GRADIENT CONJUGUE, 3=LAGRANGE)
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C     
C     APPELE PAR :            HOMERE_PIT
C     
C     SOUS-PROGRAMME APPELE : OS
C     
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
      INTEGER , INTENT(IN)             :: NPARAM,OPTID
      CHARACTER(LEN=72)                :: ESTIME
      DOUBLE PRECISION , INTENT(IN)    :: JCOUT1
      LOGICAL , INTENT(IN)             :: RSTART
      TYPE(BIEF_OBJ) , INTENT(IN)      :: GRADJ,GRADJN
      TYPE(BIEF_OBJ) , INTENT(INOUT)   :: DESC
      DOUBLE PRECISION , INTENT(INOUT) :: R02,R03,RO
C      
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C    
      INTEGER I
C     
      DOUBLE PRECISION R1,DENOM,GRAD_JN
C     
C     CALCUL DU VERITABLE GRADIENT (QUI TIENT COMPTE DU VERITABLE
C                                                  NOMBRE DE PARAMETRES)     
      DENOM=0.D0
      GRAD_JN=0.D0
      DO I = 1,NPARAM
        DENOM  = DENOM + GRADJ%R(I)**2
        GRAD_JN=GRAD_JN+GRADJN%R(I)**2 
      ENDDO
C
      IF(DENOM.LT.1.D-12) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'METGRA : GRADIENT TROP PETIT, ARRET'
        IF(LNG.EQ.2) WRITE(LU,*) 'METGRA: GRADIENT TOO SMALL, STOP'
        WRITE(LU,*) 'DENOM = ',DENOM
        CALL PLANTE(1)
        STOP
      ENDIF
C           
C-----------------------------------------------------------------------
C     RO = - JCOUT / GRADJ*GRADJ
C-----------------------------------------------------------------------
C
      IF(OPTID.EQ.1.OR.OPTID.EQ.3.OR.RSTART) THEN
C
            R02 = - JCOUT1 / DENOM
            RO = R02
            R03=0.5D0*R02
C     
C           CALCUL DE LA DIRECTION DE DESCENTE INITIALE
C     
            CALL OV('X=Y     ',DESC%R,GRADJ%R,GRADJ%R,0.D0,NPARAM)
C     
C-----------------------------------------------------------------------
C     
      ELSEIF(OPTID.EQ.2) THEN
C
            R02 = - JCOUT1 / DENOM
C      
            R1 = GRAD_JN/DENOM
C     
C           CALCUL DE LA DIRECTION DE DESCENTE
C
            CALL OV('X=Y+CZ  ',DESC%R,GRADJ%R,DESC%R,R1,NPARAM)
C     
            DENOM=0.D0
            DO I=1,NPARAM
               DENOM=DENOM+GRADJ%R(I)*DESC%R(I)
            ENDDO      
            R03 = - JCOUT1/DENOM
            RO =R03
C
      ENDIF
C     
C-----------------------------------------------------------------------
C
      RETURN
      END
