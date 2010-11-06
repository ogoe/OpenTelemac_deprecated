C                       *****************
                        SUBROUTINE PREDES
C                       *****************
C
     *(LLT,AAT)
C
C***********************************************************************
C SISYPHE VERSION 6.0                             E. PELTIER    11/09/95
C                                                 C. LENORMANT
C                                                 J.-M. HERVOUET
C 
C
C JMH 07/12/2009: KS SET TO 0 IF LLT=0
C                                               
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT   
C***********************************************************************
C
C     FONCTION  : PREPARATION DE VARIABLES QUI SERONT ECRITES SUR
C                 LE FICHIER DE RESULTATS OU SUR LE LISTING.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |      LLT       |--> | LOCAL LT (MAY BE LT-1+PERCOU) 
C |      AAT       |--> | CURRENT TIME (FOR BUILDING SOLUTIONS)
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C     - PROGRAMME APPELANT : SISYPH  
C     - SOUS-PROGRAMMES APPELES : OVD,OV
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_SISYPHE
      IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)          :: LLT
      DOUBLE PRECISION, INTENT(IN) :: AAT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C     
      INTEGER LTT,I      
      LOGICAL IMP,LEO
C
C-----------------------------------------------------------------------
C
C     THE OUTPUT VARIABLES ARE BUILT ONLY IF NECESSARY, HENCE THE
C     FOLLOWING TESTS, WHICH MUST BE THE SAME THAN IN DESIMP (BIEF LIBRARY)
C
      IMP=.FALSE.
      LEO=.FALSE.
      LTT=(LLT/LISPR)*LISPR
      IF(LLT.EQ.LTT.AND.LLT.GE.PTINIL) IMP=.TRUE.
      LTT=(LLT/LEOPR)*LEOPR
      IF(LLT.EQ.LTT.AND.LLT.GE.PTINIG) LEO=.TRUE.
C
C     PAS D'IMPRESSION, PAS DE SORTIE SUR FICHIER, ON RESSORT
      IF (.NOT.(LEO.OR.IMP)) GO TO 1000
C
C=======================================================================
C     COMPUTING SECONDARY VARIABLES
C=======================================================================
C
C     FREE SURFACE: H+ZF
C
      IF((LEO.AND.SORLEO(4)).OR.(IMP.AND.SORIMP(4))) THEN
        CALL OS('X=Y+Z   ',X=Z,Y=HN,Z=ZF)
      ENDIF
C
C     DISCHARGE 
C
      IF((LEO.AND.SORLEO(6)).OR.(IMP.AND.SORIMP(6))) THEN
        DO I=1,NPOIN
          Q%R(I)=HN%R(I)*SQRT(U2D%R(I)**2+V2D%R(I)**2)
        ENDDO
      ENDIF
C
C     DISCHARGE ALONG X
C
      IF((LEO.AND.SORLEO(7)).OR.(IMP.AND.SORIMP(7))) THEN
        CALL OS('X=YZ    ',X=QU,Y=U2D,Z=HN)
      ENDIF
C
C     DISCHARGE ALONG Y
C
      IF((LEO.AND.SORLEO(8)).OR.(IMP.AND.SORIMP(8))) THEN
        CALL OS('X=YZ    ',X=QU,Y=V2D,Z=HN)
      ENDIF
C
C=======================================================================
C
C     VARIABLES WHICH ARE NOT INITIALISED AT THE FIRST CALL OF PREDES
C
      IF(LLT.EQ.0) THEN
C       JMH ON 27/11/2009
        IF((LEO.AND.SORLEO(19)).OR.(IMP.AND.SORIMP(19))) THEN
          CALL OS('X=0     ',X=KS)
        ENDIF
      ENDIF
C
C=======================================================================
C
1000  CONTINUE
C
C=======================================================================
C
      RETURN
      END
