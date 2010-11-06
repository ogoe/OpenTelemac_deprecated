C                       *****************
                        SUBROUTINE NOEROD
C                       *****************
C
     * (H , ZF , ZR , Z , X , Y , NPOIN , CHOIX , NLISS )
C
C***********************************************************************
C SISYPHE VERSION 5.1                             C. LENORMANT
C                                                
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT   
C***********************************************************************
C
C     FONCTION  : IMPOSE LA VALEUR DE LA COTE DU FOND NON ERODABLE  ZR
C
C
C     RQ: LES METHODES DE TRAITEMENT DES FONDS NON ERODABLES PEUVENT CONDUIRE
C     A ZF < ZR A CERTAINS PAS DE TEMPS, POUR PALLIER A CELA ON PEUT CHOISIR 
C     CHOISIR DE LISSER LA SOLUTION OBTENUE i.e NLISS > 0.  
C
C     FUNCTION  : IMPOSE THE RIGID BED LEVEL  ZR
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   H            | -->| WATER DEPTH
C |   ZF           | -->| BED LEVEL
C |   ZR           |<-- | RIGID BED LEVEL
C |   Z            | -->| FREE SURFACE 
C |   X,Y          | -->| 2D COORDINATES
C |   NPOIN        | -->| NUMBER OF 2D POINTS
C |   CHOIX        | -->| SELECTED METHOD FOR THE TREATMENT OF RIGID BEDS
C |   NLISS        |<-->| NUMBER OF SMOOTHINGS
C |________________|____|______________________________________________
C MODE : -->(INPUT), <--(RESULT), <-->(MODIFIED DATA)
C-----------------------------------------------------------------------
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN):: NPOIN , CHOIX
      INTEGER, INTENT(INOUT):: NLISS 
C
      DOUBLE PRECISION, INTENT(IN)::  Z(NPOIN) , ZF(NPOIN)    
      DOUBLE PRECISION , INTENT(IN)::  X(NPOIN) , Y(NPOIN), H(NPOIN)
      DOUBLE PRECISION , INTENT(INOUT)::  ZR(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C
C--------------------
C RIGID BEDS POSITION
C---------------------
C
C     DEFAULT VALUE:       ZR=ZF-100 
C                                                              
      CALL OV( 'X=C       ',ZR,ZF,ZF,-100.D0,NPOIN)                                                    
C
C------------------
C SMOOTHING OPTION
C------------------
C
C     NLISS : NUMBER OF SMOOTHING IF  (ZF - ZR ) NEGATIVE
C             DEFAULT VALUE : NLISS = 0 (NO SMOOTHING)
C
      NLISS = 0        
C
C-----------------------------------------------------------------------
C
      RETURN
      END
