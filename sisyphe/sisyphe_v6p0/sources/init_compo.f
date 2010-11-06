C                         *********************
                          SUBROUTINE INIT_COMPO
C                         *********************
C
     *(NCOUCHES)
C
C***********************************************************************
C SISYPHE VERSION 6.0
C                             Matthieu GONZALES DE LINARES 2002
C
C                                                
C COPYRIGHT EDF-BAW-IFH   
C***********************************************************************
C
C     FONCTION  : DISTRIBUTION DES CLASSES
C                 % PAR COUCHE, STRATIFICATION 
C     SUBROUTINE A REMPLIR PAR l'UTILISATEUR
C
C 
C     FUNCTION  : INITIAL FRACTION DISTRIBUTION, STRATIFICATION, 
C                 VARIATION IN SPACE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |                |    |  
C |    AVAIL       |<-- | SEDIMENT FRACTION FOR EACH LAYER, CLASS AND NODE
C |                |    | AVAIL(NPOIN,10,NSICLA)
C |    ES          |<-- | THICKNESS FOR EACH LAYER AND NODE ES(NPOIN,10)
C |    NCOUCHES    |--> | NUMBER OF LAYER FOR EACH POINT
C |    NSICLA      |--> | NUMBER OF SIZE-CLASSES OF BED MATERIAL
C |                |    | (LESS THAN 10)
C |    NPOIN       |--> | NUMBER OF NODES
C |________________|____|______________________________________________
C MODE : -->(INPUT), <--(RESULT), <--> (MODIFIED INPUT)
C-----------------------------------------------------------------------
C PROGRAMME APPELANT : INIT_AVAI 
C PROGRAMMES APPELES : NONE
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_SISYPHE
C
      IMPLICIT NONE  
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C                                       NPOIN
      INTEGER, INTENT (INOUT)::NCOUCHES(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I , J  
C 
C-----------------------------------------------------------------------
C
      DO J=1,NPOIN
C
C       BY DEFAULT : UNIFORM BED COMPOSITION
C
          NCOUCHES(J) = 1
          DO I = 1, NSICLA
            AVAIL(J,1,I) = AVA0(I)
            AVAIL(J,2,I) = AVA0(I)
          ENDDO
C      
C  TO BE FILLED BY THE USER
!      NCOUCHES(J) = 10
!      ES(J,1) = 1.D0
!      ES(J,2) = 1.D0
!      ES(J,3) = 1.D0
!      ES(J,4) = 1.D0
!      ES(J,5) = 1.D0
!      ES(J,6) = 1.D0
!      ES(J,7) = 1.D0
!      ES(J,8) = 1.D0
!      ES(J,9) = 1.D0      
!        DO I = 1, NSICLA
!          DO K = 1, NCOUCHES(J)
!          AVAIL(J,K,I) = AVA0(I)
!          ENDDO
!        ENDDO
C          
      ENDDO 
C 
C-----------------------------------------------------------------------
C
      RETURN
      END
