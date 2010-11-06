C                       *****************
                        SUBROUTINE CORPOR
C                       *****************
C
     *(POROS)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.2          01/03/90    J-M HERVOUET
C***********************************************************************
C
C  USER SUBROUTINE CORPOR
C
C  FUNCTION  : MODIFICATION OF THE POROSITY OF ELEMENTS
C
C
C-----------------------------------------------------------------------
C  ARGUMENTS 
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|_______________________________________________
C |      POROS     |<-->| POROSITY TO BE MODIFIED.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT :
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC2D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: POROS
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DOUBLE PRECISION XSOM(4),YSOM(4),XX1,YY1
C     INTEGER NSOM,IELEM
C
C-----------------------------------------------------------------------
C
C     EXAMPLE : POROSITY IS SET TO 0.5 IN A QUADRILATERAL
C
C     NSOM = 4
C     XSOM(1) = 8020.88D0
C     XSOM(2) = 7761.81D0
C     XSOM(3) = 8679.17D0
C     XSOM(4) = 8988.75D0
C     YSOM(1) =-1547.11D0
C     YSOM(2) =-2415.26D0
C     YSOM(3) =-2604.16D0
C     YSOM(4) =-1543.75D0
C                                                                       
C-----------------------------------------------------------------------
C                                                                  
C     CALL OS( 'X=C     ' , POROS , POROS , POROS , 1.D0 )                   
C                                                                       
C--------------------------------------------------------------
C
C     DO 4 IELEM = 1 , NELEM                                         
C                                                                            
C       XX1 = (  X(IKLE%I(IELEM)          )+
C    *           X(IKLE%I(IELEM+NELMAX)   )+
C    *           X(IKLE%I(IELEM+2*NELMAX) ))/3.D0 
C       YY1 = (  Y(IKLE%I(IELEM)          )+
C    *           Y(IKLE%I(IELEM+NELMAX)   )+
C    *           Y(IKLE%I(IELEM+2*NELMAX) ))/3.D0                                     
C
C       IF(INPOLY(XX1,YY1,XSOM,YSOM,NSOM)) THEN
C         POROS%R(IELEM) = 0.5D0
C       ENDIF
C                                      
C4     CONTINUE                                                       
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
