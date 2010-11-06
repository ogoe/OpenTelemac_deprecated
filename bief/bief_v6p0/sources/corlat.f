C                       *****************
                        SUBROUTINE CORLAT
C                       *****************
C
C***********************************************************************
C BIEF VERSION 5.1          01/03/90    J-M HERVOUET
C***********************************************************************
C
C  USER SUBROUTINE CORLAT
C
C  FUNCTION  : MODIFICATION OF THE LATITUDE OF THE POINTS IN THE MESH
C
C              LINES WITH AN INITIAL CEX ARE AN EXAMPLE
C              WITH TELEMAC-2D
C
C              THIS SUBROUTINE MUST BE MODIFIED ACCORDING TO
C              THE CALLING PROGRAM AND THE NEEDED MODIFICATION
C              BY ADDING USE DECLARATIONS_"NAME OF CALLING CODE"
C              ALL THE DATA STRUCTURE OF THIS CODE IS
C              AVAILABLE
C
C-----------------------------------------------------------------------
C  ARGUMENTS USED IN THE EXAMPLE 
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|_______________________________________________
C |    X,Y         | -->|  COORDONNEES DU MAILLAGE .                   |
C |    NPOIN       | -->|  NOMBRE DE POINTS DU MAILLAGE                |
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT :
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
C APPELE PAR : INBIEF
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      USE BIEF, EX_CORLAT => CORLAT
CEX   USE DECLARATIONS_TELEMAC2D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
CEX   INTEGER I
C
C-----------------------------------------------------------------------
C
C  EXAMPLE : MULTIPLICATION BY A CONSTANT
C
CEX   DO I = 1 , NPOIN
CEX     X(I) = X(I) * 1.D0
CEX     Y(I) = Y(I) * 1.D0
CEX   ENDDO
C
C-----------------------------------------------------------------------
C
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'CORLAT : PAS DE MODIFICATION DE LA LATITUDE'
      ENDIF
      IF(LNG.EQ.2) THEN
        WRITE(LU,*) 'CORLAT :NO MODIFICATION OF LATITUDE'
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
 
