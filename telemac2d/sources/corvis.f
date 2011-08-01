C                       *****************
                        SUBROUTINE CORVIS
C                       *****************
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2    17/08/94    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C      FONCTION: CORRECTION DU COEFFICIENT DE VISCOSITE
C
C      CE SOUS-PROGRAMME EST SIMPLEMENT UN MODELE
C      IL DOIT ETRE REMPLI PAR L'UTILISATEUR
C
C      ATTENTION, CE PROGRAMME EST APPELE A CHAQUE PAS DE TEMPS
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C      FUNCTION: CORRECTING THE DIFFUSION COEFFICIENT
C
C      HERE IS JUST AN EXAMPLE, THIS SUBROUTINE MAY BE MODIFIED
C      BY A USER
C
C      BEWARE, IT IS CALLED AT EVERY TIME STEP
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |                |    |    
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : TELMAC
C
C  SOUS-PROGRAMME APPELE :
C
C**********************************************************************
C
      USE BIEF
CEX   USE DECLARATIONS_TELEMAC2D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C-----------------------------------------------------------------------
C
C     EXAMPLE : THE VALUE 0.1 IS SET INSIDE A SQUARE.
C               REMOVE THE COMMENTS CEX TO IMPLEMENT IT.
C
C
CEX   INTEGER NSOM
CEX   DOUBLE PRECISION XSOM(10),YSOM(10)
C
C-----------------------------------------------------------------------
C
C     DESCRIPTION OF THE SQUARE AS A POLYGON TO CALL ROUTINE FILPOL
C
CEX   NSOM = 4
CEX   XSOM(1) =     0.D0
CEX   YSOM(1) =     0.D0
CEX   XSOM(2) =  1000.D0
CEX   YSOM(2) =     0.D0
CEX   XSOM(3) =  1000.D0
CEX   YSOM(3) =  1000.D0
CEX   XSOM(4) =     0.D0
CEX   YSOM(4) =  1000.D0
C
CEX   CALL FILPOL( VISC , 0.1D0 , XSOM , YSOM , NSOM , MESH )
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
