C                       *****************
                        SUBROUTINE CORSTR
C                       *****************
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.6    17/08/94    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C      FONCTION: CORRECTION DU COEFFICIENT DE FROTTEMENT SUR LE FOND
C                QUAND IL EST VARIABLE EN TEMPS.
C
C      CE SOUS-PROGRAMME EST SIMPLEMENT UN MODELE
C      IL DOIT ETRE REMPLI PAR L'UTILISATEUR
C
C
C 
C
C-----------------------------------------------------------------------
C  EXAMPLE OF POSSIBLE ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    CHESTR      |<-- |  COEFFICIENT DE FROTTEMENT                   |
C |    X,Y         | -->|  COORDONNEE DU MAILLAGE .                    |
C |    NPOIN       | -->|  NOMBRE DE POINTS DU MAILLAGE                |
C |    PRIVE       | -->|  TABLEAU DE TRAVAIL DEFINI DANS PRINCI       |
C |    ZF          | -->|  COTE DU FOND                                |
C |    KFROT       | -->|  LOI DE FROTTEMENT (LINEAIRE,CHEZY,STRICKLER)|
C |    FFON        | -->|  COEFFICIENT DE FROTTEMENT ASSOCIE A LA LOI  |
C |    H           | -->|  HAUTEUR D'EAU.
C |    AT          | -->|  TIME.
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
C
C     C2D: EXAMPLE FOR TELEMAC-2D
C     C3D: EXAMPLE FOR TELEMAC-3D
C
C2D   USE DECLARATIONS_TELEMAC2D
C3D   USE DECLARATIONS_TELEMAC3D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C2D   INTEGER I
C3D   INTEGER I
C
C-----------------------------------------------------------------------
C
C2D   DO I = 1 , NPOIN
C2D     IF(AT.GT.1200.D0) THEN
C2D       CHESTR%R(I) = 40.D0
C2D     ELSE
C2D       CHESTR%R(I) = 60.D0
C2D     ENDIF
C2D   ENDDO 
C
C-----------------------------------------------------------------------
C
C3D   DO I = 1 , NPOIN2
C3D     IF(AT.GT.1200.D0) THEN
C3D       RUGOF%R(I) = 40.D0
C3D     ELSE
C3D       RUGOF%R(I) = 60.D0
C3D     ENDIF
C3D   ENDDO 
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
