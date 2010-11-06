C                       *****************
                        SUBROUTINE STRCHE
C                       *****************
C
C***********************************************************************
C  BIEF VERSION 5.2           01/10/96    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C      FONCTION : CALCUL DU COEFFICIENT DE FROTTEMENT SUR LE FOND.
C                 SI IL  EST VARIABLE EN ESPACE.
C
C      CE SOUS-PROGRAMME EST SIMPLEMENT UN MODELE.
C      IL DOIT ETRE REMPLI PAR L'UTILISATEUR.
C
C      NOTE : EN IDENTIFICATION DE PARAMETRES AVEC PLAN D'EXPERIENCE
C             LES VALEURS DONNEES ICI NE SONT PAS PRISES EN COMPTE.
C
C-----------------------------------------------------------------------
C
C      FUNCTION: SETTING THE FRICTION COEFFICIENT IF IT IS VARIABLE
C                IN SPACE.
C
C      ONLY AN EXAMPLE IS GIVEN HERE, MUST BE IMPLEMENTED BY THE USER.
C      COMMENTS CEX MUST BE REMOVED TO IMPLEMENT THE EXAMPLE.
C
C      NOTE : IN PARAMETER ESTIMATION WITH A LIST OF TESTS, VALUES
C             GIVEN HERE ARE DISCARDED.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
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
C |    MESH        | -->|  BLOC DES ENTIERS DU MAILLAGE.
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C
C**********************************************************************
C
      USE BIEF
C
C     DECLARATIONS MUST BE ADAPTED TO EVERY CODE
C     EXAMPLE OF TELEMAC2D HERE
C
CEX   USE DECLARATIONS_TELEMAC2D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
CEX   INTEGER I
C
C-----------------------------------------------------------------------
C
C     HERE A CONSTANT VALUE OF FRICTION IS GIVEN
C
CEX   DO I=1,NPOIN
CEX     CHESTR%R(I) = 60.D0
CEX   ENDDO
C
C-----------------------------------------------------------------------
C
C     COMMENTS HERE MAY BE CHANGED
C
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'STRCHE (BIEF) : PAS DE MODIFICATION DU FROTTEMENT'
        WRITE(LU,*)
      ENDIF
      IF(LNG.EQ.2) THEN
        WRITE(LU,*) 'STRCHE (BIEF): NO MODIFICATION OF FRICTION'
        WRITE(LU,*)
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END          
