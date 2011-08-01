C                       ***************************
                        SUBROUTINE SOURCE_TELEMAC2D
C                       ***************************
C
C***********************************************************************
C TELEMAC 2D VERSION 5.2       26/10/94    J-M HERVOUET LNH 30 87 80 18
C***********************************************************************
C
C  FONCTION  : SOUS-PROGRAMME UTILISATEUR POUR REDEFINIR LES
C
C              CARACTERISTIQUES DES SOURCES SANS PASSER PAR LE
C              FICHIER DES PARAMETRES.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |
C |   SOME USEFUL DATA INCLUDED IN MODULE DECLARATIONS_TELEMAC2D:
C |
C |   NREJET       |<-- | NOMBRE DE SOURCES/PUITS
C |   XSCE,YSCE    |<-- | COORDONNEES DES SOURCES.
C |   DSCE         |<-- | DEBITS DES SOURCES (ENTRANT POSITIF)
C |   TSCE         |<-- | VALEURS DU TRACEUR AUX SOURCES (FACULTATIF)
C |                |    | TSCE N'EST PRIS EN COMPTE QUE SI LE DEBIT
C |                |    | DE LA SOURCE EST ENTRANT (POSITIF).
C |   USCE,VSCE    |<-- | VITESSES AUX SOURCES.
C |   ZF           | -->| FOND.
C |   X,Y          | -->| COORDONNEES DU MAILLAGE.
C |   IKLE         | -->| NUMEROS DES NOEUDS DE CHAQUE ELEMENT.
C |   NELEM        | -->| NOMBRE D'ELEMENTS
C |   NELMAX       | -->| NOMBRE MAXIMUM D'ELEMENTS.
C |   NPOIN        | -->| NOMBRE DE POINTS DU MAILLAGE.
C |   PRIVE        | -->| TABLEAU DE TRAVAIL DE L'UTILISATEUR.
C |   MAXSCE       | -->| TAILLE MAXI DES TABLEAUX QUI FINISSENT EN SCE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : TELEMAC2D
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF
C     USE DECLARATIONS_TELEMAC2D
C
      IMPLICIT NONE
C
      INTEGER I
C
C-----------------------------------------------------------------------
C
C     EXAMPLE
C
C     NREJET = ?  (UNTIL MAXSCE)
C     NREJEU = NREJET (IF VELOCITIES GIVEN)
C
C     DO I=1,NREJET
C
C       XSCE(I) = ???
C       YSCE(I) = ???
C       DSCE(I) = ???
C       TSCE(I) = ???
C       USCE(I) = ???
C       VSCE(I) = ???
C
C     ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
