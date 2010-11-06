C                       *****************
                        SUBROUTINE MASKOB
C                       *****************
C
     *(MASKEL,X,Y,IKLE,NELEM,NELMAX,NPOIN,AT,LT)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2    17/08/94      J-M JANIN (LNH) 30 87 72 84
C
C***********************************************************************
C
C     FONCTION : SOUS-PROGRAMME UTILISATEUR
C                SERVANT A ELIMINER FORMELLEMENT DES ELEMENTS
C                DU MAILLAGE, AVEC UN TABLEAU DE MASQUE MASKEL
C
C     POUR UN ELEMENT MASQUE : MASKEL = 0.D0
C     POUR UN ELEMENT NORMAL : MASKEL = 1.D0
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   MASKEL       |<-->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |   X,Y          | -->|  COORDONNEES DES POINTS
C |   IKLE         | -->|  NUMEROS GLOBAUX DES POINTS DES ELEMENTS
C |   NELEM        | -->|  NOMBRE D'ELEMENTS.
C |   NPOIN        | -->|  NOMBRE DE POINTS
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C APPELE PAR: TELMAC
C***********************************************************************
C                                                                      *
C                                                                      *
C***********************************************************************
C
      IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NPOIN,NELMAX,LT
      INTEGER, INTENT(IN) :: IKLE(NELMAX,*)
C
      DOUBLE PRECISION, INTENT(INOUT) :: MASKEL(NELEM)
      DOUBLE PRECISION, INTENT(IN) :: X(NPOIN),Y(NPOIN),AT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
