C                       *****************
                        SUBROUTINE CLTRAC
C                       *****************
C
     *(NWEIRS,NPSING,NPSMAX,NUMDIG,ZF,ZDIG,H,T,NBOR,LITBOR,TBOR,NTRAC)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2     19/04/96     V. GUINOT   (LHF)
C              MODIFIE LE     03/10/96  J.-M. HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C      FONCTION:  GESTION DES CONDITIONS AUX LIMITES DU TRACEUR.
C      =========  SUR LES SEUILS.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   NWEIRS       | -->| NOMBRE DE SINGULARITES LINEIQUES.
C |   NPSING(N)    | -->| NOMBRE DE POINTS DE CHAQUE COTE DE LA
C |                |    | SINGULARITE N.
C |   NPSMAX       | -->| NOMBRE MAXIMUM DE POINTS POUR UN COTE D'UNE
C |                |    | SINGULARITE.
C |   NUMDIG(K,N,I)| -->| NUMERO DES POINTS DES DIGUES
C |                |    | DANS LA NUMEROTATION DES POINTS DE BORD
C |                |    | DES CONDITIONS AUX LIMITES) DU I-EME
C |                |    | POINT SUR LE COTE K DE L'OUVRAGE N
C |   ZF           | -->| COTE DU FOND.
C |   ZDIG         | -->| COTE DES POINTS DES SEUILS.
C |   H            | -->| HAUTEUR D'EAU.
C |   T            | -->| TRACEUR.
C |   NBOR         | -->| NUMEROTATION GLOBALE DES POINTS DE BORD.
C |   LITBOR(J)    |<-- | TYPE DE LA CONDITION EN TRACEUR AU
C |                |    | J-EME POINT LIMITE
C |   TBOR(J)      |<-- | VALEUR DE LA CONDITION EN TRACEUR AU
C |                |    | J-EME POINT LIMITE
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NWEIRS,NPSMAX,NTRAC
      INTEGER, INTENT(IN) :: NPSING(NWEIRS),NUMDIG(2,NWEIRS,NPSMAX)
      INTEGER, INTENT(IN) :: NBOR(*)
      DOUBLE PRECISION, INTENT(IN)    :: ZDIG(NWEIRS,NPSMAX)
      DOUBLE PRECISION, INTENT(IN)    :: ZF(*),H(*)
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: LITBOR,TBOR 
      TYPE(BIEF_OBJ), INTENT(IN)    :: T 
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,N,N1,N2,ITRAC      
C      
      DOUBLE PRECISION Z1,Z2
C
C-----------------------------------------------------------------------
C
      DO ITRAC=1,NTRAC
C
        DO 10 N=1,NWEIRS
        DO 20 I=1,NPSING(N)
C
          N1=NBOR(NUMDIG(1,N,I))
          N2=NBOR(NUMDIG(2,N,I))
          Z1=H(N1)+ZF(N1)
          Z2=H(N2)+ZF(N2)
          IF(Z1.GT.Z2.AND.Z1.GT.ZDIG(N,I)) THEN
            TBOR%ADR(ITRAC)%P%R(NUMDIG(1,N,I))=T%ADR(ITRAC)%P%R(N1)
            TBOR%ADR(ITRAC)%P%R(NUMDIG(2,N,I))=T%ADR(ITRAC)%P%R(N1)
            LITBOR%ADR(ITRAC)%P%I(NUMDIG(1,N,I))=4
            LITBOR%ADR(ITRAC)%P%I(NUMDIG(2,N,I))=5
          ELSEIF(Z2.GE.Z1.AND.Z2.GT.ZDIG(N,I)) THEN
            TBOR%ADR(ITRAC)%P%R(NUMDIG(1,N,I))=T%ADR(ITRAC)%P%R(N2)
            TBOR%ADR(ITRAC)%P%R(NUMDIG(2,N,I))=T%ADR(ITRAC)%P%R(N2)
            LITBOR%ADR(ITRAC)%P%I(NUMDIG(1,N,I))=5
            LITBOR%ADR(ITRAC)%P%I(NUMDIG(2,N,I))=4
          ELSE
            LITBOR%ADR(ITRAC)%P%I(NUMDIG(1,N,I))=2
            LITBOR%ADR(ITRAC)%P%I(NUMDIG(2,N,I))=2
          ENDIF
C
20      CONTINUE
10      CONTINUE
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
