C                       *****************
                        SUBROUTINE GESTIO
C                       *****************
C
     *(U,V,C,T,AK,EP,UTILD,VTILD,CTILD,TTILD,AKTILD,EPTILD,
     * TRAC,PROPA,CONVV,ITURB,IETAPE)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2    17/08/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C      FONCTION:    CE SOUS-PROGRAMME GERE LE CHARGEMENT DES TABLEAUX
C                   SUIVANT LES EQUATIONS CHOISIES.
C
C      ATTENTION : NE MARCHE PLUS SI ON FAIT D'AUTRES SCHEMAS DE
C                  CONVECTION QUE LES CARACTERISTIQUES.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   U,V          |<-->| VITESSE A LA FIN DE L'ETAPE TRAITEE.         |
C |   C            |<-->| CELERITE A LA FIN DE L'ETAPE TRAITEE.        |
C |   T            |<-->| TRACEUR A LA FIN DE L'ETAPE TRAITEE.         |
C |   AK           |<-->| ENERGIE TURBULENTE A LA FIN DE L'ETAPE TRAITEE
C |   EP           |<-->| DISSIPASSION A FIN DE L'ETAPE TRAITEE.       |
C |   UTILD,VTILD  | -->| VITESSE AVANT L'ETAPE TRAITEE.               |
C |   CTILD        | -->| CELERITE AVANT L'ETAPE TRAITEE.              |
C |   TTILD        | -->| TRACEUR AVANT L'ETAPE TRAITEE.               |
C |   AKTILD       | -->| ENERGIE TURBULENTE AVANT L'ETAPE TRAITEE     |
C |   EPTILD       | -->| DISSIPASSION AVANT L'ETAPE TRAITEE.          |
C |   TRAC         | -->| LOGIQUE INDIQUANT LA PRESENCE D'UN TRACEUR   |
C |   PROPA        | -->| SI PROPA=.FALSE. : PAS DE PROPAGATION.       |
C |   CONVV        | -->| LOGIQUES INDIQUANT LES VARIABLES QU'ON NE    |
C |                |    | VEUT PAS CONVECTER.                          |
C |   ITURB        | -->| MODELE DE TURBULENCE  1 : LAMINAIRE          |
C |                |    |                       2 : LONGUEUR DE MELANGE|
C |                |    |                       3 : K-EPSILON          |
C |   IETAPE       | -->| INDICATEUR D'AVANCEMENT DANS LE PROGRAMME .  |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMME APPELE : OV
C
C**********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)           :: ITURB,IETAPE
      LOGICAL, INTENT(IN)           :: TRAC,CONVV(4),PROPA
      TYPE(BIEF_OBJ), INTENT(IN)    :: T,AK,EP
      TYPE(BIEF_OBJ), INTENT(INOUT) :: U,V,C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: UTILD,VTILD,CTILD,TTILD
      TYPE(BIEF_OBJ), INTENT(INOUT) :: AKTILD,EPTILD
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C----------------------------------------------------------------------
C   CHARGEMENT DES TABLEAUX SUIVANT LES EQUATION CHOISIES
C-----------------------------------------------------------------------
C
C    CONVECTION
C
      IF(IETAPE.EQ.3) THEN
C
        IF(.NOT.CONVV(1)) THEN
          CALL OS( 'X=Y     ' , X=UTILD , Y=U )
          CALL OS( 'X=Y     ' , X=VTILD , Y=V )
        ENDIF
        IF(.NOT.CONVV(2)) THEN
          CALL OS( 'X=Y     ' , X=CTILD , Y=C )
        ENDIF
        IF(TRAC.AND.(.NOT.CONVV(3))) THEN
          CALL OS( 'X=Y     ' , X=TTILD , Y=T )
        ENDIF
        IF(ITURB.EQ.3.AND.(.NOT.CONVV(4))) THEN
          CALL OS( 'X=Y     ' , X=AKTILD , Y=AK )
          CALL OS( 'X=Y     ' , X=EPTILD , Y=EP )
        ENDIF
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C    PROPAGATION
C
      IF(IETAPE.EQ.6) THEN
C
            IF(.NOT.PROPA) THEN
C
                   CALL OS( 'X=Y     ' , X=U , Y=UTILD )
                   CALL OS( 'X=Y     ' , X=V , Y=VTILD )
                   CALL OS( 'X=Y     ' , X=C , Y=CTILD )
C
            ENDIF
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
