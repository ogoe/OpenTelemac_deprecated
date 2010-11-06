C                       *****************
                        SUBROUTINE GTSH41
C                       *****************
C
     *(U,V,WS,X,Y,SHP,SHZ,ELT,ETA,IKLE2,INDIC,NLOC,NPOIN2,NELEM2,NPLAN,
     * LV,MSK,MASKEL)
C
C  A ENLEVER : U,V,X,Y,NLOC,INDIC
C
C***********************************************************************
C BIEF VERSION 5.9           21/08/2008    J-M JANIN (LNH) 30 87 72 84
C
C***********************************************************************
C
C      FONCTION:
C
C   - FIXE, POUR LES PRISMES DE TELEMAC-3D ET,
C     AVANT LA REMONTEE DES COURBES CARACTERISTIQUES,
C     LES COORDONNEES BARYCENTRIQUES DE TOUS LES NOEUDS DU
C     MAILLAGE DANS L'ELEMENT VERS OU POINTE CETTE COURBE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    U,V         | -->| COMPOSANTES DE LA VITESSE                    |
C |    X,Y         | -->| COORDONNEES DES POINTS DU MAILLAGE.          |
C |    SHP         |<-- | COORDONNEES BARYCENTRIQUES DES NOEUDS DANS   |
C |                |    | LEURS ELEMENTS "ELT" ASSOCIES.               |
C |    ELT         |<-- | NUMEROS DES ELEMENTS CHOISIS POUR CHAQUE     |
C |                |    | NOEUD.                                       |
C |    IKLE        | -->| TRANSITION ENTRE LES NUMEROTATIONS LOCALE    |
C |                |    | ET GLOBALE.                                  |
C |    NPOIN       | -->| NOMBRE DE POINTS.                            |
C |    NELEM       | -->| NOMBRE D'ELEMENTS.                           |
C |    NELMAX      | -->| NOMBRE MAXIMAL D'ELEMENTS DANS LE MAILLAGE 2D|
C |    NPTFR       | -->| NOMBRE DE POINTS FRONTIERES.                 |
C |    MSK         | -->| SI OUI, PRESENCE D'ELEMENTS MASQUES.         |
C |    MASKEL      | -->| TABLEAU DE MASQUAGE DES ELEMENTS             |
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE.        |
C |    MASKPT      | -->| TABLEAU DE MASQUAGE DES POINTS               |
C |                |    |  =1. : NORMAL   =0. : POINT MASQUE.          |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : CARACT
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NPOIN2,NELEM2,NPLAN,LV
      INTEGER, INTENT(IN)    :: IKLE2(NELEM2,3)
      INTEGER, INTENT(INOUT) :: ELT(NPOIN2,NPLAN),ETA(NPOIN2,NPLAN)
      INTEGER, INTENT(INOUT) :: INDIC(NPOIN2),NLOC(NPOIN2)
C
      DOUBLE PRECISION, INTENT(IN) :: U(NPOIN2,NPLAN),V(NPOIN2,NPLAN)
      DOUBLE PRECISION, INTENT(IN) :: WS(NPOIN2,NPLAN)
      DOUBLE PRECISION, INTENT(IN) :: X(NPOIN2),Y(NPOIN2),MASKEL(NELEM2)
      DOUBLE PRECISION, INTENT(INOUT) :: SHP(3,NPOIN2,NPLAN)
      DOUBLE PRECISION, INTENT(INOUT) :: SHZ(NPOIN2,NPLAN)
C
      LOGICAL, INTENT(IN) :: MSK
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IPLAN,IPOIN,IELEM,I1,I2,I3
C
C***********************************************************************
C
C     LOOP ON ALL POINTS
C
      DO IPLAN = 1,NPLAN
        DO IELEM=1,NELEM2
          I1=IKLE2(IELEM,1)
          I2=IKLE2(IELEM,2)
          I3=IKLE2(IELEM,3)
          ELT(I1,IPLAN)=IELEM
          ELT(I2,IPLAN)=IELEM
          ELT(I3,IPLAN)=IELEM
          SHP(1,I1,IPLAN)=1.D0
          SHP(2,I1,IPLAN)=0.D0
          SHP(3,I1,IPLAN)=0.D0
          SHP(1,I2,IPLAN)=0.D0
          SHP(2,I2,IPLAN)=1.D0
          SHP(3,I2,IPLAN)=0.D0
          SHP(1,I3,IPLAN)=0.D0
          SHP(2,I3,IPLAN)=0.D0
          SHP(3,I3,IPLAN)=1.D0
        ENDDO
      ENDDO
C
C     ON THE VERTICAL, IT IS DONE DEPENDING ON THE VERTICAL VELOCITY
C
      DO IPLAN = 1,NPLAN
        DO IPOIN=1,NPOIN2
          IF((WS(IPOIN,IPLAN).GT.0.D0.AND.IPLAN.NE.1).OR.
     *                                              IPLAN.EQ.NPLAN) THEN
            ETA(IPOIN,IPLAN) = IPLAN-1
            SHZ(IPOIN,IPLAN) = 1.D0
          ELSE
            ETA(IPOIN,IPLAN) = IPLAN
            SHZ(IPOIN,IPLAN) = 0.D0
          ENDIF
        ENDDO
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
