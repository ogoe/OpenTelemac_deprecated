C                       *****************
                        SUBROUTINE INTERP
C                       *****************
C     
     * ( U , UTILD , SHP , NDP , SHZ , ETA , ELT , NP , NPOIN2 , NPLAN ,
     *   IELM , IKLE , NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.1           28/04/93     J-M JANIN (LNH) 30 87 72 84
C
C***********************************************************************
C
C   FONCTION : INTERPOLATION DE VALEURS D'UNE FONCTION EN CERTAINS
C              POINTS D'UN MAILLAGE EN FONCTION DES COORDONNEES
C              BARYCENTRIQUES DES POINTS ET DES VALEURS AUX NOEUDS
C              DE LA FONCTION.
C
C   ATTENTION : NE MARCHE PAS SI LES COORDONNEES BARYCENTRIQUES
C               FOURNIES NE CORRESPONDENT PAS A L'ELEMENT DE LA
C               FONCTION.
C
C               LES ELEMENTS AUTRES QUE 11 , 21 ET 41 NE SONT DONC
C               PAS PREVUS ICI.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   U...         | -->| VALEURS AUX NOEUDS.
C |   UTILD...     |<-- | VALEURS INTERPOLEES.
C |   SHP          | -->| COORDONNEES BARYCENTRIQUES 2D AU PIED DES
C |                |    | COURBES CARACTERISTIQUES.
C |   NDP          | -->| NOMBRE DE POINTS PAR ELEMENT POUR LA FONCTION
C |                |    | U
C |   SHZ          | -->| COORDONNEES BARYCENTRIQUES SUIVANT Z AU PIED
C |                |    | DES COURBES CARACTERISTIQUES (POUR TEL3D)
C |   ELT          | -->| NUMEROS DES ELEMENTS 2D AU PIED DES COURBES
C |                |    | CARACTERISTIQUES.
C |   NP           | -->| NOMBRE DE POINTS A INTERPOLER
C |   NPOIN2       | -->| NOMBRE DE POINTS EN 2D (MEME POUR LE 3D)
C |   NPLAN        | -->| NOMBRE DE PLANS EN 3D
C |   IELM         | -->| TYPE D'ELEMENT.
C |   IKLE         | -->| TABLE DE CONNECTIVITE.
C |   NELMAX       | -->| NOMBRE MAXIMUM D'ELEMENTS.
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : CARACT
C
C***********************************************************************
C
      USE BIEF, EX_INTERP => INTERP
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NP,NELMAX,NPLAN,NPOIN2,NDP,IELM
      INTEGER, INTENT(IN) :: IKLE(NELMAX,*),ELT(NP),ETA(NP)
C
      DOUBLE PRECISION, INTENT(IN)    :: U(NPOIN2,NPLAN)
      DOUBLE PRECISION, INTENT(IN)    :: SHP(NDP,NP),SHZ(NP)
      DOUBLE PRECISION, INTENT(INOUT) :: UTILD(NP)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IP
C
C-----------------------------------------------------------------------
C
      IF(IELM.EQ.11.OR.IELM.EQ.12) THEN
C
C    TRIANGLES P1
C    ============
C
      DO IP = 1 , NP
         UTILD(IP) = U(IKLE(ELT(IP),1),1) * SHP(1,IP)
     *             + U(IKLE(ELT(IP),2),1) * SHP(2,IP)
     *             + U(IKLE(ELT(IP),3),1) * SHP(3,IP)
      ENDDO
C
C-----------------------------------------------------------------------
C
      ELSEIF(IELM.EQ.13) THEN
C
C    TRIANGLES P2
C    ============
C  
      DO IP = 1 , NP
         UTILD(IP) = U(IKLE(ELT(IP),1),1) *
     *               (2.D0*SHP(1,IP)-1.D0)* SHP(1,IP)
     *             + U(IKLE(ELT(IP),2),1) * 
     *               (2.D0*SHP(2,IP)-1.D0)* SHP(2,IP)
     *             + U(IKLE(ELT(IP),3),1) *
     *               (2.D0*SHP(3,IP)-1.D0)* SHP(3,IP)
     *             + U(IKLE(ELT(IP),4),1) * 4.D0 * SHP(1,IP)*SHP(2,IP)
     *             + U(IKLE(ELT(IP),5),1) * 4.D0 * SHP(2,IP)*SHP(3,IP)
     *             + U(IKLE(ELT(IP),6),1) * 4.D0 * SHP(3,IP)*SHP(1,IP)
      ENDDO     
C
C------------------------------------------------------------------------
C
      ELSEIF(IELM.EQ.41) THEN
C
C    PRISMES DE TELEMAC-3D
C    =====================
C
      DO IP = 1 , NP
         UTILD(IP) =
     *     U(IKLE(ELT(IP),1),ETA(IP))   * SHP(1,IP) * (1.D0-SHZ(IP))
     *   + U(IKLE(ELT(IP),2),ETA(IP))   * SHP(2,IP) * (1.D0-SHZ(IP))
     *   + U(IKLE(ELT(IP),3),ETA(IP))   * SHP(3,IP) * (1.D0-SHZ(IP))
     *   + U(IKLE(ELT(IP),1),ETA(IP)+1) * SHP(1,IP) * SHZ(IP)
     *   + U(IKLE(ELT(IP),2),ETA(IP)+1) * SHP(2,IP) * SHZ(IP)
     *   + U(IKLE(ELT(IP),3),ETA(IP)+1) * SHP(3,IP) * SHZ(IP)
      ENDDO
C
C-----------------------------------------------------------------------
C
      ELSE
C
        IF(LNG.EQ.1) WRITE(LU,11) IELM
        IF(LNG.EQ.2) WRITE(LU,12) IELM
11      FORMAT(1X,'INTERP : TYPE D''ELEMENT INCONNU : ',I6)
12      FORMAT(1X,'INTERP : UNKNOWN TYPE OF ELEMENT : ',I6)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
