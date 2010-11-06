C                       *******************
                        SUBROUTINE READGEO2
C                       *******************
C
     *(NPOIN,NELEM,NPTFR,NDP,IKLES,IPOBO,IB,NFIC)
C
C***********************************************************************
C BIEF VERSION 5.5           29/04/04  J-M HERVOUET (LNH) 01 30 71 80 18
C                                      REGINA NEBAUER
C                                      LAM MINH PHUONG
C***********************************************************************
C
C   FUNCTION: 
C     
C   READS OR COMPUTES THE VALUES OF NPOIN, NELEM, NPTFR.
C   READS THE CONNECTIVITY TABLE AND THE NUMBERING FOR THE BORDER NODES
C
C   MAY BE REWRITTEN FOR ANOTHER FILE FORMAT
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   NPOIN        |--> | NOMBRE DE POINTS DU MAILLAGE.
C |   NELEM        |--> | NOMBRE D'ELEMENTS DU MAILLAGE.
C |   NPTFR        |--> | NOMBRE DE POINTS FRONTIERE DU DOMAINE.
C |   NDP          |--> | NOMBRE DE POINTS PAR ELEMENT 
C |   IKLES        |<-- | TABLE DE CONNECTIVITE 
C |   IPOBO        |<-- | TABLE NUEROS NOEUDS DE BORD 
C |   IB           | -->| 
C |   NFIC         | -->| UNITE LOGIQUE DU FICHIER GEO 
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : LIT 
C
C***********************************************************************
C
      USE BIEF, EX_READGEO2 => READGEO2
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C     
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(OUT) :: NPTFR
      INTEGER, INTENT(IN)  :: NFIC,NPOIN,NELEM,NDP,IB(10)
      INTEGER, INTENT(OUT) :: IKLES(NDP*NELEM)
      INTEGER, INTENT(OUT) :: IPOBO(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION XB(2)
      REAL RB(2)
      INTEGER ISTAT,I
      CHARACTER(LEN=1)  :: CB
C
C-----------------------------------------------------------------------
C
C     DEBUT DU FICHIER DEJA LU
C
C     REWIND NFIC
C
C     6: IKLES (LIKE IKLE BUT INDICES EXCHANGED)
C
      CALL LIT(XB,RB,IKLES,CB,NELEM*NDP,'I ',NFIC,'STD',ISTAT)
C
C     7: IPOBO (CAS DES FICHIERS SANS PARALLELISME)
C
      IF(IB(8).EQ.0.AND.IB(9).EQ.0) THEN
C
        CALL LIT(XB,RB,IPOBO,CB,NPOIN,'I ',NFIC,'STD',ISTAT)
C
        NPTFR = 0
C
        IF(NPOIN.GE.1) THEN
          DO 22 I = 1 , NPOIN
            IF(IPOBO(I).NE.0) NPTFR = NPTFR + 1
22        CONTINUE
        ENDIF
C
      ELSE
C
C       PARALLELISME
C       CAS OU IPOBO EST REMPLACE PAR KNOLG : ON NE LIT PAS IPOBO ICI
C       IL SERA LU DANS READGEO2
C       MAIS IL FAUT CALCULER NPTFR, MXPTVS ET MXELVS
        NPTFR = IB(8)
        IF(NPOIN.GE.1) THEN
          DO 122 I = 1 , NPOIN
            IPOBO(I)=1
122       CONTINUE
        ENDIF
C       IPOBO MIS A UN : MXPTVS SERA SUREVALUE DE 1
C       SINON IL FAUDRAIT RECONSTRUIRE LES VRAIS IPOBO ET DE PLUS
C       REGARDER LES POINTS INTERFACES.
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
