C                       *****************
                        SUBROUTINE SOLAUX
C                       *****************
C
     *(IPT, TB,TBB,ITB,ITBB,S)
C
C***********************************************************************
C BIEF VERSION 5.1           01/02/95  J.M. HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C PROGRAMME AUXILIAIRE POUR LE SOUS-PROGRAMME SOLVE
C
C FONCTION : TB EST UN BLOC DE VECTEURS ET TBB UN BLOC DE BLOCS.
C            SOLAUX PREPARE UN BLOC DE TBB EN LE REMPLISSANT DE
C            MAX(1,S) VECTEURS PRIS DANS TB.
C
C            L'ADRESSE RENVOYEE EST UNE ADRESSE RELATIVE A TBB
C
C            EN PRATIQUE LES OBJETS DE TB SONT DES VECTEURS
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     IPT        |<-- |  POINTEUR DE TBB VERS L'OBJET RENDU
C |     TB         | -->|  BLOC DE TABLEAUX DE TRAVAIL
C |     TBB        | -->|  BLOC DE BLOCS DE TRAVAIL
C |     ITB        | -->|  PREMIER VECTEUR LIBRE DE TB
C |     ITBB       | -->|  PREMIER BLOC LIBRE DE TBB
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES :
C
C**********************************************************************
C
      USE BIEF, EX_SOLAUX => SOLAUX
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: S
      INTEGER, INTENT(INOUT) :: ITB,IPT,ITBB
C
C-----------------------------------------------------------------------
C
C  STRUCTURES DE BLOCS DE TABLEAUX DE TRAVAIL
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: TB,TBB
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K,IAD,MAXTB,MAXTBB
C
C-----------------------------------------------------------------------
C
      MAXTBB = TBB%N
      IF(ITBB.GT.MAXTBB) THEN
        IF(LNG.EQ.1) WRITE(LU,30) ITBB
        IF(LNG.EQ.2) WRITE(LU,31) ITBB
        CALL PLANTE(1)
        STOP
      ENDIF
      IPT = ITBB
      MAXTB = TB%N
      ITBB = ITBB + 1
C     REINITIALISATION DU BLOC
      TBB%ADR(IPT)%P%N=0
      DO 20 K = 1 , MAX(S,1)
        IF(ITB.GT.MAXTB) THEN
          IF(LNG.EQ.1) WRITE(LU,10) ITB + MAX(S,1) - K + 1
          IF(LNG.EQ.2) WRITE(LU,11) ITB + MAX(S,1) - K + 1
          CALL PLANTE(1)
          STOP
        ENDIF
        IAD=ITB
        CALL ADDBLO(TBB%ADR(IPT)%P,TB%ADR(IAD)%P)
        ITB = ITB + 1
20    CONTINUE
C
C-----------------------------------------------------------------------
C
10    FORMAT(1X,'SOLAUX (BIEF) : NOMBRE DE TABLEAUX DE TRAVAIL',/,
     *       1X,'INSUFFISANT DANS LE BLOC TB. MINIMUM NECESSAIRE :',1I4)
11    FORMAT(1X,'SOLAUX (BIEF): INSUFFICIENT NUMBER OF ARRAYS',/,
     *       1X,'IN BLOCK TB. MINIMUM REQUIRED:',1I4)
30    FORMAT(1X,'SOLAUX (BIEF) : NOMBRE DE BLOCS DE TRAVAIL',/,
     *       1X,'INSUFFISANT DANS LE BLOC TBB. MINIMUM NECESSAIRE:',1I4)
31    FORMAT(1X,'SOLAUX (BIEF): INSUFFICIENT NUMBER OF BLOCKS',/,
     *       1X,'IN BLOCK TBB. MINIMUM REQUIRED:',1I4)
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
 
