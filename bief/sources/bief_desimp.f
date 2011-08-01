C                       **********************
                        SUBROUTINE BIEF_DESIMP
C                       **********************
C
     *(FORMAT_RES,VARSOR,HIST,NHIST,N,NRES,STD,AT,LT,LISPRD,LEOPRD,
     * SORLEO,SORIMP,MAXVAR,TEXTE,PTINIG,PTINIL)
C
C***********************************************************************
C BIEF VERSION 6.0     01/04/2009     J-M HERVOUET (LNHE) 01 30 71 80 18
C
C***********************************************************************
C
C     FONCTION  : ECRITURE DE RESULTATS SUR FICHIER OU LISTING
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   VARSOR       | -->| BLOC CONTENANT LES VARIABLES A METTRE DANS LES
C |                |    | RESULTATS
C |   VARLIS       | -->| BLOC CONTENANT LES VARIABLES A IMPRIMER
C |   NHIST        | -->| NOMBRE DE VALEURS DANS HIST
C |   N            | -->| NOMBRE DE POINTS DU MAILLAGE.
C |   NRES         | -->| UNITE LOGIQUE DU FICHIER DE RESULTATS.
C |   STD          | -->| BINAIRE DU FICHIER DE RESULTATS (IBM,I3E,STD)
C |   AT , LT      | -->| TEMPS , NUMERO DU PAS DE TEMPS
C |   LISPRD       | -->| PERIODE DE SORTIE SUR LISTING.
C |   LEOPRD       | -->| PERIODE DE SORTIE SUR LE FICHIER DE RESULTAT
C |   TEXTE        | -->| NOMS ET UNITES DES VARIABLES.
C |   PTINIG       | -->| 1ER PAS DE TEMPS POUR LES SORTIES GRAPHIQUES
C |   PTINIL       | -->| 1ER PAS DE TEMPS POUR LES SORTIES LISTING
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C     - FICHIERS EXTERIEURS: NRES (FICHIER DE RESULTATS).
C     - DOCUMENTION: NOTICES LEONARD ET SELAFIN.
C     - SOUS-PROGRAMMES APPELES : ECRI2 , IMPVEC
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
      TYPE(BIEF_OBJ)   , INTENT(IN) :: VARSOR
      CHARACTER(LEN=8) , INTENT(IN) :: FORMAT_RES
      INTEGER          , INTENT(IN) :: NRES,LT,LISPRD,LEOPRD
      INTEGER          , INTENT(IN) :: NHIST,PTINIG,PTINIL,N
      INTEGER          , INTENT(IN) :: MAXVAR
      DOUBLE PRECISION , INTENT(IN) :: AT,HIST(NHIST)
      CHARACTER(LEN=32), INTENT(IN) :: TEXTE(*)
      CHARACTER(LEN=3) , INTENT(IN) :: STD
      LOGICAL          , INTENT(IN) :: SORLEO(MAXVAR),SORIMP(MAXVAR)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER LTT,K
C
      LOGICAL LEO,IMP
C
C-----------------------------------------------------------------------
C
C LOGIQUES POUR DECIDER DES SORTIES
C
      IMP=.FALSE.
      LEO=.FALSE.
      LTT=(LT/LISPRD)*LISPRD
      IF(LT.EQ.LTT.AND.LT.GE.PTINIL) IMP=.TRUE.
      LTT=(LT/LEOPRD)*LEOPRD
      IF(LT.EQ.LTT.AND.LT.GE.PTINIG) LEO=.TRUE.
C
C-----------------------------------------------------------------------
C
      IF(LEO) THEN
        CALL WRITE_DATA(FORMAT_RES,       ! ID FORMAT FICHIER RESULTAT
     *                  NRES,             ! LU FICHIER RESULTAT
     *                  MAXVAR,           ! MAX NO VARIABLES DANS VARSOR
     *                  AT,               ! TEMPS
     *                  LT,               ! PAS DE TEMPS
     *                  SORLEO(1:MAXVAR), ! SI SORTIE OU NON
     *                  TEXTE(1:MAXVAR),  ! PAS DE TEMPS
     *                  VARSOR,           ! COLLECTION DES VECTEURS
     *                  N)                ! NUMBER OF VALUES
      ENDIF              
C
C SUR LISTING
C
      IF(IMP) THEN
        DO K=1,MAXVAR
          IF(SORIMP(K)) THEN 
            IF(ASSOCIATED(VARSOR%ADR(K)%P%R)) THEN
              CALL IMPVEC(VARSOR%ADR(K)%P%R,TEXTE(K),N)
            ELSE
              IF(LNG.EQ.1) THEN
                WRITE(LU,*) 'DESIMP: VARIABLE NUMERO: ',K
                WRITE(LU,*) '        PAS OU MAL ALLOUEE'
                WRITE(LU,*) '        OU POINTEUR NON ASSOCIE'
              ENDIF
              IF(LNG.EQ.2) THEN
                WRITE(LU,*) 'DESIMP: VARIABLE NUMBER: ',K
                WRITE(LU,*) '        NOT OR NOT WELL ALLOCATED'
                WRITE(LU,*) '        OR POINTER NOT ASSOCIATED '
              ENDIF
C             CALL PLANTE(1)
C             STOP
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
C=======================================================================
C
      RETURN
      END
