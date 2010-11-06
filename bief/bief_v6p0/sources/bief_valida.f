C                       **********************
                        SUBROUTINE BIEF_VALIDA
C                       **********************
C
     *(VARREF,TEXTREF,UREF,REFFORMAT,VARRES,TEXTRES,URES,RESFORMAT,
     * MAXTAB,NP,IT,MAXIT,ACOMPARER)
C
C***********************************************************************
C  BIEF VERSION 6.0     05/08/2009    J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C  FONCTION : VALIDATION DES RESULTATS PAR COMPARAISON AVEC UNE SOLUTION
C             ANALYTIQUE OU AVEC LES RESULTATS STOCKES SUR LE FICHIER DU
C             CALCUL DE REFERENCE.
C
C             CE SOUS-PROGRAMME DOIT ETRE REMPLI EN FONCTION DE CHAQUE
C             CAS PARTICULIER.
C
C             ICI ON NE COMPARE QUE LE DERNIER PAS DE TEMPS DU CALCUL
C
C             ATTENTION : A PART POUR LE FOND, ON SUPPOSE QUE LE FICHIER
C                         DE REFERENCES CONTIENT BIEN LES VARIABLES A
C                         COMPARER.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   NP           | -->| NOMBRE DE POINTS A VERIFIER.                 |
C |   UREF         | -->| UNITE LOGIQUE DU FICHIER DE REFERENCE
C |   STDPRE       | -->| STANDARD DU FICHIER DE RESULTATS DU CALCUL   |
C |                |    | PRECEDENT                                    |
C |   BINPRE       | -->| TYPE DE BINAIRE DU FICHIER DE RESULTATS DU   |
C |                |    | CALCUL PRECEDENT                             |
C |   TEXTPR       | -->| NOMS DES VARIABLES DU CALCUL PRECEDENT.      |
C |   IT           | -->| NUMERO DU PAS DE TEMPS               
C |   MAXIT        | -->| NOMBRE MAXIMUM D'ITERATIONS DU PRESENT CALCUL
C |   STDRES       | -->| STANDARD DU FICHIER DE DESSIN.               |
C |   ACOMPARER    | -->| TABLEAU DES VARIABLE A LIRE DANS LE FICHIER  |
C |                |    | DES RESULTATS DU CLCUL PRECEDENT             |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMME APPELE : SUITE
C
C***********************************************************************
C
      USE BIEF    !, EX_VALIDA => VALIDA
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER, INTENT(IN) :: NP,MAXTAB,IT,MAXIT,URES,UREF
      INTEGER, INTENT(IN) :: ACOMPARER(MAXTAB)
C
      CHARACTER(LEN=32), INTENT(IN) :: TEXTREF(MAXTAB),TEXTRES(MAXTAB)
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: VARREF,VARRES
      CHARACTER(LEN=*), INTENT(IN)  :: REFFORMAT,RESFORMAT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      INTEGER IVAR,I,IREF,IRES,IERMAX
C
      DOUBLE PRECISION TIMEREF,TIMERES,ERMAX,HIST(1),ERR
C
      INTRINSIC MAX
C
      DOUBLE PRECISION P_DMAX
      EXTERNAL         P_DMAX
C
C-----------------------------------------------------------------------
C
C     ATTENTION IL FAUT MAXTAB <300
C
      INTEGER FINDREF(300),FINDRES(300)
C
C-----------------------------------------------------------------------
C
      IF(IT.EQ.MAXIT) THEN
C
C  APPEL DE SUITE POUR RELIRE LE FICHIER DE REFERENCES
C
      IF(LNG.EQ.1) WRITE(LU,10)
      IF(LNG.EQ.2) WRITE(LU,11)                               
      CALL BIEF_SUITE(VARREF,VARREF,IREF,UREF,REFFORMAT,HIST,0,NP,
     *                TIMEREF,TEXTREF,TEXTREF,0,FINDREF,ACOMPARER,
     *                .TRUE.,.TRUE.,MAXTAB)
C
C  APPEL DE SUITE POUR RELIRE LE FICHIER DE RESULTATS
C
      IF(LNG.EQ.1) WRITE(LU,12)
      IF(LNG.EQ.2) WRITE(LU,13)    
      CALL BIEF_SUITE(VARRES,VARRES,IRES,URES,RESFORMAT,HIST,0,NP,
     *                TIMERES,TEXTRES,TEXTRES,0,FINDRES,ACOMPARER,
     *                .TRUE.,.TRUE.,MAXTAB)
C
C-----------------------------------------------------------------------
C
      IF(LNG.EQ.1) WRITE(LU,14)
      IF(LNG.EQ.2) WRITE(LU,15)
C 
      IF(ABS(TIMERES-TIMEREF).GT.1.D-4) THEN
        IF(LNG.EQ.1) WRITE(LU,16)
        IF(LNG.EQ.2) WRITE(LU,17)
      ENDIF
      IF(IRES.NE.IREF) THEN
        IF(LNG.EQ.1) WRITE(LU,18)
        IF(LNG.EQ.2) WRITE(LU,19)
      ENDIF
C
C-----------------------------------------------------------------------       
C
C     BOUCLE SUR LES VARIABLES A VERIFIER
C
      DO IVAR=1,MAXTAB
C
        IF(ACOMPARER(IVAR).EQ.1) THEN
C
C       VERIFICATION DE LA VARIABLE DE NUMERO IVAR
C
          IF(FINDREF(IVAR).EQ.1.AND.FINDRES(IVAR).EQ.1) THEN
C
            ERMAX = 0.D0
            IERMAX = 1
            DO I = 1 , NP
              ERR=ABS(VARREF%ADR(IVAR)%P%R(I)-VARRES%ADR(IVAR)%P%R(I))
              IF(ERR.GT.ERMAX) THEN
                ERMAX=ERR
                IERMAX=I
              ENDIF
            ENDDO
C
            IF(NCSIZE.GT.1) ERMAX=P_DMAX(ERMAX)
            IF(LNG.EQ.1) WRITE(LU,60) TEXTRES(IVAR)(1:16),ERMAX
            IF(LNG.EQ.2) WRITE(LU,61) TEXTRES(IVAR)(1:16),ERMAX
C
          ELSEIF(FINDREF(IVAR).EQ.1) THEN
C
            IF(LNG.EQ.1) WRITE(LU,70) TEXTRES(IVAR)(1:16)
            IF(LNG.EQ.2) WRITE(LU,71) TEXTRES(IVAR)(1:16)      
C
          ENDIF
C
        ENDIF
C
      ENDDO
C
      IF(LNG.EQ.1) WRITE(LU,50) 
      IF(LNG.EQ.2) WRITE(LU,51)
C
      ENDIF 
C
C-----------------------------------------------------------------------
C
10    FORMAT(1X,////,1X,80('='),/,
     *       25X,' PROCEDURE DE VALIDATION ',/,
     *       1X,80('-'),//,
     *       1X,' 1) RELECTURE DU FICHIER DE REFERENCE :',/,
     *       1X,' --------------------------------------',/)
11    FORMAT(1X,////,1X,80('='),/,
     *       25X,' VALIDATION PROCEDURE ',/,
     *       1X,80('-'),//,
     *       1X,' 1) READING THE REFERENCE FILE :',/,
     *       1X,' ------------------------------',/)
12    FORMAT(1X,///,
     *       1X,' 2) RELECTURE DU FICHIER DE RESULTATS :',/,
     *       1X,' --------------------------------------',/)
13    FORMAT(1X,///,
     *       1X,' 2) READING THE RESULTS FILE :',/,
     *       1X,' --------------------------------',/)
14    FORMAT(1X,///,
     *       1X,' 3) COMPARAISON :',/,
     *       1X,' ----------------',/)
15    FORMAT(1X,///,
     *       1X,' 3) COMPARISON:',/,
     *       1X,' --------------',/)
16    FORMAT(1X,///,
     *       1X,' ATTENTION : TEMPS DIFFERENTS',/,
     *       1X,' ----------------------------',/)
17    FORMAT(1X,///,
     *       1X,' BEWARE: TIMES ARE DIFFERENT',/,
     *       1X,' ---------------------------',/)
18    FORMAT(1X,///,
     *       1X,' ATTENTION : NUMEROS D''ENREGISTREMENT DIFFERENTS',/,
     *       1X,' ------------------------------------------------',/)
19    FORMAT(1X,///,
     *       1X,' BEWARE: RECORD NUMBERS ARE DIFFERENT',/,
     *       1X,' ------------------------------------',/)
C
50    FORMAT(1X,80('-'),/,23X,'FIN DU COMPTE-RENDU DE VALIDATION',/,
     *       1X,80('='),////)
51    FORMAT(1X,80('-'),/,23X,'END OF VALIDATION REPORT',/,
     *       1X,80('='),////)
C
60    FORMAT(1X,'VARIABLE : ',A16,'  DIFFERENCE : ',G16.7,/)
61    FORMAT(1X,'VARIABLE: ' ,A16,'  DIFFERENCE: ',G16.7,/)
C
70    FORMAT(1X,'VARIABLE : ',A16,'  NON TROUVEE',/)
71    FORMAT(1X,'VARIABLE: ' ,A16,'  NOT FOUND'  ,/)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
