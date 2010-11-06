C                       ***************
                        SUBROUTINE FOND
C                       ***************
C
     *(ZF  ,X,Y,NPOIN,NFON,NBOR,KP1BOR,NPTFR)
C
C***********************************************************************
C  BIEF VERSION 5.9       20/03/08    J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C   FONCTION : CE SOUS-PROGRAMME INITIALISE LA COTE DU FOND
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    ZF          |<-- |  COTE DU FOND AUX NOEUDS DU MAILLAGE.
C |    X,Y         | -->|  COORDONNEES DU MAILLAGE.
C |    NPOIN       | -->|  NOMBRE DE POINTS DU MAILLAGE.
C |    NFON        | -->|  NUMERO D'UNITE LOGIQUE DU FICHIER DES FONDS.
C |    NBOR        | -->|  NUMEROS GLOBAUX DES NOEUDS DE BORD.
C |    NPTFR       | -->|  NOMBRE DE POINTS FRONTIERES.
C |    PRIVE       | -->|  TABLEAU PRIVE POUR L'UTILISATEUR.
C |    COTE        |<-- |  TABLEAU DE TRAVAIL DE DIMENSION NPMAX.
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PREDON
C
C SOUS-PROGRAMME APPELE : FASP
C
C***********************************************************************
C
      USE BIEF, EX_FOND => FOND
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NFON,NPOIN,NPTFR
      DOUBLE PRECISION, INTENT(OUT) :: ZF(NPOIN)
      DOUBLE PRECISION, INTENT(IN)  :: X(NPOIN),Y(NPOIN)
      INTEGER, INTENT(IN) :: NBOR(NPTFR),KP1BOR(NPTFR)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NP,ERR
C
      DOUBLE PRECISION BID
C
      CHARACTER*1 C
C
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: XRELV,YRELV,COTE
C
C-----------------------------------------------------------------------
C        LECTURE DES POINTS RELEVES SUR LA TABLE A DIGITALISER
C                       SUR UNITE LOGIQUE NFON
C-----------------------------------------------------------------------
C
C  EVALUATION DU NOMBRE DE DONNEES
C
      NP = 0
20    READ(NFON,120,END=24,ERR=124) C
120   FORMAT(A1)
      IF(C(1:1).NE.'C'.AND.C(1:1).NE.'B') THEN
        BACKSPACE ( UNIT = NFON )
        NP = NP + 1
        READ(NFON,*) BID,BID,BID
      ENDIF
      GO TO 20
124   CONTINUE
      IF(LNG.EQ.1) WRITE(LU,18) NP
      IF(LNG.EQ.2) WRITE(LU,19) NP
18    FORMAT(1X,'FOND (BIEF)'
     *      ,/,1X,'ERREUR DANS LE FICHIER DES FONDS'
     *      ,/,1X,'A LA LIGNE ',I7)
19    FORMAT(1X,'FOND (BIEF)'
     *      ,/,1X,'ERROR IN THE BOTTOM FILE'
     *      ,/,1X,'AT LINE ',I7)
      STOP 
24    CONTINUE
C
C  ALLOCATION DYNAMIQUE DES TABLEAUX
C
      ALLOCATE(XRELV(NP),STAT=ERR)
      ALLOCATE(YRELV(NP),STAT=ERR)
      ALLOCATE(COTE(NP) ,STAT=ERR)
C
      IF(ERR.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,10) NP
        IF(LNG.EQ.2) WRITE(LU,11) NP
10      FORMAT(1X,'FOND (BIEF)'
     *      ,/,1X,'ERREUR A L''ALLOCATION DE 3 TABLEAUX'
     *      ,/,1X,'DE TAILLE ',I7)
11      FORMAT(1X,'FOND (BIEF)'
     *      ,/,1X,'ERROR DURING ALLOCATION OF 3 ARRAYS'
     *      ,/,1X,'OF SIZE ',I7)
        STOP
      ENDIF
C
C  LECTURE REELLE
C
      REWIND(NFON)
      NP = 0
23    READ(NFON,120,END=22,ERR=122) C
      IF(C(1:1).NE.'C'.AND.C(1:1).NE.'B') THEN
        BACKSPACE ( UNIT = NFON )
        NP = NP + 1
        READ(NFON,*) XRELV(NP) , YRELV(NP) , COTE(NP)
      ENDIF
      GO TO 23
C
122   CONTINUE
      IF(LNG.EQ.1) WRITE(LU,12) NP
      IF(LNG.EQ.2) WRITE(LU,13) NP
12    FORMAT(1X,'FOND (BIEF)'
     *      ,/,1X,'ERREUR DANS LE FICHIER DES FONDS'
     *      ,/,1X,'A LA LIGNE ',I7)
13    FORMAT(1X,'FOND (BIEF)'
     *      ,/,1X,'ERROR IN THE BOTTOM FILE'
     *      ,/,1X,'AT LINE ',I7)
      STOP 
C
22    CONTINUE
C
      IF(LNG.EQ.1) WRITE(LU,112) NP
      IF(LNG.EQ.2) WRITE(LU,113) NP
112   FORMAT(1X,'FOND (BIEF) :'
     *      ,/,1X,'NOMBRE DE POINTS DANS LE FICHIER DES FONDS : ',I7)
113   FORMAT(1X,'FOND (BIEF):'
     *      ,/,1X,'NUMBER OF POINTS IN THE BOTTOM FILE: ',I7)
C
C-----------------------------------------------------------------------
C   LE FOND EST CALCULE PAR INTERPOLATION SUR LES POINTS INTERIEURS
C                         AU DOMAINE
C-----------------------------------------------------------------------
C
      CALL FASP(X,Y,ZF,NPOIN,XRELV,YRELV,COTE,NP,NBOR,KP1BOR,NPTFR,0.D0)
C
C-----------------------------------------------------------------------
C
      DEALLOCATE(XRELV)
      DEALLOCATE(YRELV)
      DEALLOCATE(COTE)
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
