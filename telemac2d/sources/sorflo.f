C                       *****************
                        SUBROUTINE SORFLO
C                       *****************
C
     *(XFLOT,YFLOT,IKLFLO,DEBFLO,FINFLO,NFLOT,NITFLO,FLOPRD,
     * NRBI,TITCAS,BINRES,NOMRBI,NIT,MAXVAR,DATE,TIME,MESH,
     * I_ORIG,J_ORIG)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.7    08/03/07       J-M JANIN (LNH) 30 87 72 84
C
C***********************************************************************
C
C      FONCTION: REMPLIT AU FORMAT SELAFIN LE FICHIER DE RESULTATS
C                BINAIRE AVEC LES INFORMATIONS SUR LES TRAJECTOIRES
C                SUIVIES PAR LES FLOTTEURS.
C
C-----------------------------------------------------------------------
C
C      FUNCTION: FILLS THE BINARY RESULTS FILE AT SELAFIN FORMAT WITH
C                THE DATA ON TRAJECTORIES OF FLOATING BODIES.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  XFLOT,YFLOT   | -->| POSITIONS SUCCESSIVES DES FLOTTEURS.         |
C |    IKLFLO      | -- | TABLE DE CONNECTIVITE BIDON UTILISEE POUR LA |
C |                |    | SORTIE DES TRAJECTOIRES SOUS FORME DE MAILLAGE
C |    DEBFLO      | -->| NUMEROS DES PAS DE TEMPS DE LARGAGE DE       |
C |                |    | CHAQUE FLOTTEUR.                             |
C |    FINFLO      | -->| NUMEROS DES PAS DE TEMPS DE FIN DE CALCUL DE |
C |                |    | DERIVE POUR CHAQUE FLOTTEUR.                 |
C |    NFLOT       | -->| NOMBRE DE FLOTTEURS.                         |
C |    NITFLO      | -->| NOMBRE MAXIMAL D'ENREGISTREMENTS DES         |
C |                |    | POSITIONS SUCCESSIVES DES FLOTTEURS.         |
C |    FLOPRD      | -->| NOMBRE DE PAS DE TEMPS ENTRE 2 ENREGITREMENTS|
C |                |    | DES POSITIONS SUCCESSIVES DES FLOTTEURS.     |
C |    NRBI        | -->| FICHIER DE RESULTATS BINAIRE SUPPLEMENTAIRE  |
C |                |    | POUR STOCKER LES TRAJECTOIRES DE FLOTTEURS   |
C |    TITCAS      | -->| TITRE DU FICHIER CAS                         |
C |    BINRES      | -->| TYPE DE BINAIRE DU FICHIER DE RESULTATS      |
C |    NOMRBI      | -->| NOM DU FICHIER DE RESULTATS BINAIRE SUP.     |
C |    NIT         | -->| NOMBRE DE PAS DE TEMPS                       |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : TELMAC
C
C  SOUS-PROGRAMME APPELE : ECRGEO,ECRIT
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
      INTEGER, INTENT(IN)             :: I_ORIG,J_ORIG,FLOPRD,NITFLO
      INTEGER, INTENT(IN)             :: NFLOT,NRBI,NIT,MAXVAR
      INTEGER, INTENT(IN)             :: DATE(3),TIME(3)
      INTEGER, INTENT(INOUT)          :: IKLFLO(3*NITFLO*NFLOT)
      INTEGER, INTENT(IN)             :: DEBFLO(NFLOT),FINFLO(NFLOT)
      DOUBLE PRECISION, INTENT(INOUT) :: XFLOT(NITFLO*NFLOT)
      DOUBLE PRECISION, INTENT(INOUT) :: YFLOT(NITFLO*NFLOT)
      CHARACTER(LEN=72), INTENT(IN)   :: TITCAS,NOMRBI
      CHARACTER(LEN=3), INTENT(IN)    :: BINRES
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IIT,NVAR,ISTAT,IDEB,IFIN,IFLOT,NELFLO,IELFLO
C
      INTEGER BIDI(1)
C
      DOUBLE PRECISION A(2)
C
      LOGICAL DRAPO(26)
C
      CHARACTER*32 TEXTE(26),BID32(26)
      CHARACTER*3  BIDC
C
      DATA DRAPO /.TRUE.,25*.FALSE./
      DATA TEXTE(1) /'MAILLAGE                        '/
C
C-----------------------------------------------------------------------
C
      IF(NOMRBI(1:1).NE.' ') THEN
C
         NELFLO= 0
C
         DO 10 IFLOT=1,NFLOT
            IF (DEBFLO(IFLOT).LT.NIT) THEN
               IDEB=(DEBFLO(IFLOT)-1)/FLOPRD + 1
               IFIN=(MIN0(FINFLO(IFLOT),NIT)-1)/FLOPRD + 1
               DO 20 IIT=IDEB,IFIN
                  NELFLO = NELFLO + 1
                  XFLOT (NELFLO) = XFLOT (NITFLO*(IFLOT-1)+IIT)
                  YFLOT (NELFLO) = YFLOT (NITFLO*(IFLOT-1)+IIT)
                  IKLFLO(NELFLO) = NELFLO + 1
20             CONTINUE
               IKLFLO(NELFLO) = NELFLO - 1
            ENDIF
10       CONTINUE
C
         IF(NELFLO.NE.0) THEN
C
            DO 30 IELFLO=1,NELFLO
               IKLFLO(  NELFLO+IELFLO) = IKLFLO(IELFLO)
               IKLFLO(2*NELFLO+IELFLO) = IKLFLO(IELFLO)
               IKLFLO(IELFLO) = IELFLO
30          CONTINUE
C
C           STANDARD SELAFIN
C
            CALL ECRGEO(XFLOT,YFLOT,NELFLO,IKLFLO,
     *                  NRBI,NVAR,TEXTE,BID32,0,TITCAS,DRAPO,26,
     *                  IKLFLO,NELFLO,NELFLO,3,DATE,TIME,
     *                  NCSIZE,NPTIR,MESH%KNOLG%I,
     *                  I3=I_ORIG,I4=J_ORIG)
C
C           ECRITURE DU TEMPS
C
            A(1) = 0.D0
            CALL ECRI2(A,BIDI,BIDC,1,'R4',NRBI,BINRES,ISTAT)
C
C ECRITURE D'ENREGISTREMENTS BIDON POUR SIMULER DES RESULTATS
C DE CALCUL (ON SE SERT DE XFLOT CAR IL FAUT UN TABLEAU DE REELS)
C
            CALL ECRI2 (XFLOT,BIDI,BIDC,NELFLO,'R4',NRBI,BINRES,ISTAT)
C
         ELSE
C
            IF(LNG.EQ.1) WRITE(LU,11) NFLOT
            IF(LNG.EQ.2) WRITE(LU,12) NFLOT
11          FORMAT(1X,'ATTENTION : VOUS AVEZ PREVU',I4,' FLOTTEURS',/,
     *             1X,'MAIS AUCUN N''A ETE LARGUE',/,
     *             1X,'AVANT LE DERNIER PAS DE TEMPS')
12          FORMAT(1X,'ATTENTION : YOU ASKED FOR',I4,' FLOATS',/,
     *             1X,'BUT NONE OF THEM HAS BEEN RELEASED',/,
     *             1X,'SINCE THE LAST TIME STEP')
C
         ENDIF
C
      ELSE
C
         IF(LNG.EQ.1) WRITE(LU,21)
         IF(LNG.EQ.2) WRITE(LU,22)
21       FORMAT(1X,'LE FICHIER DE DERIVES DE FLOTTEURS',/,
     *          1X,'N''A PU ETRE CONSTITUE,',/,
     *          1X,'IL FAUT FOURNIR DANS LE FICHIER DES PARAMETRES',/,
     *          1X,'UN NOM DE FICHIER DE RESULTATS BINAIRE')
22       FORMAT(1X,'THE FILE OF FLOAT DRIFTS COULD NOT BE FIELD,',/,
     *          1X,'YOU SHOULD INCLUDE IN THE STEERING FILE',/,
     *          1X,'A NAME OF BINARY RESULTS FILE')
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
