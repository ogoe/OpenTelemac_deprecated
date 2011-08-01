C                       ******************
                        SUBROUTINE P_INIT1
C                       ******************
C
     *(REP,IREP,IPID,NCSIZE,ITID,NPTIR)
C
C***********************************************************************
C  PARA VERSION 5.1       /06/96           HANS HERRMANN (HANOVRE)
C             MODIFIED        08/06/96     REINHARD HINKELMANN (HANOVRE)
C             MODIFIED        17/12/96            J-M HERVOUET (LNH)
C             MODIFIED        11/12/97                F MARCOS (LNH)
C***********************************************************************
C
C      FONCTIONS: INITIALISATIONS POUR TOMAWAC
C      ==========
C
C      MELDET PROGRAMM BEI PARASTATION AN.
C
C      VARIABLEN
C      ---------
C      IPID     : NUMMER DES KNOTENS / NUMERO DU PROCESSEUR
C      NCSIZE   : PROZESSORENANZAHL  / NOMBRE DE PROCESSEURS
C      ITID     : PROZESSNUMMERNFELD (WIRD AUF NCUBE NICHT BENOETIGT)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     REP        |--> |  NOM DU REPERTOIRE DE TRAVAIL                |
C |     IREP       |--> |  NOMBRE MAX DE CARACTERE DE LA CHAINE REP    |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PRINCI
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER NPTIR,IPID,NCSIZE,ITID(0:255),NPID1,NDEMM
      LOGICAL PREPA
      COMMON/TOMAPA/PREPA,NPID1,NDEMM
C
      CHARACTER*144 REP
      INTEGER IREP
C
      REP = ' '
      IREP = 0
      NCSIZE = 1
      NPTIR=0
      IPID=0
      ITID(0)=0
C
C     IF(LNG.EQ.1) WRITE(LU,*) 'APPEL DE P_INIT1 VERSION VIDE'
C     IF(LNG.EQ.2) WRITE(LU,*) 'CALLING P_INIT1 IN ITS VOID VERSION'
C
C-----------------------------------------------------------------------
C
       RETURN
       END           
