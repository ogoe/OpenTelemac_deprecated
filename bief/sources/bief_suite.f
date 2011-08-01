C                       *********************
                        SUBROUTINE BIEF_SUITE
C                       *********************
C
     *(VARSOR,CLAND,NUMDEB,
     * NPRE,STD,HIST,NHIST,NPOIN,AT,TEXTPR,VARCLA,NVARCL,
     * TROUVE,ALIRE,LISTIN,FIN,MAXVAR,NPLAN,DT,NDT)
C
C***********************************************************************
C BIEF VERSION 6.0      09/04/2009    J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C   FONCTION  : LECTURE DES RESULTATS INSCRITS SUR UN FICHIER
C               DE RESULTATS.
C
C   09/12/2008 : STD IS NOW A STRING OF ANY SIZE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   VARSOR       |<-- | BLOC DES TABLEAUX CONTENANT LES VARIABLES
C |   CLAND        |<-- | BLOC DES VARIABLES CLANDESTI-NES
C |   NUMDEB       |<-->| FIN = .TRUE. NUMERO DU DERNIER ENREGISTREMENT
C |                |    | FIN = .FALSE. : NUMERO DE L'ENREGISTREMENT
C |                |    |                 QUE L'ON VEUT LIRE.
C |   NPRE         | -->| NUMERO DE CANAL DU FICHIER
C |   STD          | -->| BINAIRE DU FICHIER : STD, IBM OU I3E
C |   HIST         | -->| TABLEAU DE VALEURS MISES DANS L'ENREGISTREMENT
C |                |    | DU TEMPS.
C |   NHIST        | -->| NOMBRE DE VALEURS DANS LE TABLEAU HIST.
C |   NPOIN        | -->| NOMBRE DE POINTS DANS LE MAILLAGE
C |   AT           | -->| TEMPS
C |   TEXTPR       | -->| NOMS ET UNITES DES VARIABLES.
C |   VARCLA       | -->| TABLEAU OU L'ON RANGE LES VARIABLES
C |                |    | CLANDESTiINES.
C |   NVARCL       | -->| NOMBRE DE VARIABLES CLANDESTi-NES.
C |   TROUVE       |<-- | INDIQUE (TROUVE(K)=1) LES VARIABLES TROUVEES
C |                |    | DANS LE FICHIER.
C |                |    | DE K =  1 A 26 VARIABLES NORMALES
C |                |    | DE K = 27 A 36 VARIABLES CLANDESTi-NES.
C |   ALIRE        | -->| VARIABLES QU'IL FAUT LIRE (POUR LES AUTRES ON
C |                |    | SAUTE L'ENREGISTREMENT CORRESPONDANT)
C |                |    | LES VARIABLES CLANDESTi-NES SONT LUES
C |                |    | SYSTEMATIQUEMENT.
C |   LISTIN       | -->| SI OUI, IMPRESSION D'INFORMATIONS SUR LISTING
C |   FIN          | -->| VOIR LE TROISIEME ARGUMENT NUMDEB
C |   MAXVAR       | -->| DIMENSION DES TABLEAUX DES VARIABLES : ALIRE, ETC
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : LIT
C
C***********************************************************************
C
      USE BIEF, EX_BIEF_SUITE => BIEF_SUITE
      USE M_MED
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: VARSOR,CLAND
      INTEGER, INTENT(IN), OPTIONAL :: NPLAN
      INTEGER, INTENT(IN)           :: NHIST,NVARCL,MAXVAR
      INTEGER                       :: NUMDEB,NPRE,NPOIN,TROUVE(MAXVAR)
      INTEGER                       :: ALIRE(MAXVAR)        
      CHARACTER(LEN=*)              :: STD
      CHARACTER(LEN=32)             :: TEXTPR(MAXVAR),VARCLA(NVARCL)
      DOUBLE PRECISION              :: HIST(*),AT
      LOGICAL                       :: FIN,LISTIN 
      DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: DT
      INTEGER, INTENT(OUT), OPTIONAL :: NDT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NNPLAN, NNDT
      DOUBLE PRECISION DDT
C
C-----------------------------------------------------------------------
C
      IF(PRESENT(NPLAN)) THEN
        NNPLAN=NPLAN
      ELSE
        NNPLAN=1
      ENDIF
C
      SELECT CASE(STD)
C
        CASE ('SERAFIN ','SERAFIND')
C
          IF(PRESENT(DT)) THEN
            CALL SUITE_SERAFIN(VARSOR,CLAND,NUMDEB,NPRE,STD,HIST,
     *                         NHIST,NPOIN,AT,TEXTPR,VARCLA,NVARCL,
     *                         TROUVE,ALIRE,LISTIN,FIN,MAXVAR,NNPLAN,
     *                         DT)
          ELSE
            CALL SUITE_SERAFIN(VARSOR,CLAND,NUMDEB,NPRE,STD,HIST,
     *                         NHIST,NPOIN,AT,TEXTPR,VARCLA,NVARCL,
     *                         TROUVE,ALIRE,LISTIN,FIN,MAXVAR,NNPLAN)
          ENDIF
C
        CASE ('MED     ')
C
          CALL SUITE_MED(VARSOR,CLAND,NUMDEB,NPRE,STD,HIST,NHIST,
     *                   NPOIN,AT,TEXTPR,VARCLA,NVARCL,TROUVE,ALIRE,
     *                   LISTIN,FIN,MAXVAR,NNPLAN,DDT,NNDT)
          IF(PRESENT(DT)) DT=DDT
          IF(PRESENT(NDT)) NDT=NNDT
          
C
        CASE DEFAULT
C
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'BIEF_SUITE : MAUVAIS FORMAT : ',STD
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'BIEF_SUITE: BAD FILE FORMAT : ',STD
          ENDIF          
          CALL PLANTE(1)
          STOP
C
      END SELECT
C
C-----------------------------------------------------------------------
C
      RETURN
      END
