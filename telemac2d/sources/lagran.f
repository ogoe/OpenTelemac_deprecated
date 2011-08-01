C                       *****************
                        SUBROUTINE LAGRAN
C                       *****************
C
     *(NLAG,DEBLAG,FINLAG)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2    17/08/94    J-M JANIN    (LNH) 30 87 72 84
C
C***********************************************************************
C
C      FONCTION :
C
C - INITIALISATION DES PAS DE TEMPS DE DEBUT DE CALCUL DE CHAQUE DERIVE
C - INITIALISATION DES PAS DE TEMPS DE FIN DE CALCUL DE CHAQUE DERIVE
C
C      ATTENTION :
C
C - ENTRE 2 PAS DE TEMPS DE FIN DE CALCUL IL DOIT Y AVOIR AU MOINS UNE
C   ECRITURE DES RESULTATS SUR FICHIER POUR NE PAS PERDRE D'INFORMATION
C - SI 2 CALCULS FINISSENT AU MEME PAS, UN DES 2 (LE PREMIER DANS LA
C   NUMEROTATION) NE FIGURERA PAS DANS LES RESULTATS.
C
C
C                     - A REMPLIR PAR L'UTILISATEUR -
C
C-----------------------------------------------------------------------
C
C      FUNCTION :
C
C - INITIALISING FIRST AND FINAL TIME STEP OF THE LAGRANGIAN DRIFTS
C
C      BEWARE :
C
C - TWO DRIFTS CANNOT BE COMPLETED IN THE SAME TIME-STEP.
C - BETWEEN TWO DRIFT COMPUTATION ENDINGS, THE RESULTS MUST BE SAVED.
C
C
C                     - TO BE FILLED BY USER -
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    NLAG        | -->| NUMBER OF LAGRANGIAN DRIFTS
C |    DEBLAG      |<-- | TIME STEP AT THE BEGINNING                               |
C |    FINLAG      |<-- | TIME STEP AT THE END                                   |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
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
      INTEGER, INTENT(IN)    :: NLAG
      INTEGER, INTENT(INOUT) :: DEBLAG(NLAG) , FINLAG(NLAG)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER ILAG
C
C-----------------------------------------------------------------------
C
C   INITIALISATION DU PAS DE TEMPS DE DEBUT ET DE FIN DE CALCUL DE
C   LE DERIVE LAGRANGIENNE - PAR DEFAUT RIEN N'EST FAIT
C
C-----------------------------------------------------------------------
C
C     THIS WARNING AND THE CALL PLANTE MUST BE REMOVED IF
C     SOMETHING IS IMPLEMENTED BELOW  
C
      IF(LNG.EQ.1) WRITE(LU,20)
      IF(LNG.EQ.2) WRITE(LU,21)
20    FORMAT(1X,'ATTENTION, VOUS APPELEZ LE SOUS-PROGRAMME LAGRAN',/,1X,
     *          'DE LA BIBLIOTHEQUE.   COMME VOUS CALCULEZ UN OU',/,1X,
     *          'PLUSIEURS CHAMPS DE DERIVES LAGRANGIENNES, VOUS',/,1X,
     *          'DEVEZ RAPATRIER "LAGRAN" DANS VOTRE FORTRAN, ET',/,1X,
     *          'LE COMPLETER',/////)
21    FORMAT(1X,'ATTENTION, YOU CALL SUBROUTINE LAGRAN OF THE LIBRARY.',
     *     /,1X,'AS YOU COMPUTE ONE OR MORE FIELDS OF LAGRANGIAN',/,1X,
     *          'DRIFTS, YOU NEED TO COPY THIS SUBROUTINE IN YOUR',/,1X,
     *          'OWN FORTRAN FILE AND COMPLETE IT.',/////)
C
      CALL PLANTE(1)
C
C-----------------------------------------------------------------------
C
C  EXAMPLE :
C
C      DO 10 ILAG=1,NLAG
C         DEBLAG(ILAG) = 1
C         FINLAG(ILAG) = 299
C 10   CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
