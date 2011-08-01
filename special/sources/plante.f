C                       *****************
                        SUBROUTINE PLANTE
C                       *****************
C
     *(IVAL)
C
C***********************************************************************
C BIEF VERSION 5.5           17/08/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : PROVOQUE UNE DIVISION PAR ZERO SI IVAL = 0 DE FACON
C            A FOURNIR L'ARBRE DES APPELS LORS D'UN ARRET DE PROGRAMME
C            APRES LA DETECTION D'UNE ERREUR.
C
C            A EMPLOYER A LA PLACE D'UN ORDRE "STOP"
C
C            EN CAS DE PROBLEME DE COMPILATEUR AVEC CE SOUS-PROGRAMME
C            EFFACER SES DEUX LIGNES.
C
C            CALL PLANTE EST TOUJOURS SUIVI D'UN ORDRE STOP
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  IVAL          | -->| VALEUR ENTIERE
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C  PRECAUTIONS D'EMPLOI : FAIRE SUIVRE L'APPEL DE PLANTE D'UN ORDRE
C                         STOP QUI AIDE LES COMPILATEURS A COMPRENDRE
C                         QU'ON NE RESSORT JAMAIS DE PLANTE.
C
C  ATTENTION : EXISTE AUSSI DANS LA BIBLIOTHEQUE BIEF
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER IVAL,N,ICODE
C
C-----------------------------------------------------------------------
C
      IF(LNG.EQ.1) WRITE(LU,10)
      IF(LNG.EQ.2) WRITE(LU,20)
10    FORMAT(1X,///,1X,'PLANTE : ARRET DU PROGRAMME APRES ERREUR')
20    FORMAT(1X,///,1X,'PLANTE: PROGRAM STOPPED AFTER AN ERROR')
C
      IF(NCSIZE.GT.1) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'SORTIE DE PVM : APPEL DE P_EXIT'
        IF(LNG.EQ.2) WRITE(LU,*) 'EXITING PVM: CALLING P_EXIT'
        CALL P_EXIT
      ENDIF
C
C     ON CRAY, THE 2 FOLLOWING LINES GIVE A LIST OF CALLING ROUTINES
C     UP TO THE ERROR. THEY MAY BE REMOVED IF IT MAKES PROBLEMS
C     ON SOME COMPILER.
C
Cjaj      N = 1/IVAL
Cjaj      WRITE(LU,*) N
C
C-----------------------------------------------------------------------
C parallelism
C
C     IF (NCSIZE > 0) CALL P_EXIT()
C
Cjaj setting exit values according to the ival value 
C    in code ival=0 or ival=1 are used non-consequently
C
      if(ival < 0) then 
        icode = 0      ! just assumed for non-error stop
      else if ((ival==0) .or. (ival==1)) then 
        icode = 2      ! exit code 1 indicating a "controlled" error
      else             
        icode = 1     ! something else? but an error!
      endif 
      write(lu,*) 'Returning exit code: ', icode
C
C     THIS IS NOT STANDARD FORTRAN
C     CALL EXIT(ICODE)
C
      STOP    ! which is usually equivalent to call EXIT(0)
C
C-----------------------------------------------------------------------
C
      END
 
 
