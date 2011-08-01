C                       *****************
                        SUBROUTINE INTANG
C                       *****************
C
     *( LAVANT, LAPRES, IPLAN , NPLAN , DELTAD)
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 1.2 ----                         C
C                                                                     C
C  INTANG :  CALCUL DES INDICES ANGULAIRES ENTOURANT UNE DIRECTION    C
C  ********  DONNEE POUR LE TERME D'INTERACTIONS NON-LINEAIRES AVEC   C
C            DE LA METHODE "DISCRETE INTERACTION APPROXIMATION (DIA)" C
C            PROPOSEE PAR HASSELMANN ET HASSELMANN (1985).            C
C            PROCEDURE SPECIFIQUE AU CAS OU LES DIRECTIONS SONT       C
C            REPARTIES REGULIEREMENT SUR [0;2.PI].                    C
C                                                                     C
C   - CREE POUR VERSION 1.2  LE 26/06/96 PAR M. BENOIT                C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  ARGUMENTS D'APPEL :                                                C
C  *******************                                                C
C  +-------------+----+--------------------------------------------+  C
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! LAVANT      !<-- ! INDICE ANGULAIRE PRECEDANT LA DIR. ARRIVEE !  C
C  ! LAPRES      !<-- ! INDICE ANGULAIRE SUIVANT   LA DIR. ARRIVEE !  C
C  ! IPLAN       ! -->! INDICE DE LA DIRECTION DE DEPART           !  C
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !  C
C  ! DELTAD      ! -->! DEVIATION PAR RAPPORT A LA DIRECTION DEPART!  C
C  +-------------+----+--------------------------------------------+  C
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !  C
C  +---------------------------------------------------------------+  C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  PRENL1                     C
C  ********    - PROGRAMME(S) APPELE(S) :    -                        C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  REMARQUES :                                                        C
C  ***********                                                        C
C                                                                     C
C  - LA DEVIATION DELTAD EST A FOURNIR EN DEGRES.                     C
C                                                                     C
C  - LAVANT ET LAPRES SONT TOUS DEUX COMPRIS ENTRE 1 ET NPLAN.        C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES.
C     """""""""""""""""""""
      INTEGER  LAVANT, LAPRES, NPLAN , IPLAN
      DOUBLE PRECISION DELTAD
C
C.....VARIABLES LOCALES.
C     """"""""""""""""""
      DOUBLE PRECISION TETA  , DTETAD
C
C
      DTETAD=360.D0/DBLE(NPLAN)
      TETA=DBLE(IPLAN-1)*DTETAD+DELTAD
C
C.....TETA EST AJUSTE POUR SE TROUVER COMPRIS ENTRE 0 ET 360 DEG.
C     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  100 IF (TETA.GE.360.D0) THEN
        TETA=TETA-360.D0
        GOTO 100
      ENDIF
  110 IF (TETA.LT.0.D0) THEN
        TETA=TETA+360.D0
        GOTO 110
      ENDIF
C
C.....CALCUL DES INDICES ANGULAIRES PRECEDANT ET SUIVANT TETA.
C     """"""""""""""""""""""""""""""""""""""""""""""""""""""""
      LAVANT=INT(TETA/DTETAD)+1
      LAPRES=LAVANT+1
      IF (LAPRES.GT.NPLAN) LAPRES=1
C
      RETURN
      END
