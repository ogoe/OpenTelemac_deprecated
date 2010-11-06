C                       *****************
                        SUBROUTINE LECSNG
C                       *****************
C
     *(NWEIRS,NWRMAX,NPSING,NUMDIG,ZDIG,PHIDIG,IOPTAN,NPSMAX,NPOIN,IFIC)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2     19/04/96     V. GUINOT   (LHF)
C              MODIFIE LE     03/10/96  J.-M. HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C      FONCTIONS:
C      ==========
C
C           LECTURE DES DONNEES DEFINISSANT LES SINGULARITES SUR LE
C           FICHIER FORMATTE 1.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   NWEIRS       | -->| NOMBRE DE SINGULARITES LINEIQUES
C |   NWRMAX       | -->| NOMBRE MAXIMUM DE SINGULARITES PREVUES
C |   NPSING       |<-- | NOMBRE DE POINTS POUR 1 COTE DE SINGULARITE
C |   NUMDIG(K,N,I)|<-- | NUMERO DES POINTS DES DIGUES
C |                |    | DANS LA NUMEROTATION DES POINTS DE BORD
C |   ZDIG         |<-- | COTE DES POINTS DES DIGUES
C |   PHIDIG       |<-- | COEFFICIENT DE DEBIT DES POINTS DES DIGUES
C |   IOPTAN       |<-- | OPTION DE TRAITEMENT DES VITESSES TANGENTES
C |   NPSMAX       | -- | NOMBRE MAXIMUM DE POINTS POUR UN COTE
C |                |    | D'UNE SINGULARITE.
C |   NPOIN        | -->| NOMBRE DE NOEUDS DU MAILLAGE
C |                |    | LES POINTEURS POUR ZDIG, NUMDIG ET PHIDIG
C |                |    | SONT FAITS DANS POINT EN SUPPOSANT QUE :
C |                |    | NPOIN >= NWEIRS * NPSMAX.
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C       
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NWRMAX,NPOIN,IFIC,NWEIRS
      INTEGER, INTENT(INOUT) :: NPSMAX,IOPTAN
C                                                              NPSMAX 
      INTEGER, INTENT(INOUT) :: NPSING(NWEIRS),NUMDIG(2,NWEIRS,*     )
C                                                    NPSMAX
      DOUBLE PRECISION, INTENT(INOUT) :: ZDIG(NWEIRS,*     )
C                                                      NPSMAX
      DOUBLE PRECISION, INTENT(INOUT) :: PHIDIG(NWEIRS,*     )
C         
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER N,I,NNWEIRS
C
C-----------------------------------------------------------------------
C
      NNWEIRS=NWEIRS
C
      READ(IFIC,*,END=900)
      READ(IFIC,*,ERR=998) N,IOPTAN
C
C     CONTROLE DU DIMENSIONNEMENT
C
      IF(N.GT.NWRMAX) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'LECSNG : NOMBRE DE SEUILS : ',N
          WRITE(LU,*) '         TROP GRAND'
          WRITE(LU,*) '         LE MAXIMUM EST :',NWRMAX
        ELSEIF(LNG.EQ.2) THEN
          WRITE(LU,*) 'LECSNG : NUMBER OF WEIRS:',N
          WRITE(LU,*) '         EXCEEDIND THE MAXIMUM OF: ',NWRMAX
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
C
C     COHERENCE AVEC LE FICHIER CAS
C
      IF(N.NE.NWEIRS) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'LECSNG : NOMBRE DE SEUILS : ',N
          WRITE(LU,*) '         DIFFERENT DE LA VALEUR DONNEE DANS LE'
          WRITE(LU,*) '         FICHIER DES PARAMETRES :',NWEIRS
        ELSEIF(LNG.EQ.2) THEN
          WRITE(LU,*) 'LECSNG : NUMBER OF WEIRS:',N
          WRITE(LU,*) '         DIFFERENT FROM THE ONE GIVEN IN THE'
          WRITE(LU,*) '         PARAMETER FILE: ',NWEIRS
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
C
      IF(LNG.EQ.1) THEN
        WRITE(LU,*)'LECSNG : NOMBRE DE DIGUES :',NWEIRS
      ELSEIF(LNG.EQ.2) THEN
        WRITE(LU,*)'LECSNG : NUMBER OF WEIRS :',NWEIRS
      ENDIF
C
      NPSMAX = 0
C
      DO 10 N=1,NWEIRS
        READ(IFIC,*,END=900)
        READ(IFIC,*,END=900)
        READ(IFIC,*,ERR=997) NPSING(N)
        NPSMAX=MAX(NPSING(N),NPSMAX)
        READ(IFIC,*,END=900)
        READ(IFIC,*,ERR=996) (NUMDIG(1,N,I),I=1,NPSING(N))
        READ(IFIC,*,END=900)
        READ(IFIC,*,ERR=994) (NUMDIG(2,N,I),I=1,NPSING(N))
        READ(IFIC,*,END=900)
        READ(IFIC,*,ERR=992) (ZDIG(N,I),I=1,NPSING(N))
        READ(IFIC,*,END=900)
        READ(IFIC,*,ERR=991) (PHIDIG(N,I),I=1,NPSING(N))
10    CONTINUE
C
      GO TO 1000
C
C-----------------------------------------------------------------------
C     MESSAGES D'ERREURS
C-----------------------------------------------------------------------
C
998   CONTINUE
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'LECSNG : ERREUR DE LECTURE SUR LE'
        WRITE(LU,*) '         FICHIER DE DONNEES FORMATE 1'
        WRITE(LU,*) '         2EME LIGNE DU FICHIER NON CONFORME.'
      ELSEIF(LNG.EQ.2) THEN
        WRITE(LU,*) 'LECSNG : READ ERROR ON THE'
        WRITE(LU,*) '         FORMATTED DATA FILE 1'
        WRITE(LU,*) '         AT LINE 2'
      ENDIF
      GO TO 2000
C
997   CONTINUE
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'LECSNG : ERREUR DE LECTURE SUR LE'
        WRITE(LU,*) '         FICHIER DE DONNEES FORMATE 1'
        WRITE(LU,*) '         POUR LA SINGULARITE ',N
        WRITE(LU,*) '         NOMBRE DE POINTS ILLISIBLE'
      ELSEIF(LNG.EQ.2) THEN
        WRITE(LU,*) 'LECSNG : READ ERROR ON THE'
        WRITE(LU,*) '         FORMATTED DATA FILE 1'
        WRITE(LU,*) '         FOR SINGULARITY NUMBER ',N
        WRITE(LU,*) '         THE NUMBER OF POINTS CANNOT BE READ'
      ENDIF
      GO TO 2000
C
996   CONTINUE
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'LECSNG : ERREUR DE LECTURE SUR LE'
        WRITE(LU,*) '         FICHIER DE DONNEES FORMATE 1'
        WRITE(LU,*) '         POUR LA SINGULARITE ',N
        WRITE(LU,*) '         NUMDIGS DES POINTS ILLISIBLE'
        WRITE(LU,*) '         POUR LE COTE 1'
      ELSEIF(LNG.EQ.2) THEN
        WRITE(LU,*) 'LECSNG : READ ERROR ON THE'
        WRITE(LU,*) '         FORMATTED DATA FILE 1'
        WRITE(LU,*) '         FOR SINGULARITY NUMBER ',N
        WRITE(LU,*) '         THE NUMBER OF THE POINTS CANNOT BE READ'
        WRITE(LU,*) '         FOR SIDE NUMBER 1'
      ENDIF
      GO TO 2000
C
994   CONTINUE
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'LECSNG : ERREUR DE LECTURE SUR LE'
        WRITE(LU,*) '         FICHIER DE DONNEES FORMATE 1'
        WRITE(LU,*) '         POUR LA SINGULARITE ',N
        WRITE(LU,*) '         NUMDIGS DES POINTS ILLISIBLE'
        WRITE(LU,*) '         POUR LE COTE 2'
      ELSEIF(LNG.EQ.2) THEN
        WRITE(LU,*) 'LECSNG : READ ERROR ON THE'
        WRITE(LU,*) '         FORMATTED DATA FILE 1'
        WRITE(LU,*) '         FOR SINGULARITY NUMBER ',N
        WRITE(LU,*) '         THE NUMBER OF THE POINTS CANNOT BE READ'
        WRITE(LU,*) '         FOR SIDE NUMBER 2'
      ENDIF
      GO TO 2000
C
992   CONTINUE
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'LECSNG : ERREUR DE LECTURE SUR LE'
        WRITE(LU,*) '         FICHIER DE DONNEES FORMATE 1'
        WRITE(LU,*) '         POUR LA SINGULARITE ',N
        WRITE(LU,*) '         COTES SUR LA DIGUE ILLISIBLES'
      ELSEIF(LNG.EQ.2) THEN
        WRITE(LU,*) 'LECSNG : READ ERROR ON THE'
        WRITE(LU,*) '         FORMATTED DATA FILE 1'
        WRITE(LU,*) '         FOR SINGULARITY NUMBER ',N
        WRITE(LU,*) '         ELEVATIONS ON THE WEIR CANNOT BE READ'
      ENDIF
      GO TO 1000
C
991   CONTINUE
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'LECSNG : ERREUR DE LECTURE SUR LE'
        WRITE(LU,*) '         FICHIER DE DONNEES FORMATE 1'
        WRITE(LU,*) '         POUR LA SINGULARITE ',N
        WRITE(LU,*) '         COEFFICIENTS DE DEBIT ILLISIBLES'
      ELSEIF(LNG.EQ.2) THEN
        WRITE(LU,*) 'LECSNG : READ ERROR ON THE'
        WRITE(LU,*) '         FORMATTED DATA FILE 1'
        WRITE(LU,*) '         FOR SINGULARITY NUMBER ',N
        WRITE(LU,*) '         DISCHARGE COEFFICIENTS CANNOT BE READ'
      ENDIF
      GO TO 2000
C
900   CONTINUE
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'LECSNG : ERREUR DE LECTURE SUR LE'
        WRITE(LU,*) '         FICHIER DE DONNEES FORMATE 1'
        WRITE(LU,*) '         FIN DE FICHIER PREMATUREE'
      ELSEIF(LNG.EQ.2) THEN
        WRITE(LU,*) 'LECSNG : READ ERROR ON THE'
        WRITE(LU,*) '         FORMATTED DATA FILE 1'
        WRITE(LU,*) '         UNEXPECTED END OF FILE'
      ENDIF
C
2000  CONTINUE
C
      NNWEIRS = 0
C
1000  CONTINUE
C
      IF(NNWEIRS.EQ.0) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*)
          WRITE(LU,*)'LECSNG : ERREUR DE LECTURE'
          WRITE(LU,*)'         AUCUNE SINGULARITE NE SERA'
          WRITE(LU,*)'         PRISE EN COMPTE.'
          WRITE(LU,*)
        ELSEIF(LNG.EQ.2) THEN
          WRITE(LU,*)
          WRITE(LU,*)'LECSNG : READ ERROR'
          WRITE(LU,*)'         NO SINGULARITY WILL BE TAKEN'
          WRITE(LU,*)'         INTO ACCOUNT'
          WRITE(LU,*)
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C
C VERIFICATION DU DIMENSIONNEMENT DES TABLEAUX NUMDIG, PHIDIG, ZDIG
C
      IF(NPOIN.LT.NWEIRS*NPSMAX) THEN
C
        IF(LNG.EQ.1) THEN
          WRITE(LU,*)
          WRITE(LU,*)'LECSNG : TROP DE POINTS SUR LES SINGULARITES'
          WRITE(LU,*)'         UN CHANGEMENT DES POINTEURS EST'
          WRITE(LU,*)'         NECESSAIRE'
          WRITE(LU,*)
        ELSEIF(LNG.EQ.2) THEN
          WRITE(LU,*)
          WRITE(LU,*)'LECSNG : TOO MANY POINTS IN THE SINGULARITIES'
          WRITE(LU,*)'         IT IS NECESSARY TO CHANGE THE MEMORY'
          WRITE(LU,*)'         POINTERS'
          WRITE(LU,*)
        ENDIF
C
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
