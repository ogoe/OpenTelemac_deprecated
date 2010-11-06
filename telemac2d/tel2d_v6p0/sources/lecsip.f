C                       *****************
                        SUBROUTINE LECSIP
C                       *****************
C
     * (RELAXS,NSIPH,ENTSIP,SORSIP,SECSCE,
     *  ALTSCE,CSSCE,CESCE,DELSCE,ANGSCE,LSCE,MAXSCE,IFIC)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2     19/04/96     V. GUINOT   (LHF)
C              MODIFIE LE     03/10/96  J.-M. HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C      FONCTIONS: LECTURE DES DONNEES SUR LES SIPHONS.
C      ==========
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   RELAXS       |<-- | COEFFICIENT DE RELAXATION.
C |   NPSIPH       |<-- | NOMBRE DE SIPHONS.
C |   ENTSIP       |<-- | INDICES DANS LA NUMEROTATION DES SOURCES
C |   SORSIP       |<-- | INDICES DANS LA NUMEROTATION DES SOURCES
C |   SECSCE       |<-- | SECTION DES SIPHONS (NUMEROTATION DES SOURCES)
C |   ALTSCE       |<-- | COTE DES ENTREES ET SORTIES DE BUSES
C |   CSSCE        |<-- | COEFFICIENTS DE PERTE DE CHARGE
C |                |    | LORS D'UN FONCTIONNEMENT EN SORTIE.
C |   CESCE        |<-- | COEFFICIENTS DE PERTE DE CHARGE
C |                |    | LORS D'UN FONCTIONNEMENT EN ENTREE.
C |   DELSCE       |<-- | ANGLE DES BUSES AVEC LA VERTICALE
C |   ANGSCE       |<-- | ANGLE DES BUSES AVEC L'AXE OX.
C |   LSCE         |<-- | PERTE DE CHARGE LINEAIRE DE LA CONDUITE.
C |   MAXSCE       | -->| NOMBRE MAXIMUM DE POINTS SOURCES.
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
      INTEGER, INTENT(IN) :: MAXSCE,IFIC
      INTEGER, INTENT(INOUT) :: ENTSIP(*),SORSIP(*),NSIPH
      DOUBLE PRECISION, INTENT(INOUT) :: RELAXS
      DOUBLE PRECISION, INTENT(INOUT) :: SECSCE(*),ALTSCE(*),CSSCE(*)
      DOUBLE PRECISION, INTENT(INOUT) :: DELSCE(*),ANGSCE(*)
      DOUBLE PRECISION, INTENT(INOUT) :: CESCE(*),LSCE(*)
C      
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER N      
C
      DOUBLE PRECISION DELTA1,DELTA2,S12,ALT1,ALT2
      DOUBLE PRECISION ANG1,ANG2,CS1,CS2,CE1,CE2,L12      
C
      DOUBLE PRECISION PI
      PARAMETER(PI=3.141592653589D0)
C
C-----------------------------------------------------------------------
C
      READ(IFIC,*,END=900)
      READ(IFIC,*,ERR=998) RELAXS,N
      READ(IFIC,*,END=900)
C
C     CONTROLE DU DIMENSIONNEMENT
C
      IF(N.GT.MAXSCE/2) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'LECSIP : NOMBRE DE SIPHONS : ',N
          WRITE(LU,*) '         TROP GRAND'
          WRITE(LU,*) '         LE MAXIMUM EST :',MAXSCE/2
        ELSEIF(LNG.EQ.2) THEN
          WRITE(LU,*) 'LECSIP : NUMBER OF CULVERTS:',N
          WRITE(LU,*) '         EXCEEDIND THE MAXIMUM OF: ',MAXSCE/2
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
C
C     COHERENCE AVEC LE FICHIER CAS
C
      IF(N.NE.NSIPH) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'LECSIP : NOMBRE DE SIPHONS : ',N
          WRITE(LU,*) '         DIFFERENT DE LA VALEUR DONNEE DANS LE'
          WRITE(LU,*) '         FICHIER DES PARAMETRES :',NSIPH
        ELSEIF(LNG.EQ.2) THEN
          WRITE(LU,*) 'LECSIP : NUMBER OF CULVERTS:',N
          WRITE(LU,*) '         DIFFERENT FROM THE ONE GIVEN IN THE'
          WRITE(LU,*) '         PARAMETER FILE: ',NSIPH
        ENDIF
        CALL PLANTE(1)
        STOP
      ENDIF
C
      DO 10 N=1,NSIPH
        READ(IFIC,*,ERR=997) ENTSIP(N),SORSIP(N),DELTA1,DELTA2,
     *                       CE1,CE2,CS1,CS2,S12,L12,ALT1,ALT2,ANG1,ANG2
        DELSCE(ENTSIP(N)) = DELTA1*PI/180.D0
        DELSCE(SORSIP(N)) = DELTA2*PI/180.D0
        CESCE(ENTSIP(N))  = CE1
        CESCE(SORSIP(N))  = CE2
        CSSCE(ENTSIP(N))  = CS1
        CSSCE(SORSIP(N))  = CS2
        SECSCE(ENTSIP(N)) = S12
        SECSCE(SORSIP(N)) = S12
        LSCE(ENTSIP(N))   = L12
        LSCE(SORSIP(N))   = L12
        ALTSCE(ENTSIP(N)) = ALT1
        ALTSCE(SORSIP(N)) = ALT2
        ANGSCE(ENTSIP(N)) = ANG1*PI/180.D0
        ANGSCE(SORSIP(N)) = ANG2*PI/180.D0
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
        WRITE(LU,*) 'LECSIP : ERREUR DE LECTURE SUR LE'
        WRITE(LU,*) '         FICHIER DE DONNEES FORMATE 1'
        WRITE(LU,*) '         2EME LIGNE DU FICHIER NON CONFORME.'
      ELSEIF(LNG.EQ.2) THEN
        WRITE(LU,*) 'LECSIP : READ ERROR ON THE'
        WRITE(LU,*) '         FORMATTED DATA FILE 1'
        WRITE(LU,*) '         AT LINE 2'
      ENDIF
      GO TO 2000
C
997   CONTINUE
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'LECSIP : ERREUR DE LECTURE SUR LE'
        WRITE(LU,*) '         FICHIER DE DONNEES FORMATE 1'
        WRITE(LU,*) '         POUR LE SIPHON ',N
        WRITE(LU,*) '         DONNEES ILLISIBLES'
      ELSEIF(LNG.EQ.2) THEN
        WRITE(LU,*) 'LECSIP : READ ERROR ON THE'
        WRITE(LU,*) '         FORMATTED DATA FILE 1'
        WRITE(LU,*) '         FOR CULVERT NUMBER ',N
        WRITE(LU,*) '         THE DATA CANNOT BE READ'
      ENDIF
      GO TO 2000
C
900   CONTINUE
      IF(LNG.EQ.1) THEN
        WRITE(LU,*) 'LECSIP : ERREUR DE LECTURE SUR LE'
        WRITE(LU,*) '         FICHIER DE DONNEES FORMATE 1'
        WRITE(LU,*) '         FIN DE FICHIER PREMATUREE'
      ELSEIF(LNG.EQ.2) THEN
        WRITE(LU,*) 'LECSIP : READ ERROR ON THE'
        WRITE(LU,*) '         FORMATTED DATA FILE 1'
        WRITE(LU,*) '         UNEXPECTED END OF FILE'
      ENDIF
C
2000  CONTINUE
C
      NSIPH = 0
C
1000  CONTINUE
C
      IF(NSIPH.EQ.0) THEN
        IF(LNG.EQ.1) THEN
          WRITE(LU,*)
          WRITE(LU,*)'LECSIP : ERREUR DE LECTURE'
          WRITE(LU,*)'         AUCUN SIPHON NE SERA'
          WRITE(LU,*)'         PRIS EN COMPTE.'
          WRITE(LU,*)
        ELSEIF(LNG.EQ.2) THEN
          WRITE(LU,*)
          WRITE(LU,*)'LECSIP : READ ERROR'
          WRITE(LU,*)'         NO CULVERT WILL BE TAKEN'
          WRITE(LU,*)'         INTO ACCOUNT'
          WRITE(LU,*)
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
