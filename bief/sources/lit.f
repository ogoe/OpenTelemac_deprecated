C                       **************
                        SUBROUTINE LIT
C                       **************
C
     *( X , W , I , C , NVAL , TYPE , CANAL , STD2 , ISTAT )
C
C***********************************************************************
C BIEF VERSION 6.0      01/04/2009    J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C FONCTION  :  LECTURE DE VALEURS SUIVANT DIFFERENTS STANDARDS
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |      X         |<-- | TABLEAU A LIRE S'IL EST REEL
C |      W         |<-- | TABLEAU DE TRAVAIL (UTILISE EN CAS DE
C |                |    | CONVERSION DE SIMPLE EN DOUBLE PRECISION)
C |      I         |<-- | TABLEAU A LIRE S'IL EST ENTIER
C |      C         |<-- | CHAINE DE CARACTERES A LIRE
C |      NVAL      | -->| NOMBRE DE VALEURS DANS LE TABLEAU
C |                |    | OU NOMBRE DE CARACTERES DE LA CHAINE
C |      TYPE      | -->| TYPE DES DONNEES A LIRE
C |      CANAL     | -->| UNITE LOGIQUE POUR L'ECRITURE
C |      STD       | -->| STANDARD DE LECTURE : STD , IBM OU I3E
C |      ISTAT     |<-- | ENTIER EN CAS D'ERREUR
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C----------------------------------------------------------------------
C  PRECAUTIONS D'EMPLOI : SI LA CHAINE DE CARACTERES STD VAUT
C                         IBM OU I3E, ON APPELLE LES SOUS-PROGRAMMES
C                         LECIBM OU LECI3E QUI SONT DEPENDANTS DE LA
C                         MACHINE UTILISEE.
C----------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NVAL,CANAL
      INTEGER, INTENT(INOUT)          :: ISTAT
      CHARACTER*(*), INTENT(IN)       :: TYPE,STD2
      INTEGER, INTENT(INOUT)          :: I(NVAL)
      DOUBLE PRECISION, INTENT(INOUT) :: X(NVAL)
      REAL, INTENT(INOUT)             :: W(NVAL)
      CHARACTER*(*), INTENT(INOUT)    :: C
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER J 
      CHARACTER(LEN=8) STD  
C
      INTRINSIC DBLE,MIN,LEN
C
C-----------------------------------------------------------------------
C
      ISTAT = 0
C
C-----------------------------------------------------------------------
C
C     STD2 MAY BE SHORTER THAN 8 CHARACTERS
      STD='        '
      STD(1:MIN(8,LEN(STD2)))=STD2(1:MIN(8,LEN(STD2)))
C
C-----------------------------------------------------------------------
C
      IF(STD(1:3).EQ.'STD'.OR.STD(1:7).EQ.'SERAFIN') THEN
C
         IF(TYPE(1:2).EQ.'R4') THEN
            IF(STD(1:8).EQ.'SERAFIND') THEN
!             IF SERAFIN DOUBLE, R4 SHOULD BE R8
              READ(CANAL,END=100,ERR=101)(X(J),J=1,NVAL)
            ELSE
              READ(CANAL,END=100,ERR=101)(W(J),J=1,NVAL)
              DO J=1,NVAL
                X(J) = DBLE(W(J))
              ENDDO
            ENDIF
         ELSEIF(TYPE(1:2).EQ.'R8') THEN
            READ(CANAL,END=100,ERR=101)(X(J),J=1,NVAL)
         ELSEIF (TYPE(1:1).EQ.'I') THEN
            READ(CANAL,END=100,ERR=101)(I(J),J=1,NVAL)
         ELSEIF(TYPE(1:2).EQ.'CH') THEN
            READ(CANAL,END=100,ERR=101) C(1:NVAL)
         ELSE
            IF(LNG.EQ.1) WRITE(LU,20) TYPE
            IF(LNG.EQ.2) WRITE(LU,21) TYPE
20          FORMAT(1X,'LIT : TYPE INCONNU :',A2)
21          FORMAT(1X,'LIT : UNKNOWN TYPE :',A2)
            CALL PLANTE(1)
            STOP
         ENDIF
C
         GO TO 102
C
100      CONTINUE
         IF(LNG.EQ.1) THEN
          WRITE(LU,'(1X,A)')       'LIT : FIN DE FICHIER ANORMALE'
          WRITE(LU,'(1X,A)')       'ON VOULAIT LIRE UN'
          WRITE(LU,'(1X,A,1I6,A)') 'ENREGISTREMENT DE ',NVAL,' VALEURS'
          WRITE(LU,'(1X,A,A)')     'DE TYPE : ',TYPE
          WRITE(LU,'(1X,A,1I6)')   'SUR LE CANAL : ',CANAL
         ENDIF
         IF(LNG.EQ.2) THEN
          WRITE(LU,'(1X,A)')       'LIT : ABNORMAL END OF FILE'
          WRITE(LU,'(1X,A)')       'ONE INTENDED TO READ'
          WRITE(LU,'(1X,A,1I6,A)') 'A RECORD OF ',NVAL,' VALUES'
          WRITE(LU,'(1X,A,A)')     'OF TYPE : ',TYPE
          WRITE(LU,'(1X,A,1I6)')   'ON LOGICAL UNIT : ',CANAL
         ENDIF
C        ISTAT = -6
         CALL PLANTE(1)
         STOP
C
101      CONTINUE
         IF(LNG.EQ.1) THEN
          WRITE(LU,'(1X,A)')       'LIT : ERREUR DE LECTURE'
          WRITE(LU,'(1X,A)')       'ON VOULAIT LIRE UN'
          WRITE(LU,'(1X,A,1I6,A)') 'ENREGISTREMENT DE ',NVAL,' VALEURS'
          WRITE(LU,'(1X,A,A)')     'DE TYPE : ',TYPE
          WRITE(LU,'(1X,A,1I6)')   'SUR LE CANAL : ',CANAL
         ENDIF
         IF(LNG.EQ.2) THEN
          WRITE(LU,'(1X,A)')       'LIT : READ ERROR'
          WRITE(LU,'(1X,A)')       'ONE INTENDED TO READ'
          WRITE(LU,'(1X,A,1I6,A)') 'A RECORD OF ',NVAL,' VALUES'
          WRITE(LU,'(1X,A,A)')     'OF TYPE : ',TYPE
          WRITE(LU,'(1X,A,1I6)')   'ON LOGICAL UNIT : ',CANAL
         ENDIF
C        ISTAT = -6
         CALL PLANTE(1)
         STOP
C
102      CONTINUE
C
C-----------------------------------------------------------------------
C
C     ELSEIF(STD(1:3).EQ.'IBM') THEN
C
C        IF (TYPE(1:2).EQ.'R4') THEN
C           CALL LECIBM( W , NVAL , TYPE , CANAL )
C           DO 77 J=1,NVAL
C             X(J)=DBLE(W(J))
C77          CONTINUE
C        ELSEIF (TYPE(1:2).EQ.'R8') THEN
C           CALL LECIBM( X , NVAL , TYPE , CANAL )
C        ELSEIF (TYPE(1:1).EQ.'I') THEN
C           CALL LECIBM( I , NVAL , TYPE , CANAL )
C        ELSEIF (TYPE(1:2).EQ.'CH') THEN
C           CALL LECIBM( C , NVAL , TYPE , CANAL )
C        ELSE
C           IF(LNG.EQ.1) WRITE(LU,20) TYPE
C           IF(LNG.EQ.2) WRITE(LU,21) TYPE
C           CALL PLANTE(0)
C           STOP
C        ENDIF
C
C-----------------------------------------------------------------------
C
C     ELSEIF(STD(1:3).EQ.'I3E') THEN
C  LECTURE R4 ET R8 A VERIFIER
C        IF (TYPE(1:2).EQ.'R4') THEN
C           CALL LECI3E( W , NVAL , 'F' , CANAL , ISTAT )
C           DO 78 J=1,NVAL
C             X(J)=DBLE(W(J))
C78          CONTINUE
C        ELSEIF (TYPE(1:2).EQ.'R8') THEN
C           CALL LECI3E( X , NVAL , 'F' , CANAL , ISTAT )
C        ELSEIF (TYPE(1:1).EQ.'I') THEN
C           CALL LECI3E( I , NVAL , 'I' , CANAL , ISTAT )
C        ELSEIF (TYPE(1:2).EQ.'CH') THEN
C           CALL LECI3E( C , NVAL , 'C' , CANAL , ISTAT )
C        ELSE
C           IF(LNG.EQ.1) WRITE(LU,20) TYPE
C           IF(LNG.EQ.2) WRITE(LU,21) TYPE
C           CALL PLANTE(0)
C           STOP
C        ENDIF
C
C
C-----------------------------------------------------------------------
C
      ELSE
C
        IF(LNG.EQ.1) WRITE(LU,10) STD
        IF(LNG.EQ.2) WRITE(LU,11) STD
10      FORMAT(1X,'LIT : STANDARD DE LECTURE INCONNU :',A8)
11      FORMAT(1X,'LIT : UNKNOWN STANDARD:',A8)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
