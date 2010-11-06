C                       *****************
                        SUBROUTINE SORTIE
C                       *****************
C
     *( CHAINE , MNEMO , NBRE , SORLEO )
C
C***********************************************************************
C  BIEF VERSION 6.0    03/11/2009                    J-M HERVOUET (LNHE)
C
C  03/11/2009 : JOKER '*' ALLOWED IN NAMES.
C
C***********************************************************************
C
C     FONCTION  : AFFECTATION DES VARIABLES SORLEO ET SORIMP
C
C----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !   CHAINE       !<-->! CHAINE DES VARIABLES DE SORTIE               !
C !   MNEMO        !<-->! MNEMONIQUE DES VARIABLES                     !
C !   NBRE         !<-->! NOMBRE MAXIMAL DE VARIABLES A IMPRIMER       !
C !   SORLEO       !<-->! TABLEAU DE LOGIQUES POUR IMPRESSION          !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C SOUS-PROGRAMME APPELANT :
C SOUS-PROGRAMMES APPELES :
C **********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NBRE
C
      CHARACTER*(*), INTENT(INOUT) :: CHAINE
      CHARACTER(LEN=8), INTENT(IN) :: MNEMO(NBRE)
C
      LOGICAL, INTENT(INOUT) :: SORLEO(NBRE)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C VARIABLES INTERNES :
C
      CHARACTER C(2)
      CHARACTER(LEN=8) MOT(100)
      INTEGER I,J,LONG,I1,I2,NMOT,L
      LOGICAL OK
C
      INTRINSIC LEN
C
C-----------------------------------------------------------------------
C
C  SEPARATEURS ADMIS DANS LA CHAINE
C
      C(1) = ','
      C(2) = ';'
      LONG = LEN(CHAINE)
      IF (LONG.EQ.0) THEN
        IF(LNG.EQ.1) WRITE(LU,1002)
        IF(LNG.EQ.2) WRITE(LU,1003)
1002    FORMAT(1X,'SORTIE (BIEF) : CHAINE VIDE')
1003    FORMAT(1X,'SORTIE (BIEF): EMPTY STRING')
        CALL PLANTE(1)
        STOP
      ENDIF
C
      DO I=1,LONG
        DO J=1,2
          IF(CHAINE(I:I).EQ.C(J)) CHAINE(I:I) = ' '
        ENDDO
      ENDDO
C
C CHAINE EST MAINTENANT COMPOSEE DE MOTS SEPARES PAR DES BLANCS
C
      I1 = 0
      NMOT=0
C
 10   CONTINUE
      IF (I1.GE.LONG) GOTO 30
      I1=I1+1
      IF (CHAINE(I1:I1).EQ.' ') GOTO 10
C
      I2=0
C
 20   CONTINUE
      I2=I2+1
      IF (CHAINE(I1+I2:I1+I2).NE.' ') GOTO 20
C
      NMOT=NMOT+1
      IF (I2.GT.8) THEN
        IF(LNG.EQ.1) WRITE(LU,1004) CHAINE
        IF(LNG.EQ.2) WRITE(LU,1005) CHAINE
1004    FORMAT(1X,'SORTIE (BIEF) : PLUS DE 8 CARACTERES PAR MOT',/,1X,
     *            '                 DANS LA CHAINE :',A)
1005    FORMAT(1X,'SORTIE (BIEF): MORE THAN 8 LETTERS PER WORD',/,1X,
     *            '                 IN THE CHAIN: ',A)
        CALL PLANTE(1)
        STOP
      ENDIF
      MOT(NMOT)=CHAINE(I1:I1+I2)
      I1=I1+I2
      GOTO 10
C
30    CONTINUE
C
C     COMPARAISON DE MOT ET DE MNEMO
C
      DO I=1,NBRE
        DO J=1,NMOT
          OK=.TRUE.
          DO L=1,8
C           A JOKER '*' IS ALLOWED
            IF(MOT(J)(L:L).NE.MNEMO(I)(L:L).AND.MOT(J)(L:L).NE.'*') THEN
              OK=.FALSE.
              EXIT
            ENDIF
          ENDDO
          SORLEO(I)=OK
          IF(SORLEO(I)) EXIT
        ENDDO
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
