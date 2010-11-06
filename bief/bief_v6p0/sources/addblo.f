C                       *****************
                        SUBROUTINE ADDBLO
C                       *****************
C
     *( BLOC , OBJ )
C
C***********************************************************************
C BIEF VERSION 5.5            01/03/95    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C  FONCTION  : AJOUT D'UN OBJET A UNE STRUCTURE DE BLOC
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   BLOC         |<-->| NOM FORTRAN DU BLOC
C |   OBJ          | -->| NOUVEL OBJET.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : AUCUN
C
C***********************************************************************
C
      USE BIEF, EX_ADDBLO => ADDBLO
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ), INTENT(INOUT)      :: BLOC
      TYPE(BIEF_OBJ), INTENT(IN), TARGET :: OBJ
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     ACCROISSEMENT DU NOMBRE D'OBJETS DANS LE BLOC
C
      BLOC%N = BLOC%N + 1
      IF(BLOC%N.GT.BLOC%MAXBLOCK) THEN
       IF(LNG.EQ.1) THEN
         WRITE(LU,*) 'ADDBLO : ',OBJ%NAME,' TROP PETIT'
         WRITE(LU,*) '         AUGMENTER MAXBLOCK DANS ALLBLO'
         WRITE(LU,*) '         (ACTUELLEMENT : ',BLOC%MAXBLOCK,')'
       ENDIF
       IF(LNG.EQ.2) THEN
         WRITE(LU,*) 'ADDBLO : ',OBJ%NAME,' TOO SMALL'
         WRITE(LU,*) '         INCREASE MAXBLOCK IN ALLBLO'
         WRITE(LU,*) '         (CURRENTLY : ',BLOC%MAXBLOCK,')'
       ENDIF
       STOP
      ENDIF
C
C     AFFECTATION DE LA CIBLE OBJ AU POINTEUR DE RANG BLOC%N
C
      BLOC%ADR(BLOC%N)%P => OBJ
C
C-----------------------------------------------------------------------
C
C      IF(LNG.EQ.1) THEN
C        WRITE(LU,*) 'ADDBLO : ',OBJ%NAME,' AJOUTE A ',BLOC%NAME
C      ENDIF
C      IF(LNG.EQ.2) THEN
C        WRITE(LU,*) 'ADDBLO : ',OBJ%NAME,'  ADDED TO ',BLOC%NAME
C      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
