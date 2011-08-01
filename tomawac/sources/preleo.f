C                       *****************
                        SUBROUTINE PRELEO
C                       *****************
C
     *(XLEO,YLEO,NLEO,X,Y,NPOIN2,NOLEO)
C
C***********************************************************************
C  COWADIS  VERSION 1.0     1/02/95    F. MARCOS (LNH) 30 87 72 66
C***********************************************************************
C
C     FONCTION  : CHOISI LES POINTS DE CALCUL LES PLUS PROCHES
C                 DES POINTS DE SORTIE DEMANDES
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !   XLEO         ! -->! TABLEAU DES ABSCISSES DES POINTS DE SORTIE   !
C !   YLEO         ! -->! TABLEAU DES ORDONNEES DES POINTS DE SORTIE   !
C !   NLEO         ! -->! NOMBRE DE POINTS DE SORTIE                   !
C !   X            ! -->! ABSCISSES DES POINTS                         !
C !   Y            ! -->! ORDONNEES DES POINTS                         !
C !   NPOIN2       ! -->! NOMBRE DE POINTS 2D                          !
C !   NOLEO        !<-- ! TABLEAU DES NUMERO DES POINTS CHOISIS        !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : COW
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER I,ILEO,NLEO,NPOIN2
C
      DOUBLE PRECISION X(NPOIN2)  , Y(NPOIN2)
      DOUBLE PRECISION XLEO(NLEO)  , YLEO(NLEO)
      DOUBLE PRECISION DIST,DIST2
C
      INTEGER NOLEO(NLEO)
C
C-----------------------------------------------------------------------
C
      DO 10 ILEO=1,NLEO
        DIST=1.D99
        DO 20 I=1,NPOIN2
         DIST2=(XLEO(ILEO)-X(I))**2+(YLEO(ILEO)-Y(I))**2
         IF (DIST2.LT.DIST) THEN
             DIST=DIST2
             NOLEO(ILEO)=I
         ENDIF
20      CONTINUE
10    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
