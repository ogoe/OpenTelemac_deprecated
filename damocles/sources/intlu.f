!                    **********************
                     INTEGER FUNCTION INTLU
!                    **********************
!
     &( ICOL , LIGNE )
!
!***********************************************************************
! DAMOCLES   V6P0                                   21/08/2010
!***********************************************************************
!
!brief    DECODES AN INTEGER, FROM COLUMN ICOL+1 OF THE LINE.
!+             MOVES THE POINTER ICOL TO THE LAST DECODED CHARACTER.
!+             IF THE STRING IS NOT COMPLETE, GOES TO THE NEXT LINE
!+             IF NEED BE.
!+             MOVES THE POINTER ICOL TO THE LAST DECODED CHARACTER.
!+             OR TO ICOL=0 IF THE NEXT LINE WAS READ.
!
!note     PORTABILITY : IBM,CRAY,HP,SUN
!
!warning  IF THE VALUE READ IS NOT AN INTEGER, COULD YIELD A
!+            NON-CONTROLLED ERROR BY THE PROGRAM
!
!history  J.M. HERVOUET (LNH); A. YESSAYAN
!+        30/09/1993
!+        V5P1
!+
!
!history  O. QUIQUEMPOIX (LNH)
!+        15/12/1993
!+
!+
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| ICOL           |<->| POSITION COURANTE DU POINTEUR DANS LA LIGNE
!| LIGNE          |<->| LIGNE EN COURS DE DECODAGE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IMPLICIT NONE
!
      INTEGER          ICOL
      CHARACTER*(*)    LIGNE
!
      INTEGER          NEXT,PREVAL
      EXTERNAL NEXT,PREVAL
!
      INTEGER          LNG,LU
      INTEGER          NLIGN,LONGLI
      INTEGER          NFIC
      LOGICAL          ERREUR , RETOUR
!
!-----------------------------------------------------------------------
!
      INTRINSIC DLOG10,DBLE,INT,CHAR
!
      INTEGER          I1,I2,ILONG,ISIGNE,IVAL,JD1,I3
      LOGICAL          LUFIC,LISUIV
      CHARACTER*1      CDEB,TABUL
      CHARACTER*3      LLONG
      CHARACTER*72     LIGNE2,FORMA
!
!-----------------------------------------------------------------------
!
      COMMON / DCINFO / LNG,LU
      COMMON / DCRARE / ERREUR , RETOUR
      COMMON / DCMLIG / NLIGN , LONGLI
      COMMON / DCCHIE / NFIC
!
!***********************************************************************
!
      LUFIC = .FALSE.
      LISUIV = .FALSE.
      TABUL = CHAR(9)
!
      I1     = NEXT( ICOL+1 , LIGNE )
!
!     //// DECODES THE SIGN IF NEED BE ////
!
      IF ( LIGNE(I1:I1).EQ.'-' ) THEN
           ISIGNE = -1
           I1     =   NEXT ( I1+1      , LIGNE )
      ELSE IF ( LIGNE(I1:I1).EQ.'+' ) THEN
           ISIGNE = +1
           I1     =   NEXT ( I1+1      , LIGNE )
      ELSE
           ISIGNE = +1
      ENDIF
!
!     //// SEEKS THE FIRST WHITE CHARACTER FOLLOWING THE NUMBER ////
!                       OR A SEPARATOR ';'
!
      I2 = PREVAL (  I1  , LIGNE ,  ' ' , ';' , TABUL)
!
!     CASE WHERE THE INTEGER DOES NOT FINISH ON THE LINE                                                                                                                                                                                                                                                                                                                                                                                                                                              LINE
!
      IF (I2.GT.LONGLI) THEN
         LUFIC=.TRUE.
         READ(NFIC,END=900,ERR=998,FMT='(A)') LIGNE2
         CDEB = LIGNE2(1:1)
         IF (CDEB.EQ.'0'.OR.CDEB.EQ.'1'.OR.CDEB.EQ.'2'.OR.
     &       CDEB.EQ.'3'.OR.CDEB.EQ.'4'.OR.CDEB.EQ.'5'.OR.
     &       CDEB.EQ.'6'.OR.CDEB.EQ.'7'.OR.CDEB.EQ.'8'.OR.
     &       CDEB.EQ.'9'.OR.CDEB.EQ.'.') THEN
            LISUIV = .TRUE.
            I3=1
            I3=PREVAL(I3,LIGNE2 , ' ' , ';', TABUL)
            IF (I1.LE.LONGLI) THEN
              LIGNE = LIGNE(I1:LONGLI)//LIGNE2(1:I3)
            ELSE
              LIGNE =LIGNE2(1:I3)
            ENDIF
            I2 = LONGLI-I1+1+I3
            I1 = 1
         ENDIF
       ENDIF
       GOTO 910
!
 900  CONTINUE
      RETOUR = .TRUE.
 910  CONTINUE
!     ACCEPTS THE CASE WHERE A USER WRITES AN INTEGER IN
!     REAL FORM WITH A POINT AT THE END
      IF(LIGNE(I2-1:I2-1).EQ.'.') THEN
        LIGNE(I2-1:I2-1)=' '
        I2 = I2 - 1
      ENDIF
!
!     ILONG: LENGTH OF THE INTEGER
      ILONG  = I2 - I1
!
!     //// DECODING FORMAT ////
!
      JD1 = 3 - INT(DLOG10(DBLE(ILONG)))
      WRITE ( LLONG , '(I3)' ) ILONG
!
      IF(I1.EQ.1) THEN
         WRITE (FORMA , 1101 )  LLONG(JD1:3)
      ELSE
         WRITE (FORMA , 1100 )  I1-1 , LLONG(JD1:3)
      ENDIF
!
!     ////  DECODES ////
!
      READ  ( LIGNE , FORMA , ERR=995 ) IVAL
      INTLU = ISIGNE * IVAL
!
!     //// UPDATES THE POINTER ////
!
      IF (LUFIC) THEN
        NLIGN = NLIGN + 1
        LIGNE = LIGNE2
        IF (LISUIV) THEN
          ICOL = I3-1
        ELSE
          ICOL = 0
        ENDIF
      ELSE
        ICOL = I2 - 1
      ENDIF
!
1100  FORMAT('(',I3,'X,I',A,')')
1101  FORMAT('(I',A,')')
!
!-----------------------------------------------------------------------
!
      RETURN
!
! TREATS THE ERRORS DUE TO THE INTERNAL READ FOR CONVERSION
!
995   CONTINUE
      IF(LNG.EQ.1) WRITE(6,996) NLIGN
      IF(LNG.EQ.2) WRITE(6,1996) NLIGN
      WRITE(6,*) LIGNE
996   FORMAT(1X,'ERREUR LIGNE ',1I6,', UN ENTIER EST ATTENDU : ',/)
1996  FORMAT(1X,'ERREUR LINE ',1I6,', INTEGER EXPECTED : ',/)
      ERREUR=.TRUE.
      RETURN
!
! TREATS THE ERRORS DUE TO FILE MISREADING
!
998   CONTINUE
      IF(LNG.EQ.1) WRITE(6,999) NFIC,NLIGN+1
      IF(LNG.EQ.2) WRITE(6,1999) NFIC,NLIGN+1
999   FORMAT(1X,'UNITE LOGIQUE ',1I2,'   ERREUR LIGNE ',1I6)
1999  FORMAT(1X,'LOGICAL UNIT ',1I2,'   ERROR LINE ',1I6)
      RETOUR = .TRUE.
      RETURN
!
      END