!                    ************************
                     SUBROUTINE FRICTION_SCAN
!                    ************************
!
     &(NCOF,NOMCOF,TYP,LINE)
!
!***********************************************************************
! TELEMAC2D   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    READS FRICTION FILE.
!
!history  F. HUVELIN
!+        20/04/2004
!+
!+
!
!history  J-M HERVOUET (LNHE)
!+
!+        V5P5
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
!| LINE           |<--| LINE READ IN THE FRICTION FILE
!| NCOF           |-->| LOGICAL UNIT OF FRICTION FILE
!| NOMCOF         |-->| NAME OF FRICTION FILE
!| TYP            |<--| 2: STARTING LINE
!|                |   | 3: ENDING LINE
!|                |   | 1: CURRENT LINE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER,       INTENT(IN)    :: NCOF
      CHARACTER(LEN=144), INTENT(IN)    :: NOMCOF
      INTEGER,       INTENT(INOUT) :: LINE
      INTEGER,       INTENT(OUT)   :: TYP
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      CHARACTER*20                 :: C
!
!=======================================================================!
!=======================================================================!
!                               PROGRAMME                               !
!=======================================================================!
!=======================================================================!
!
      DO
        READ(NCOF,*,END=989,ERR=988) C
        LINE = LINE + 1
        IF (C(1:1).NE.'*') EXIT
      ENDDO
!
      CALL MAJUS(C)
!
      IF(C(1:4).EQ.'FROM'.OR.C(1:3).EQ.'VON'.OR.C(1:2).EQ.'DE') THEN
        TYP = 2
      ELSEIF(C(1:3).EQ.'END'.OR.C(1:3).EQ.'FIN') THEN
        TYP = 3
      ELSE
        TYP = 1
      ENDIF
!
      BACKSPACE(NCOF)
!
      GOTO 987
!
! -------------------------------------------------------------- !
!                         WARNING MESSAGE                        !
! -------------------------------------------------------------- !
!
989   CONTINUE
      IF (LNG.EQ.1) THEN
         WRITE(LU,*) 'FICHIER DE DONNEES POUR LE FROTTEMENT : ',NOMCOF
         WRITE(LU,*) 'FIN DE FICHIER ANORMALE'
      ENDIF
      IF (LNG.EQ.2) THEN
         WRITE(LU,*) 'FRICTION DATA FILE : ',NOMCOF
         WRITE(LU,*) 'ABNORMAL END OF FILE'
      ENDIF
      CALL PLANTE(1)
      STOP
!
988   CONTINUE
      IF (LNG.EQ.1) THEN
         WRITE(LU,*) 'FICHIER DE DONNEES POUR LE FROTTEMENT : ',NOMCOF
         WRITE(LU,*) 'ERREUR DE LECTURE'
         WRITE(LU,*) 'ERREUR LORS DE LA LECTURE DE : ',C
      ENDIF
      IF (LNG.EQ.2) THEN
         WRITE(LU,*) 'FRICTION DATA FILE : ',NOMCOF
         WRITE(LU,*) 'READ ERROR'
         WRITE(LU,*) 'ERROR WHEN READING: ',C
      ENDIF
      CALL PLANTE(1)
      STOP
!
! -------------------------------------------------------------- !
! -------------------------------------------------------------- !
!
987   CONTINUE
!
!=======================================================================!
!=======================================================================!
!
      RETURN
      END
