C                     ************************
                      SUBROUTINE FRICTION_SCAN
C                     ************************
C
     &(NCOF,NOMCOF,TYP,LINE)

C***********************************************************************
C  TELEMAC-2D VERSION 5.5              J-M HERVOUET (LNH) 01 30 87 80 18
C***********************************************************************
C
C 20/04/04 : subroutine written by F. Huvelin
C
C
C           -----------------------------------------------  
C                       READING FRICTION FILE                  
C           -----------------------------------------------  
C
C
C----------------------------------------------------------------------C
C                             ARGUMENTS                                C
C .________________.____.______________________________________________C
C |      NOM       |MODE|                   ROLE                       C
C |________________|____|______________________________________________C
C |                | => |                                              C
C |________________|____|______________________________________________C
C                    <=  input value                                   C
C                    =>  output value                                  C 
C ---------------------------------------------------------------------C
!
!=======================================================================!
!=======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS                !
!=======================================================================!
!=======================================================================!
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER,       INTENT(IN)    :: NCOF
      CHARACTER(LEN=144), INTENT(IN)    :: NOMCOF
      INTEGER,       INTENT(INOUT) :: LINE
      INTEGER,       INTENT(OUT)   :: TYP
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
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
         IF (C(1:1) /= '*') EXIT
      ENDDO
!
      CALL MAJUS(C)
!
      IF ((C(1:4) == 'FROM').OR.(C(1:3) == 'VON').OR.
     &    (C(1:2) == 'DE')) THEN
         TYP = 2
      ELSE IF ((C(1:3) == 'END').OR.(C(1:4) == 'ENDE').OR.
     &    (C(1:3) == 'FIN')) THEN
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
         WRITE(LU,*) 'ERROR FOR THE READING OF : ',C
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

