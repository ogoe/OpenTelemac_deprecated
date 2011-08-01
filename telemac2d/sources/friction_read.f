C                       ************************
                        SUBROUTINE FRICTION_READ
C                       ************************
C
     & (NCOF, NZONMX, ITURB, LISRUG, LINDNER, NOMCOF, NZONES, FRTAB)
C
C***********************************************************************
C  TELEMAC-2D VERSION 5.5                 J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C 20/04/04 : subroutine written by F. Huvelin
C
C
! ----------------------------------------------- !
!             FRICTION FILE READ                  !
! ----------------------------------------------- !
C
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
      USE FRICTION_DEF
      USE INTERFACE_TELEMAC2D, EX_FRICTION_READ => FRICTION_READ
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER,            INTENT(IN)    :: NCOF, NZONMX
      INTEGER,            INTENT(IN)    :: ITURB, LISRUG
      LOGICAL,            INTENT(IN)    :: LINDNER
      CHARACTER(LEN=144), INTENT(IN)    :: NOMCOF
      INTEGER,            INTENT(OUT)   :: NZONES
      TYPE(FRICTION_OBJ), INTENT(INOUT) :: FRTAB
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER          :: I, J, LOOP, NLOOP, N, IZ1, IZ2
      INTEGER          :: IZONE, TYP, LINE
      DOUBLE PRECISION :: R1, R2
      CHARACTER*4      :: LAW
      CHARACTER*20     :: CHAINE(10)
!
!=======================================================================!
!=======================================================================!
!                               PROGRAMME                               !
!=======================================================================!
!=======================================================================!
!
      ! Check that there is a friction file
      ! -----------------------------------
      IF (NOMCOF(1:1) == ' ') THEN
         IF (LNG == 1) WRITE(LU,1)
         IF (LNG == 2) WRITE(LU,2)
         CALL PLANTE(1)
         STOP
      ENDIF
!
      REWIND NCOF
!
      DO I=1,NZONMX
         DO J = 1, 2
            FRTAB%ADR(I)%P%GNUMB(J) = 0
            FRTAB%ADR(I)%P%RTYPE(J) = 0
            FRTAB%ADR(I)%P%RCOEF(J) = 0.D0
            FRTAB%ADR(I)%P%NDEF (J) = 0.D0
         ENDDO
         FRTAB%ADR(I)%P%DP          = 0.D0
         FRTAB%ADR(I)%P%SP          = 0.D0
      ENDDO
!
      ! K-epsilon : read parameters for 2 laws (bottom and boundary conditions)
      ! Else      : read parameters for 1 law  (bottom)
      ! ----------------------------------------------------------------------
      IF ((ITURB == 3).AND.(LISRUG == 2)) THEN
         NLOOP=2
      ELSE
         NLOOP=1
      ENDIF
!
      ! Listing
      ! -------
      WRITE(LU,3)
      WRITE(LU,4)
!
      ! Read data for each zones
      ! ------------------------
      IZONE = 0
      LINE = 0
      DO I = 1, NZONMX
!
         ! Find the next zone index and law
         ! --------------------------------
         CALL FRICTION_SCAN(NCOF,NOMCOF,TYP,LINE)
!
         ! End of file => exit
         ! -------------------
         IF (TYP == 3) THEN
            NZONES = IZONE
            EXIT
         ENDIF
!
         ! Increment number of zones
         ! -------------------------
         IZONE = IZONE + 1
!
         ! Not enough zones allocated
         ! --------------------------
         IF (IZONE > NZONMX) THEN
            IF (LNG == 1) WRITE(LU,5)
            IF (LNG == 2) WRITE(LU,6)
            CALL PLANTE(1)
            STOP
         ENDIF
!
         ! Read and save parameters of the zone
         ! ------------------------------------
         IF (TYP == 1) N = 1
         IF (TYP == 2) N = 4
!
         DO LOOP = 1, NLOOP
!
            IF (LOOP == 2) BACKSPACE(NCOF)

            ! Read the name of the law and numbering of the zone
            ! --------------------------------------------------
            CHAINE(1) = ' '
            IF (TYP == 1) THEN
               READ(NCOF,*,END=999,ERR=998) IZ1,(CHAINE(J),J=2,N),LAW
               IZ2 = IZ1
            ELSE IF (TYP == 2) THEN
               READ(NCOF,*,END=999,ERR=998) CHAINE(1),IZ1,CHAINE(3),IZ2,
     &                                      (CHAINE(J),J=5,N),LAW
            ENDIF
!
            ! Local-Global number of the zone
            ! -------------------------------
            IF (LOOP == 1) FRTAB%ADR(I)%P%GNUMB(1) = IZ1
            IF (LOOP == 1) FRTAB%ADR(I)%P%GNUMB(2) = IZ2
!
            BACKSPACE(NCOF)
            CALL MAJUS(LAW)
!
            ! Find the law and the number of parameters to read
            ! -------------------------------------------------
            SELECT CASE (LAW)
!
            CASE('NOFR')
               READ(NCOF,*,END=999,ERR=900) (CHAINE(J),J=1,N),LAW
               FRTAB%ADR(I)%P%RTYPE(LOOP) = 0
               FRTAB%ADR(I)%P%RCOEF(LOOP) = 0.D0
               FRTAB%ADR(I)%P%NDEF (LOOP) = 0.D0
               N = N + 1
            CASE('HAAL')
               READ(NCOF,*,END=999,ERR=901) (CHAINE(J),J=1,N),LAW,R1
               FRTAB%ADR(I)%P%RTYPE(LOOP) = 1
               FRTAB%ADR(I)%P%RCOEF(LOOP) = R1
               FRTAB%ADR(I)%P%NDEF (LOOP) = 0.D0
!
               ! The boundary coefficient coef must be the same as the bottom
               ! -----------------------------------------------------------    
               IF ( ( ((N>1).AND.(TYP==1)).OR.((N>4).AND.(TYP==2)) )
     &             .AND.
     &              (FRTAB%ADR(I)%P%RCOEF(1)/=FRTAB%ADR(I)%P%RCOEF(2))
     &            ) THEN
                  IF (LNG==1) WRITE(LU,15) NOMCOF,I,
     &                                     FRTAB%ADR(I)%P%RCOEF(1)
                  IF (LNG==2) WRITE(LU,16) NOMCOF,I,
     &                                     FRTAB%ADR(I)%P%RCOEF(1)
                  CALL PLANTE(1)
                  STOP
               ENDIF
               N = N + 2
            CASE('CHEZ')
               READ(NCOF,*,END=999,ERR=901) (CHAINE(J),J=1,N),LAW,R1
               FRTAB%ADR(I)%P%RTYPE(LOOP) = 2
               FRTAB%ADR(I)%P%RCOEF(LOOP) = R1
               FRTAB%ADR(I)%P%NDEF (LOOP) = 0.D0
               N = N + 2
            CASE('STRI')
               READ(NCOF,*,END=999,ERR=901) (CHAINE(J),J=1,N),LAW,R1
               FRTAB%ADR(I)%P%RTYPE(LOOP) = 3
               FRTAB%ADR(I)%P%RCOEF(LOOP) = R1
               FRTAB%ADR(I)%P%NDEF (LOOP) = 0.D0
               N = N + 2
            CASE('MANN')
               READ(NCOF,*,END=999,ERR=901) (CHAINE(J),J=1,N),LAW,R1
               FRTAB%ADR(I)%P%RTYPE(LOOP) = 4
               FRTAB%ADR(I)%P%RCOEF(LOOP) = R1
               FRTAB%ADR(I)%P%NDEF (LOOP) = 0.D0
               N = N + 2
            CASE('NIKU')
               READ(NCOF,*,END=999,ERR=901) (CHAINE(J),J=1,N),LAW,R1
               FRTAB%ADR(I)%P%RTYPE(LOOP) = 5
               FRTAB%ADR(I)%P%RCOEF(LOOP) = R1
               FRTAB%ADR(I)%P%NDEF (LOOP) = 0.D0
               N = N + 2
            CASE('LOGW')
               IF (((N == 1).AND.(TYP == 1))  .OR.
     &             ((N == 4).AND.(TYP == 2))) THEN
                  IF (LNG == 1) WRITE(LU,7) NOMCOF, I
                  IF (LNG == 2) WRITE(LU,8) NOMCOF, I
                  CALL PLANTE(1)
                  STOP
               ENDIF
               READ(NCOF,*,END=999,ERR=901) (CHAINE(J),J=1,N),LAW,R1
               FRTAB%ADR(I)%P%RTYPE(LOOP) = 6
               FRTAB%ADR(I)%P%RCOEF(LOOP) = R1
               FRTAB%ADR(I)%P%NDEF (LOOP) = 0.D0
               N = N + 2
            CASE('COWH')
               READ(NCOF,*,END=999,ERR=907) (CHAINE(J),J=1,N),LAW,R1,R2
               FRTAB%ADR(I)%P%RTYPE(LOOP) = 7
               FRTAB%ADR(I)%P%RCOEF(LOOP) = R1
               FRTAB%ADR(I)%P%NDEF (LOOP) = R2
               N = N + 3
            CASE DEFAULT
               IF (LNG==1) WRITE(LU, 9) LINE, IZ1, IZ2, LAW
               IF (LNG==2) WRITE(LU,10) LINE, IZ1, IZ2, LAW
               CALL PLANTE(0)
               STOP
            END SELECT
         ENDDO
!
         ! Read non-submerged coefficient if needed
         ! ----------------------------------------  
         IF (LINDNER) THEN
            BACKSPACE(NCOF)
            READ(NCOF,*,END=999,ERR=888) (CHAINE(J),J=1,N),R1,R2
            FRTAB%ADR(I)%P%DP = R1
            FRTAB%ADR(I)%P%SP = R2
         ELSE
            FRTAB%ADR(I)%P%DP = 0.D0
            FRTAB%ADR(I)%P%SP = 0.D0
         ENDIF

         WRITE(LU,11) FRTAB%ADR(I)%P%GNUMB(1),
     &                FRTAB%ADR(I)%P%GNUMB(2),
     &                FRTAB%ADR(I)%P%RTYPE(1),
     &                FRTAB%ADR(I)%P%RCOEF(1),
     &                FRTAB%ADR(I)%P%NDEF (1),
     &                FRTAB%ADR(I)%P%RTYPE(2),
     &                FRTAB%ADR(I)%P%RCOEF(2),
     &                FRTAB%ADR(I)%P%NDEF (2),
     &                FRTAB%ADR(I)%P%DP,
     &                FRTAB%ADR(I)%P%SP
      ENDDO
!
!
      ! End
      ! ---
      WRITE(LU,3)
!
      IF (LNG == 1) WRITE(LU,13) NZONES
      IF (LNG == 2) WRITE(LU,14) NZONES
!
      GOTO 997
!
      ! ============ !
      ! Error format !
      ! ============ !
!
1     FORMAT('PAS DE FICHIER DE DONNEES POUR LE FROTTEMENT')
2     FORMAT('NO FRICTION DATA FILE')

3     FORMAT('-------------------------------------------------------'
     &     , '------------------------------------------------------')
4     FORMAT('                      Bottom                         '
     &     , 'Boundary Condition              Non-submerged vegetation'
     &     ,/'no                    Law   RCoef        NDef        '
     &     , 'Law   RCoef        NDef         dp           sp')

5     FORMAT('NOMBRE DE ZONES DE FROTTEMENT DEFINI TROP NOMBREUSES'
     &     ,/'AUGMENTER LE NOMBRE DE ZONES MAXIMALES AVEC LE MOT-CLE :'
     &     ,/'NOMBRE MAXIMALE DE ZONES POUR LE FROTTEMENT')
6     FORMAT('TOO MANY NUMBER OF FRICTION ZONES DEFINED'
     &     ,/'INCREASED THE NUMBER OF MAXIMAL ZONES WITH THE KEYWORD :'
     &     ,/'MAXIMUM NUMBER OF ZONES FOR THE FRICTION')

7     FORMAT('FICHIER DE DONNEES POUR LE FROTTEMENT : ',a144
     &      ,/'ZONE : ',i9
     &      ,/'LA LOI LOG NE PEUT ETRE UTILISEE SUR LE FOND')
8     FORMAT('FRICTION DATA FILE : ',a144
     &      ,/'ZONE : ',i9
     &      ,/'LOG LAW CAN''T BE USED FOR THE BOTTOM')

9     FORMAT('FICHIER DE DONNEES POUR LE FROTTEMENT'
     &     ,/'ERREUR DE LECTURE LIGNE : ',i10
     &     ,/'ZONE DE ',i10,' A ',i10
     &     ,/'LOI ',a4)

10    FORMAT('FRICTION DATA FILE'
     &     ,/'READ ERROR LINE',i10
     &     ,/'ZONE FROM ',i10,' TO ',i10
     &     ,/'LAW ',a4)

11    FORMAT(2(1x,i9),1x,i4,2(1x,e12.4),1x,i4,4(1x,e12.4))
!
13    FORMAT(i5,' types de zones definies')
14    FORMAT(i5,' zones type specifications') 

15    FORMAT('FICHIER DE DONNEES POUR LE FROTTEMENT : ',a144
     &      ,/'ZONE : ',i9
     &      ,/'LE COEFFICIENT DE FROTTEMENT DE LA LOI DE HAALAND POUR'
     &      ,/'LES CONDITIONS LIMITES DOIT ETER EGALE A CELUI DU FOND :'
     &      , e12.4)
16    FORMAT('FRICTION DATA FILE : ',a144
     &      ,/'ZONE : ',i9
     &      ,/'FRICTION COEFFICIENT OF HAALAND LAW FOR'
     &      ,/'BOUNDARY CONDITION MUST BE THE SAME AS THE BOTTOM :'
     &      , e12.4)
!
      ! End of File
      ! -----------
999   CONTINUE
      IF (LNG.EQ.1) THEN
         WRITE(LU,*) 'FICHIER DE DONNEES POUR LE FROTTEMENT : ',NOMCOF
         WRITE(LU,*) 'FIN DE FICHIER ANORMALE'
         WRITE(LU,*) 'VERIFIER QUE TOUTES LES VALEURS SONT ENTREES'
      ENDIF
      IF (LNG.EQ.2) THEN
         WRITE(LU,*) 'FRICTION DATA FILE : ',NOMCOF
         WRITE(LU,*) 'ABNORMAL END OF FILE'
         WRITE(LU,*) 'CHECK ALL VALUE HAVE BEEN WRITTEN'
      ENDIF
      CALL PLANTE(1)
      STOP
!
      ! Index and first law of the zone
      ! -------------------------------
998   CONTINUE
      IF (LNG.EQ.1) THEN
         WRITE(LU,*)'FICHIER DE DONNEES POUR LE FROTTEMENT : ',NOMCOF
         WRITE(LU,*)'ERREUR DE LECTURE ZONE : ',CHAINE(1)
      ENDIF
      IF (LNG.EQ.2) THEN
         WRITE(LU,*) 'FRICTION DATA FILE : ',NOMCOF
         WRITE(LU,*) 'READ ERROR ZONE : ',CHAINE(1)
      ENDIF
      CALL PLANTE(0)
      STOP

      ! Non-submerged vegetation parameter
      ! ----------------------------------
888   CONTINUE
      IF (LNG.EQ.1) THEN
         WRITE(LU,*) 'FICHIER DE DONNEES POUR LE FROTTEMENT'
         WRITE(LU,*)'ERREUR DE LECTURE POUR DP ET SP, ZONE : ',CHAINE(1)
      ENDIF
      IF (LNG.EQ.2) THEN
         WRITE(LU,*) 'FRICTION DATA FILE'
         WRITE(LU,*) 'READ ERROR FOR DP AND SP, ZONE : ',CHAINE(1)
      ENDIF
      CALL PLANTE(0)
      STOP
!
      ! No friction Law
      ! ---------------
900   CONTINUE
      IF (LNG.EQ.1) THEN
         WRITE(LU,*)'FICHIER DE DONNEES POUR LE FROTTEMENT'
         WRITE(LU,*) 'ERREUR DE LECTURE ZONE  : ',CHAINE(1)
         IF ((ITURB==3).AND.(LISRUG==2)) THEN
            IF (LOOP==1) WRITE(LU,*) 'Pour la 1ere loi definir '//
     &                               'seulement le nom de la loi : NOFR'
            IF (LOOP==2) WRITE(LU,*) 'Pour la 2nde loi definir'//
     &                               'seulement le nom de la loi : NOFR'
         ELSE
            WRITE(LU,*) 'Definir seulement le nom de la loi'
         ENDIF
      ENDIF
      IF (LNG.EQ.2) THEN
         WRITE(LU,*) 'FRICTION DATA FILE'
         WRITE(LU,*) 'READ ERROR ZONE : ',CHAINE(1)
         IF ((ITURB==3).AND.(LISRUG==2)) THEN
            IF (LOOP==1) WRITE(LU,*) 'For the 1st law define '//
     &                               'only the name of the law : NOFR'
            IF (LOOP==2) WRITE(LU,*) 'For the 2nd law define '//
     &                               'only the name of the law : NOFR'
         ELSE
            WRITE(LU,*) 'Define only the name of the law : NOFR'
         ENDIF
      ENDIF
      CALL PLANTE(1)
      STOP
!
      ! Haaland-Chezy-Strickler-Manning-Nikuradse-Log Wall Laws
      ! -------------------------------------------------------
901   CONTINUE
      IF (LNG.EQ.1) THEN
         WRITE(LU,*)'FICHIER DE DONNEES POUR LE FROTTEMENT : ',NOMCOF
         WRITE(LU,*) 'ERREUR DE LECTURE ZONE: ',CHAINE(1)
         IF ((ITURB==3).AND.(LISRUG==2)) THEN
            IF (LOOP==1) WRITE(LU,*) 'Pour la 1ere loi definir ' //
     &                               'le nom de la loi : ',LAW,' et'//
     &                               ' le coefficient de frottement'
            IF (LOOP==1) WRITE(LU,*) 'Pour la 2nde loi definir ' //
     &                               'le nom de la loi : ',LAW,' et'//
     &                               ' le coefficient de frottement'
         ELSE
            WRITE(LU,*) 'Definir le nom de la loi : ',LAW,' et'//
     &                  ' le coefficient de frottement'
         ENDIF
      ENDIF
!
      IF (LNG.EQ.2) THEN
         WRITE(LU,*) 'FRICTION DATA FILE : ',NOMCOF
         WRITE(LU,*) 'READ ERROR ZONE : ',CHAINE(1)
         IF ((ITURB==3).AND.(LISRUG==2)) THEN
            IF (LOOP==1) WRITE(LU,*) 'For the 1st law define '//
     &                               'the name of the law : ',LAW//
     &                               ' and the friction coefficient'
            IF (LOOP==2) WRITE(LU,*) 'For the 2nd law define '//
     &                               'the name of the law : ',LAW//
     &                               ' and the friction coefficient'
         ELSE
            WRITE(LU,*) 'Define the name of the law : ',LAW//
     &                  ' and the friction coefficient'
         ENDIF
      ENDIF
      CALL PLANTE(1)
      STOP
!
      ! Colebrook White Law
      ! -------------------
907   CONTINUE
      IF (LNG.EQ.1) THEN
         WRITE(LU,*)'FICHIER DE DONNEES POUR LE FROTTEMENT : ',NOMCOF
         WRITE(LU,*) 'ERREUR DE LECTURE ZONE : ',CHAINE(1)
         IF ((ITURB==3).AND.(LISRUG==2)) THEN
            IF (LOOP==1) WRITE(LU,*) 'Pour la 1ere loi definir ' //
     &                               'le nom de la loi : ',LAW,' et'//
     &                               ' le coefficient de frottement'//
     &                               ' et le Manning'
            IF (LOOP==1) WRITE(LU,*) 'Pour la 2nde loi definir ' //
     &                               'le nom de la loi : ',LAW,' et'//
     &                               ' le coefficient de frottement'//
     &                               ' et le Manning'
         ELSE
            WRITE(LU,*) 'Definir le nom de la loi : ',LAW,' et'//
     &                  ' le coefficient de frottement et le Manning'
         ENDIF
      ENDIF

      IF (LNG.EQ.2) THEN
         WRITE(LU,*) 'FRICTION DATA FILE : ',NOMCOF
         WRITE(LU,*) 'READ ERROR ZONE : ',CHAINE(1)
         IF ((ITURB==3).AND.(LISRUG==2)) THEN
            IF (LOOP==1) WRITE(LU,*) 'For the 1st law define '//
     &                               'the name of the law : ',LAW//
     &                               ' and the friction coefficient'//
     &                               ' and default Manning'
            IF (LOOP==2) WRITE(LU,*) 'For the 2nd law define '//
     &                               'the name of the law : ',LAW//
     &                               ' and the friction coefficient'//
     &                               ' and default Manning'
         ELSE
            WRITE(LU,*) 'Define the name of the law : ',LAW//
     &                  ' and the friction coefficient'//
     &                  ' and default Manning'
         ENDIF
      ENDIF
      CALL PLANTE(1)
      STOP     
!
997   CONTINUE
!
!=======================================================================!
!=======================================================================!
!
      RETURN
      END
