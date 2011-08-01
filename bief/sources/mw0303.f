C                       *****************
                        SUBROUTINE MW0303
C                       *****************
C
     *(OP, X , DA,TYPDIA,XAS,TYPEXT, Y,C,
     * IKLEM1,DIMIKM,LIMVOI,MXPTVS,NPMAX,NPOIN,TRAV)
C
C***********************************************************************
C BIEF VERSION 5.1           05/02/91    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C  FONCTION : OPERATIONS MATRICE VECTEUR POUR TRIANGLES P1
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES VECTEURS X,Y ET LA MATRICE M. LE RESULTAT
C   EST LE VECTEUR X.
C
C   CES OPERATIONS SONT DIFFERENTES SUIVANT LE TYPE DE DIAGONALE
C   ET LE TYPE DES TERMES EXTRADIAGONAUX.
C
C  OPERATIONS PROGRAMMEES :
C
C      OP = 'X=AY    '  : X = AY
C      OP = 'X=-AY   '  : X = -AY
C      OP = 'X=X+AY  '  : X = X + AY
C      OP = 'X=X-AY  '  : X = X - AY
C      OP = 'X=X+CAY '  : X = X + C AY
C      OP = 'X=TAY   '  : X = TA Y (TRANSPOSEE DE A)
C      OP = 'X=-TAY  '  : X = - TA Y (- TRANSPOSEE DE A)
C      OP = 'X=X+TAY '  : X = X + TA Y
C      OP = 'X=X-TAY '  : X = X - TA Y
C      OP = 'X=X+CTAY'  : X = X + C TA Y
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!_______________________________________________
C !      OP        ! -->!OPERATION A EFFECTUER
C !      X         !<-- !VECTEUR IMAGE
C !      DA        ! -->! DIAGONALE DE LA MATRICE
C !      TYPDIA    ! -->! TYPE DE LA DIAGONALE (CHAINE DE CARACTERES)
C !                !    ! TYPDIA = 'Q' : DIAGONALE QUELCONQUE
C !                !    ! TYPDIA = 'I' : DIAGONALE IDENTITE.
C !                !    ! TYPDIA = '0' : DIAGONALE NULLE.
C !      XAS       ! -->! TERMES EXTRA-DIAGONAUX DE LA MATRICE
C !      TYPEXT    ! -->! TYPE DES TERMES EXTRADIAGONAUX
C !                !    ! TYPEXT = 'Q' : QUELCONQUES.
C !                !    ! TYPEXT = 'S' : SYMETRIQUES.
C !                !    ! TYPEXT = '0' : NULS.
C !      Y         ! -->! VECTEUR OPERANDE
C !      C         ! -->! CONSTANTE DONNEE
C !      IKLEM1    ! -->!
C !      DIMIKM    ! -->! PREMIERE DIMENSION DE IKLEM1
C !      LIMVOI    ! -->!
C !      MXPTVS    ! -->!
C !      NPMAX     ! -->! NOMBRE MAXIMUM DE POINTS.
C !      NPOIN     ! -->! NOMBRE DE POINTS.
C !      TRAV      ! -->! TABLEAU DE TRAVAIL.
C !________________!____!_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   : ASSVEC , OV
C
C***********************************************************************
C
      USE BIEF, EX_MW0303 => MW0303
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: DIMIKM,MXPTVS,NPMAX,NPOIN
      INTEGER, INTENT(IN) :: IKLEM1(DIMIKM,4,2),LIMVOI(MXPTVS,2)
C
      DOUBLE PRECISION, INTENT(INOUT) :: X(*),TRAV(*)
      DOUBLE PRECISION, INTENT(IN)    :: DA(*),Y(*)
      DOUBLE PRECISION, INTENT(IN)    :: XAS(*),C
C
      CHARACTER(LEN=8), INTENT(IN)    :: OP
      CHARACTER(LEN=1), INTENT(IN)    :: TYPDIA,TYPEXT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
      DOUBLE PRECISION Z(1)
C
C-----------------------------------------------------------------------
C
C   TRAITEMENT SPECIFIQUE A LA TRANSPOSITION :
C
      I = 1
      IF(OP(3:3).EQ.'T'.OR.OP(4:4).EQ.'T'.OR.
     *   OP(5:5).EQ.'T'.OR.OP(6:6).EQ.'T') I = 3
C
C-----------------------------------------------------------------------
C
C   PRODUIT MATRICE VECTEUR SIMPLE FONCTION DE LA FORME DE LA MATRICE :
C
      IF(TYPEXT(1:1).EQ.'S'.OR.TYPEXT(1:1).EQ.'Q') THEN
C
        IF(TYPEXT(1:1).EQ.'Q') THEN
        CALL OPASS('X=WY    ',TRAV,XAS,IKLEM1(1,I,1),
     *             Y,IKLEM1(1,I+1,1),LIMVOI,MXPTVS,NPMAX)
        ELSEIF(TYPEXT(1:1).EQ.'S') THEN
        CALL OPASS('X=WY    ',TRAV,XAS,IKLEM1(1,I,2),
     *             Y,IKLEM1(1,I+1,2),LIMVOI,MXPTVS,NPMAX)
        ENDIF
C
        IF(TYPDIA(1:1).EQ.'Q') THEN
          CALL OV ('X=X+YZ  ', TRAV , Y , DA , C , NPOIN )
        ELSEIF(TYPDIA(1:1).EQ.'I') THEN
          CALL OV ('X=X+Y   ', TRAV , Y , Z , C , NPOIN )
        ELSEIF(TYPDIA(1:1).NE.'0') THEN
          IF (LNG.EQ.1) WRITE(LU,2000) TYPDIA
          IF (LNG.EQ.2) WRITE(LU,2001) TYPDIA
          CALL PLANTE(0)
          STOP
        ENDIF
C
      ELSEIF(TYPEXT(1:1).EQ.'0') THEN
C
         IF(TYPDIA(1:1).EQ.'Q') THEN
           CALL OV ('X=YZ    ', TRAV , Y , DA , C , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'I') THEN
           CALL OV ('X=Y     ', TRAV , Y , Z , C , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'0') THEN
           CALL OV ('X=C     ', TRAV , Y , Z , 0.D0 , NPOIN )
         ELSE
            IF (LNG.EQ.1) WRITE(LU,2000) TYPDIA
            IF (LNG.EQ.2) WRITE(LU,2001) TYPDIA
            CALL PLANTE(0)
            STOP
         ENDIF
C
      ELSE
C
         IF (LNG.EQ.1) WRITE(LU,1000) TYPEXT
         IF (LNG.EQ.2) WRITE(LU,1001) TYPEXT
         CALL PLANTE(0)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C   OPTION D'OPERATION :
C
      IF(OP(1:8).EQ.'X=AY    '.OR.OP(1:8).EQ.'X=TAY   ') THEN
         CALL OV ('X=Y     ', X , TRAV , Z , C , NPOIN )
      ELSEIF(OP(1:8).EQ.'X=-AY   '.OR.OP(1:8).EQ.'X=-TAY  ') THEN
         CALL OV ('X=-Y    ', X , TRAV , Z , C , NPOIN )
      ELSEIF(OP(1:8).EQ.'X=X+AY  '.OR.OP(1:8).EQ.'X=X+TAY ') THEN
         CALL OV ('X=X+Y   ', X , TRAV , Z , C , NPOIN )
      ELSEIF(OP(1:8).EQ.'X=X-AY  '.OR.OP(1:8).EQ.'X=X-TAY ') THEN
         CALL OV ('X=X-Y   ', X , TRAV , Z , C , NPOIN )
      ELSEIF(OP(1:8).EQ.'X=X+CAY '.OR.OP(1:8).EQ.'X=X+CTAY') THEN
         CALL OV ('X=X+CY  ', X , TRAV , Z , C , NPOIN )
      ELSEIF(OP(1:8).EQ.'X=CAY   ') THEN
         CALL OV ('X=CY    ', X , TRAV , Z , C , NPOIN ) 
      ELSE
         IF (LNG.EQ.1) WRITE(LU,3000) OP
         IF (LNG.EQ.2) WRITE(LU,3001) OP
         CALL PLANTE(0)
         STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
C
1000  FORMAT(1X,'MW0303 (BIEF) : TERMES EXTRADIAG. TYPE INCONNU: ',A1)
1001  FORMAT(1X,'MW0303 (BIEF) : EXTRADIAG. TERMS  UNKNOWN TYPE : ',A1)
2000  FORMAT(1X,'MW0303 (BIEF) : DIAGONALE : TYPE INCONNU: ',A1)
2001  FORMAT(1X,'MW0303 (BIEF) : DIAGONAL : UNKNOWN TYPE : ',A1)
3000  FORMAT(1X,'MW0303 (BIEF) : OPERATION INCONNUE : ',A8)
3001  FORMAT(1X,'MW0303 (BIEF) : UNKNOWN OPERATION : ',A8)
C
      END
