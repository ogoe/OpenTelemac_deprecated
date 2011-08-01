C                        ****************
                         SUBROUTINE SOLVE
C                        ****************
C
     *(X, A,B,TB,CFG,INFOGR,MESH,AUX)
C
C***********************************************************************
C BIEF VERSION 5.9        18/02/08    J-M HERVOUET (LNHE) 01 30 87 80 18
C BIEF VERSION 6.0        19/03/10    C. DENIS (SINETICS)
C                                                        
C***********************************************************************
C
C  FONCTION : RESOLUTION D'UN SYSTEME LINEAIRE DE LA FORME A X = B
C
C  PAR DES METHODES ITERATIVES
C
C  AVEC PRECONDITIONNEMENT EVENTUEL
C
C  ATTENTION : POUR CERTAINS PRECONDITIONNEMENTS LA MATRICE A 
C              PEUT ETRE MODIFIEE.
C
C-----------------------------------------------------------------------
C                        CHOIX DE LA METHODE
C-----------------------------------------------------------------------
C  VALEUR DE METHOD   I    SIGNIFICATION      I
C-----------------------------------------------------------------------
C                     I                       I
C        1            I   GRADIENT CONJUGUE   I
C                     I                       I
C        2            I   RESIDU   CONJUGUE   I
C                     I                       I
C        3            I   GRADIENT CONJUGUE   I
C                     I  SUR EQUATION NORMALE I
C                     I                       I
C        4            I   ERREUR MINIMALE     I
C                     I                       I
C        5            I  GRADIENT CONJUGUE    I
C                     I       CARRE           I
C                     I                       I
C        6            I  GRADIENT CONJUGUE    I
C                     I       CARRE           I
C                     I   (VERSION CGSTAB)    I
C                     I                       I
C        7            I       GMRES           I
C                     I                       I
C        8            I   direct : YSMP       I
C                     I                       I
C        9            I   direct : MUMPS      I
C                     I                       I
C-----------------------------------------------------------------------
C
C                        PRECONDITIONNEMENT
C
C     NOTES     : CERTAINS PRECONDITIONNEMENTS SONT CUMULABLES
C                 (LES DIAGONAUX 2, 3 OU 5 AVEC LES AUTRES)
C                 POUR CETTE RAISON ON NE RETIENT QUE LES NOMBRES
C                 PREMIERS POUR DESIGNER LES PRECONDITIONNEMENTS.
C                 SI L'ON SOUHAITE EN CUMULER PLUSIEURS ON MET DANS
C                 LE PRODUIT DES NOMBRES CORRESPONDANTS.
C
C-----------------------------------------------------------------------
C  VALEUR DE PRECON   I                  SIGNIFICATION
C-----------------------------------------------------------------------
C        0 OU 1       I  RIEN
C                     I
C        2            I  PRECONDITIONNEMENT DIAGONAL AVEC LA DIAGONALE
C                     I  DE LA MATRICE.
C                     I
C        3            I  PRECONDITIONNEMENT BLOC-DIAGONAL.
C                     I
C        5            I  PRECONDITIONNEMENT DIAGONAL AVEC LA VALEUR
C                     I  ABSOLUE DE LA DIAGONALE DE LA MATRICE.
C                     I
C        7            I  PRECONDITIONNEMENT DE CROUT PAR ELEMENT
C                     I
C       11            I  PRECONDITIONNEMENT DE GAUSS-SEIDEL PAR ELEMENT
C                     I
C       13            I  MATRICE DE PRECONDITIONNEMENT DONNEE
C                     I  PAR LE PROGRAMME APPELANT.
C                     I
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     X          |<-- |  VALEUR INITIALE, PUIS SOLUTION
C |     A          | -->|  MATRICE DU SYSTEME
C |     B          | -->|  SECOND MEMBRE DU SYSTEME.
C |     TB         | -->|  BLOC DE VECTEURS DE TRAVAIL AVEC AU MOINS
C |                |    |  MAX(7,2+2*CFG%KRYLOV)*S VECTEURS, S VALANT 1
C |                |    |  SI A EST UNE MATRICE, 2 SI C'EST UN BLOC DE 4
C |                |    |  ET 3 SI C'EST UN BLOC DE 9.
C |                |    |  CFG%KRYLOV N'EST PRIS EN COMPTE QUE SI
C |                |    |  CFG%SLV = 7 (GMRES)
C |     INFOGR     | -->|  SI OUI, ON IMPRIME UN COMPTE-RENDU
C |     MESH       | -->|  MAILLAGE.
C |     AUX        | -->|  MATRICE DE TRAVAIL DE MEME STRUCTURE QUE A
C |                |    |  (UTILISEE SEULEMENT POUR CERTAINS PRECONDI-
C |                |    |   TIONNEMENTS : 7 , 11 , 13)
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : SOLAUX , ADDBLO , GRACJG
C
C**********************************************************************
C
      USE BIEF, EX_SOLVE => SOLVE
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(SLVCFG), INTENT(INOUT) :: CFG
C
C     STRUCTURES DE VECTEURS OU DE BLOCS DE VECTEURS
C
      TYPE(BIEF_OBJ), TARGET, INTENT(INOUT) :: X,B
      TYPE(BIEF_OBJ), INTENT(INOUT)         :: TB
C
C     STRUCTURES DE MATRICE OU DE BLOC DE MATRICES
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: A,AUX
C
c      INTEGER, INTENT(IN)            :: LT
      LOGICAL, INTENT(IN) :: INFOGR
C
C     STRUCTURE DE MAILLAGE
C
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER PRESTO,IG,LV,S,NBL,I
      INTEGER IT1,IT2,IT3,IT4,IT5,IT6,IT7,IBL1,IBL2,K,IAD,ITB,ITBB
C
c      EXTERNAL TIME_IN_SECONDS2
c      DOUBLE PRECISION TIME_IN_SECONDS2
      DOUBLE PRECISION C
      DOUBLE PRECISION TDEB1, TFIN1
C
      LOGICAL DIADON,PREXSM
      INTEGER NPOIN_TOT,DIFF1,DIFF2,NPOIN_TEMP,NPLAN
      EXTERNAL P_ISUM
      INTEGER P_ISUM
      INTEGER NPOIN2
C-----------------------------------------------------------------------
C
C     STRUCTURES DE BLOCS DE TABLEAUX DE TRAVAIL
C
      TYPE(BIEF_OBJ)          :: TBB
      TYPE(BIEF_OBJ), TARGET  :: BB,BX
      TYPE(BIEF_OBJ), POINTER :: PB,PX
C
C-----------------------------------------------------------------------
C
      LOGICAL FIRST
      DATA FIRST/.TRUE./
C
      SAVE TBB,BB,BX
C-----------------------------------------------------------------------
C
C  ALLOCATION DU BLOC DE BLOCS TBB ET ALLOCATION DE BLOCS DANS TBB
C
c      IF (IPID .EQ. 0) TDEB1=TIME_IN_SECONDS2()
      
       IF(FIRST) THEN
        CALL ALLBLO(TBB,'TBB   ')
        CALL ALLBLO(BB ,'BB    ')
        CALL ALLBLO(BX ,'BX    ')
        FIRST=.FALSE.
      ENDIF
      NBL = 7
      IF(CFG%SLV.EQ.7) NBL = MAX(NBL,4+2*CFG%KRYLOV)
      IF(NBL.GT.TBB%N) THEN
        TBB%N=0
        CALL ALLBLO_IN_BLOCK(TBB,NBL,'BL    ')
      ENDIF
C
C-----------------------------------------------------------------------
C
      LV = MESH%LV
C
C-----------------------------------------------------------------------
C
C  DIFFERENTS TYPES DE SYSTEMES RESOLUS S = 0 : MATRICE NORMALE
C                                       S = 1 : 1 MATRICE  DANS UN BLOC
C                                       S = 2 : 4 MATRICES DANS UN BLOC
C                                       S = 3 : 9 MATRICES DANS UN BLOC
C
       IF(A%TYPE.EQ.3) THEN
        S = 0
        BB%N = 0
        BX%N = 0
        CALL ADDBLO(BB,B)
        CALL ADDBLO(BX,X)
        PX => BX
        PB => BB
      ELSEIF(A%TYPE.EQ.4) THEN
        IF(A%N.EQ.1) THEN
          S = 1
        ELSEIF(A%N.EQ.4) THEN
          S = 2
        ELSEIF(A%N.EQ.9) THEN
          S = 3
        ENDIF
        PX => X
        PB => B
      ENDIF
     
C----------------
C  DIRECT SOLVERS
C  ---> YSMP
C----------------
      IF(CFG%SLV.EQ.8) THEN
C
        IF(NCSIZE.GT.1) THEN
          IF(LNG.EQ.1) WRITE(LU,2018)
          IF(LNG.EQ.2) WRITE(LU,2019)
2018      FORMAT(1X,'UTILISER LE SOLVEUR DIRECT PARALLEL MUMPS,',/,1X,
     *              'SOLVEUR = 9',///)
2019      FORMAT(1X,'USE THE PARALLEL DIRECT SOLVER MUMPS,',/,1X,
     *              'SOLVER = 9',///)
          CALL PLANTE(1)
          STOP
        ENDIF
C
        IF(S.EQ.0) THEN
          IF(A%TYPEXT.NE.'S'.AND.A%TYPEXT.NE.'Q') THEN
            IF(LNG.EQ.1) THEN
              WRITE(LU,*) 'SOLVE (BIEF) : TERMES EXTRA-DIAGONAUX'
              WRITE(LU,*) '               DE TYPE ',A%TYPEXT
              WRITE(LU,*) '               NON TRAITE'
            ENDIF
            IF(LNG.EQ.2) THEN
              WRITE(LU,*) 'SOLVE (BIEF): OFF-DIAGONAL TERMS'
              WRITE(LU,*) '              OF TYPE ',A%TYPEXT
              WRITE(LU,*) '              NOT IMPLEMENTED'
            ENDIF
            CALL PLANTE(1)
            STOP
          ENDIF
          CALL SD_SOLVE_1(A%D%DIM1,MESH%NSEG,MESH%GLOSEG%I,
     *                    MESH%GLOSEG%DIM1,
     *                    A%D%R,A%X%R,X%R,B%R,INFOGR,A%TYPEXT)	
        ELSEIF(S.EQ.1) THEN
          IF(A%ADR(1)%P%TYPEXT.NE.'S'.AND.A%ADR(1)%P%TYPEXT.NE.'Q') THEN
            IF(LNG.EQ.1) THEN
             WRITE(LU,*) 'SOLVE (BIEF) : SOLVEUR DIRECT POUR LES'
             WRITE(LU,*) '               SYSTEMES SYMETRIQUES SEULEMENT'
            ENDIF
            IF(LNG.EQ.2) THEN
             WRITE(LU,*) 'SOLVE (BIEF): DIRECT SOLVER FOR SYMMETRIC'
             WRITE(LU,*) '              SYSTEMS ONLY'
            ENDIF
            CALL PLANTE(1)
            STOP
          ENDIF
          CALL SD_SOLVE_1(A%ADR(1)%P%D%DIM1,MESH%NSEG,MESH%GLOSEG%I,
     *                    MESH%GLOSEG%DIM1,
     *                    A%ADR(1)%P%D%R,A%ADR(1)%P%X%R,X%ADR(1)%P%R,
     *                    B%ADR(1)%P%R,INFOGR,A%ADR(1)%P%TYPEXT)
        ELSEIF(S.EQ.2) THEN
C
          CALL SD_SOLVE_4(MESH%NPOIN,MESH%NSEG,MESH%GLOSEG%I,
     *                    A%ADR(1)%P%D%R,A%ADR(2)%P%D%R,
     *                    A%ADR(3)%P%D%R,A%ADR(4)%P%D%R,
     *                    A%ADR(1)%P%X%R,A%ADR(2)%P%X%R,
     *                    A%ADR(3)%P%X%R,A%ADR(4)%P%X%R,
     *                    X%ADR(1)%P%R,X%ADR(2)%P%R,
     *                    B%ADR(1)%P%R,B%ADR(2)%P%R,INFOGR,
     *                    A%ADR(1)%P%TYPEXT)
!       ELSEIF(S.EQ.3) THEN
        ELSE
          IF(LNG.EQ.1) WRITE(LU,301) S
          IF(LNG.EQ.2) WRITE(LU,401) S
301       FORMAT(1X,'SOLVE (BIEF) : S=',1I6,' CAS NON PREVU')
401       FORMAT(1X,'SOLVE (BIEF): S=',1I6,' CASE NOT IMPLEMENTED')
          CALL PLANTE(1)
          STOP
        ENDIF
        RETURN
      ELSEIF(CFG%SLV.EQ.9) THEN
C----------------
C  DIRECT SOLVERS
C  ---> MUMPS
C----------------
         
         IF(NCSIZE.LT. 1) THEN
           
         
         IF(LNG.EQ.1) WRITE(LU,3018)
         IF(LNG.EQ.2) WRITE(LU,3019)
3018      FORMAT(1X,'MUMPS NON DISPONIBLE POUR DES TESTS SEQUENTIELS',
     *         /,1X,
     *              'UTILISER LE SOLVEUR SEQUENTIEL (SOLVEUR =8)',///)
3019     FORMAT(1X,'MUMPS ARE NOT AVAILABLE FOR SEQUENTIAL RUNS,',/,1X,
     *         'USE SEQUENITAL DIRECT SOLVER (SOLVER = 8) ',///)
          CALL PLANTE(1)
          STOP
       ENDIF 
       OPEN(UNIT=25,FILE='front_glob.dat')
       READ(25,*) NPOIN_TOT
       CLOSE(25)
       IF(S.EQ.0) THEN
            IF(LNG.EQ.1) WRITE(LU,30110) S
            IF(LNG.EQ.2) WRITE(LU,40110) S
30110       FORMAT(1X,'SOLVE (BIEF) : S=',1I6,' CAS NON ENCORE PREVU 
     c           POUR MUMPS')
40110       FORMAT(1X,'SOLVE (BIEF): S=',1I6,' CASE NOT YET 
     c           IMPLEMENTED FOR MUMPS')
            CALL PLANTE(1)
            STOP

         ELSEIF(S.EQ.1) THEN
            
            IF(LNG.EQ.1) WRITE(LU,3011) S
            IF(LNG.EQ.2) WRITE(LU,4011) S
 3011       FORMAT(1X,'SOLVE (BIEF) : S=',1I6,' CAS NON ENCORE PREVU 
     c           POUR MUMPS')
 4011       FORMAT(1X,'SOLVE (BIEF): S=',1I6,' CASE NOT YET 
     c           IMPLEMENTED FOR MUMPS')
            CALL PLANTE(1)
            STOP
         ELSEIF(S.EQ.2) THEN
             CALL PRE4_MUMPS(MESH%NPOIN,MESH%NSEG,MESH%GLOSEG%I,
     *            A%ADR(1)%P%D%R,A%ADR(2)%P%D%R,
     *            A%ADR(3)%P%D%R,A%ADR(4)%P%D%R,
     *            A%ADR(1)%P%X%R,A%ADR(2)%P%X%R,
     *            A%ADR(3)%P%X%R,A%ADR(4)%P%X%R,
     *            X%ADR(1)%P%R,X%ADR(2)%P%R,
     *            B%ADR(1)%P%R,B%ADR(2)%P%R,INFOGR,
     *            A%ADR(1)%P%TYPEXT,MESH%KNOLG%I,NPOIN_TOT,IPID)
          ELSE
             IF(LNG.EQ.1) WRITE(LU,301) S
             IF(LNG.EQ.2) WRITE(LU,401) S
             CALL PLANTE(1)
             STOP
          END IF
             RETURN
      ENDIF
C-----------------------------------------------------------------------
C
      PRESTO = CFG%PRECON
      IF(CFG%PRECON.EQ.0) CFG%PRECON = 1
C
C-----------------------------------------------------------------------
C
C  GESTION DES TABLEAUX DE TRAVAIL : ITB --> PROCHAIN VECTEUR DISPONIBLE
C
C  ITB  --> PROCHAIN VECTEUR DISPONIBLE DANS TB
C  ITBB --> PROCHAIN BLOC    DISPONIBLE DANS TBB
C
      ITB  = 1
      ITBB = 1
C
C  ALLOCATION DE DEUX BLOCS DE TRAVAIL CONTENANT, SOIT UN VECTEUR,
C  SOIT UN BLOC DE VECTEURS (CAS OU S EST <>0)
C  CES DEUX BLOCS SONT COMMUNS A TOUTES LES METHODES.
C
C     POUR LES MATRICES DE PRECONDITIONNEMENT
      IF(3*(CFG%PRECON/3).EQ.CFG%PRECON) THEN
C       PRECONDITIONNEMENT BLOC-DIAGONAL : S**2 DIAGONALES
        CALL SOLAUX(IT1, TB,TBB,ITB,ITBB,S**2)
      ELSE
C       AUTRES : S DIAGONALES
        CALL SOLAUX(IT1, TB,TBB,ITB,ITBB,S)
      ENDIF
C
      CALL SOLAUX(IT2, TB,TBB,ITB,ITBB,S)
C
      IF(CFG%SLV.EQ.7) THEN
C       SPECIAL  GMRES : TABLEAUX DEPENDANT DE LA DIMENSION DE KRYLOV
C                       TBB(IBL1) : BLOC DE CFG%KRYLOV VECTEURS
C                       OU BLOC DE CFG%KRYLOV BLOCS DE S VECTEURS
C                       TBB(IBL2) : IDEM
C
        IBL1=ITBB
        ITBB = ITBB + 1
        IBL2=ITBB
        ITBB = ITBB + 1
        TBB%ADR(IBL1)%P%N=0
        TBB%ADR(IBL2)%P%N=0
        DO 10 K=1,CFG%KRYLOV
          CALL SOLAUX(IAD, TB,TBB,ITB,ITBB,S)
          CALL ADDBLO(TBB%ADR(IBL1)%P,TBB%ADR(IAD)%P)
          CALL SOLAUX(IAD, TB,TBB,ITB,ITBB,S)
          CALL ADDBLO(TBB%ADR(IBL2)%P,TBB%ADR(IAD)%P)
10      CONTINUE
C       POUR EVITER UN WARNING DU COMPILATEUR INTEL
        IT3=-1
        IT4=-1
        IT5=-1
        IT6=-1
        IT7=-1
      ELSE
C       AUTRES METHODES (ON POURRAIT PARFOIS NE PAS ALLOUER IT6 OU IT7)
        CALL SOLAUX(IT3, TB,TBB,ITB,ITBB,S)
        CALL SOLAUX(IT4, TB,TBB,ITB,ITBB,S)
        CALL SOLAUX(IT5, TB,TBB,ITB,ITBB,S)
        CALL SOLAUX(IT6, TB,TBB,ITB,ITBB,S)
        CALL SOLAUX(IT7, TB,TBB,ITB,ITBB,S)
C       POUR EVITER UN WARNING DU COMPILATEUR CRAY
        IBL1=1
        IBL2=1
C
      ENDIF
C
C     PRECONDITIONNEMENT DE CROUT : IL FAUT AU PREALABLE UN PRECONDITION
C     NEMENT QUI METTE DES DIAGONALES A 1.
C     DE PLUS LE GRADIENT SERA DISTINGUE DU RESIDU (POINTEUR IG)
C
      IF(  7*(CFG%PRECON/ 7).EQ.CFG%PRECON.OR.
     *    11*(CFG%PRECON/11).EQ.CFG%PRECON.OR.
     *    13*(CFG%PRECON/13).EQ.CFG%PRECON.OR.
     *    17*(CFG%PRECON/17).EQ.CFG%PRECON     ) THEN
        IG=IT6
        IF(2*(CFG%PRECON/2).NE.CFG%PRECON.AND.
     *     3*(CFG%PRECON/3).NE.CFG%PRECON) THEN
C         ON CHOISIT DIAGONAL
          CFG%PRECON=2*CFG%PRECON
        ENDIF
      ELSE
C       NOTE IT5 =-1 SI CFG%SLV.EQ.7 MAIS DANS CE CAS IG EST INUTILE
        IG=IT5
      ENDIF
C
C  FIN DE LA GESTION DES TABLEAUX DE TRAVAIL
C
C-----------------------------------------------------------------------
C                                   -1/2      -1/2  1/2        -1/2
C  PRECONDITIONNEMENTS DIAGONAUX : D     A  D      D     X  = D      B
C
      DIADON = .FALSE.
      PREXSM = .TRUE.
C
      IF(3*(CFG%PRECON/3).EQ.CFG%PRECON.AND.(S.EQ.2.OR.S.EQ.3)) THEN
C       PRECONDITIONNEMENT DIAGONAL-BLOC (4 OU 9 MATRICES)
        CALL PREBDT(X,A,B,TBB%ADR(IT1)%P,MESH,PREXSM,DIADON,S)
C       NE PAS MOFIFIER D11,D22,D33 A L'APPEL DE PRECDT
        DIADON = .TRUE.
      ENDIF
C
      IF(2*(CFG%PRECON/2).EQ.CFG%PRECON.OR.
     *   3*(CFG%PRECON/3).EQ.CFG%PRECON.OR.
     *   5*(CFG%PRECON/5).EQ.CFG%PRECON) THEN
        CALL PRECDT(X,A,B,TBB%ADR(IT1)%P,MESH,
     *              CFG%PRECON,PREXSM,DIADON,S)
      ENDIF
C
C-----------------------------------------------------------------------
C
C  CONSTRUCTION DES MATRICES DE PRECONDITIONNEMENT :
C
      IF(7*(CFG%PRECON/7).EQ.CFG%PRECON) THEN
        CALL DCPLDU(AUX,A,MESH,.TRUE.,LV)
      ELSEIF(11*(CFG%PRECON/11).EQ.CFG%PRECON) THEN
        CALL GSEBE(AUX,A,MESH)
      ELSEIF(13*(CFG%PRECON/13).EQ.CFG%PRECON) THEN
C       ON NE FAIT RIEN, AUX EST DONNEE PAR LE PROGRAMME APPELANT
      ELSEIF(17*(CFG%PRECON/17).EQ.CFG%PRECON) THEN
        IF(CFG%SLV.NE.1.AND.CFG%SLV.NE.2) THEN
          WRITE(LU,*) 'PRECONDITIONING 17'
          WRITE(LU,*) 'NOT IMPLEMENTED FOR SOLVER ',CFG%SLV
          CALL PLANTE(1)
          STOP
        ENDIF
        IF(MESH%TYPELM.NE.40) THEN
          WRITE(LU,*) 'PRECONDITIONING 17'
          WRITE(LU,*) 'IMPLEMENTED ONLY FOR PRISMS'
          CALL PLANTE(1)
          STOP
        ENDIF
        IF(AUX%TYPE.NE.3) THEN
          WRITE(LU,*) 'PRECONDITIONING 17'
          WRITE(LU,*) 'NOT IMPLEMENTED FOR BLOCKS OF MATRICES'
          CALL PLANTE(1)
          STOP
        ENDIF
C
        IF(AUX%STO.EQ.1) THEN      
        CALL PREVEREBE(AUX%X%R,A%D%R,A%X%R,A%TYPDIA,A%TYPEXT,
     *              MESH%IKLE%I,MESH%NPOIN,MESH%NELEM,MESH%NELMAX,MESH)
        ELSE
        CALL PREVERSEG(AUX%X%R,A%D%R,A%X%R,A%TYPDIA,A%TYPEXT,
     *                 MESH%NPOIN,MESH,MESH%NSEG)
        ENDIF
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  PARALLELISME : SECOND MEMBRE
C 
      IF(NCSIZE.GT.1) THEN
        CALL PARCOM(B,2,MESH)
      ENDIF
C
C-----------------------------------------------------------------------
C
C  RESOLUTION DU SYSTEME LINEAIRE :
C
      IF(CFG%SLV.EQ.1) THEN
C
C       GRADIENT CONJUGUE
C
        CALL GRACJG(PX, A,PB, MESH,
     *              TBB%ADR(IT2)%P,TBB%ADR(IT3)%P,
     *              TBB%ADR(IT5)%P,TBB%ADR(IG)%P,
     *              CFG,INFOGR,AUX)
C     
      ELSEIF(CFG%SLV.EQ.2) THEN
C
C       RESIDU CONJUGUE
C
        CALL RESCJG(PX, A,PB, MESH,
     *              TBB%ADR(IT2)%P,TBB%ADR(IT3)%P,
     *              TBB%ADR(IT4)%P,TBB%ADR(IT5)%P,
     *              TBB%ADR(IG)%P,
     *              CFG,INFOGR,AUX)
C
      ELSEIF(CFG%SLV.EQ.3) THEN
C
C       EQUATION NORMALE
C
        CALL EQUNOR(PX, A,PB, MESH,
     *              TBB%ADR(IT2)%P,TBB%ADR(IT3)%P,
     *              TBB%ADR(IT4)%P,TBB%ADR(IT5)%P,
     *              TBB%ADR(IG)%P,
     *              CFG,INFOGR,AUX)
C
      ELSEIF(CFG%SLV.EQ.4) THEN
C
C       ERREUR MINIMALE
C
        CALL ERRMIN(PX, A,PB, MESH,
     *              TBB%ADR(IT2)%P,TBB%ADR(IT3)%P,TBB%ADR(IT5)%P,
     *              TBB%ADR(IG)%P,
     *              CFG,INFOGR,AUX)
C
      ELSEIF(CFG%SLV.EQ.5) THEN
C
C       GRADIENT CONJUGUE CARRE
C
        CALL CGSQUA(PX, A,PB, MESH,
     *              TBB%ADR(IT2)%P,TBB%ADR(IT3)%P,TBB%ADR(IT4)%P,
     *              TBB%ADR(IT5)%P,TBB%ADR(IT6)%P,TBB%ADR(IT7)%P,
     *              CFG,INFOGR)
C
      ELSEIF(CFG%SLV.EQ.6) THEN
C
C       GRADIENT CONJUGUE CARRE STABILISE
C
        CALL CGSTAB(PX, A,PB, MESH,
     *              TBB%ADR(IT2)%P,TBB%ADR(IT3)%P,TBB%ADR(IT4)%P,
     *              TBB%ADR(IT5)%P,TBB%ADR(IT6)%P,TBB%ADR(IT7)%P,
     *              CFG,INFOGR,AUX)
C
      ELSEIF(CFG%SLV.EQ.7) THEN
C
C       GENERALIZED MINIMUM RESIDUAL
C
        CALL GMRES(PX, A,PB,MESH,
     *             TBB%ADR(IT2)%P,TBB%ADR(IBL1)%P,TBB%ADR(IBL2)%P,
     *             CFG,INFOGR,AUX)
C
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,300) CFG%SLV
        IF (LNG.EQ.2) WRITE(LU,400) CFG%SLV
300     FORMAT(1X,'SOLVE (BIEF) :',1I6,' METHODE NON PREVUE :')
400     FORMAT(1X,'SOLVE (BIEF) :',1I6,' METHOD NOT AVAILABLE :')
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  INVERSION DU CHANGEMENT DE VARIABLE EN CAS DE PRECONDITIONNEMENT
C                                                DIAGONAL
C                                                DIAGONAL-BLOC
C
      IF(2*(CFG%PRECON/2).EQ.CFG%PRECON.OR.
     *   3*(CFG%PRECON/3).EQ.CFG%PRECON.OR.
     *   5*(CFG%PRECON/5).EQ.CFG%PRECON    ) THEN
        CALL OS( 'X=XY    ' , PX , TBB%ADR(IT1)%P , PX , C )
      ENDIF
C
      IF(3*(CFG%PRECON/3).EQ.CFG%PRECON.AND.(S.EQ.2.OR.S.EQ.3)) THEN
        CALL UM1X(X,TBB%ADR(IT1)%P,S)
      ENDIF
C
C-----------------------------------------------------------------------
C
      CFG%PRECON = PRESTO
C
C-----------------------------------------------------------------------
C

   
      RETURN
      END
      
