C                       *****************
                        SUBROUTINE LIMWAC
C                       *****************
C
     *(F     , FBOR  , LIFBOR, NPTFR , NPLAN , NF    ,  TETA , FREQ  ,
     * NPOIN2, NBOR  , AT    , LT    , DDC   , LIMSPE, FPMAXL, FETCHL,
     * SIGMAL, SIGMBL, GAMMAL, FPICL , HM0L  , APHILL, TETA1L, SPRE1L,
     * TETA2L, SPRE2L, XLAMDL, X ,Y  , KENT  , KSORT , NFO1  , NBI1  ,
     * BINBI1, UV    , VV    , SPEULI, VENT  , VENSTA, GRAVIT, DEUPI , 
     * PRIVE , NPRIV , SPEC  , FRA   , DEPTH , FRABL )
C
C***********************************************************************
C TOMAWAC   V1.0            01/02/95        F. MARCOS  (LNH) 30 87 72 66
C***********************************************************************
C
C      FONCTION:
C      =========
C
C    CONDITIONS AUX LIMITES
C
C    ATTENTION
C    PAR DEFAUT, ON DUPLIQUE SUR L'ENSEMBLE DES DIRECTIONS ET DES
C    FREQUENCES LA CONDITION A LA LIMITE DONNEE DANS LE FICHIER DYNAM
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    F           ! -->!  DENSITE SPECTRALE                           !
C !    FBOR        !<-->!  DENSITE SPECTRALE AU BORD                   !
C !    LIFBOR      ! -->!  TYPE DE CONDITION LIMITE SUR F              !
C !    NPTFR       ! -->!  NOMBRE DE POINTS FRONTIERE 2D               !
C !    NPLAN       ! -->!  NOMBRE DE DIRECTIONS                        !
C !    NF          ! -->!  NOMBRE DE FREQUENCES                        !
C !    TETA        ! -->! DIRECTIONS DE PROPAGATION                    !
C !    FREQ        ! -->! FREQUENCES DISCRETISEES                      !
C !    NPOIN2      ! -->!  NOMBRE DE POINTS 2D                         !
C !    NBOR        ! -->!  NUMEROTATION DES POINTS DE BORD 2D          !
C !    AT          ! -->!  TEMPS                                       !
C !    LT          ! -->!  NUMERO DU PAS DE TEMPS                      !
C !    DDC         ! -->!  DATE DU DEBUT DU CALCUL                     !
C !    X           ! -->!  ABSCISSES DES POINTS 2D                     !
C !    Y           ! -->!  ORDONNEES DES POINTS 2D                     !
C !    KENT        ! -->!  C.L. INDIQUANT UNE FRONTIERE MARITIME       !
C !    KSORT       ! -->!  C.L. INDIQUANT UNE FRONTIERE SOLIDE         !
C !    NFO1        ! -->!  NUMERO DU FICHIER FORMATE UTILISATEUR       !
C !    NBI1        ! -->!  NUMERO DU FICHIER BINAIRE UTILISATEUR       !
C !    BINBI1      ! -->!  BINAIRE DU FICHIER BINAIRE UTILISATEUR      !
C !    PRIVE       ! -->!  TABLEAU DE L'UTILISATEUR                    !
C !    NPRIV       ! -->!  DIMENSION DU TABLEAU PRIVE                  !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE PAR : WAC
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
      INTEGER NPLAN,NF,NPOIN2,NPTFR,LT,NPRIV
C
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF),X(NPOIN2),Y(NPOIN2)
      DOUBLE PRECISION FBOR(NPTFR,NPLAN,NF),TETA(NPLAN),FREQ(NF)
      DOUBLE PRECISION UV(NPOIN2),VV(NPOIN2), SPEC(NF), FRA(NPLAN)
      DOUBLE PRECISION PRIVE(NPOIN2,NPRIV),DDC, DEPTH(NPOIN2)
      DOUBLE PRECISION HM0L,FPICL,GAMMAL,SIGMAL,SIGMBL,APHILL,FETCHL
      DOUBLE PRECISION FPMAXL,TETA1L,SPRE1L,TETA2L,SPRE2L,XLAMDL
      DOUBLE PRECISION GRAVIT,DEUPI,E2FMIN
C
      DOUBLE PRECISION AT
C
      LOGICAL SPEULI, VENT, VENSTA
C
      INTEGER NBOR(NPTFR),LIFBOR(NPTFR),NFO1,NBI1,NPB
      INTEGER KENT,KSORT,IFF,IPLAN,IPTFR,LIMSPE,FRABL
C
      DOUBLE PRECISION, ALLOCATABLE :: TRAV(:)
      DOUBLE PRECISION, ALLOCATABLE :: UV2D(:),VV2D(:),PROF(:)
      DOUBLE PRECISION, ALLOCATABLE :: FB_CTE(:,:)
      LOGICAL FLAG
C
      CHARACTER*3 BINBI1
C
      SAVE NPB,UV2D,VV2D,PROF,FB_CTE
C
C***********************************************************************
C
C   MODIFICATION EVENTUELLE DU TYPE DE CONDITION A LA LIMITE
C
C   A REMPLIR PAR L'UTILISATEUR
C
C   LIFBOR(IPTFR)=KENT OU KSORT
C
      IF (LIMSPE.EQ.0 .AND. .NOT.SPEULI) RETURN
C
      FLAG=.FALSE.
      IF (VENT .AND. (LIMSPE.EQ.1 .OR. LIMSPE.EQ.2 .OR. LIMSPE.EQ.3
     * .OR. LIMSPE.EQ.5)) FLAG=.TRUE.
C
C     AU PREMIER PASSAGE, ON ALLOUE DE LA MEMOIRE AUX TABLEAUX UTILES
C     ---------------------------------------------------------------
      IF (LT.LT.1) THEN
        NPB=1
        IF (FLAG) THEN
           ALLOCATE(UV2D(1:NPTFR),VV2D(1:NPTFR))
           NPB=NPTFR
        ENDIF
        IF (LIMSPE.EQ.7 .OR. SPEULI) THEN
           ALLOCATE(PROF(1:NPTFR))
           NPB=NPTFR
        ENDIF
        IF (NPB.EQ.1) THEN
           ALLOCATE(FB_CTE(1:NPLAN,1:NF))
        ENDIF
      ENDIF
      IF (.NOT.allocated(UV2D)) ALLOCATE(UV2D(NPTFR))
      IF (.NOT.allocated(VV2D)) ALLOCATE(VV2D(NPTFR)) 
      IF (.NOT.allocated(PROF)) ALLOCATE(PROF(NPTFR))     
      IF (.NOT.allocated(FB_CTE)) ALLOCATE(FB_CTE(1:NPLAN,1:NF))
C
C     AU PREMIER PASSAGE (ET EVENTUELLEMENT AUX AUTRES SI LE VENT EST 
C     INSTATIONNAIRE ET QUE LE SPECTRE A LA LIMITE EN DEPEND),
C     ON CALCULE LE SPECTRE AUX LIMITES
C     ----------------------------------------------------------------
      IF (LT.LT.1 .OR. (.NOT.VENSTA.AND.FLAG) .OR. SPEULI) THEN
        IF (FLAG) THEN
          DO IPTFR=1,NPTFR
            UV2D(IPTFR)=UV(NBOR(IPTFR))
            VV2D(IPTFR)=VV(NBOR(IPTFR))
          ENDDO
        ENDIF
        IF(LIMSPE.EQ.7 .OR. SPEULI) THEN
          DO IPTFR=1,NPTFR
            PROF(IPTFR)=DEPTH(NBOR(IPTFR))
          ENDDO
        ENDIF
C
C       APPEL A SPEINI
C     ----------------------------------------------------------------
        E2FMIN = 1.D-30
C
        IF (NPB.EQ.NPTFR) THEN
          CALL SPEINI
     *( FBOR  , SPEC  , FRA    , UV2D  , VV2D  , FREQ ,
     *  TETA  , GRAVIT, FPMAXL , FETCHL, SIGMAL, SIGMBL, GAMMAL, FPICL,
     *  HM0L  , APHILL, TETA1L , SPRE1L, TETA2L, SPRE2L, XLAMDL,
     *  NPTFR , NPLAN , NF     , LIMSPE, E2FMIN, PROF  , FRABL )
        ELSE
          CALL SPEINI
     *( FB_CTE, SPEC  , FRA    , UV2D  , VV2D  , FREQ ,
     *  TETA  , GRAVIT, FPMAXL , FETCHL, SIGMAL, SIGMBL, GAMMAL, FPICL,
     *  HM0L  , APHILL, TETA1L , SPRE1L, TETA2L, SPRE2L, XLAMDL,
     *  NPB   , NPLAN , NF     , LIMSPE, E2FMIN, PROF  , FRABL )
	ENDIF
C
C     ===========================================================
C     ZONE UTILISATEUR - ON PEUT Y MODIFIER RESU
C     ===========================================================
        IF (SPEULI) THEN
C
C        EXEMPLE DE MODIFICATION DE FRA - A MODIFIER SUIVANT VOTRE CAS
C        EXAMPLE OF MODIFICATION OF FRA - TO BE MODIFIED DEPENDING 
C        ON YOUR CASE
C        ALLOCATE(TRAV(1:NF))
C
C        DO IFREQ=1,NF
C             IF (FREQ(IFF).LT.FPIC) THEN
C              TRAV(IFF)=0.4538D0*(FREQ(IFF)/FPIC)**(-2.03D0)
C           ELSE
C              TRAV(IFF)=0.4538D0*(FREQ(IFF)/FPIC)**(1.04D0)
C           ENDIF
C        ENDDO
C
C        DO IPLAN=1,NPLAN
C             DTETA=TETA(IPLAN)-TETA1
C           IF ((TETA(IPLAN)-TETA1).GT.DEUPI/2) THEN
C              DTETA=DEUPI-DTETA
C           ENDIF
C           DO IFF=1,NF
C              FRA(IPLAN)=1.D0/SQRT(DEUPI)*TRAV(IFF)*
C     *                       EXP(-DTETA**2/(2.D0*TRAV(IFF)**2))
C              DO IPTFR=1,NPTFR
C                FBOR(IPTFR,IPLAN,IFF)= SPEC(IFF)*FRA(IPLAN)
C              ENDDO
C           ENDDO
C        ENDDO
C        DEALLOCATE(TRAV)
C
C        PARTIE A SUPPRIMER SI ON FAIT DES MODIFICATIONS
C        LINES TO ERASE IF YOU DO MODIFICATIONS 
C
        IF (LNG.EQ.1) THEN
          WRITE(LU,*)'*****  ERREUR LIMWAC  ******'
          WRITE(LU,*)
     *      ' VOUS NE MODIFIEZ PAS LE SPECTRE AUX LIMITES ALORS QUE'
          WRITE(LU,*)' VOUS EN DEMANDEZ LA POSSIBILITE'
        ELSE
          WRITE(LU,*)'*****  ERROR LIMWAC  ******'
          WRITE(LU,*)
     *      ' YOU DID NOT MODIFY THE BOUNDARY SPECTRUM WHEREAS '
          WRITE(LU,*)' YOU ASK FOR THAT '
        ENDIF
        STOP
      ENDIF
C
C     ===========================================================
C     FIN DE LA ZONE UTILISATEUR 
C     ===========================================================
      ENDIF
C
C
C     -----------------------------------------------------------------
C     DUPLICATION SUR TOUTES LES DIRECTIONS ET TOUTES LES FREQUENCES
C     DE LA C.L. DE DYNAM SI ON EST EN CONDITION DE FRONT. LIQUIDE
C     -----------------------------------------------------------------
      IF (FLAG .OR. LIMSPE.EQ.7 .OR. SPEULI) THEN
        DO IPTFR=1,NPTFR
          IF (LIFBOR(IPTFR).EQ.KENT) THEN
            DO IFF=1,NF
              DO IPLAN=1,NPLAN
                F(NBOR(IPTFR),IPLAN,IFF)=FBOR(IPTFR,IPLAN,IFF)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ELSE
        DO IPTFR=1,NPTFR
          IF (LIFBOR(IPTFR).EQ.KENT) THEN
            DO IFF=1,NF
              DO IPLAN=1,NPLAN
                F(NBOR(IPTFR),IPLAN,IFF)=FB_CTE(IPLAN,IFF)
              ENDDO
            ENDDO     
          ENDIF
        ENDDO
      ENDIF
C-----------------------------------------------------------------------
C
      RETURN
      END
