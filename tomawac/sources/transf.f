C                       *****************
                        SUBROUTINE TRANSF
C                       *****************
C
     *( FA    , FR    , FREQ  , DFREQ , COSTET, SINTET, UC    , VC    ,
     *  XK    , KNEW  , NEWF  , NEWF1 , TAUX1 , TAUX2 , NPOIN2, NPLAN ,
     *  NF    , RAISF , LT    , GRADEB, GRAPRD)
C
C***********************************************************************
C TOMAWAC     V5P6         12/01//2006    M. BENOIT (LNHE)
C***********************************************************************
C
C      FONCTION:
C      =========
C      CONEVRSION D'UN SPECTRE DONNE EN FREQUENCE RELATIVE FR(-,-,-) EN
C      SPECTRE EN FREQUENCE ABSOLUE FA(-,-,-).
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    FA(-,-,-)   !<-- ! SPECTRE EN FREQUENCE ABSOLUE                 !
C !    FR(-,-,-)   ! -->! SPECTRE EN FREQUENCE RELATIVE                !
C !    FREQ(-)     ! -->! FREQUENCES DISCRETISEES                      !
C !    DFREQ(-)    ! -->! ECARTS DE FREQUENCES DISCRETISES             !
C !    COSTET(-)   ! -->! COSINUS DES DIRECTIONS DE PROPAGATION        !
C !    SINTET(-)   ! -->! SINUS DES DIRECTIONS DE PROPAGATION          !
C !    UC(-)       ! -->! COMPOSANTE OUEST-EST DU CHAMP DE COURANT     !
C !    VC(-)       ! -->! COMPOSANTE SUD-NORD  DU CHAMP DE COURANT     !
C !    XK(-,-)     ! -->! NOMBRE D'ONDE                                !
C !    KNEW(-)     !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)        !
C !    NEWF(-)     !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)        !
C !    NEWF1(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)        !
C !    TAUX1(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)        !
C !    TAUX2(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)        !
C !    NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE 2D              !
C !    NPLAN       ! -->! NOMBRE DE DIRECTIONS DE PROPAGATION          !
C !    NF          ! -->! NOMBRE DE FREQUENCES                         !
C !    RAISF       ! -->! RAISON FREQUENTIELLE                         !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE PAR : WAC
C SOUS-PROGRAMMES APPELES : NEANT
C
C***********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER          NPOIN2, NPLAN , NF    , LT    , GRADEB, GRAPRD
      INTEGER          KNEW(NPOIN2)  , NEWF(NPOIN2)  , NEWF1(NPOIN2)
      DOUBLE PRECISION RAISF
      DOUBLE PRECISION FA(NPOIN2,NPLAN,NF),FR(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION FREQ(NF),DFREQ(NF),COSTET(NPLAN),SINTET(NPLAN)
      DOUBLE PRECISION UC(NPOIN2),VC(NPOIN2),TAUX1(NPOIN2),TAUX2(NPOIN2)
      DOUBLE PRECISION XK(NPOIN2,NF)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER          IP    , JP    , JF    , NEWM  , NEWM1 , KH
      DOUBLE PRECISION F0    , UK    , DEUPI , AUXI  , Y     , Z
      DOUBLE PRECISION FNEW  , UNSLRF
      LOGICAL          IMP
C
C
C-----------------------------------------------------------------------
C     ON NE FAIT LA TRANSFORMATION QUE POUR LES DATES DE SORTIE
C-----------------------------------------------------------------------
      IMP=.FALSE.
      IF ((LT.GE.GRADEB).AND.(MOD(LT-GRADEB,GRAPRD).EQ.0)) IMP=.TRUE.
      IF (.NOT.(IMP)) RETURN
C-----------------------------------------------------------------------
C
C
      DEUPI=2.D0*3.141592654D0
      F0=FREQ(1)
      UNSLRF=1.0D0/DLOG(RAISF)
C
      CALL OV( 'X=C     ' , FA , Y , Z , 0.D0 , NPOIN2*NPLAN*NF)
C
      DO JF=1,NF
C
        DO JP=1,NPLAN
C
          DO IP=1,NPOIN2
C
C           ---------------------------------------------------------
C           CALCUL DE DIFFERENCE ENTRE FREQENCES ABSOLUE ET RELATIVE
C                                            -> ->
C                 Z = Freq_abs - Freq_rel = (k .U)/(2.pi)
C           ON NE REPROJETTE LE SPECTRE SUR LES FREQUENCES ABSOLUES
C           QUE SI VARIATION RELATIVE Z/Freq_rel EST SIGNIFICATIVE.
C           ---------------------------------------------------------
            UK=SINTET(JP)*UC(IP)+COSTET(JP)*VC(IP)
            Z=UK*XK(IP,JF)/DEUPI
C
            IF (DABS(Z)/FREQ(JF).LT.1.0D-3) THEN
              KNEW (IP)=JP
              NEWF (IP)=JF
              NEWF1(IP)=-1
              TAUX1(IP)=FR(IP,JP,JF)
              TAUX2(IP)=0.0D0
            ELSE
C
C             -------------------------------------------------------
C             CALCUL DE FNEW ET KNEW
C             -------------------------------------------------------
              FNEW = FREQ(JF)+Z
              IF (FNEW.GT.0.D0) THEN 
                KNEW(IP)=JP
              ELSE
                KNEW(IP)=1+MOD(JP+NPLAN/2-1,NPLAN)
                FNEW=-FNEW
              ENDIF
C
C             -------------------------------------------------------
C             CALCUL DE NEWF INDICE FREQUENTIEL DE LA FREQUENCE DE
C             DISCRETISATION IMMEDIATEMENT INFERIEURE A FNEW
C             -------------------------------------------------------
              IF (FNEW.LT.F0/RAISF) THEN
                NEWF(IP)=-1
              ELSE
                NEWF(IP)=INT(1.0D0+DLOG(FNEW/F0)*UNSLRF)
              ENDIF
C
C             -------------------------------------------------------
C             CALCUL DES COEFFICIENTS ET INDICES POUR LA PROJECTION
C             -------------------------------------------------------
              IF ((NEWF(IP).LT.NF).AND.(NEWF(IP).GE.1)) THEN  
                NEWF1(IP)=NEWF(IP)+1
                AUXI=FR(IP,JP,JF)*DFREQ(JF)
     *               /(FREQ(NEWF1(IP))-FREQ(NEWF(IP)))
                TAUX1(IP)=AUXI*(FREQ(NEWF1(IP))-FNEW)/DFREQ(NEWF(IP))
                TAUX2(IP)=AUXI*(FNEW-FREQ(NEWF(IP)))/DFREQ(NEWF1(IP))
              ELSEIF (NEWF(IP).EQ.0) THEN
                AUXI=FR(IP,JP,JF)*DFREQ(JF)/(F0*(1.D0-1.D0/RAISF))
                TAUX2(IP)=AUXI*(FNEW-F0/RAISF)/DFREQ(1)
                NEWF (IP)=-1
                NEWF1(IP)= 1
              ELSEIF (NEWF(IP).EQ.NF) THEN
                AUXI=FR(IP,JP,JF)*DFREQ(JF)/(FREQ(NF)*(RAISF-1.D0))
                TAUX1(IP)=AUXI*(FREQ(NF)*RAISF-FNEW)/DFREQ(NF)
                NEWF1(IP)=-1
              ELSE
                NEWF (IP)=-1
                NEWF1(IP)=-1
              ENDIF
C
            ENDIF
C
          ENDDO
C
C
C         -------------------------------------------------------
C         PROJECTION DU SPECTRE PROPREMENT DITE
C         -------------------------------------------------------
          DO IP=1,NPOIN2
            NEWM =NEWF (IP)
            NEWM1=NEWF1(IP)
            KH=KNEW(IP)
            IF (NEWM .NE.-1) FA(IP,KH,NEWM )=FA(IP,KH,NEWM )+TAUX1(IP)
            IF (NEWM1.NE.-1) FA(IP,KH,NEWM1)=FA(IP,KH,NEWM1)+TAUX2(IP)
          ENDDO
C
        ENDDO
C
      ENDDO
C
      RETURN
      END
