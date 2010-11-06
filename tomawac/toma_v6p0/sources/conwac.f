C                       *****************
                        SUBROUTINE CONWAC
C                       *****************
C
     *( CX    , CY    , CT    , XK    , CG    , COSF  , TGF   , DEPTH ,
     *  DZX   , DZY   , FREQ  , COSTET, SINTET, NPOIN2, NPLAN , JF    ,
     *  NF    , PROINF, SPHE  , PROMIN, TRA01 , TRA02 )
C
C***********************************************************************
C TOMAWAC  V5.4      19/01/2004     M. BENOIT (EDF LNHE) 01 30 87 83 51
C***********************************************************************
C
C      FONCTION:
C      =========
C
C    CALCULE LE CHAMP CONVECTEUR (3D SANS COURANT)
C                 -->                            -->
C    ATTENTION ICI X EST VERTICAL VERS LE HAUT ET Y EST HORIZONTAL
C      VERS LA DROITE ET TETA EST LA DIRECTION / NORD COMPTE DANS LE
C      SENS INVERSE DU SENS TRIGO
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    CX,CY,CT    !<-- !  CHAMP CONVECTEUR SELON X(OU PHI),           !
C !                !    !  Y(OU LAMBDA) ET TETA                        !
C !    XK          ! -->!  NOMBRE D'ONDE DISCRETISE                    !
C !    CG          ! -->!  VITESSE DE GROUPE DISCRETISEE               !
C !    COSF        ! -->!  COSINUS DES LATITUDES DES POINTS 2D         !
C !    TGF         ! -->!  TANGENTES DES LATITUDES DES POINTS 2D       !
C !    DEPTH       ! -->!  PROFONDEUR                                  !
C !    DZX         ! -->!  GRADIENT DE FOND SELON X                    !
C !    DZY         ! -->!  GRADIENT DE FOND SELON Y                    !
C !    FREQ        ! -->!  FREQUENCES DISCRETISEES                     !
C !    COSTET      ! -->!  COSINUS TETA                                !
C !    SINTET      ! -->!  SINUS TETA                                  !
C !    NPOIN2      ! -->!  NOMBRE DE POINTS DU MAILLAGE 2D             !
C !    NPLAN       ! -->!  NOMBRE DE PLANS OU DE DIRECTIONS            !
C !    JF          ! -->!  FREQUENCES COURANTE                         !
C !    NF          ! -->!  NOMBRE DE FREQUENCES                        !
C !    PROINF      ! -->!  LOGIQUE INDIQUANT SI ON EST EN PROF INFINIE !
C !    SPHE        ! -->!  LOGIQUE INDIQUANT SI ON EST EN COORD. SPHER.!
C !    PROMIN      ! -->!  VALEUR MINIMALE DE LA PROFONDEUR D'EAU      !
C !    TRA01       !<-->!  TABLEAU DE TRAVAIL                          !
C !    TRA02       !<-->!  TABLEAU DE TRAVAIL                          !
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
C.....VARIABLES TRANSMISES.
C     """""""""""""""""""""
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
      INTEGER          NF    , NPLAN , NPOIN2, JF
      DOUBLE PRECISION PROMIN
      DOUBLE PRECISION DEPTH(NPOIN2) , DZX(NPOIN2)   , DZY(NPOIN2)
      DOUBLE PRECISION COSF(NPOIN2)  , TGF(NPOIN2)   , FREQ(NF)
      DOUBLE PRECISION COSTET(NPLAN) , SINTET(NPLAN)
      DOUBLE PRECISION TRA01(NPLAN)  , TRA02(NPLAN)
      DOUBLE PRECISION CG(NPOIN2,NF) , XK(NPOIN2,NF)
      DOUBLE PRECISION CX(NPOIN2,NPLAN),CY(NPOIN2,NPLAN)
      DOUBLE PRECISION CT(NPOIN2,NPLAN)
      LOGICAL          PROINF, SPHE
C
C.....VARIABLES LOCALES.
C     """"""""""""""""""
      INTEGER          JP    , IP
      DOUBLE PRECISION GSQP  , SR    , R     , SRCF  , TFSR
      DOUBLE PRECISION DDDN  , DSDNSK, GRADEG, DEUKD , DEUPI
C
C
      GSQP=0.780654996D0
      R=6400.D3
      DEUPI=6.283185307D0
C
      IF (PROINF) THEN
C-----------------------------------------------------------------------
C     EN PROFONDEUR INFINIE ...
C-----------------------------------------------------------------------
C
        DO JP=1,NPLAN
          TRA01(JP)=GSQP/FREQ(JF)*COSTET(JP)
          TRA02(JP)=GSQP/FREQ(JF)*SINTET(JP)
        ENDDO
C
        IF (.NOT.SPHE) THEN
C       ----------------------------------------------------------------
C       ... ET COORDONNEES CARTESIENNES
C       ----------------------------------------------------------------
          DO IP=1,NPOIN2
            DO JP=1,NPLAN
              CX(IP,JP)=TRA01(JP)
              CY(IP,JP)=TRA02(JP)
              CT(IP,JP)=0.D0
            ENDDO
          ENDDO
C
        ELSE
C       ----------------------------------------------------------------
C       ... ET COORDONNEES SPHERIQUES
C       ----------------------------------------------------------------
          SR=1.D0/R
          GRADEG=180.D0/3.1415926D0
          DO IP=1,NPOIN2
            SRCF=SR/COSF(IP)
            TFSR=TGF(IP)*SR
            DO JP=1,NPLAN
              CX(IP,JP)=TRA01(JP)*SR*GRADEG
              CY(IP,JP)=TRA02(JP)*SRCF*GRADEG
              CT(IP,JP)=TRA02(JP)*TFSR
            ENDDO
          ENDDO
C
        ENDIF
C
C
      ELSE
C-----------------------------------------------------------------------
C     EN PROFONDEUR FINIE ....
C-----------------------------------------------------------------------
C
        IF (.NOT.SPHE) THEN
C       ----------------------------------------------------------------
C       ... ET COORDONNEES CARTESIENNES
C       ----------------------------------------------------------------
          DO IP=1,NPOIN2
            IF (DEPTH(IP).GT.PROMIN) THEN
              DO JP=1,NPLAN
                DDDN=-SINTET(JP)*DZX(IP)+COSTET(JP)*DZY(IP)
                CX(IP,JP)=CG(IP,JF)*COSTET(JP)
                CY(IP,JP)=CG(IP,JF)*SINTET(JP)
                DEUKD=2.0D0*XK(IP,JF)*DEPTH(IP)
                IF (DEUKD.GT.7.0D2) THEN
                  DSDNSK=0.0D0
                ELSE
                  DSDNSK=DEUPI*FREQ(JF)/SINH(DEUKD)
                ENDIF
                CT(IP,JP)=-DSDNSK*DDDN
              ENDDO
            ELSE
              DO JP=1,NPLAN
                CX(IP,JP)=0.0D0
                CY(IP,JP)=0.0D0
                CT(IP,JP)=0.0D0
              ENDDO
            ENDIF
          ENDDO
C
        ELSE
C       ----------------------------------------------------------------
C       ... ET COORDONNEES SPHERIQUES
C       ----------------------------------------------------------------
          SR=1.D0/R
          GRADEG=180.D0/3.1415926D0
          DO IP=1,NPOIN2
            IF (DEPTH(IP).GT.PROMIN) THEN
              SRCF=SR/COSF(IP)
              TFSR=TGF(IP)*SR
              DO JP=1,NPLAN
                DDDN=-SINTET(JP)*DZX(IP)*SR+COSTET(JP)*DZY(IP)*SRCF
                CX(IP,JP)=(CG(IP,JF)*COSTET(JP))*SR*GRADEG
                CY(IP,JP)=(CG(IP,JF)*SINTET(JP))*SRCF*GRADEG
                DEUKD=2.D0*XK(IP,JF)*DEPTH(IP)
                IF (DEUKD.GT.7.D2) THEN
                  DSDNSK=0.D0
                ELSE
                  DSDNSK=DEUPI*FREQ(JF)/SINH(DEUKD)
                ENDIF
                CT(IP,JP)=CG(IP,JF)*SINTET(JP)*TFSR-DSDNSK*DDDN*GRADEG
              ENDDO
            ELSE
              DO JP=1,NPLAN
                CX(IP,JP)=0.0D0
                CY(IP,JP)=0.0D0
                CT(IP,JP)=0.0D0
              ENDDO
            ENDIF
          ENDDO
C
        ENDIF
C
      ENDIF
C-----------------------------------------------------------------------
C
      RETURN
      END
