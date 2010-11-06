C
C  DECLARATIONS COMMON TO ALL PROGRAMMES
C
C  VERSION 6.0
C
      MODULE DECLARATIONS_TELEMAC
C
C----------------------------------------------------------------------
C
C 1./ INTEGER VALUES FOR DESCRIBING BOUNDARY CONDITIONS:
C
C
C     FOR THE BOUNDARY CONDITIONS FILE:
C
C     ENTRANCE: PRESCRIBED VALUES (SAVE VELOCITIES)
      INTEGER, PARAMETER :: KENT  =  5
C
C     VELOCITY IMPOSED (INSTEAD OF DISCHARGE)
      INTEGER, PARAMETER :: KENTU =  6
C
C     FREE OUTPUT
      INTEGER, PARAMETER :: KSORT =  4
C
C     NO-SLIP CONDITION
      INTEGER, PARAMETER :: KADH  =  0
C
C     WALL WITH OR WITHOUT FRICTION
      INTEGER, PARAMETER :: KLOG  =  2
C
C     OPEN BOUNDARY WITH INCIDENT WAVE
      INTEGER, PARAMETER :: KINC  =  1
C
C     ESTEL-2D : DRAINAGE LIBRE
      INTEGER, PARAMETER :: KDRAIN  =  3
C
C     ESTEL-2D : CONDITION MIXTE
      INTEGER, PARAMETER :: KMIX  =  4
C
C     DEPENDING ON ALGORITHMS AND CASES, THESE VALUES WILL BE TRANSFORMED
C     INTO:
C
C     TECHNICAL BOUNDARY CONDITIONS
C
C     NEUMANN
      INTEGER, PARAMETER :: KNEU  =  1
C
C     DIRICHLET
      INTEGER, PARAMETER :: KDIR  =  2
C
C     DEGREE OF FREEDOM
      INTEGER, PARAMETER :: KDDL  =  3
C
C     INCIDENT WAVE
      INTEGER, PARAMETER :: KOND  =  4
C
C----------------------------------------------------------------------
C
C 2./ INTEGER VALUES FOR DESCRIBING ADVECTION SCHEMES:
C
C     CHARACTERISTICS
      INTEGER, PARAMETER :: ADV_CAR     =  1
C     SUPG
      INTEGER, PARAMETER :: ADV_SUP     =  2
C     LEO POSTMA
      INTEGER, PARAMETER :: ADV_LPO     =  3
C     DISTRIBUTIVE SCHEME N
      INTEGER, PARAMETER :: ADV_NSC     =  4
C     DISTRIBUTIVE SCHEME PSI
      INTEGER, PARAMETER :: ADV_PSI     =  5
C     NON CONSERVATIVE EQUATION, DISTRIBUTIVE SCHEME PSI
      INTEGER, PARAMETER :: ADV_PSI_NC  =  6
C     NON CONSERVATIVE EQUATION, DISTRIBUTIVE SCHEME N
      INTEGER, PARAMETER :: ADV_NSC_NC  =  7
C     LEO POSTMA, EDGE-BASED FOR TIDAL FLATS
      INTEGER, PARAMETER :: ADV_LPO_TF  = 13
C     DISTRIBUTIVE SCHEME N, EDGE-BASED FOR TIDAL FLATS
      INTEGER, PARAMETER :: ADV_NSC_TF  = 14
C     DISTRIBUTIVE SCHEME PSI, EDGE-BASED FOR TIDAL FLATS
      INTEGER, PARAMETER :: ADV_PSI_TF  = 15
C
C-----------------------------------------------------------------------
C
C 3./ CODE COUPLING
C       
      CHARACTER*144 COUPLING
C
C 4./ NAME OF CURRENT CODE (SEE BIEF_OPEN_FILES AND CONFIG_CODE)
C       
      CHARACTER(LEN=24) :: NAMECODE,NNAMECODE(3)
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C
C 2b./ Number for each code in TELEMAC System (this was for protection)
C
      INTEGER, PARAMETER :: TMCOD_ARTEMIS      = 1
      INTEGER, PARAMETER :: TMCOD_TOMAWAC      = 2
      INTEGER, PARAMETER :: TMCOD_ESTEL2D      = 3
      INTEGER, PARAMETER :: TMCOD_ESTEL3D      = 4
      INTEGER, PARAMETER :: TMCOD_MATISSE      = 5
      INTEGER, PARAMETER :: TMCOD_RUBENS       = 6
      INTEGER, PARAMETER :: TMCOD_SISYPHE      = 7
      INTEGER, PARAMETER :: TMCOD_STBTEL       = 8
      INTEGER, PARAMETER :: TMCOD_SUBIEF2D     = 9
      INTEGER, PARAMETER :: TMCOD_SUBIEF3D     =10
      INTEGER, PARAMETER :: TMCOD_TELEMAC2D    =11
      INTEGER, PARAMETER :: TMCOD_TELEMAC3D    =12
      INTEGER, PARAMETER :: TMCOD_POSTEL3D     =13
      INTEGER, PARAMETER :: TMCOD_TEL2DSIS     =13
      INTEGER, PARAMETER :: TMCOD_TEL3DSIS     =14
      INTEGER, PARAMETER :: TMCOD_TEL2DE3D     =15
      INTEGER, PARAMETER :: TMCOD_SPARTACUS2D  =16
C
C 2c./ Global Data (this was for protection)
C
      INTEGER IFVP,ICVP,IVVP
      DOUBLE PRECISION RVVP
      COMMON/INFOS/RVVP,ICVP,IFVP,IVVP
C
C Other functions (this was for protection)
C
      INTERFACE
        INTEGER FUNCTION ISQRT0 (IVAL)
          INTEGER  , INTENT(IN) :: IVAL
        END FUNCTION
      END INTERFACE
C
      INTERFACE
        INTEGER FUNCTION ISQRT (IVAL)
          INTEGER  , INTENT(IN) :: IVAL
        END FUNCTION
      END INTERFACE
C
      INTERFACE
        DOUBLE PRECISION FUNCTION RSQRT (RVAL)
          DOUBLE PRECISION , INTENT(IN) :: RVAL
        END FUNCTION
      END INTERFACE
C
      INTERFACE
        INTEGER FUNCTION ISQRTF (IVAL)
          INTEGER , INTENT(IN) :: IVAL
        END FUNCTION
      END INTERFACE
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
      CONTAINS    
C                       ****************************
                        CHARACTER*11 FUNCTION EXTENS
C                       ****************************
C
     *(N,IPID)
C
C***********************************************************************
C  BIEF VERSION 5.9       26/05/2008 J-M HERVOUET (LNHE)  01 30 87 80 18
C
C***********************************************************************
C
C      FONCTIONS: EXTENSION DES FICHIERS SUR CHAQUE PROCESSEUR.
C      ==========
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________|
C |     N          | -->| NOMBRE DE PROCESSEURS MOINS UN = NCSIZE-1
C |     IPID       | -->| NUMERO DU PROCESSEUR
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR :
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER, INTENT(IN) :: IPID,N
C
C-----------------------------------------------------------------------
C
      IF(N.GT.0) THEN
C
        EXTENS='00000-00000'
C
        IF(N.LT.10) THEN
          WRITE(EXTENS(05:05),'(I1)') N
        ELSEIF(N.LT.100) THEN
          WRITE(EXTENS(04:05),'(I2)') N
        ELSEIF(N.LT.1000) THEN
          WRITE(EXTENS(03:05),'(I3)') N
        ELSEIF(N.LT.10000) THEN
          WRITE(EXTENS(02:05),'(I4)') N
        ELSE
          WRITE(EXTENS(01:05),'(I5)') N
        ENDIF
C
        IF(IPID.LT.10) THEN
          WRITE(EXTENS(11:11),'(I1)') IPID
        ELSEIF(IPID.LT.100) THEN
          WRITE(EXTENS(10:11),'(I2)') IPID
        ELSEIF(IPID.LT.1000) THEN
          WRITE(EXTENS(09:11),'(I3)') IPID
        ELSEIF(IPID.LT.10000) THEN
          WRITE(EXTENS(08:11),'(I4)') IPID
        ELSE
          WRITE(EXTENS(07:11),'(I5)') IPID
        ENDIF
C
      ELSE
C
        EXTENS='       '
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END FUNCTION
C
C-----------------------------------------------------------------------
C
      END MODULE DECLARATIONS_TELEMAC

