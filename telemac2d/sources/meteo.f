C                       ****************
                        SUBROUTINE METEO
C                       ****************
C
     *(PATMOS,WINDX,WINDY,FUAIR,FVAIR,X,Y,AT,LT,NPOIN,VENT,ATMOS,
     * HN,TRA01,GRAV,ROEAU,NORD,PRIVE)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.4    02/01/04  J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C     FONCTION  :  CALCUL DES CHAMPS DE VENT ET DE PRESSION
C                  EN GENERAL A PARTIR DE FICHIERS DE DONNEES
C
C     CE SOUS-PROGRAMME PEUT ETRE COMPLETE PAR L'UTILISATEUR
C
C-----------------------------------------------------------------------
C
C     FUNCTION:  SETTING ATMOSPHERIC PRESSURE AND WIND VELOCITIES
C
C     MUST BE ADAPTED BY USER.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    PATMOS      |<-- | ATMOSPHERIC PRESSURE
C |    WINDX,Y     |<-- | TWO COMPONENTS OF WIND VELOCITY
C | FUAIR,FVAIR    | -->| IDEM IF WIND CONSTANT.
C |   X , Y        | -->| COORDINATES OF POINTS IN THE MESH
C |     AT,LT      | -->| TIME, ITERATION NUMBER
C |    NPOIN       | -->| NUMBER OF POINTS IN THE MESH
C |    VENT        | -->| YES IF WIND TAKEN INTO ACCOUNT
C |    ATMOS       | -->| YES IF PRESSURE TAKEN INTO ACCOUNT
C |     HN         | -->| DEPTH
C |    TRA01       | -->| WORKING ARRAY
C |    GRAV        | -->| GRAVITY ACCELERATION
C |   ROEAU        | -->| WATER DENSITY
C |   NORD         | -->| DIRECTION OF NORTH, COUNTER-CLOCK-WISE
C |                |    | STARTING FROM VERTICAL AXIS
C |   PRIVE        | -->| USER WORKING ARRAYS (BIEF_OBJ BLOCK)
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS PROGRAMME APPELE : OV
C
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C      
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: LT,NPOIN
      LOGICAL, INTENT(IN)             :: ATMOS,VENT
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN),Y(NPOIN),HN(NPOIN)  
      DOUBLE PRECISION, INTENT(INOUT) :: WINDX(NPOIN),WINDY(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: PATMOS(NPOIN),TRA01(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: FUAIR,FVAIR,AT,GRAV,ROEAU,NORD
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: PRIVE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+      
C
      DOUBLE PRECISION P0,Z(1)
C
C-----------------------------------------------------------------------
C
C     BEWARE, HERE ONLY ONE COMPUTATION AT FIRST TIME-STEP
C
      IF(LT.EQ.0) THEN
C
C-----------------------------------------------------------------------
C
C     ATMOSPHERIC PRESSURE
C
      IF(ATMOS) THEN
        P0 = 100000.D0
        CALL OV( 'X=C     ' , PATMOS , Y , Z , P0 , NPOIN )
      ENDIF
C
C-----------------------------------------------------------------------
C
C VENT : ICI ON LE PREND CONSTANT AVEC LES VALEURS DONNEES DANS LE
C        FICHIER CAS.
C
C     SUIVANT LE REPERE DANS LEQUEL LA VITESSE DU VENT EST FOURNIE
C     IL PEUT Y AVOIR UNE ROTATION A FAIRE.
C
      IF(VENT) THEN
        CALL OV( 'X=C     ' , WINDX , Y , Z , FUAIR , NPOIN )
        CALL OV( 'X=C     ' , WINDY , Y , Z , FVAIR , NPOIN )
      ENDIF
C
C     FIN DE IF(LT.EQ.0)
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
