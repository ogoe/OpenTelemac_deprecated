C                      *************************
                       SUBROUTINE COEFRO_SISYPHE
C                      *************************
C
     *(CF,H,KFROT,CHESTR,GRAV,NPOIN,HMIN,KARMAN)
C
C***********************************************************************
C  SISYPHE VERSION 5.4                   C. VILLARET (LNHE)   01/10/2003
C
C***********************************************************************
C
C     FONCTION  : CALCUL DU COEFFICIENT DE FROTTEMENT QUADRATIQUE CF
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|   
C |       CF       | <--|  COEFFICIENT DE FROTTEMENT QUADRATIQUE       |   
C |       H        | -->| HAUTEUR D'EAU
C |     KFROT      | -->| LOI DE FROTTEMENT SUR LE FOND                |
C |     CHESTR     | -->| COEFFICIENTS DE FROTTEMENT SUR LE  FOND.
C |     GRAV       | -->| ACCELERATION DE LA PESANTEUR        
C |    NPOIN       | -->| NOMBRE DE POINTS
C |    HMIN        | -->| HAUTEUR D'EAU MINIMALE
C |    KARMAN      | -->| CONSTANTE DE KARMAN                          |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : TOB_SISYPHE
C
C  SOUS-PROGRAMME APPELES : OV
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
      INTEGER, INTENT(IN):: NPOIN,KFROT
      DOUBLE PRECISION,INTENT(IN):: GRAV,KARMAN,HMIN
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: CF
      TYPE(BIEF_OBJ),INTENT(IN) :: CHESTR,H
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER N
      DOUBLE PRECISION HC, AUX, TIERS
      INTRINSIC MAX,LOG
C
C-----------------------------------------------------------------------
C
      TIERS  = 1.D0/3.D0
C
C  CONSTRUCTION DU COEFFICIENT DE FROTTEMENT
C
C     LOIS DE FROTTEMENT :
C
C     KFROT = 0 :  FOND PLAT  (KS=3D50)
C     KFROT = 1 :  RIDES A L'EQUILIBRE (HOULE SEULEMENT) KS=(MAX 3D50,ETA)
C     KFROT = 2 :  LOI DE CHEZY
C     KFROT = 3 :  LOI DE STRICKLER
C     KFROT = 4 :  LOI DE MANNING
C     KFROT = 5 :  LOI DE NIKURADSE
C
      DO N=1,NPOIN
        IF(CHESTR%R(N).LE.0.D0) THEN
          WRITE(LU,*) 'FROTTEMENT NON DEFINI DANS COEFRO AU POINT ',N
          CALL PLANTE(1)
          STOP
        ENDIF
      ENDDO
C
C     ***********************
      IF(KFROT.LE.1.OR.KFROT.EQ.5) THEN
C     *********************** 
        AUX=30.D0/EXP(1.D0)
        DO N=1,NPOIN    
           HC = MAX(H%R(N),HMIN)
           CF%R(N) = 2.D0 / (LOG( AUX*HC/CHESTR%R(N))/KARMAN )**2
        ENDDO 
C SUGGESTION:
C       DO N=1,NPOIN 
C          AUX=MAX( 1.001D0 , H%R(N)/EXP(1.D0)/(CHESTR%R(N)/30.D0) )  
C          CF%R(N) = 2.D0 / (  LOG(AUX) / KARMAN  )**2
C       ENDDO      
C     ***********************
      ELSEIF(KFROT.EQ.2) THEN
C     ***********************
C
        DO N=1,NPOIN
           CF%R(N) = 2.D0 * GRAV / CHESTR%R(N)**2
        ENDDO
C
C     ***********************
      ELSEIF(KFROT.EQ.3) THEN
C     ***********************
C
        DO N=1,NPOIN
           HC = MAX(H%R(N),HMIN)
           CF%R(N) = 2.D0 * GRAV / CHESTR%R(N)**2 / HC**TIERS
        ENDDO
C
C     ***********************
      ELSEIF(KFROT.EQ.4) THEN
C     ***********************
C
        DO N=1,NPOIN
           HC = MAX(H%R(N),HMIN)
           CF%R(N) = 2.D0 * CHESTR%R(N)**2 * GRAV / HC**TIERS
        ENDDO
C
C     ****
      ELSE
C     ****
C
        IF(LNG.EQ.1) WRITE(LU,300) KFROT
        IF(LNG.EQ.2) WRITE(LU,301) KFROT
300     FORMAT(1X,'COEFRO : LOI DE FROTTEMENT INCONNUE :',1I6)
301     FORMAT(1X,'COEFRO: UNKNOWN LAW OF BOTTOM FRICTION: ',1I6)
        CALL PLANTE(1)
        STOP
C
C     *****
      ENDIF
C     *****
C
C-----------------------------------------------------------------------
C
      RETURN
      END
