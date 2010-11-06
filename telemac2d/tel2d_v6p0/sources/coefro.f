C                       *****************
                        SUBROUTINE COEFRO
C                       *****************
C
     *(CF,H,U,V,KARMAN,KFROT,CHESTR,GRAV,MESH,T1)
C
C***********************************************************************
C  TELEMAC 2D VERSION 6.0  27/07/09   J-M HERVOUET (LNHE) 01 30 87 80 18
C                            
C***********************************************************************
C
C     FONCTION  : CALCUL DU COEFFICIENT DE FROTTEMENT CF
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      CF        |<-- | COEFFICIENT DE FROTTEMENT POUR K-EPSILON     |
C |       H        | -->| HAUTEUR D'EAU
C |     KARMAN     | -->| CONSTANTE DE KARMAN                          |
C |     KFROT      | -->| LOI DE FROTTEMENT SUR LE FOND                |
C |     CHESTR     | -->| TABLEAU DES COEFFICIENTS DE FROTTEMENT SUR LE|
C |                |    | FOND.
C |     GRAV       | -->| ACCELERATION DE LA PESANTEUR                 |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : TELMAC
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
      INTEGER, INTENT(IN)            :: KFROT
      DOUBLE PRECISION, INTENT(IN)   :: GRAV,KARMAN
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: CF,T1
      TYPE(BIEF_OBJ), INTENT(IN)     :: CHESTR,H,U,V
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NPOIN,N,IELMC,IELMH     
C
      DOUBLE PRECISION TIERS,HC,UNORM,AUX,INLOG,C
      DOUBLE PRECISION, POINTER :: HH(:)
C
      INTRINSIC SQRT,MAX,LOG
C
C-----------------------------------------------------------------------
C
      IELMC = CF%ELM
      IELMH = H%ELM
C
C  CONSTRUCTION D'UNE HAUTEUR AVEC LA MEME DISCRETISATION QUE CF
C  DANS LES CAS OU ON S'EN SERT.
C
      IF(KFROT.NE.0.AND.KFROT.NE.2) THEN
C
        IF(IELMC.EQ.IELMH) THEN
          HH=>H%R
        ELSE
          CALL OS( 'X=Y     ' , X=T1 , Y=H )
          CALL CHGDIS( T1 , IELMH , IELMC , MESH )
          HH=T1%R
        ENDIF
C
      ENDIF
C
      NPOIN = CF%DIM1
C
C-----------------------------------------------------------------------
C
      TIERS  = 1.D0/3.D0
C
C  CONSTRUCTION DU COEFFICIENT DE FROTTEMENT
C
C     LOIS DE FROTTEMENT :
C
C     KFROT = 0 :  PAS DE FROTTEMENT
C     KFROT = 1 :  LOI DE HAALAND
C     KFROT = 2 :  LOI DE CHEZY
C     KFROT = 3 :  LOI DE STRICKLER
C     KFROT = 4 :  LOI DE MANNING
C     KFROT = 5 :  LOI DE NIKURADSE
C
C     *******************
      IF(KFROT.EQ.0) THEN
C     *******************
C
        DO N=1,NPOIN
          CF%R(N) = 0.D0
        ENDDO
C
C     ***********************
      ELSEIF(KFROT.EQ.1) THEN
C     ***********************
C
        DO N=1,NPOIN
          HC = MAX(HH(N),1.D-4)
          UNORM = MAX(SQRT(U%R(N)**2+V%R(N)**2),1.D-6)
C                       1.D-6 : VISCOSITE LAMINAIRE DE L'EAU
          INLOG =(6.9D0*1.D-6/4.D0/HC/UNORM)**3+
     *                                  (CHESTR%R(N)/14.8D0/HC)**3.33
          INLOG = MIN(1.D0-1.D-6,INLOG)
          AUX   = -0.6D0*LOG(INLOG)/LOG(10.D0)
          CF%R(N) = 0.25D0 / AUX**2
        ENDDO
C
C     ***********************
      ELSEIF(KFROT.EQ.2) THEN
C     ***********************
C
        DO N=1,NPOIN
          CF%R(N) = 2 * GRAV / CHESTR%R(N)**2
        ENDDO
C
C     ***********************
      ELSEIF(KFROT.EQ.3) THEN
C     ***********************
C
        DO N=1,NPOIN
          HC = MAX(HH(N),1.D-4)
          CF%R(N) = 2 * GRAV / CHESTR%R(N)**2 / HC**TIERS
        ENDDO
C
C     ***********************
      ELSEIF(KFROT.EQ.4) THEN
C     ***********************
C
        DO N=1,NPOIN
          HC = MAX(HH(N),1.D-4)
          CF%R(N) = 2 * CHESTR%R(N)**2 * GRAV / HC**TIERS
        ENDDO
C
C     ***********************
      ELSEIF(KFROT.EQ.5) THEN
C     ***********************
C
        DO N=1,NPOIN
          HC = MAX(HH(N),1.D-4)
          CF%R(N) = 2.D0 / (LOG( 11.D0*HC/CHESTR%R(N))/KARMAN )**2
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
