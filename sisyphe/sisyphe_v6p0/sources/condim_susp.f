C                       **********************
                        SUBROUTINE CONDIM_SUSP
C                       **********************
C
C      
     *(CS,CS0,NSICLA,X,Y,AT,NPOIN)
C
C***********************************************************************
C SISYPHE VERSION 5.9                   M. GONZALES DE LINARES    2004
C                                                
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT      
C***********************************************************************
C
C     FONCTION  : VALEURS IMPOSEES INITIALES
C                         - DE LA CONCENTRATION POUR CHAQUE TRACEUR
C (la routine condim_sisyphe.f est lue meme si CHARR=NO)
C
C     FUNCTION  : INITIALIZATION OF SUSPENDED SEDIMENT CONCENTRATION
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   CSi          |<-- | SUSPENDED SEDIMENT CONCENTRATION FOR CLASS I
C |________________|____|______________________________________________
C MODE : -->(INPUT), <--(RESULT), <-->(MODIFIED DATA)
C-----------------------------------------------------------------------
C CALLED BY : SISYPHE,CALCUL_SUSP
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
      INTEGER, INTENT(IN)           :: NPOIN,NSICLA 
      DOUBLE PRECISION,INTENT(IN)   :: AT,CS0(NSICLA)
      DOUBLE PRECISION,INTENT(IN)   :: X(NPOIN),Y(NPOIN)
C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: CS
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C
C-----------------------------------------------------------------------
C   
C  --------------------------------------------------------------
C  INITIALISATION DES TABLEAUX NON LUS DANS LE FICHIER RESULTATS:
C  --------------------------------------------------------------
C
      IF(NSICLA.GT.0) THEN
        DO I=1,NSICLA
          CALL OS('X=C     ',X=CS%ADR(I)%P,C=CS0(I))
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
     
