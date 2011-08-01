C                       *****************
                        SUBROUTINE CHPCON
C                       *****************
C
     *(UCONV,VCONV,U,V,UN,VN,TETAU)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2    17/08/94    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C     FONCTION  :  CALCUL DU CHAMP CONVECTEUR UCONV,VCONV
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   UCONV,VCONV  | -->| COMPOSANTES DU CHAMP CONVECTEUR.             |
C |   U,V          | -->| COMPOSANTES DE LA VITESSE.                   |
C |   UN,VN        | -->| COMPOSANTES DE LA VITESSE A L'ETAPE N.       |
C |   TETAU        | -->| IMPLICITATION SUR U.                         |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMME APPELE : OV
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
      DOUBLE PRECISION, INTENT(IN)  :: TETAU
      TYPE(BIEF_OBJ), INTENT(IN)    :: U,UN,V,VN
      TYPE(BIEF_OBJ), INTENT(INOUT) :: UCONV,VCONV
C    
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C  FABRICATION DES CONVECTEURS UCONV ET VCONV
C
      CALL OS( 'X=CY    ' , UCONV , UN , U , 1.D0-TETAU )
      CALL OS( 'X=X+CY  ' , UCONV , U  , U ,      TETAU )
      CALL OS( 'X=CY    ' , VCONV , VN , U , 1.D0-TETAU )
      CALL OS( 'X=X+CY  ' , VCONV , V  , U ,      TETAU )
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
