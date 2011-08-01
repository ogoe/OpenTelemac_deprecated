!                       ******************
                        SUBROUTINE ELEBD31
!                       ******************
!
     *(NELBOR,NULONE,IKLBOR,IFABOR,NBOR,IKLE,
     * NELEM,NELEB,NELMAX,NPOIN,NPTFR,IELM)     
!    
!***********************************************************************
! BIEF VERSION 5.5           09/04/04    J-M HERVOUET (LNH) 30 87 80 18
! 
!***********************************************************************
!              
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________.
! |      NOM       |MODE|                   ROLE                       |
! |________________|____|______________________________________________|
! |    NELBOR      |<-- | NUMERO DE L'ELEMENT ADJACENT AU KIEME SEGMENT|
! |    NULONE      |<-- | NUMERO LOCAL D'UN POINT DE BORD DANS         |
! |                |    | L'ELEMENT ADJACENT DONNE PAR NELBOR          
! |    IKLBOR      |<-- | NUMERO LOCAL DES NOEUDS A PARTIR D'UN ELEMENT
! |                |    |  DE BORD
! |    IFABOR      | -->| TABLEAU DES VOISINS DES FACES.
! |    NBOR        | -->| NUMERO GLOBAL D'UN NOEUD A PARTIR DU N° LOCAL
! |    IKLE        | -->| NUMEROS GLOBAUX DES POINTS DE CHAQUE ELEMENT.
! |    NELEM       | -->| NOMBRE TOTAL D'ELEMENTS DANS LE MAILLAGE.
! |    NPOIN       | -->| NOMBRE TOTAL DE POINTS DU DOMAINE.
! |    NPTFR       | -->| NOMBRE DE POINTS FRONTIERES.
C |    NELEB       | -->| NOMBRE D'ELEMENTS DE BORD.
! |    IELM        | -->| TYPE D'ELEMENT.
! |________________|____|______________________________________________|
!  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
!
!-----------------------------------------------------------------------
!
! Subroutine written by Lam Minh-Phuong
!
!-----------------------------------------------------------------------
!
! Subroutine: elebd31
!
! Function: Construction de NELBOR, NULONE, IKLBORD 
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      USE BIEF, EX_ELEBD31 => ELEBD31
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)    :: NELEM,NELEB,NELMAX
      INTEGER, INTENT(IN)    :: NPOIN,NPTFR,IELM
      INTEGER, INTENT(IN)    :: NBOR(NPTFR)
      INTEGER, INTENT(IN)    :: IFABOR(NELMAX,4)
      INTEGER, INTENT(IN)    :: IKLE(NELEM,4)
      INTEGER, INTENT(OUT)   :: NELBOR(NELEB),NULONE(NELEB,3)
      INTEGER, INTENT(OUT)   :: IKLBOR(NELEB,3)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER   :: IELEM, IELEB, J,K,IPOIN
      INTEGER   :: IPOBO(NPOIN)
!
      INTEGER SOMFAC(3,4)
      DATA SOMFAC / 1,2,3 , 4,1,2 , 2,3,4 , 3,4,1  /
!     face numero:    1       2       3       4
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!

      IF(IELM /= 31) THEN

        IF(LNG.EQ.1) WRITE(LU,98) IELM
        IF(LNG.EQ.2) WRITE(LU,99) IELM
98      FORMAT(1X,'VOISIN: IELM=',1I6,' TYPE D''ELEMENT NON PREVU')
99      FORMAT(1X,'VOISIN: IELM=',1I6,' TYPE OF ELEMENT NOT AVAILABLE')
        CALL PLANTE(1)
        STOP
      ENDIF
      
      
! Creation de IPOBO qui permet de passer du numero global au numero local
      
      DO IPOIN=1,NPOIN
        IPOBO(IPOIN) = 0
      ENDDO
      
      DO K = 1, NPTFR
         IPOBO(NBOR(K)) = K
      ENDDO
           
       
! Construction de NELBOR, NULONE, IKLBORD 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IELEB = 0
      
      DO IELEM = 1,NELEM
         DO J = 1,4
            IF (IFABOR(IELEM,J)== 0) THEN
               IELEB           = IELEB + 1
               IF ( IELEB .GT. NELEB ) THEN
                 IF(LNG.EQ.1) WRITE(LU,101)
                 IF(LNG.EQ.2) WRITE(LU,102)
101              FORMAT(1X,'ELEBD31 : Erreur dans le Maillage')
102              FORMAT(1X,'ELEBD31 : Error in Mesh')
                 CALL PLANTE(1)
                 STOP
               END IF
               NELBOR(IELEB)   = IELEM
               NULONE(IELEB,1) = SOMFAC(1,J)
               NULONE(IELEB,2) = SOMFAC(2,J)
               NULONE(IELEB,3) = SOMFAC(3,J)
               IKLBOR(IELEB,1) = IPOBO(IKLE(NELBOR(IELEB),SOMFAC(1,J)))
               IKLBOR(IELEB,2) = IPOBO(IKLE(NELBOR(IELEB),SOMFAC(2,J)))
               IKLBOR(IELEB,3) = IPOBO(IKLE(NELBOR(IELEB),SOMFAC(3,J)))               
            END IF
         END DO
      END DO
 
!
!-----------------------------------------------------------------------
!
      RETURN
      END SUBROUTINE ELEBD31    
