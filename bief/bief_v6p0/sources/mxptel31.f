!                         *******************
                          SUBROUTINE MXPTEL31
!                         *******************
!      
     & (NELEM,NPOIN,MXELVS,IKLES,MXPTVS)
!          
!***********************************************************************
! BIEF VERSION 5.5           31/08/2009                  Lam Minh-Phuong
!                                                         F. DECUNG
!***********************************************************************          
!
! Function: Donne le nombre maximum de noeuds voisins pour un noeud  
!           du maillage.
!     
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________
! |      NOM       |MODE|                   ROLE
! |________________|____|______________________________________________
! |   MXPTVS       |<-- | NOMBRE MAXIMUM DE POINTS VOISINS.
! |   MXELVS       | -->| NOMBRE MAXIMUM D'ELEMENTS VOISINS.
! |   IKLES        | -->| TABLE DE CONNECTIVITE (DU FORMAT SELAFIN)
! |   NPOIN        | -->| NOMBRE DE POINTS DU MAILLAGE.
! |   NELEM        | -->| NOMBRE D'ELEMENTS DU MAILLAGE.
! |________________|____|______________________________________________
!
!  MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
!-----------------------------------------------------------------------
!    
      USE BIEF, EX_MXPTEL31 => MXPTEL31
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)                       :: NELEM
      INTEGER, INTENT(IN)                       :: NPOIN
      INTEGER, INTENT(IN)                       :: MXELVS
      INTEGER, INTENT(IN), DIMENSION(4,NELEM)   :: IKLES
      INTEGER, INTENT(OUT)                      :: MXPTVS
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C        
      INTEGER                            :: IPOIN,I,J,K,IKLJ,IKL
      INTEGER                            :: NVOIS
      INTEGER,DIMENSION(:),ALLOCATABLE   :: VOIS
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: IND_ELEM
!
      ALLOCATE(VOIS(3*MXELVS))
      ALLOCATE(IND_ELEM(NPOIN,MXELVS+1))
!
!-----------------------------------------------------------------------
!
! IND_ELEM donne le nombre d'elements autour d'un noeud et leurs numeros      
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      
      DO I = 1, NPOIN
        IND_ELEM(I,1) = 0
      ENDDO
!      
      DO J=1, 4      
        DO I=1,NELEM
          IKL = IKLES(J,I)
          IND_ELEM(IKL,1)=IND_ELEM(IKL,1)+1
          IND_ELEM(IKL,IND_ELEM(IKL,1)+1)=I                       
        ENDDO  
      ENDDO 
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!      
      MXPTVS = 0
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Boucle sur tous les noeuds du maillage !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!      
      DO IPOIN = 1, NPOIN
!      
        NVOIS   = 1
        VOIS(1) = IPOIN
!             
!       initialiation de VOIS a 0
!
        DO I = 1, 3*MXELVS
          VOIS(I) = 0
        ENDDO
!        
!       remplissage du tableau VOIS contenant tous les numeros des noeuds 
!       voisins de IPOIN. 
!        
        DO J = 1,4
          DO I = 2, IND_ELEM(IPOIN,1)+1         
            IKLJ = IKLES(J,IND_ELEM(IPOIN,I))           
            DO K = 1, NVOIS
              IF ( VOIS(K) == IKLJ ) EXIT
            ENDDO                      
            IF( K > NVOIS ) THEN
              NVOIS       = NVOIS + 1
              VOIS(NVOIS) = IKLJ
            ENDIF           
          ENDDO
        ENDDO
!        
        NVOIS = NVOIS - 1       
        IF( MXPTVS < NVOIS) MXPTVS = NVOIS
!              
      ENDDO          
!    
!-----------------------------------------------------------------------
!
      DEALLOCATE (VOIS)
      DEALLOCATE (IND_ELEM)
!
      RETURN
      END
