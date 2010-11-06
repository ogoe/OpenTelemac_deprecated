C                       *****************
                        SUBROUTINE SMAGOR
C                       *****************
C
     *(VISC,CF,U,V,MESH,T1,T2,T3,T4,MSK,MASKEL,PROPNU)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.6    06/10/97    ADRIAN KLINGS (ENPC)
C
C***********************************************************************
C
C    FONCTION:
C
C    CE SOUS-PROGRAMME CALCULE LES VISCOSITES GRACE AU MODELE DE 
C    SMAGORINSKY
C                                                                     
C                                       (1/2)                     
C    NU    =   CS2 * ( 2.0 * Sij * Sij )  * (MESH SIZE)**2         
C                                                                      
C                                                                      
C                                                                      
C                         2        2            2                          
C                      DU       DV     DU   DV                          
C     2*Sij*Sij = ( 2*(--) + 2*(--) + (-- + --)                            
C                      DX       DY     DY   DX 
C       
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     VISC       | -->| DIFFUSION TURBULENTE                         |
C |      CF        | -->| COEFFICIENT DE FROTTEMENT POUR K-EPSILON     |
C |     U , V      | -->| COMPOSANTES DE LA VITESSE                    |
C |   T1,2,3,4     | -- | TABLEAUX DE TRAVAIL                          |
C |     MSK        | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |     MASKEL     | -->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMMES APPELES : CLIP , GRADF , LUMPIN , MATDIF , MATMAS ,
C                           MATMAT , MATVEC , OV , PRIDIR , SOLV01
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
      LOGICAL, INTENT(IN) :: MSK
      DOUBLE PRECISION, INTENT(IN)   :: PROPNU
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: VISC,T1,T2,T3,T4
      TYPE(BIEF_OBJ), INTENT(IN)     :: MASKEL,CF,U,V
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER N,NPOIN,IELMU,IELMC
      DOUBLE PRECISION CS,CS2
C
C-----------------------------------------------------------------------
C
      INTRINSIC SQRT
C
C-----------------------------------------------------------------------
C
      CS = 0.1D0
      CS2 = CS**2
C
C-----------------------------------------------------------------------
C
      IELMU = U%ELM
      IELMC = VISC%ELM
C
C     COMPUTATION OF GRADIENTS (IN FACT AVERAGED GRADIENT MULTIPLIED BY
C     A SURFACE WHICH IS THE INTEGRAL OF TEST FUNCTIONS ON THE DOMAIN,
C     THIS SURFACE IS CONSIDERED TO BE (MESH SIZE)**2
C
      CALL VECTOR(T1,'=','GRADF          X',IELMU,
     *            1.D0,U,U,U,U,U,U,MESH,MSK,MASKEL)
      CALL VECTOR(T2,'=','GRADF          Y',IELMU,
     *            1.D0,U,U,U,U,U,U,MESH,MSK,MASKEL)
      CALL VECTOR(T3,'=','GRADF          X',IELMU,
     *            1.D0,V,V,V,V,V,V,MESH,MSK,MASKEL)
      CALL VECTOR(T4,'=','GRADF          Y',IELMU,
     *            1.D0,V,V,V,V,V,V,MESH,MSK,MASKEL)
C
      NPOIN = VISC%DIM1
C
      DO N = 1,NPOIN
        VISC%R(N)=PROPNU+SQRT((2*T1%R(N)**2+2*T4%R(N)**2
     *                                       +(T2%R(N)+T3%R(N))**2))*CS2
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
