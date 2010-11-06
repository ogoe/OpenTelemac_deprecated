C                       ****************                                                        
                        SUBROUTINE CALDT                                    
C                       ****************                                    
C     
     *(NS,G,H,U,V,DTHAUT,DT,CFL)
C
C***********************************************************************
C TELEMAC-2D VERSION 5.2         22/10/01  ?????? TEL: ??????
C                                                                                                      
C***********************************************************************
C     
C     FUNCTION : CALCUL DU PAS DE TEMPS SATISFAISANT LA CONDITION CFL
C          
C-----------------------------------------------------------------------
C     ARGUMENTS                                         
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |                | <--|  
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C     
C     APPELE PAR :  HOMERE_TELEMAC2D
C     
C     SOUS-PROGRAMME APPELE : OS
C     
C**********************************************************************
C
      IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NS
      DOUBLE PRECISION, INTENT(INOUT) :: DT
      DOUBLE PRECISION, INTENT(IN) :: H(NS),U(NS),V(NS),DTHAUT(NS)
      DOUBLE PRECISION, INTENT(IN) :: G,CFL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C  
      INTEGER IS
C
      DOUBLE PRECISION RA3,EPSL 
      DOUBLE PRECISION SIGMAX ,UA2, UA3, UNORM  
C
C-----------------------------------------------------------------------
C
      RA3 = SQRT(1.5D0*G)
C
      DT = 1.E+12
      EPSL = 0.01D0
C
      DO IS=1,NS
        SIGMAX = H(IS)
        UA2    = U(IS)
        UA3    = V(IS)                                       
        UNORM=SQRT(UA2*UA2 + UA3*UA3)
        SIGMAX= MAX(EPSL, RA3*SQRT(SIGMAX) +UNORM )
        DT = MIN(DT, CFL*DTHAUT(IS)/SIGMAX)
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
