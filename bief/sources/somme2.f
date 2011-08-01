C                       ******************************** 
                        DOUBLE PRECISION FUNCTION SOMME2                 
C                       ********************************                 
C                                                                       
     *( X , NPX )                                                       
C                                                                       
C***********************************************************************
C BIEF VERSION 5.1          08/12/98                   A. DESITTER (NAG)
c                                       
C***********************************************************************
C                                                                       
C   FONCTION : SOMME DES COMPOSANTES D'UN VECTEUR EN MINIMISANT
C              LES ERREURS DE TRONCATURE.                       
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C |      X         | -->| TABLEAU FORTRAN                               
C |      NPX       | -->| NOMBRE DE VALEURS A ADDITIONNER               
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE) 
C-----------------------------------------------------------------------
C                                                                       
C PROGRAMMES APPELES   :                                                
C                                                                       
C***********************************************************************
C                                                                       
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                       
      INTEGER, INTENT(IN) :: NPX                                                     
C                                                                       
      DOUBLE PRECISION, INTENT(IN) :: X(*)                                             
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I      
      DOUBLE PRECISION C,Y,T
C
C-----------------------------------------------------------------------
c Kahan Summation Formula
c
c Goldberg, David  1991. What Every Computer Scientist Should Know 
c About Floating-Point Arithmetic, ACM  Computing Surveys, 23(1),
c pp5-48 (Corrigendum, Computing Surveys, 1991, 23(3))
c (article reproduced in the Numerical Computation Guide of Sun 
c  Microsystem, http://docs.sun.com )
c
c Knuth, D.E. 1981. The Art of Programming. Addison-Wesley, Reading,
c Mass., vol. II, 2nd ed.
c Proof p 572
C                                                                       
      SOMME2 = X(1)
      C = 0.D0
      DO I = 2 , NPX 
         Y = X(I) - C
         T = SOMME2 + Y
         C = (T - SOMME2) - Y
         SOMME2 = T
      ENDDO
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END                                                               


 
