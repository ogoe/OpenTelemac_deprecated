C                       ********************
                        SUBROUTINE SD_STRSSD
C                       ********************
C
     *(NPBLK,NSEGBLK,GLOSEG1,GLOSEG2,IN,IP,ISEGIP,IW)   
C
C***********************************************************************
C BIEF VERSION 5.9     18/02/08    E.RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C    FONCTION : CONSTRUCTION D'UN STOCKAGE COMPACT 
C               (IN,IP) = (XADJ, ADJNCY) DES TERMES EXTRADIAGONAUX 
C               VIA LE STOCKAGE SEGMENT
C 
C
C    ATTENTION : NSEG ET IW SONT LE MEME TABLEAU DANS LA MEMOIRE               
C
C***********************************************************************
C
      USE BIEF,EX_SD_STRSSD => SD_STRSSD
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NSEGBLK,NPBLK
      INTEGER, INTENT(IN)    :: GLOSEG1(NSEGBLK),GLOSEG2(NSEGBLK)
      INTEGER, INTENT(INOUT) :: IN(NPBLK+1),IP(NSEGBLK*2+1)
      INTEGER, INTENT(INOUT) :: ISEGIP(NSEGBLK*2+1)
      INTEGER, INTENT(INOUT) :: IW(NPBLK)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,J1,J2,ISEG      
C      
C  -----------------------
C
C---> DEGRE D'UN POINT : NOMBRE DE POINTS VOISIN
C      
      DO I=1,NPBLK 
        IW(I)=0
      ENDDO 
C
C     IW : NOMBRE DE VOISINS DE CHAQUE POINT
C      
      DO ISEG=1,NSEGBLK
        J1 = GLOSEG1(ISEG)	    
        J2 = GLOSEG2(ISEG)
        IW(J1) = IW(J1) + 1
        IW(J2) = IW(J2) + 1	   
      ENDDO
C
C---> STOCKAGE COMPACT SANS LA DIAGONALE : (XADJ,ADJNCY) = (IN,IP)
C
C     COEFFICIENTS DU POINT I : DE IN(I) A IN(I+1)-1
C
      IN(1)=1
      DO I=1,NPBLK
        IN(I+1)=IN(I)+IW(I)
      ENDDO
C    
C     MAINTENANT IW N'EST PLUS LE NOMBRE DE VOISINS
C
      DO I=1,NPBLK
        IW(I)=IN(I)
      ENDDO
C
      DO ISEG=1,NSEGBLK
        J1 = GLOSEG1(ISEG)	    
        J2 = GLOSEG2(ISEG)		    
C--> TABLEAU DE CONNECTIVITE : SEGMENT ---> POINT	    	    	    
        IP(IW(J1))=J2
        IP(IW(J2))=J1	
C--> TABLEAU DE CONNECTIVITE INVERSE : POINT ---> SEGMENT
C    NOTATION POUR COEF. TRIANGULAIRE SUPERIEUR            	    
        ISEGIP(IW(J1))=-ISEG
C    NOTATION POUR COEF. TRIANGULAIRE INFERIEUR            	    	  
        ISEGIP(IW(J2))= ISEG
C		 
        IW(J1) = IW(J1) + 1
        IW(J2) = IW(J2) + 1    	       
      ENDDO
C
C-----------------------------------------------------------------------
C          
      RETURN
      END
