C                         ********************
                          SUBROUTINE SD_FABSG4
C                         ********************
C
     * (NPOIN,NSEG,DAB1,DAB2,DAB3,DAB4,XAB1,XAB2,XAB3,XAB4,
     *  NPBLK,NSEGBLK,DA,XA)
C
C***********************************************************************
C BIEF VERSION 5.7     20/11/06   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  FONCTION: CONSTRUCTION DE LA MATRICE EN UN SEUL BLOC
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   NPOIN        | -->| NOMBRE D'INCONNUES D'UNE MATRICE DU BLOC
C |   NSEG         | -->| NOMBRE DE SEGMENTS 
C |   DA,XA        | -->| DIAGONALES ET TERMES EXTRA-DIAGONAUX DES
C |                |    | MATRICES
C |   XX1,XX2      |<-- | SOLUTIONS
C |   CVB1,CVB2    | -->| SECONDS MEMBRES
C |   INFOGR       | -->| IF, YES INFORMATIONS ON LISTING
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_SD_FABSG4 => SD_FABSG4
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NSEGBLK,NPBLK,NSEG,NPOIN
      DOUBLE PRECISION, INTENT(IN)    :: XAB1(NSEG),XAB2(NSEG)
      DOUBLE PRECISION, INTENT(IN)    :: XAB3(NSEG),XAB4(NSEG)
      DOUBLE PRECISION, INTENT(IN)    :: DAB1(NPOIN),DAB2(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: DAB3(NPOIN),DAB4(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: XA(2*NSEGBLK),DA(NPBLK)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,ISEG,JSEG
C      
C----------------------------------------------
C     INFO :      NPBLK   = NPOIN*NBLOC
C                 NSEGBLK = NSEG*4 + 2*NPOIN
C----------------------------------------------
C
C
C
C-------------------
C 1.  Diagonale bloc
C------------------- 
C          
      DO I=1,NPOIN 	     
         DA(I) = DAB1(I)
         DA(I+NPOIN) = DAB4(I)
      ENDDO 
C
C---------------------------      
C 2.   Termes Extradiagonaux  
C---------------------------
C      
C
C     BLOC 1 
C     ------
C  
      JSEG=0
      DO ISEG=1,NSEG
         JSEG=JSEG+1
         XA(JSEG)        =XAB1(ISEG)
         XA(JSEG+NSEGBLK)=XAB1(ISEG)
      ENDDO	    
C      
C     BLOC 2 ET 3 (EXTRA-DIAG)
C     ------------------------  
      DO I=1,NPOIN
         JSEG=JSEG+1
         XA(JSEG)        =DAB2(I)
         XA(JSEG+NSEGBLK)=DAB3(I)	    	    
      ENDDO	    
C
      DO ISEG=1,NSEG
         JSEG=JSEG+1
         XA(JSEG)        =XAB2(ISEG)
         XA(JSEG+NSEGBLK)=XAB3(ISEG)
         JSEG=JSEG+1
         XA(JSEG)        =XAB2(ISEG)
         XA(JSEG+NSEGBLK)=XAB3(ISEG)
      ENDDO
C	    
C     BLOC 4 (EXTRA) 
C     --------------     
      DO ISEG=1,NSEG
         JSEG=JSEG+1
         XA(JSEG)        =XAB4(ISEG)
         XA(JSEG+NSEGBLK)=XAB4(ISEG)
      ENDDO
C
C-----------------------------------------------------------------------
C	    
      RETURN
      END
