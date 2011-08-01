C                         ********************
                          SUBROUTINE SD_STRSG4
C                         ********************
C
     *(NPOIN,NSEG,GLOSEGB,NPBLK,NSEGBLK,GLOSEG4)
C
C***********************************************************************
C BIEF VERSION 5.7     20/11/06   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  FONCTION: CONSTRUCTION DES SEGMENTS DE LA MATRICE EN UN SEUL BLOC
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   NPOIN        | -->| NOMBRE D'INCONNUES D'UNE MATRICE DU BLOC
C |   NSEG         | -->| NOMBRE DE SEGMENTS 
C |   GLOSEGB      | -->| NUMEROS GLOBAUX DES POINTS DES SEGMENTS
C |   NPBLK        | -->| COMME NPOIN MAIS POUR LE BLOC
C |   NSEGBLK      | -->| COMME NSEG MAIS POUR LE BLOC
C |   GLOSEG4      |<-- | IF, YES INFORMATIONS ON LISTING
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_SD_STRSG4 => SD_STRSG4
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NSEGBLK,NPBLK,NSEG,NPOIN
      INTEGER, INTENT(IN)    :: GLOSEGB(NSEG,2)
      INTEGER, INTENT(INOUT) :: GLOSEG4(2*NSEGBLK)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C            
      INTEGER I,ISEG,JSEG
C      
C----------------------------------------------
C     INFO :      NPBLK = NPOIN*NBLOC
C                 NSEGBLK=NSEG*4 + 2*NPOIN 
C----------------------------------------------
C
C     MATRICE ASSEMBLE TOTAL BLOCS :  
C 
C----------------------------------------------
C
      JSEG=0
C
C     BLOC 1
C     ------   
C
      DO ISEG=1,NSEG
         JSEG=JSEG+1
         GLOSEG4(JSEG)= GLOSEGB(ISEG,1)
         GLOSEG4(JSEG+NSEGBLK)= GLOSEGB(ISEG,2)
      ENDDO	    
C
C     BLOC 2 ET 3 (EXTRA-DIAG)
C     ------------------------
C  
      DO I=1,NPOIN
         JSEG=JSEG+1
         GLOSEG4(JSEG)= I
         GLOSEG4(JSEG+NSEGBLK)= I+NPOIN
      ENDDO 	    
C     
      DO ISEG=1,NSEG
         JSEG=JSEG+1
         GLOSEG4(JSEG)= GLOSEGB(ISEG,1)
         GLOSEG4(JSEG+NSEGBLK)= GLOSEGB(ISEG,2)+ NPOIN
         JSEG=JSEG+1
         GLOSEG4(JSEG)= GLOSEGB(ISEG,2)
         GLOSEG4(JSEG+NSEGBLK)= GLOSEGB(ISEG,1)+ NPOIN
      ENDDO
C	    
C     BLOC 4 (EXTRA) 
C     --------------
C     
      DO ISEG=1,NSEG
         JSEG=JSEG+1
         GLOSEG4(JSEG)         = GLOSEGB(ISEG,1)+NPOIN
         GLOSEG4(JSEG+NSEGBLK) = GLOSEGB(ISEG,2)+NPOIN
      ENDDO	    
C
C-----------------------------------------------------------------------
C
      RETURN
      END
