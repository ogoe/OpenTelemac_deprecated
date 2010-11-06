C                       ********************
                        SUBROUTINE SD_STRTRI
C                       ********************
C
     *(IS,N,IND)
C
C***********************************************************************
C BIEF VERSION 5.7     20/11/06   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  IMPORTANT NOTICE: INSPIRED FROM N3S 3.3  22/04/92  B.THOMAS 
C
C  FONCTION : TRI EN ORDRE CROISSANT DU TABLEAU D'ENTIERS IS :
C             EN SORTIE  IS(IND(I+1) >= IS(IND(I)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |  IS(I=1,N)     | -->| TABLEAU A TRIER  
C |  N             | -->| LONGUEUR DE IS                        
C |  IND(I=1,N)    |<-- | CF CI-DESSUS                 
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_SD_STRTRI => SD_STRTRI
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: N
      INTEGER, INTENT(IN)    :: IS(N)
      INTEGER, INTENT(INOUT) :: IND(N)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,I1,K
C
C-----------------------------------------------------------------------
C
      IND(1) = 1
C
      DO 1 I = 2 , N
C
C--->    IS(1:I-1) EST TRIE
C
         I1 = I-1
         DO 5 K = I1 , 1 , -1
C
C--->       POUR L > K+1  IS(IND(L)) > IS(I)
C
            IF(IS(IND(K)).GT.IS(I)) THEN
              IND(K+1) = IND(K)
            ELSE
              GO TO 2
            ENDIF
C
5        CONTINUE
C
C--->    ASSERTION : IS(IND(K)) <= IS(I)
C
2        IND(K+1)=I
C
1     CONTINUE       
C
C-----------------------------------------------------------------------
C
      RETURN
      END
