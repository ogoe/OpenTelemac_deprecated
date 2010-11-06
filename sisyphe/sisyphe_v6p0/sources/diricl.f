C                       *****************
                        SUBROUTINE DIRICL
C                       *****************
C
     *( ZF1 , ZF , EBOR , LIEBOR , NBOR , NPOIN  , NPTFR  , KENT )
C
C***********************************************************************
C SISYPHE VERSION 5.1            10/97      C.LE NORMANT 01-30-87-78-54
C                                                
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT                                       
C***********************************************************************
C
C     FONCTION  : FIXE POUR LES POINTS DE DIRICHLET LES
C                 VALEURS AUX LIMITES SUR E 
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |     ZF1        |<-->| COTE DU FOND
C |     ZF         | -->| COTE DU FOND
C |     EBOR       |<-->| EVOLUTION AUX POINTS DE BORD
C |     LIEBOR     |<-->| TYPES DE CONDITIONS AUX LIMITES SUR E
C |     NBOR       | -->| TABLEAU DES NUMEROS GLOBAUX DES POINTS
C |                |    | DE BORD
C |     NPOIN      | -->| NOMBRE DE POINTS DU MAILLAGE
C |     NPTFR      | -->| NOMBRE DE POINTS FRONTIERES
C |     KENT       | -->| TYPE DE CONDITION LIMITE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C PROGRAMME APPELANT : RESOLU
C PROGRAMMES APPELES :
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN):: KENT,NPOIN,NPTFR
      INTEGER, INTENT(IN):: NBOR(NPTFR)
      INTEGER, INTENT(IN):: LIEBOR(NPTFR) 
C
      DOUBLE PRECISION, INTENT(IN)::  ZF(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT):: ZF1(NPOIN)
      DOUBLE PRECISION, INTENT(IN) :: EBOR(NPTFR)  
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K, N
C
C-----------------------------------------------------------------------
C
      DO 10 K=1,NPTFR
C
        N = NBOR(K)
C
        IF (LIEBOR(K).EQ.KENT) THEN
          ZF1(N)   = EBOR(K)+ZF(N)
        ENDIF
C
10    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE DIRICL

