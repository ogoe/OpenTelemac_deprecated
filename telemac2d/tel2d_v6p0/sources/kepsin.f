C                       *****************
                        SUBROUTINE KEPSIN
C                       *****************
C
     *(LIMKEP,LIUBOR,NPTFR,
     * KENT,KENTU,KSORT,KADH,KLOG,KINC,KNEU,KDIR)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.2    17/08/94    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C      FONCTION:
C      ========:
C
C     CE SOUS-PROGRAMME INITIALISE LES CONDITIONS AUX LIMITES POUR
C     L'ETAPE DE DIFFUSION-TERMES SOURCES DU MODELE K-EPSILON
C
C     ATTENTION: LIMKEP EST CONSTRUIT A PARTIR DE LIUBOR
C                  (LIKBOR ET LIEBOR N'EXISTENT PAS)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   LIMKEP(*,1)  |<-- | TYPES DE CONDITIONS AUX LIMITES POUR AK      |
C |   LIMKEP(*,2)  |<-- | TYPES DE CONDITIONS AUX LIMITES POUR EP      |
C |   LIUBOR       | -->| TYPES DE CONDITIONS AUX LIMITES SUR U.       |
C |   NPTFR        | -->| DIMENSION DES TABLEAUX .                     |
C |   KENT         | -->| INDICATEUR DE POINT D'ENTREE FLUIDE .        |
C |   KENTU        | -->| INDICATEUR DE VITESSES IMPOSEES.             |
C |   KSORT        | -->| INDICATEUR DE POINT DE SORTIE FLUIDE .       |
C |   KADH         | -->| INDICATEUR DE POINT DIRICHLET.               |
C |   KLOG         | -->| INDICATEUR DE PAROI SOLIDE .                 |
C |   KINC         | -->| INDICATEUR POUR L'ONDE INCIDENTE.            |
C |   KNEU         | -->| INDICATEUR DE POINT DE NEUMANN               |
C |   KDIR         | -->| INDICATEUR DE POINT DE DIRICHLET             |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
C
C SOUS-PROGRAMME APPELE : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NPTFR,KENT,KSORT,KADH,KLOG
      INTEGER, INTENT(IN)    :: KINC,KNEU,KDIR,KENTU
      INTEGER, INTENT(INOUT) :: LIMKEP(NPTFR,2)
      INTEGER, INTENT(IN)    :: LIUBOR(NPTFR)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K      
C
C-----------------------------------------------------------------------
C
C  CONSTRUCTION DU TABLEAU LIMKEP
C
      DO 1 K=1,NPTFR
C
        IF(LIUBOR(K).EQ.KENT) THEN
C
          LIMKEP(K,1) = KDIR
          LIMKEP(K,2) = KDIR
C
        ELSEIF(LIUBOR(K).EQ.KENTU) THEN
C
          LIMKEP(K,1) = KDIR
          LIMKEP(K,2) = KDIR
C
        ELSEIF(LIUBOR(K).EQ.KADH) THEN
C
          LIMKEP(K,1) = KDIR
          LIMKEP(K,2) = KDIR
C
        ELSEIF(LIUBOR(K).EQ.KSORT) THEN
C
          LIMKEP(K,1) = KNEU
          LIMKEP(K,2) = KNEU
C
        ELSEIF(LIUBOR(K).EQ.KINC) THEN
C
          LIMKEP(K,1) = KNEU
          LIMKEP(K,2) = KNEU
C
        ELSEIF(LIUBOR(K).EQ.KLOG ) THEN
C
          LIMKEP(K,1) = KDIR
          LIMKEP(K,2) = KDIR
C
        ELSE
C
          IF(LNG.EQ.1) WRITE(LU,100) K,LIUBOR(K)
          IF(LNG.EQ.2) WRITE(LU,101) K,LIUBOR(K)
100       FORMAT(1X,'KEPSIN: K=',1I6,' LIUBOR=',1I6,' CAS NON PREVU')
101       FORMAT(1X,'KEPSIN: K=',1I6,' LIUBOR=',1I6,' UNKNOWN CASE')
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
1     CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
