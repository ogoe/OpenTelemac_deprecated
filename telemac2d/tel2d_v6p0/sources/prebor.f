C                       *****************
                        SUBROUTINE PREBOR
C                       *****************
C
     *(HBOR,UBOR,VBOR,TBOR,U,V,H,HN,T,NBOR,KP1BOR,NPOIN,NPTFR,
     * NTRAC,DEBLIQ,FINLIQ,NFRLIQ)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.8    03/09/07    E DAVID (LHF) 04 76 33 42 36
C
C***********************************************************************
C
C FONCTION :
C
C PREPARATION DES CONDITIONS AUX LIMITES POUR TRAITEMENT PAR THOMPSON
C UBOR,VBOR,HBOR,TBOR SONT INITIALISES ICI AUX VALEURS AU TEMPS N
C APRES BORD, CES TABLEAUX CONTIENDRONT DONC SOIT LA VALEUR AU TEMPS N
C SOIT LA VALEUR IMPOSEE
C
C CHARGEMENT DE H DANS UN TABLEAU INTERMEDIAIRE POUR SAUVEGARDER
C SA VALEUR AU BORD AU TEMPS N (MODIFIE DANS BORD).
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   HBOR         |<-- |  HAUTEUR IMPOSEE.                           
C |   UBOR         |<-- |  VITESSE U IMPOSEE.                    
C |   VBOR         |<-- |  VITESSE V IMPOSEE.        
C |   TBOR         |<-- |  TRACEUR IMPOSE AU BORD     
C |   U,V          | -->|  COMPOSANTES DE LA VITESSE AU TEMPS N 
C |   H            | -->|  HAUTEUR AU TEMPS N                  
C |   HN           | -->|  HAUTEUR DE PROPAGATION (OPTION H-U) 
C |   T            | -->|  TRACEUR AU TEMPS N      
C |   NBOR         | -->|  ADRESSES DES POINTS DE BORD 
C |   KP1BOR       | -->|  NUMERO DU POINT FRONTIERE SUIVANT 
C |   NPOIN        | -->|  NOMBRE DE POINTS DU MAILLAGE
C |   NPTFR        | -->|  NOMBRE DE POINTS FRONTIERE
C |   NTRAC        | -->|  NOMBRE DE TRACEURS
C |   DEBLIQ       | -->|  TABLEAU D'INDICES DE DEBUT DE FRONTIERE LIQ.
C |   FINLIQ       | -->|  TABLEAU D'INDICES DE FIN DE FRONTIERE LIQUI.
C |   NFRLIQ       | -->|  NOMBRE DE FRONTIERES LIQUIDES
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : TELMAC
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
      INTEGER, INTENT(IN)             :: NPOIN,NPTFR,NFRLIQ,NTRAC
      INTEGER, INTENT(IN)             :: DEBLIQ(NFRLIQ),FINLIQ(NFRLIQ)
      INTEGER, INTENT(IN)             :: NBOR(NPTFR),KP1BOR(NPTFR,2)
      DOUBLE PRECISION, INTENT(IN)    :: U(NPOIN),V(NPOIN),H(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: HBOR(NPTFR),UBOR(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: VBOR(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: HN(NPOIN)
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: T
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: TBOR
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K,IFRLIQ,ITRAC
C           
      DOUBLE PRECISION C
C
      LOGICAL DEP
C
C-----------------------------------------------------------------------
C
C  BOUCLE SUR LES FRONTIERES LIQUIDES
C
      IF(NFRLIQ.NE.0) THEN
C
      DO 10 IFRLIQ = 1 , NFRLIQ
C
        DEP = .FALSE.
        K = DEBLIQ(IFRLIQ)
11      CONTINUE
        UBOR(K)=U(NBOR(K))
        VBOR(K)=V(NBOR(K))
        HBOR(K)=H(NBOR(K))
        IF(NTRAC.GT.0) THEN
          DO ITRAC=1,NTRAC
            TBOR%ADR(ITRAC)%P%R(K)=T%ADR(ITRAC)%P%R(NBOR(K))
          ENDDO
        ENDIF
        IF(K.EQ.FINLIQ(IFRLIQ).AND.DEP) THEN
          GO TO 12
        ELSE
          DEP=.TRUE.
          K = KP1BOR(K,1)
          GO TO 11
        ENDIF
12      CONTINUE
C
10    CONTINUE
C
      ENDIF
C
      CALL OV( 'X=Y     ' , HN , H , H , C , NPOIN )
C
C-----------------------------------------------------------------------
C
      RETURN
      END
