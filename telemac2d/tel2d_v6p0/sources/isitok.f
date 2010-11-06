C                       *****************
                        SUBROUTINE ISITOK
C                       *****************
C
     *(H,NPH,U,NPU,V,NPV,NTRAC,T,NPT,X,Y,BORNES,ARRET)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.8    05/09/07  J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C     FONCTION  : CONTROLE DE VRAISEMBLANCE SUR LES GRANDEURS PHYSIQUES
C
C     NOTE : ARRET N'EST PAS INITIALISE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    H,NPH       | -->| HAUTEUR ET NOMBRE DE POINTS DE HAUTEUR.
C |    U,NPU       | -->| VITESSE U ET NOMBRE DE POINTS DE VITESSE U.
C |    V,NPV       | -->| VITESSE V ET NOMBRE DE POINTS DE VITESSE V.
C |    NTRAC       | -->| NOMBRE DE TRACEURS.
C |    T,NPT       | -->| TRACEUR ET NOMBRE DE POINTS DE TRACEUR.
C |    X,Y         | -->| COORDONNEES
C |    BORNES      | -->| VALEURS LIMITES DES TABLEAUX H,U,V,T
C |                |    | DANS L'ORDRE SUIVANT : HMIN,HMAX,UMIN,UMAX,...
C |    ARRET       |<-- | LOGIQUE MIS A TRUE SI BORNES SONT DEPASSEES
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
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)          :: NPH,NPU,NPV,NPT,NTRAC
      LOGICAL, INTENT(INOUT)       :: ARRET
      DOUBLE PRECISION, INTENT(IN) :: H(NPH),U(NPU),V(NPV)
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*),BORNES(8)
      TYPE(BIEF_OBJ)  , INTENT(IN) :: T
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,ITRAC
C
C-----------------------------------------------------------------------
C
C  CONTROLES DE LA HAUTEUR
C
      DO 10 I = 1 , NPH
        IF(H(I).LT.BORNES(1)) THEN
          ARRET = .TRUE.
          IF(LNG.EQ.1) THEN
           WRITE(LU,100) 'INFERIEURE','H',I,X(I),Y(I),'H',H(I),BORNES(1)
          ENDIF
          IF(LNG.EQ.2) THEN
           WRITE(LU,101) 'LOWER','H',I,X(I),Y(I),'H',H(I),BORNES(1)
          ENDIF
        ENDIF
        IF(H(I).GT.BORNES(2)) THEN
          ARRET = .TRUE.
          IF(LNG.EQ.1) THEN
           WRITE(LU,100) 'SUPERIEURE','H',I,X(I),Y(I),'H',H(I),BORNES(2)
          ENDIF
          IF(LNG.EQ.2) THEN
           WRITE(LU,101) 'UPPER','H',I,X(I),Y(I),'H',H(I),BORNES(2)
          ENDIF
        ENDIF
10    CONTINUE
C
C-----------------------------------------------------------------------
C
C  CONTROLES DE LA VITESSE U
C
      DO 20 I = 1 , NPU
        IF(U(I).LT.BORNES(3)) THEN
          ARRET = .TRUE.
          IF(LNG.EQ.1) THEN
           WRITE(LU,100) 'INFERIEURE','U',I,X(I),Y(I),'U',U(I),BORNES(3)
          ENDIF
          IF(LNG.EQ.2) THEN
           WRITE(LU,101) 'LOWER','U',I,X(I),Y(I),'U',U(I),BORNES(3)
          ENDIF
        ENDIF
        IF(U(I).GT.BORNES(4)) THEN
          ARRET = .TRUE.
          IF(LNG.EQ.1) THEN
           WRITE(LU,100) 'SUPERIEURE','U',I,X(I),Y(I),'U',U(I),BORNES(4)
          ENDIF
          IF(LNG.EQ.2) THEN
           WRITE(LU,101) 'UPPER','U',I,X(I),Y(I),'U',U(I),BORNES(4)
          ENDIF
        ENDIF
20    CONTINUE
C
C-----------------------------------------------------------------------
C
C  CONTROLES DE LA VITESSE V
C
      DO 30 I = 1 , NPV
        IF(V(I).LT.BORNES(5)) THEN
          ARRET = .TRUE.
          IF(LNG.EQ.1) THEN
           WRITE(LU,100) 'INFERIEURE','V',I,X(I),Y(I),'V',V(I),BORNES(5)
          ENDIF
          IF(LNG.EQ.2) THEN
           WRITE(LU,101) 'LOWER','V',I,X(I),Y(I),'V',V(I),BORNES(5)
          ENDIF
        ENDIF
        IF(V(I).GT.BORNES(6)) THEN
          ARRET = .TRUE.
          IF(LNG.EQ.1) THEN
           WRITE(LU,100) 'SUPERIEURE','V',I,X(I),Y(I),'V',V(I),BORNES(6)
          ENDIF
          IF(LNG.EQ.2) THEN
           WRITE(LU,101) 'UPPER','V',I,X(I),Y(I),'V',V(I),BORNES(6)
          ENDIF
        ENDIF
30    CONTINUE
C
C-----------------------------------------------------------------------
C
C  CONTROLES DU TRACEUR
C
      IF(NTRAC.GT.0) THEN
C
      DO ITRAC=1,NTRAC
C
      DO I = 1 , NPT
        IF(T%ADR(ITRAC)%P%R(I).LT.BORNES(7)) THEN
          ARRET = .TRUE.
          IF(LNG.EQ.1) THEN
           WRITE(LU,100) 'INFERIEURE','T',I,X(I),Y(I),
     *                   'T',T%ADR(ITRAC)%P%R(I),BORNES(7)
          ENDIF
          IF(LNG.EQ.2) THEN
           WRITE(LU,101) 'LOWER','T',I,X(I),Y(I),
     *                   'T',T%ADR(ITRAC)%P%R(I),BORNES(7)
          ENDIF
        ENDIF
        IF(T%ADR(ITRAC)%P%R(I).GT.BORNES(8)) THEN
          ARRET = .TRUE.
          IF(LNG.EQ.1) THEN
           WRITE(LU,100) 'SUPERIEURE','T',I,X(I),Y(I),
     *                   'T',T%ADR(ITRAC)%P%R(I),BORNES(8)
          ENDIF
          IF(LNG.EQ.2) THEN
           WRITE(LU,101) 'UPPER','T',I,X(I),Y(I),
     *                   'T',T%ADR(ITRAC)%P%R(I),BORNES(8)
          ENDIF
        ENDIF
      ENDDO
C
      ENDDO
C
      ENDIF
C
C-----------------------------------------------------------------------
C
100   FORMAT(/,1X,'LIMITE ',A10,' SUR ',A1,' ATTEINTE AU POINT ',I6,/,
     *         1X,'DE COORDONNEES ',G16.7,' ET ',G16.7,/,1X,
     *         A1,' VAUT : ',G16.7,' ET LA LIMITE EST :',G16.7)
101   FORMAT(/,1X,A5,' LIMIT ON ',A1,' REACHED AT POINT ',I6,/,1X,
     *         'WITH COORDINATES',G16.7,' AND ',G16.7,/,1X,
     *         'THE VALUE OF ',A1,' IS ',G16.7,' THE LIMIT IS: ',G16.7)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
