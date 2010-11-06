C                       *****************
                        SUBROUTINE NORMAB
C                       *****************
C
     *(XNEBOR,YNEBOR,XSGBOR,YSGBOR,DISBOR,SURFAC,NELEM,
     * NBOR,KP1BOR,NELBOR,LGSEG,NPTFR,X,Y,MESH,T1)
C
C***********************************************************************
C BIEF VERSION 5.9      26/06/2008    J-M HERVOUET (LNHE) 01 30 87 80 18
C                         
C***********************************************************************
C
C     FONCTION  : 1) CALCUL DES COMPOSANTES DU VECTEUR NORMAL SORTANT
C
C                    - POUR LES POINTS DE BORD      (XNEBOR,YNEBOR)
C                    - POUR LES SEGMENTS DE BORD    (XSGBOR,YSGBOR)
C
C                    
C     ATTENTION :  XSGBOR ET YSGBOR SONT DE DIMENSION (NPTFR,4) :
C
C                    (K,1) : VERSION NORMEE    , SEGMENT SUIVANT K
C                    (K,2) : VERSION NORMEE    , SEGMENT PRECEDENT K
C                    (K,3) : VERSION NON NORMEE, SEGMENT SUIVANT K
C                    (K,4) : VERSION NON NORMEE, SEGMENT PRECEDENT K
C
C                    XSGBOR(K,1) ET YSGBOR(K,1) SONT LES COMPOSANTES
C                    DU SEGMENT SUIVANT LE POINT K.
C
C                    XSGBOR(K,2) ET YSGBOR(K,2) SONT LES COMPOSANTES
C                    DU SEGMENT PRECEDENT LE POINT K.
C
C                 2) DISTANCE AU BORD DES PREMIERS POINTS DE L'ELEMENT.
C
C                 3) LONGUEUR DES SEGMENTS DE BORD.
C
C                 4) DISTANCE DU BORD AUX PREMIERS POINTS INTERIEURS
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   XNEBOR       |<-- | COMPOSANTE SUIVANT X DE LA NORMALE AU POINT. |
C |   YNEBOR       |<-- | COMPOSANTE SUIVANT Y DE LA NORMALE AU POINT. |
C |   XSGBOR       |<-- | COMPOSANTE SUIVANT X DE LA NORMALE DU SEGMENT|
C |   YSGBOR       |<-- | COMPOSANTE SUIVANT Y DE LA NORMALE DU SEGMENT|
C |                |    | (SEGMENT COMPRIS ENTRE LES POINTS K ET K+1)  |
C |   DISBOR       |<-- | DISTANCE DU BORD AU PREMIER POINT INTERIEUR. |
C |   SURFAC       | -->| SURFACE DES TRIANGLES                        |
C |   NELEM        | -->| NOMBRE D'ELEMENTS                            |
C |   NBOR         | -->| ADRESSES DES POINTS FRONTIERES.              |
C |   KP1BOR       | -->| NUMERO DE L'EXTREMITE DES SEGMENTS DE BORD   |
C |   NELBOR       | -->| NUMERO DE L'ELEMENT DE BORD                  |
C |   LGSEG        |<-- | LONGUEUR DES SEGMENTS DE BORD                |
C |   NPTFR        | -->| NOMBRE DE POINTS FRONTIERES                  |
C |   X,Y          | -->| COORDONNEES X,Y DONNEES PAR POINT.           |
C |   MESH         |<-->| STRUCTURE DE MAILLAGE                        |
C |   T1           |<-->| STRUCTURE BIEF_OBJ POUR TABLEAU DE TRAVAIL   |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PREDON
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      USE BIEF, EX_NORMAB => NORMAB
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NPTFR,NELEM
      INTEGER, INTENT(IN) :: NBOR(NPTFR),KP1BOR(NPTFR),NELBOR(NPTFR)
C
      DOUBLE PRECISION, INTENT(INOUT) :: XNEBOR(NPTFR,2),YNEBOR(NPTFR,2)
      DOUBLE PRECISION, INTENT(INOUT) :: XSGBOR(NPTFR,4),YSGBOR(NPTFR,4)
      DOUBLE PRECISION, INTENT(INOUT) :: DISBOR(NPTFR),LGSEG(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: SURFAC(NELEM),X(*),Y(*)
C
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: T1
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K1,K2,N1,N2,I,IELEM
      DOUBLE PRECISION X12,Y12,XNORM,X1,X2,Y1,Y2,Z(1)
C
      INTRINSIC SQRT
C
C-----------------------------------------------------------------------
C
C     CALCUL DES VECTEURS NORMAUX ET DES LONGUEURS DES SEGMENTS
C
C     0) INITIALISATION DE LGSEG, XSGBOR ET YSGBOR A ZERO
C
      IF(NPTFR.GT.0) THEN
C
        DO K1=1,NPTFR 
          LGSEG(K1)    = 0.D0      
          XSGBOR(K1,1) = 0.D0
          YSGBOR(K1,1) = 0.D0
          XSGBOR(K1,2) = 0.D0
          YSGBOR(K1,2) = 0.D0  
        ENDDO
C
C       1) NORMALES PAR SEGMENT ET LONGUEUR DU SEGMENT DE BORD
C          VERSION COMMUNE SCALAIRE/PARALLELE
C
        DO K1=1,NPTFR
C       
          K2=KP1BOR(K1)
          IF(K2.NE.K1) THEN
            N1=NBOR(K1)
            N2=NBOR(K2)
            X1 = X(N1)
            Y1 = Y(N1)
            X2 = X(N2)
            Y2 = Y(N2)
            X12 = X2 - X1
            Y12 = Y2 - Y1
C           LONGUEUR DU SEGMENT DE BORD
            LGSEG(K1) = SQRT( X12**2 + Y12**2 )
C           NORMALE AU SEGMENT SUIVANT K1 :
            XSGBOR(K1,1) =  Y12
            YSGBOR(K1,1) = -X12
C           NORMALE AU SEGMENT PRECEDENT LE SUIVANT DE K1
            XSGBOR(K2,2) =  Y12
            YSGBOR(K2,2) = -X12   
          ENDIF
C
        ENDDO
C
      ENDIF
C
C     2) COMPLEMENT EN PARALLELE, AVEC PARCOM OPTION 1
C        (VALEUR DE PLUS GRANDE VALEUR ABSOLUE)
C
      IF(NCSIZE.GT.1) THEN
        IF(NPTFR.GT.0) THEN
          CALL PARCOM_BORD(LGSEG            ,1,MESH)
          CALL PARCOM_BORD(XSGBOR(1:NPTFR,1),1,MESH)
          CALL PARCOM_BORD(XSGBOR(1:NPTFR,2),1,MESH)
          CALL PARCOM_BORD(YSGBOR(1:NPTFR,1),1,MESH)
          CALL PARCOM_BORD(YSGBOR(1:NPTFR,2),1,MESH)
        ELSE
          CALL PARCOM_BORD(Z,1,MESH)
          CALL PARCOM_BORD(Z,1,MESH)
          CALL PARCOM_BORD(Z,1,MESH)
          CALL PARCOM_BORD(Z,1,MESH)
          CALL PARCOM_BORD(Z,1,MESH)
        ENDIF
      ENDIF      
C
C     3) NORMALES PAR NOEUDS, DISTANCE APPROXIMATIVE DES BORDS
C        PUIS ON NORME LES VECTEURS.
C
      IF(NPTFR.GT.0) THEN
C
      DO K1=1,NPTFR
C
C       NORMALE AU POINT : MOYENNE DES NORMALES NON NORMEES AUX
C       DEUX SEGMENTS ADJACENTS
C
C       VERSION NON NORMEE XNEBOR(*,2) ET YNEBOR(*,2)
        XNEBOR(K1,2)=(XSGBOR(K1,1)+XSGBOR(K1,2))*0.5D0
        YNEBOR(K1,2)=(YSGBOR(K1,1)+YSGBOR(K1,2))*0.5D0
C
C       VERSION NON NORMEE XSGBOR(*,3) ET XSGBOR(*,4)
C                          YSGBOR(*,3) ET YSGBOR(*,4)
        XSGBOR(K1,3)=XSGBOR(K1,1)
        XSGBOR(K1,4)=XSGBOR(K1,2)
        YSGBOR(K1,3)=YSGBOR(K1,1)
        YSGBOR(K1,4)=YSGBOR(K1,2)
C
C       VERSION NORMEE XNEBOR(*,1) ET YNEBOR(*,1)
        XNORM=SQRT(XNEBOR(K1,2)**2+YNEBOR(K1,2)**2)
        XNEBOR(K1,1)=XNEBOR(K1,2)/XNORM
        YNEBOR(K1,1)=YNEBOR(K1,2)/XNORM
C
C       VERSION NORMEE DE XSGBOR ET YSGBOR POUR SEGMENT SUIVANT
        XNORM=SQRT(XSGBOR(K1,1)**2+YSGBOR(K1,1)**2)
        XSGBOR(K1,1)=XSGBOR(K1,1)/XNORM
        YSGBOR(K1,1)=YSGBOR(K1,1)/XNORM
C
C       VERSION NORMEE DE XSGBOR ET YSGBOR POUR SEGMENT PRECEDENT
        XNORM=SQRT(XSGBOR(K1,2)**2+YSGBOR(K1,2)**2)
        XSGBOR(K1,2)=XSGBOR(K1,2)/XNORM
        YSGBOR(K1,2)=YSGBOR(K1,2)/XNORM
C
C       THIS CAN BE APPROXIMATION OF THE MESH SIZE AT THE BOUNDARY
C       AND IS USED FOR LOG LAW AT THE BOUNDARIES
        IELEM=NELBOR(K1)
        IF(IELEM.GT.0) THEN
          DISBOR(K1) = 2.D0*SURFAC(NELBOR(K1))/LGSEG(K1)
        ELSE
          DISBOR(K1) = 0.D0
        ENDIF
C
      ENDDO
C
      ENDIF
C
C     DISBOR EST POSITIF, ON PEUT PRENDRE LE MAX
      IF(NCSIZE.GT.1) THEN
        CALL PARCOM_BORD(DISBOR,3,MESH)
      ENDIF   
C 
C-----------------------------------------------------------------------
C
      RETURN
      END
