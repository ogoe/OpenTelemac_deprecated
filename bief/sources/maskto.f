C                       *****************
                        SUBROUTINE MASKTO
C                       *****************
C
     *(MASKEL,MASKPT,IFAMAS,IKLE,IFABOR,ELTSEG,NSEG,
     * NELEM,NPOIN,IELM,MESH)
C
C***********************************************************************
C BIEF VERSION 5.9       21/10/08     J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C     FONCTION :
C
C     SOUS PROGRAMME APPELE A L'ISSUE DU REMPLISSAGE DE MASKEL
C     EFFECTUE DANS MASKBD ET MASKOB
C
C   - REMPLISSAGE DE MASKPT (UN POINT EST MASQUE SI IL EST ENTOURE
C                            D'ELEMENTS MASQUES)
C
C   - REMPLISSAGE DE IFAMAS (TABLEAU ANALOGUE A IFABOR MAIS EN
C           CONSIDERANT TOUTE FACE SEPARANT UN ELEMENT MASQUE
C           D'UN ELEMENT NON MASQUE COMME UNE FRONTIERE SOLIDE)
C
C   - POUR LES PRISMES, CALCUL DES NOUVELLES NORMALES COMPATIBLES
C           XNEBOR ET YNEBOR
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   MASKEL       | -->| TABLEAU DE MASQUAGE DES ELEMENTS             |
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE.        |
C |   MASKPT       |<-- | TABLEAU DE MASQUAGE DES POINTS               |
C |                |    |  =1. : NORMAL   =0. : POINT MASQUE.          |
C |   IFAMAS       |<-- | NUMERO DES ELEMENTS VOISINS AVEC MASQUE.     |
C |   IKLE         | -->| TABLE DE CONNECTIVITE.                       |
C |   IFABOR       | -->| NUMERO DES ELEMENTS VOISINS SANS MASQUE.     |
C |   NELEM        | -->| NOMBRE D'ELEMENTS.                           |
C |   NPOIN        | -->| NOMBRE DE POINTS.                            |
C |   IELM         | -->| TYPE D'ELEMENTS.                             |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C APPELE PAR: TELEMAC-2D, TELEMAC-3D, SISYPHE, ETC.
C***********************************************************************
C                                                                      *
C                                                                      *
C***********************************************************************
C
      USE BIEF, EX_MASKTO => MASKTO
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER, INTENT(IN)    :: NELEM,NPOIN,IELM,NSEG
      INTEGER, INTENT(IN)    :: IKLE(NELEM,3),IFABOR(NELEM,3)
      INTEGER, INTENT(IN)    :: ELTSEG(NELEM,3)
      INTEGER, INTENT(INOUT) :: IFAMAS(NELEM,3)
C
      DOUBLE PRECISION, INTENT(IN)    :: MASKEL(NELEM)
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: MASKPT
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER IELEM,N,I1,I2,I3
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: WSEG
C
      INTRINSIC MAX
C
C-----------------------------------------------------------------------
C
      WSEG => MESH%MSEG%X%R
C
      CALL OS('X=0     ',X=MASKPT)
C
      IF(IELM.EQ.11.OR.IELM.EQ.41) THEN
C
C 1) ON MASQUE LES POINTS QUI N'APPARTIENNENT A AUCUN ELEMENT NON GELE
C    (CEUX QUI APPARTIENNENT A UN ELEMENT NORMAL SONT REMIS A 1.)
C
        DO IELEM = 1,NELEM
          I1 = IKLE(IELEM,1)
          I2 = IKLE(IELEM,2)
          I3 = IKLE(IELEM,3)
          MASKPT%R(I1) = MAX(MASKPT%R(I1),MASKEL(IELEM))
          MASKPT%R(I2) = MAX(MASKPT%R(I2),MASKEL(IELEM))
          MASKPT%R(I3) = MAX(MASKPT%R(I3),MASKEL(IELEM))
        ENDDO
C
C       IN PARALLEL FOR INTERFACE POINTS, MAXIMUM RETAINED
C
        IF(NCSIZE.GT.1) CALL PARCOM(MASKPT,3,MESH)      
C
C 2) ON RECOPIE IFABOR DANS IFAMAS
C
        DO IELEM = 1,NELEM
          IFAMAS(IELEM,1) = IFABOR(IELEM,1)
          IFAMAS(IELEM,2) = IFABOR(IELEM,2)
          IFAMAS(IELEM,3) = IFABOR(IELEM,3)
        ENDDO
C
C 3) ON MARQUE LES BORDS DES ELEMENTS GELES PAR 0 (FRONTIERE LIQUIDE)
C    POUR PROVOQUER UN ARRET DES CARACTERISTIQUES
C    ON PASSE PAR UN TABLEAU DEFINI PAR SEGMENT POUR POUVOIR 
C    LE COMMUNIQUER EN PARALLELE
C
C       WSEG SET TO 1
C
        DO N=1,NSEG
          WSEG(N)=1.D0
        ENDDO
C
C       THEN WSEG PUT TO 0 FOR DRY ELEMENTS
C
        DO IELEM=1,NELEM
          IF(MASKEL(IELEM).LT.0.5D0) THEN
            WSEG(ELTSEG(IELEM,1))=0.D0
            WSEG(ELTSEG(IELEM,2))=0.D0
            WSEG(ELTSEG(IELEM,3))=0.D0
          ENDIF          
        ENDDO
C
C       IN PARALLEL FOR INTERFACE EDGES, MINIMUM RETAINED
C
        IF(NCSIZE.GT.1) THEN
          CALL PARCOM2_SEG(WSEG,WSEG,WSEG,NSEG,1,4,1,MESH,1)
        ENDIF
C
C       WSEG = 0.D0 TRANSLATED INTO IFAMAS = 0
C
        DO IELEM=1,NELEM
          IF(WSEG(ELTSEG(IELEM,1)).LT.0.5D0) IFAMAS(IELEM,1)=0
          IF(WSEG(ELTSEG(IELEM,2)).LT.0.5D0) IFAMAS(IELEM,2)=0
          IF(WSEG(ELTSEG(IELEM,3)).LT.0.5D0) IFAMAS(IELEM,3)=0
        ENDDO
C
      ELSE
C
        IF(LNG.EQ.1) WRITE(LU,1000) IELM
        IF(LNG.EQ.2) WRITE(LU,1100) IELM
1000    FORMAT(1X,'MASKTO: TYPE D''ELEMENT INCONNU :',1I6)
1100    FORMAT(1X,'MASKTO: UNKNOWN TYPE OF ELEMENT :',1I6)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
