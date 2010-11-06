C                       *****************
                        SUBROUTINE FROPRO
C                       *****************
C
     *(NBOR,IKLE,NELEM,NELMAX,NPOIN,NPMAX,NPTFR,IELM,
     * IKLEM1,LIMVOI,OPTASS,PRODUC,MXPTVS,T1,
     * GLOSEG,SIZGLO,NSEG)
C
C***********************************************************************
C BIEF VERSION 5.9        20/03/08    J-M HERVOUET (LNHE) 01 30 87 80 18
C COPYRIGHT 1999
C***********************************************************************
C
C    FONCTION : COMPUTATION OF ARRAYS GIVING ADRESSES FOR 
C               FRONTAL MATRIX-VECTOR PRODUCT
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    NBOR        | -->| NUMERO GLOBAL DU POINT DE BORD K.
C |    IKLE        | -->| NUMEROS GLOBAUX DES POINTS DE CHAQUE ELEMENT.
C |    NELEM       | -->| NOMBRE TOTAL D'ELEMENTS DANS LE MAILLAGE.
C |    NPOIN       | -->| NOMBRE TOTAL DE POINTS DU DOMAINE.
C |    T1          | -->| TABLEAUX DE TRAVAIL ENTIERS.
C |    NPTFR       | -->| NOMBRE DE POINTS FRONTIERES.
C |    IELM        | -->| TYPE D'ELEMENT.
C |                |    | 11 : TRIANGLES.
C |                |    | 21 : QUADRILATERES.
C |    MXPTVS      | -->| NOMBRE MAXIMUM DE VOISINS D'UN POINT
C |________________|____|______________________________________________|
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : 
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      USE BIEF, EX_FROPRO => FROPRO
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NELMAX,NPMAX,MXPTVS,NELEM
      INTEGER, INTENT(IN)    :: NPOIN,NPTFR,IELM,OPTASS,PRODUC
      INTEGER, INTENT(IN)    :: NSEG,SIZGLO,NBOR(*)
      INTEGER, INTENT(IN)    :: IKLE(NELMAX,*),GLOSEG(SIZGLO,2)
      INTEGER, INTENT(OUT)   :: IKLEM1(NPMAX,MXPTVS,4,2)
      INTEGER, INTENT(OUT)   :: LIMVOI(MXPTVS,2)
      INTEGER, INTENT(OUT)   :: T1(NPOIN) 
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IPTFR,IPOIN,ISG,K,I,I1,I2,NBVOIS
C
C-----------------------------------------------------------------------
C
      IF(IELM.NE.11) THEN
        IF(LNG.EQ.1) WRITE(LU,900) IELM
        IF(LNG.EQ.2) WRITE(LU,901) IELM
900     FORMAT(1X,'FROPRO : IELM=',1I6,' TYPE D''ELEMENT INCONNU')
901     FORMAT(1X,'FROPRO: IELM=',1I6,' UNKNOWN TYPE OF ELEMENT')
        STOP
      ENDIF
C
      IF(PRODUC.EQ.2) THEN
C
C=======================================================================
C CALCUL DU NOMBRE DE POINTS ET ELEMENTS VOISINS
C=======================================================================
C
         DO 100 IPOIN = 1,NPOIN
            T1(IPOIN) = 0
100      CONTINUE
C
C        NUMBER OF ELEMENTS NEIGHBOUR OF A POINT
C
         DO 110 IELEM = 1,NELEM
            T1(IKLE(IELEM,1)) = T1(IKLE(IELEM,1)) + 1
            T1(IKLE(IELEM,2)) = T1(IKLE(IELEM,2)) + 1
            T1(IKLE(IELEM,3)) = T1(IKLE(IELEM,3)) + 1
110      CONTINUE
C
C        NUMBER OF POINTS NEIGHBOUR OF A POINT
C     =  NUMBER OF ELEMENTS NEIGHBOUR OF A POINT
C     +  1 ON BOUNDARIES
C
         DO 112 IPTFR = 1,NPTFR
            T1(NBOR(IPTFR)) = T1(NBOR(IPTFR)) + 1
112      CONTINUE
C
C=======================================================================
C VERIFICATION QUE LA RENUMEROTATION A BIEN ETE FAITE DANS STBTEL
C REMPLISSAGE DE LIMVOI ET DE IKLEM1
C=======================================================================
C
         IF(T1(1).EQ.0) THEN
           IF(LNG.EQ.1) WRITE(LU,96)
           IF(LNG.EQ.2) WRITE(LU,97)
96         FORMAT(1X,'FROPRO : LE POINT 1 N''A PAS DE VOISIN')
97         FORMAT(1X,'FROPRO: POINT NUMBER 1 HAS NO NEIGHBOUR')
           STOP
         ENDIF
C
         DO 120 IPOIN = 2,NPOIN
          IF(T1(IPOIN).LT.T1(IPOIN-1)) THEN
            IF(LNG.EQ.1) WRITE(LU,98)
            IF(LNG.EQ.2) WRITE(LU,99)
98          FORMAT(1X,'FROPRO : PRODUIT FRONTAL, IL FAUT UNE',/,1X,
     *                'RENUMEROTATION DES POINTS AVEC STBTEL')
99          FORMAT(1X,'FROPRO: FRONTAL PRODUCT REQUIRES A',/,1X,
     *                'RENUMBERING OF POINTS WITH STBTEL')
            STOP
          ELSEIF(T1(IPOIN).GT.MXPTVS) THEN
            IF(LNG.EQ.1) WRITE(LU,94) IPOIN
            IF(LNG.EQ.2) WRITE(LU,95) IPOIN
94          FORMAT(1X,'FROPRO : LE POINT ',1I6,' A TROP DE VOISINS')
95          FORMAT(1X,'FROPRO: POINT ',1I6,' HAS TOO MANY NEIGHBOURS')
            STOP
          ENDIF
120      CONTINUE
C
C  BUILDING ARRAY LIMVOI
C  LIMVOI(K,1) : BEGINNING OF SERIES WITH K NEIGHBOURS
C  LIMVOI(K,2) : END       OF SERIES WITH K NEIGHBOURS
C
         DO K=1,MXPTVS
           LIMVOI(K,1) = 0
           LIMVOI(K,2) = 0           
         ENDDO
C        POINT 1 IS THE BEGINNING OF SERIES WITH T1(1) NEIGHBOURS
         NBVOIS = T1(1)
         LIMVOI(NBVOIS,1) = 1
         DO I=2,NPOIN
           IF(T1(I).NE.NBVOIS) THEN
C          PREVIOUS POINT WAS AN END OF A SERIES
           LIMVOI(NBVOIS,2) = I-1
C          CURRENT POINT IS THE BEGINNING OF A SERIES
           NBVOIS = T1(I)
           LIMVOI(NBVOIS,1) = I
           ENDIF
         ENDDO
C        POINT NPOIN IS THE END OF ITS SERIES
         LIMVOI(NBVOIS,2) = NPOIN           
C
C   ARRAYS FOR FRONTAL MATRIX-VECTOR PRODUCT :
C
         DO 130 IPOIN = 1,NPOIN
           T1(IPOIN) = 1
130      CONTINUE
C
         IF(OPTASS.EQ.3) THEN
C
C        ICI IY NE DEPEND PAS DU CARACTERE DIRECT OU TRANSPOSE
         DO ISG = 1,NSEG
            I1 = GLOSEG(ISG,1)
            I2 = GLOSEG(ISG,2)
C
C           MATRICE QUELCONQUE
C           IXM EN PRODUIT DIRECT
            IKLEM1(I1,T1(I1),1,1) = ISG
            IKLEM1(I2,T1(I2),1,1) = ISG + NSEG
C           IY EN PRODUIT DIRECT
            IKLEM1(I1,T1(I1),2,1) = I2
            IKLEM1(I2,T1(I2),2,1) = I1
C           IXM EN PRODUIT TRANSPOSE
            IKLEM1(I1,T1(I1),3,1) = ISG + NSEG
            IKLEM1(I2,T1(I2),3,1) = ISG
C           IY EN PRODUIT TRANSPOSE
            IKLEM1(I1,T1(I1),4,1) = I2
            IKLEM1(I2,T1(I2),4,1) = I1
C
C           MATRICE SYMETRIQUE
C           IXM EN PRODUIT DIRECT
            IKLEM1(I1,T1(I1),1,2) = ISG
            IKLEM1(I2,T1(I2),1,2) = ISG 
C           IY EN PRODUIT DIRECT
            IKLEM1(I1,T1(I1),2,2) = I2
            IKLEM1(I2,T1(I2),2,2) = I1
C           IXM EN PRODUIT TRANSPOSE
            IKLEM1(I1,T1(I1),3,2) = ISG
            IKLEM1(I2,T1(I2),3,2) = ISG
C           IY EN PRODUIT TRANSPOSE
            IKLEM1(I1,T1(I1),4,2) = I2
            IKLEM1(I2,T1(I2),4,2) = I1
C
C           MISE A JOUR DU NOMBRE DE VOISINS
            T1(I1) = T1(I1) + 1
            T1(I2) = T1(I2) + 1
         ENDDO
C
         ELSE
            IF(LNG.EQ.1) THEN 
              WRITE(LU,*) 'STOCKAGE INCONNU DANS FROPRO :',OPTASS
            ENDIF
            IF(LNG.EQ.2) THEN 
              WRITE(LU,*) 'UNKNOWN STORAGE IN FROPRO :',OPTASS
            ENDIF
            CALL PLANTE(1)
            STOP
         ENDIF
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
