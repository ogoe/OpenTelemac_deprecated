!                    *****************
                     SUBROUTINE FROPRO
!                    *****************
!
     &(NBOR,IKLE,NELEM,NELMAX,NPOIN,NPMAX,NPTFR,IELM,
     & IKLEM1,LIMVOI,OPTASS,PRODUC,MXPTVS,T1,
     & GLOSEG,SIZGLO,NSEG)
!
!***********************************************************************
! BIEF   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    COMPUTES THE ARRAYS GIVING ADRESSES FOR
!+                FRONTAL MATRIX-VECTOR PRODUCT.
!
!history  J-M HERVOUET (LNHE)
!+        20/03/08
!+        V5P9
!+
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| GLOSEG         |-->| FIRST AND SECOND POINT OF SEGMENTS
!| IELM           |-->| TYPE OF ELEMENT.
!|                |   | 11 : TRIANGLES.
!|                |   | 21 : QUADRANGLES.
!| IKLE           |-->| CONNECTIVITY TABLE.
!| IKLEM1         |<--| VARIOUS ADDRESSES IN OFF-DIAGONAL PART
!|                |   | OF A MATRIX
!| LIMVOI         |<--| LIMVOI(K,1) : BEGINNING OF SERIES WITH K NEIGHBOURS
!|                |   | LIMVOI(K,2) : END       OF SERIES WITH K NEIGHBOURS
!| MXPTVS         |-->| MAXIMUM NUMBER OF NEIGHBOURS OF A POINT
!| NBOR           |-->| GLOBAL NUMBER OF BOUNDARY POINTS
!| NELEM          |-->| NUMBER OF ELEMENTS
!| NELMAX         |-->| MAXIMUM NUMBER OF ELEMENTS
!| NPMAX          |-->| MAXIMUM NUMBER OF POINTS
!| NPOIN          |-->| NUMBER OF POINTS
!| NPTFR          |-->| NUMBER OF BOUNDARY POINTS
!| NSEG           |-->| NUMBER OF SEGMENTS
!| OPTASS         |-->| OPTION OF MATRIX STORAGE 
!|                |   | 1: ELEMENT PER ELEMENT 3: EDGE-BASED 
!| PRODUC         |-->| CHOICE OF MATRIX-VECTOR PRODUCT
!|                |   | 1: NORMAL 2: FRONTAL
!| SIZGLO         |-->| FIRST DIMENSION OF GLOSEG
!| T1             |-->| INTEGER WORK ARRAY
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF, EX_FROPRO => FROPRO
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)    :: NELMAX,NPMAX,MXPTVS,NELEM
      INTEGER, INTENT(IN)    :: NPOIN,NPTFR,IELM,OPTASS,PRODUC
      INTEGER, INTENT(IN)    :: NSEG,SIZGLO,NBOR(*)
      INTEGER, INTENT(IN)    :: IKLE(NELMAX,*),GLOSEG(SIZGLO,2)
      INTEGER, INTENT(OUT)   :: IKLEM1(NPMAX,MXPTVS,4,2)
      INTEGER, INTENT(OUT)   :: LIMVOI(MXPTVS,2)
      INTEGER, INTENT(OUT)   :: T1(NPOIN)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER IELEM,IPTFR,IPOIN,ISG,K,I,I1,I2,NBVOIS
!
!-----------------------------------------------------------------------
!
      IF(IELM.NE.11) THEN
        IF(LNG.EQ.1) WRITE(LU,900) IELM
        IF(LNG.EQ.2) WRITE(LU,901) IELM
900     FORMAT(1X,'FROPRO : IELM=',1I6,' TYPE D''ELEMENT INCONNU')
901     FORMAT(1X,'FROPRO: IELM=',1I6,' UNKNOWN TYPE OF ELEMENT')
        STOP
      ENDIF
!
      IF(PRODUC.EQ.2) THEN
!
!=======================================================================
! COMPUTES THE NUMBER OF NEIGHBOURING POINTS AND ELEMENTS
!=======================================================================
!
         DO 100 IPOIN = 1,NPOIN
            T1(IPOIN) = 0
100      CONTINUE
!
!        NUMBER OF ELEMENTS NEIGHBOURING A POINT
!
         DO 110 IELEM = 1,NELEM
            T1(IKLE(IELEM,1)) = T1(IKLE(IELEM,1)) + 1
            T1(IKLE(IELEM,2)) = T1(IKLE(IELEM,2)) + 1
            T1(IKLE(IELEM,3)) = T1(IKLE(IELEM,3)) + 1
110      CONTINUE
!
!        NUMBER OF POINTS NEIGHBOURING A POINT
!     =  NUMBER OF ELEMENTS NEIGHBOURING A POINT
!     +  1 ON BOUNDARIES
!
         DO 112 IPTFR = 1,NPTFR
            T1(NBOR(IPTFR)) = T1(NBOR(IPTFR)) + 1
112      CONTINUE
!
!=======================================================================
! CHECKS THAT THE RENUMBERING WAS MADE CORRECTLY IN STBTEL
! FILLS IN LIMVOI AND IKLEM1
!=======================================================================
!
         IF(T1(1).EQ.0) THEN
           IF(LNG.EQ.1) WRITE(LU,96)
           IF(LNG.EQ.2) WRITE(LU,97)
96         FORMAT(1X,'FROPRO : LE POINT 1 N''A PAS DE VOISIN')
97         FORMAT(1X,'FROPRO: POINT NUMBER 1 HAS NO NEIGHBOUR')
           STOP
         ENDIF
!
         DO 120 IPOIN = 2,NPOIN
          IF(T1(IPOIN).LT.T1(IPOIN-1)) THEN
            IF(LNG.EQ.1) WRITE(LU,98)
            IF(LNG.EQ.2) WRITE(LU,99)
98          FORMAT(1X,'FROPRO : PRODUIT FRONTAL, IL FAUT UNE',/,1X,
     &                'RENUMEROTATION DES POINTS AVEC STBTEL')
99          FORMAT(1X,'FROPRO: FRONTAL PRODUCT REQUIRES A',/,1X,
     &                'RENUMBERING OF POINTS WITH STBTEL')
            STOP
          ELSEIF(T1(IPOIN).GT.MXPTVS) THEN
            IF(LNG.EQ.1) WRITE(LU,94) IPOIN
            IF(LNG.EQ.2) WRITE(LU,95) IPOIN
94          FORMAT(1X,'FROPRO : LE POINT ',1I6,' A TROP DE VOISINS')
95          FORMAT(1X,'FROPRO: POINT ',1I6,' HAS TOO MANY NEIGHBOURS')
            STOP
          ENDIF
120      CONTINUE
!
!  BUILDS ARRAY LIMVOI
!  LIMVOI(K,1) : BEGINNING OF SERIES WITH K NEIGHBOURS
!  LIMVOI(K,2) : END       OF SERIES WITH K NEIGHBOURS
!
         DO K=1,MXPTVS
           LIMVOI(K,1) = 0
           LIMVOI(K,2) = 0
         ENDDO
!        POINT 1 IS THE BEGINNING OF SERIES WITH T1(1) NEIGHBOURS
         NBVOIS = T1(1)
         LIMVOI(NBVOIS,1) = 1
         DO I=2,NPOIN
           IF(T1(I).NE.NBVOIS) THEN
!          PREVIOUS POINT WAS AN END OF A SERIES
           LIMVOI(NBVOIS,2) = I-1
!          CURRENT POINT IS THE BEGINNING OF A SERIES
           NBVOIS = T1(I)
           LIMVOI(NBVOIS,1) = I
           ENDIF
         ENDDO
!        POINT NPOIN IS THE END OF ITS SERIES
         LIMVOI(NBVOIS,2) = NPOIN
!
!   ARRAYS FOR FRONTAL MATRIX-VECTOR PRODUCT :
!
         DO 130 IPOIN = 1,NPOIN
           T1(IPOIN) = 1
130      CONTINUE
!
         IF(OPTASS.EQ.3) THEN
!
!        IY DOES NOT DEPEND HERE ON THE DIRECT OR TRANSPOSE CHARACTER
         DO ISG = 1,NSEG
            I1 = GLOSEG(ISG,1)
            I2 = GLOSEG(ISG,2)
!
!           ANY MATRIX
!           IXM IN DIRECT PRODUCT
            IKLEM1(I1,T1(I1),1,1) = ISG
            IKLEM1(I2,T1(I2),1,1) = ISG + NSEG
!           IY IN DIRECT PRODUCT
            IKLEM1(I1,T1(I1),2,1) = I2
            IKLEM1(I2,T1(I2),2,1) = I1
!           IXM IN TRANSPOSE PRODUCT
            IKLEM1(I1,T1(I1),3,1) = ISG + NSEG
            IKLEM1(I2,T1(I2),3,1) = ISG
!           IY IN TRANSPOSE PRODUCT
            IKLEM1(I1,T1(I1),4,1) = I2
            IKLEM1(I2,T1(I2),4,1) = I1
!
!           SYMMETRICAL MATRIX
!           IXM IN DIRECT PRODUCT
            IKLEM1(I1,T1(I1),1,2) = ISG
            IKLEM1(I2,T1(I2),1,2) = ISG
!           IY IN DIRECT PRODUCT
            IKLEM1(I1,T1(I1),2,2) = I2
            IKLEM1(I2,T1(I2),2,2) = I1
!           IXM IN TRANSPOSE PRODUCT
            IKLEM1(I1,T1(I1),3,2) = ISG
            IKLEM1(I2,T1(I2),3,2) = ISG
!           IY IN TRANSPOSE PRODUCT
            IKLEM1(I1,T1(I1),4,2) = I2
            IKLEM1(I2,T1(I2),4,2) = I1
!
!           UPDATES THE NUMBER OF NEIGHBOURS
            T1(I1) = T1(I1) + 1
            T1(I2) = T1(I2) + 1
         ENDDO
!
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
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
