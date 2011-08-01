C                       *******************
                        SUBROUTINE COMP_FAC
C                       *******************
C
     *(ELTSEG,IFABOR,NELEM,NPOIN,FAC)
C
C***********************************************************************
C BIEF VERSION 5.9        24/10/08    J-M HERVOUET (LNHE) 01 30 87 80 18
C COPYRIGHT EDF 2008                 
C                                     
C***********************************************************************
C
C    FONCTION : COMPLETING THE ARRAY FAC FOR QUADRATIC POINTS
C               AT INTERFACE BETWEEN 2 SUBDOMAINS
C
C               FAC%R(1:NPOIN) IS FILLED IN PARINI
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  ELTSEG        | -->| GIVES THE SEGMENT NUMBER OF EDGES OF ELEMENTS
C |  IFABOR        | -->| -2 MEANS INTERFACE WITH ANOTHER SUB-DOMAIN
C |  NELEM         | -->| NUMBER OF ELEMENTS
C |  NPOIN         | -->| NUMBER OF POINTS
C |  FAC           |<-->| COEFFICIENT FOR COMPUTING DOT PRODUCTS IN //
C |________________|____|______________________________________________|
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : INBIEF
C
C SOUS-PROGRAMME APPELE :
C
C***********************************************************************
C 
      USE BIEF, EX_COMP_FAC => COMP_FAC
C    
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NELEM,NPOIN
      INTEGER, INTENT(IN)    :: IFABOR(NELEM,3),ELTSEG(NELEM,3)
      TYPE(BIEF_OBJ), INTENT(INOUT) :: FAC
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM
C
C-----------------------------------------------------------------------
C
      DO IELEM=1,NELEM
C
        IF(IFABOR(IELEM,1).EQ.-2) FAC%R(NPOIN+ELTSEG(IELEM,1))=0.5D0
        IF(IFABOR(IELEM,2).EQ.-2) FAC%R(NPOIN+ELTSEG(IELEM,2))=0.5D0
        IF(IFABOR(IELEM,3).EQ.-2) FAC%R(NPOIN+ELTSEG(IELEM,3))=0.5D0
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
