C                       *****************
                        SUBROUTINE ININDS
C                       *****************
C
     *(NPOIN,NPTFR,NELEM,NPMAX,NPTFX,NELMAX,NPLAN,NSEGBOR,NELEB)
C
C***********************************************************************
C BIEF VERSION 5.9        24/10/2008   J-M HERVOUET (LNH) 01 30 87 80 18
C                                      REGINA NEBAUER
C                                      LAM MINH PHUONG
C***********************************************************************
C
C
C FONCTION : 1) INITIALISATION DU COMMON NODES QUI EST UTILISE PAR LES
C               FONCTIONS NBPTS, NBMPTS, NBFEL, ET NBPEL.
C
C
C
C INDEX 1 OF NDS : NUMBER OF DEGREES OF FREEDOM IN THE MESH FOR THIS ELEMENT
C INDEX 2 OF NDS : NUMBER OF SEGMENTS IN THE MESH FOR THIS ELEMENT
C INDEX 3 OF NDS : NUMBER OF POINTS PER ELEMENT
C INDEX 4 OF NDS : NUMBER OF FACES PER ELEMENT (ex 3 FOR A TRIANGLE) 
C INDEX 5 OF NDS : MAXIMUM NUMBER OF DEGREES OF FREEDOM IN THE MESH
C                  FOR THIS ELEMENT
C INDEX 6 OF NDS : NUMBER OF SEGMENTS FOR THIS ELEMENT
C
C FUNCTIONS IN BIEF:
C
C INDEX 1 OF NDS : NBPTS
C INDEX 2 OF NDS : NBSEG
C INDEX 3 OF NDS : NBPEL
C INDEX 4 OF NDS : NBFEL
C INDEX 5 OF NDS : NBMPTS
C INDEX 6 OF NDS : NBSEGEL
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  NPOIN         | -->| NOMBRE DE SOMMETS DU MAILLAGE 2D
C |                |    | (SAUF TETRAEDRES OU C'EST LE NOMBRE TOTAL)
C |  NPTFR         | -->| NOMBRE DE SOMMETS DE FRONTIERE
C |  NELEM         | -->| NOMBRE D'ELEMENTS DU MAILLAGE
C |  NPMAX         | -->| NOMBRE MAXIMUM DE SOMMETS DU MAILLAGE
C |  NPTFX         | -->| NOMBRE MAXIMUM DE SOMMETS DE FRONTIERE
C |  NELMAX        | -->| NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |  NELEB         | -->| NOMBRE D'ELEMENTS DE BORD 3D NON STRUCTURE
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER , INTENT(IN) :: NPOIN
      INTEGER , INTENT(IN) :: NPTFR
      INTEGER , INTENT(IN) :: NELEM
      INTEGER , INTENT(IN) :: NPMAX
      INTEGER , INTENT(IN) :: NPTFX
      INTEGER , INTENT(IN) :: NELMAX
      INTEGER , INTENT(IN) :: NPLAN
      INTEGER , INTENT(IN) :: NSEGBOR
      INTEGER , INTENT(IN), OPTIONAL :: NELEB
C
      INTEGER NDS(0:81,7)
C
      COMMON/NODES/NDS
C
C-----------------------------------------------------------------------
C
C     INITIALISATION DU COMMON NODES
C
C     LES VALEURS QUI NE DOIVENT PAS SERVIR NE SONT PAS INITIALISEES
C
C     0) ELEMENT P0 EN DIMENSION 1
C
C                 NSEGBOR (ICI DIMENSIONNEMENT PAR EXCES)
      NDS(00,1) = NPTFR
C     NDS(00,2) = ?????
      NDS(00,3) = 1
C     NDS(00,4) = ?????
      NDS(00,5) = NPTFX
C     NDS(00,6) = ?????
C
C-----------------------------------------------------------------------
C
C     1) ELEMENT P1 EN DIMENSION 1
C
      NDS(01,1) = NPTFR
C     NDS(01,2) = ?????
      NDS(01,3) = 2
C     NDS(01,4) = ?????
      NDS(01,5) = NPTFX
C     NDS(01,6) = ?????
C
C-----------------------------------------------------------------------
C
C     2) ELEMENT QUADRATIQUE EN DIMENSION 1
C
C                 NPTFR+NSEGBOR (ICI DIMENSIONNEMENT PAR EXCES)
      NDS(02,1) = 2*NPTFR
C     NDS(02,2) = ?????
      NDS(02,3) = 3
C     NDS(02,4) = ?????
      NDS(02,5) = 2*NPTFX
C     NDS(02,6) = ?????
C
C     NOTE : THE QUADRATIC BOUNDARY POINT OF SEGMENT K WILL HAVE
C            THE NUMBER K+NPTFR, THUS WE HAVE HERE TO DECLARE
C            NDS(02,1) AS 2*NPTFR, THOUGH IN PARALLEL THERE MAY
C            EXISTING POINT K WITHOUT SEGMENT K IN THE DOMAIN
C
C
C-----------------------------------------------------------------------
C
C     10) ELEMENT P0 SUR DES TRIANGLES
C
      NDS(10,1) = NELEM
C     NDS(10,2) = ?????
      NDS(10,3) = 1
      NDS(10,4) = 3
      NDS(10,5) = NELMAX
C     NDS(10,6) = ??????
C
C-----------------------------------------------------------------------
C
C     11) ELEMENT P1 SUR DES TRIANGLES
C
      NDS(11,1) = NPOIN
      NDS(11,2) = (3*NELEM+NSEGBOR)/2
      NDS(11,3) = 3
      NDS(11,4) = 3
      NDS(11,5) = NPMAX
      NDS(11,6) = 3
C
C-----------------------------------------------------------------------
C
C     12) ELEMENT QUASI-BULLE SUR DES TRIANGLES
C
      NDS(12,1) = NPOIN+NELEM
      NDS(12,2) = (9*NELEM+NSEGBOR)/2
      NDS(12,3) = 4
      NDS(12,4) = 3
      NDS(12,5) = NPMAX+NELMAX
      NDS(12,6) = 6
C
C-----------------------------------------------------------------------
C
C     13) ELEMENT QUADRATIQUE SUR DES TRIANGLES
C
C     NOMBRE DE DDL = NOMBRE DE NOEUDS + NOMBRE DE SEGMENTS DU P1
      NDS(13,1) = NDS(11,1)+NDS(11,2)
C     NOMBRE TOTAL DE SEGMENTS (3 PAR SEGMENT LINEAIRE +
C                               6 SEGMENTS INTERIEURS A CHAQUE ELEMENT) 
      NDS(13,2) = 6*NELEM+3*NDS(11,2)
      NDS(13,3) = 6
      NDS(13,4) = 3
      NDS(13,5) = NPMAX+(3*NELMAX+NSEGBOR)/2
      NDS(13,6) = 15
C
C-----------------------------------------------------------------------
C
C     14) ELEMENT P1-ISO P1 SUR DES TRIANGLES
C
C     NOMBRE DE DDL = NOMBRE DE NOEUDS + NOMBRE DE SEGMENTS DU P1
      NDS(14,1) = NDS(11,1)+NDS(11,2)
C     NOMBRE TOTAL DE SEGMENTS (3 PAR SEGMENT LINEAIRE +
C                               3 SEGMENTS INTERIEURS A CHAQUE ELEMENT) 
      NDS(14,2) = 3*NELEM+3*NDS(11,2)
      NDS(14,3) = 6
      NDS(14,4) = 3
      NDS(14,5) = NPMAX+(3*NELMAX+NSEGBOR)/2
C     3 SEGMENTS LINEAIRES SUIVIS DES 9 P1-ISO P1
      NDS(14,6) = 12
C
C-----------------------------------------------------------------------
C
C     20) ELEMENT P0 SUR DES QUADRILATERES (AVEC CAS PARTICULIER DU 3D)
C
      NDS(20,1) = NELEM
C     NDS(20,2) = ?????
      NDS(20,3) = 1
      NDS(20,4) = 4
      NDS(20,5) = NELMAX
C     NDS(20,6) = ??????
C     BORDS LATERAUX DES MAILLAGES 3D DE PRISMES
C     IF(NPLAN.GE.2) THEN 
C       NDS(20,1) = NPTFR*(NPLAN-1)
C       NDS(20,5) = NPTFX*(NPLAN-1)
C     ENDIF
C
C-----------------------------------------------------------------------
C
C     21) ELEMENT P1 SUR DES QUADRILATERES (AVEC CAS PARTICULIER DU 3D)
C
      NDS(21,1) = NPOIN
C     NDS(21,2) = ?????
      NDS(21,3) = 4
      NDS(21,4) = 4
      NDS(21,5) = NPMAX
C     NDS(21,6) = ?????
C     BORDS LATERAUX DES MAILLAGES 3D DE PRISMES
C     IF(NPLAN.GE.2) THEN
C       NDS(21,1) = NPTFR*NPLAN
C       NDS(21,5) = NPTFX*NPLAN
C     ENDIF      
C
C-----------------------------------------------------------------------
C
C     ELEMENTS TRIDIMENSIONNELS (SEULEMENT EN MAILLAGE 3D)
C
C
C     30) ELEMENT T0 SUR DES TETRAEDRES 
C
      NDS(30,1) = NELEM
C     NDS(30,2) = ?????
      NDS(30,3) = 1
      NDS(30,4) = 4
      NDS(30,5) = NELMAX
C     NDS(30,6) = ?????
C
C-----------------------------------------------------------------------
C
C     31) ELEMENT T1 SUR DES TETRAEDRES 
C
CC MAILLAGE-3D
      NDS(31,1) = NPOIN
C     NDS(31,2) = ?????
      NDS(31,3) = 4
      NDS(31,4) = 4
      NDS(31,5) = NPMAX
C     NDS(31,6) = ?????
C
C-----------------------------------------------------------------------
C
C     IF(NPLAN.GT.1) : TO AVOID ERASING WHAT HAS BEEN DONE BY A PREVIOUS
C                      CALL BY TELEMAC-3D WHEN COUPLING WITH SISYPHE
C
C
      IF(NPLAN.GT.1) THEN
C
C     40) ELEMENT P0 SUR DES PRISMES
C
      NDS(40,1) = NELEM*(NPLAN-1)
C     NDS(40,2) = ??????????
      NDS(40,3) = 1
      NDS(40,4) = 5
      NDS(40,5) = NELMAX*(NPLAN-1)
      NDS(40,6) = 15
C
C-----------------------------------------------------------------------
C
C     41) ELEMENT P1 SUR DES PRISMES
C
      NDS(41,1) = NPOIN*NPLAN
      NDS(41,2) = NDS(11,2)*(3*NPLAN-2)+NPOIN*(NPLAN-1)
      NDS(41,3) = 6
      NDS(41,4) = 5
      NDS(41,5) = NPMAX*NPLAN
      NDS(41,6) = 15
C
C-----------------------------------------------------------------------
C
C     50) PRISMES DECOUPES EN TETRAEDRES T0
C
      NDS(50,1) = NELEM*(NPLAN-1)*3
C     NDS(50,2) = ?????
      NDS(50,3) = 1
      NDS(50,4) = 4
      NDS(50,5) = NELMAX*(NPLAN-1)*3
C     NDS(50,6) = ?????
C
C-----------------------------------------------------------------------
C
C     51) PRISMES DECOUPES EN TETRAEDRES  T1
C
      NDS(51,1) = NPOIN*NPLAN
C     NDS(51,2) = ?????
      NDS(51,3) = 4
      NDS(51,4) = 4
      NDS(51,5) = NPMAX*NPLAN
C     NDS(51,6) = ?????
C
C-----------------------------------------------------------------------
C
C     60) TRIANGLES P0 SUR FACE LATERALE DE MAILLAGE DE PRISMES 3D
C         (MAILLAGE DE PRISMES DECOUPE EN TETRAEDRES)
C
      NDS(60,1) = 2*NPTFR*(NPLAN-1)
C     NDS(60,2) = ?????
      NDS(60,3) = 1
      NDS(60,4) = 3
      NDS(60,5) = 2*NPTFX*(NPLAN-1)
C     NDS(60,6) = ??????
C
C-----------------------------------------------------------------------
C
C     61) TRIANGLES P1 SUR FACE LATERALE DE MAILLAGE DE PRISMES 3D
C         (MAILLAGE DE PRISMES DECOUPE EN TETRAEDRES)
C
      NDS(61,1) = NPTFR*NPLAN
C     NDS(61,2) = ?????????
      NDS(61,3) = 3
      NDS(61,4) = 3
      NDS(61,5) = NPTFX*NPLAN
      NDS(61,6) = 3
C
C-----------------------------------------------------------------------
C
C     70) QUADRILATERES Q0 SUR FACE LATERALE DE MAILLAGE DE PRISMES 3D
C
      NDS(70,1) = NPTFR*(NPLAN-1)
C     NDS(70,2) = ?????
      NDS(70,3) = 1
      NDS(70,4) = 4
      NDS(70,5) = NPTFX*(NPLAN-1)
C     NDS(70,6) = ??????
C
C-----------------------------------------------------------------------
C
C     71) QUADRILATERES Q1 SUR FACE LATERALE DE MAILLAGE DE PRISMES 3D
C
      NDS(71,1) = NPTFR*NPLAN
C     NDS(71,2) = ?????????
      NDS(71,3) = 4
      NDS(71,4) = 4
      NDS(71,5) = NPTFX*NPLAN
C     NDS(71,6) = ???????????
C
C     CORRESPOND A : IF(NPLAN.GT.1) THEN
      ENDIF
C
C-----------------------------------------------------------------------
C
C     80) TRIANGLES DE BORD P0 DES TETRAEDRES DANS UN MAILLAGE 3D 
C         NON STRUCTURE
C
      NDS(80,1) = NPTFR
C     NDS(80,2) = ?????
      NDS(80,3) = 1
      NDS(80,4) = 3
      NDS(80,5) = 3*NELEB
C     NDS(80,6) = ??????
C
C-----------------------------------------------------------------------
C
C     81) TRIANGLES DE BORDS P1 DES TETRAEDRES DANS UN MAILLAGE 3D
C         NON STRUCTURE
C
      NDS(81,1) = NPTFR
C     NDS(81,2) = ?????????
      NDS(81,3) = 3
      NDS(81,4) = 3
      NDS(81,5) = NPTFR
      NDS(81,6) = 3
C
C-----------------------------------------------------------------------
      
      RETURN
      END
