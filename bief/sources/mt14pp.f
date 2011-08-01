C                       *****************
                        SUBROUTINE MT14PP
C                       *****************
C
     *( T,XM,PPQ,LEGO,XMUL,SU,SV,SW,U,V,W,SF,SG,SH,F,G,H,
     *  X,Y,Z,SURFAC,IKLE,NELEM,NELMAX,SIGMAG,SPECAD)
C
C***********************************************************************
C BIEF VERSION 6.0     27/04/2010    J-M HERVOUET / A. DECOENE    (LNHE) 
C                                         28/11/94    J-M JANIN    (LNH)                          
C***********************************************************************
C
C FONCTION : CONSTRUCTION DES COEFFICIENTS LAMBDA(I,J) DU SCHEMA MURD
C            DE TYPE N (VOIR THESE OU LIVRE DE J-M HERVOUET : SCHEMA
C            MURD EN DIMENSION 3 DANS LE CHAPITRE 6
C
C
C  MODIF JMH LE 16/08/99 : EPS AJOUTE POUR LES TESTS DE DIVISION
C  PAR ZERO (AVANT IL Y AVAIT IF(ALFA.GT.0.D0), ETC J'AI MIS EPS
C  POUR L'INSTANT EPS=1.D-10, A DEBATTRE
C
C  MODIF JMH LE 21/10/04 : MASS-LUMPING POUR COMPATIBILITE AVEC NOUVELLE
C  VERSION DE TRIDW2, TOUTES LES OCCURRENCES DE F123+W(...)
C  REMPLACEES PAR 4*W(...)
C
C  MODIF JMH LE 04/08/08 : INVERSION DES DIMENSIONS DE XM (VOIR AUSSI
C                          MURD3D)
C
C-----------------------------------------------------------------------
C
C     L'ELEMENT EST LE PRISME P1
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     A11,A12... |<-- |  ELEMENTS DE LA MATRICE
C |     XMUL       | -->|  FACTEUR MULTIPLICATIF
C |     SF,SG,SH   | -->|  STRUCTURES DE F,G ET H.
C |     SU,SV,SW   | -->|  STRUCTURES DE U,V ET W.
C |     F,G,H      | -->|  FONCTIONS INTERVENANT DANS LE CALCUL DE LA
C |                |    |  MATRICE.
C |     U,V,W      | -->|  COMPOSANTES D'UN VECTEUR INTERVENANT DANS LE
C |                |    |  CALCUL DE LA MATRICE.
C |     X,Y,Z      | -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |     SURFAC     | -->|  SURFACE DES TRIANGLES.
C |     IKLE1..6   | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
C |     NELEM      | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |     NELMAX     | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : AUCUN
C
C**********************************************************************
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF, EX_MT14PP => MT14PP
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE(NELMAX,6),PPQ(6,6)
C
      DOUBLE PRECISION, INTENT(INOUT) :: T(NELMAX,6),XM(30,NELMAX)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
      DOUBLE PRECISION, INTENT(IN) :: U(*),V(*),W(*),F(*),G(*),H(*)
C
      LOGICAL, INTENT(IN) :: LEGO,SIGMAG,SPECAD
C
C     STRUCTURES DE U,V,W
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SU,SV,SW,SF,SG,SH
C
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*),Z(*)
      DOUBLE PRECISION, INTENT(IN) :: SURFAC(NELMAX)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DECLARATIONS SPECIFIQUES A CE SOUS-PROGRAMME
C
      INTEGER IELEM,IELMW
      DOUBLE PRECISION XSUR3
      DOUBLE PRECISION LI0J0,LI0K0,LJ0I0,LJ0K0,LK0I0,LK0J0
      DOUBLE PRECISION LI3J3,LI3K3,LJ3I3,LJ3K3,LK3I3,LK3J3
      DOUBLE PRECISION LI3J0,LI3K0,LJ3I0,LJ3K0,LK3I0,LK3J0
      DOUBLE PRECISION LI3I0,LJ3J0,LK3K0,LI0I3
      DOUBLE PRECISION ALFA,ALFAJ,ALFAK,ALFA1,ALFA2,BETA,SOM0,SOM3
      DOUBLE PRECISION ALFAI0,ALFAJ0,ALFAK0,ALFAI3,ALFAJ3,ALFAK3
      DOUBLE PRECISION ALFAII,ALFAIJ,ALFAIK,ALFAJI,ALFAJJ,ALFAJK
      DOUBLE PRECISION ALFAKI,ALFAKJ,ALFAKK,EPS,DELTA
      DATA EPS /1.D-10/
C
      INTEGER IPLUS1(6),IPLUS3(6),INDIC(0:7)
      INTEGER IXM,I0,J0,K0,I3,J3,K3
C
      CHARACTER(LEN=16) :: FORMUL
C
C-----------------------------------------------------------------------
C
C     DONNEES
C
      DATA IPLUS1 / 2 , 3 , 1 , 5 , 6 , 4 /
      DATA IPLUS3 / 4 , 5 , 6 , 1 , 2 , 3 /
      DATA INDIC  / 4 , 4 , 5 , 3 , 6 , 2 , 1 , 1 /
C
C=======================================================================
C
      IELMW = SW%ELM
C                                                         *   *
C     CALCUL DES COEFFICIENTS A(I) (INTEGRALES DE H U.GRAD(PSI), EN PRENANT
C     SEULEMENT LA PARTIE HORIZONTALE DE U). VOIR THESE JMH.
C
C     NOTE: THIS CALL VC04PP IS ALREADY DONE BY FLUX3D IN TELEMAC-3D 
C           (THROUGH A CALL TO VECTOR, BUT T IS AT THAT TIME MESH3D%W 
C            AND IS NOT KEPT). MAYBE OPTIMIZATION POSSIBLE.
C
      FORMUL='             HOR'
      CALL VC04PP(XMUL,SU,SV,SW,U,V,W,SF,SG,SH,F,G,H,X,Y,Z,
     * IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),IKLE(1,5),IKLE(1,6),
     * NELEM,NELMAX,T(1,1),T(1,2),T(1,3),T(1,4),T(1,5),T(1,6),
     * SPECAD,FORMUL)
C     FORMUL NORMALLY STARTS WITH VGRADP OR VGRADP2, OR VGRADP22 TO KNOW 
C     SIGMAG AND SPECAD, USELESS HERE. 
      XSUR3 = XMUL/3.D0
C
C-----------------------------------------------------------------------
C
      IF(IELMW.NE.41) THEN
        STOP 'MT14PP DISCRETISATION NON PREVUE'
      ENDIF
C
C     CONSTRUCTION DES COEFFICIENTS B(I) PAR ELEMENTS
C     ET CONSTRUCTION DES LAMBDA(I,J)
C
      DO 1 IELEM = 1 , NELEM
C
         IXM = 0
         IF (W(IKLE(IELEM,1)).GT.0.D0) IXM = 1
         IF (W(IKLE(IELEM,2)).GT.0.D0) IXM = IXM + 2
         IF (W(IKLE(IELEM,3)).GT.0.D0) IXM = IXM + 4
C
C        NOTE JMH:
C        IT SEEMS THAT INDIC IS THE POINT THAT RECEIVES THE VERTICAL
C        VELOCITY WHICH IS THE ONLY ONE OF ITS SIGN
C
C        ALL W < 0         INDIC = 4 (RANDOM, WHY NOT 1 ?)
C        ONLY W1 > 0       INDIC = 4
C        ONLY W2 > 0       INCIC = 5
C        ONLY W3 > 0       INDIC = 6
C        ALL W > 0         INDIC = 1 (RANDOM, WHY NOT 4 ?)
C        ONLY W1 < 0       INDIC = 1
C        ONLY W2 < 0       INDIC = 2
C        ONLY W3 < 0       INDIC = 3
C
         I0 = INDIC(IXM)
         J0 = IPLUS1(I0)
         K0 = IPLUS1(J0)
         I3 = IPLUS3(I0)
         J3 = IPLUS1(I3)
         K3 = IPLUS1(J3)
C
         ALFA = XSUR3 * SURFAC(IELEM)
C
         LI3I0 = ALFA * ABS(W(IKLE(IELEM,MIN(I0,I3))))
         LI0I3 = 0.D0
         IF (IXM.GE.1.AND.IXM.LE.6) THEN
            LI0I3 = LI3I0
            LI3I0 = 0.D0
         ENDIF
         LJ3J0 = ALFA * ABS(W(IKLE(IELEM,MIN(J0,J3))))
         XM(PPQ(J0,J3),IELEM) = 0.D0
         LK3K0 = ALFA * ABS(W(IKLE(IELEM,MIN(K0,K3))))
         XM(PPQ(K0,K3),IELEM) = 0.D0
C
         ALFAI0 = -T(IELEM,I0)
         ALFAJ0 = -T(IELEM,J0)
         ALFAK0 = -T(IELEM,K0)
         ALFAI3 = -T(IELEM,I3)
         ALFAJ3 = -T(IELEM,J3)
         ALFAK3 = -T(IELEM,K3)
C
         LJ0I0 = MAX(0.D0,MIN(ALFAI0,-ALFAJ0))
         LK0I0 = MAX(0.D0,MIN(ALFAI0,-ALFAK0))
         LI3J3 = MAX(0.D0,MIN(ALFAJ3,-ALFAI3))
         LI3K3 = MAX(0.D0,MIN(ALFAK3,-ALFAI3))
C
         ALFAJ = MIN(LJ3J0,MAX(LI3J3,LJ0I0))
         ALFAK = MIN(LK3K0,MAX(LI3K3,LK0I0))
         ALFA  = MIN(LI0I3,ALFAJ+ALFAK)
         IF (ALFA.GT.EPS) THEN
!
!           PROGRAMMATION JMH
!
C           IF(ALFA.LT.ALFAJ+ALFAK) THEN
!             LIMITATION OF ALFAJ AND ALFAK WITH THE IDEA THAT THE
!             TWO CORRECTED FLUXES WILL BE AS EQUAL AS POSSIBLE
!             MOREOVER WE MUST HAVE ALFAJ+ALFAK=ALFA 
C             IF(ALFAK.GT.ALFAJ) THEN
C               DELTA=MIN(ALFA,ALFAK-ALFAJ)
C               ALFAK=0.5D0*(ALFA+DELTA)
C               ALFAJ=0.5D0*(ALFA-DELTA)
C             ELSE
C               DELTA=MIN(ALFA,ALFAJ-ALFAK)
C               ALFAK=0.5D0*(ALFA-DELTA)
C               ALFAJ=0.5D0*(ALFA+DELTA)              
C             ENDIF
C           ENDIF
!
!           PROGRAMMATION JMJ
!
            BETA = ALFA / (ALFAJ+ALFAK)
            ALFAJ = ALFAJ * BETA
            ALFAK = ALFAK * BETA
C
            LI0I3 = LI0I3-ALFA
            LJ3J0 = LJ3J0-ALFAJ
            LK3K0 = LK3K0-ALFAK
C
            ALFAI3 = ALFAI3+ALFA
            ALFAI0 = ALFAI0-ALFA
            ALFAJ0 = ALFAJ0+ALFAJ
            ALFAJ3 = ALFAJ3-ALFAJ
            ALFAK0 = ALFAK0+ALFAK
            ALFAK3 = ALFAK3-ALFAK
C
            LJ0I0 = MAX(0.D0,MIN(ALFAI0,-ALFAJ0))
            LK0I0 = MAX(0.D0,MIN(ALFAI0,-ALFAK0))
            LI3J3 = MAX(0.D0,MIN(ALFAJ3,-ALFAI3))
            LI3K3 = MAX(0.D0,MIN(ALFAK3,-ALFAI3))
C
         ENDIF
C
         LI0J0 = MAX(0.D0,MIN(ALFAJ0,-ALFAI0))
         LK0J0 = MAX(0.D0,MIN(ALFAJ0,-ALFAK0))
         LI0K0 = MAX(0.D0,MIN(ALFAK0,-ALFAI0))
         LJ0K0 = MAX(0.D0,MIN(ALFAK0,-ALFAJ0))
         LJ3I3 = MAX(0.D0,MIN(ALFAI3,-ALFAJ3))
         LJ3K3 = MAX(0.D0,MIN(ALFAK3,-ALFAJ3))
         LK3I3 = MAX(0.D0,MIN(ALFAI3,-ALFAK3))
         LK3J3 = MAX(0.D0,MIN(ALFAJ3,-ALFAK3))
C
         XM(PPQ(J0,I3),IELEM) = 0.D0
         XM(PPQ(K0,I3),IELEM) = 0.D0
         XM(PPQ(I0,J3),IELEM) = 0.D0
         XM(PPQ(K0,J3),IELEM) = 0.D0
         XM(PPQ(I0,K3),IELEM) = 0.D0
         XM(PPQ(J0,K3),IELEM) = 0.D0
C
         ALFAI0 = 0.D0
         ALFAJ0 = 0.D0
         ALFAK0 = 0.D0
         ALFAI3 = 0.D0
         ALFAJ3 = 0.D0
         ALFAK3 = 0.D0
C
         SOM0 = LJ0I0+LK0I0
         SOM3 = LI3J3+LI3K3
C
         ALFA1 = MIN(LI0I3,SOM0)
         IF (ALFA1.GT.EPS) THEN
C
            BETA = 1.D0 / SOM0
            ALFAJ0 = LJ0I0 * BETA
            ALFAK0 = LK0I0 * BETA
            ALFA = MAX(0.D0,ALFA1-SOM3)
C
            LJ0I0 = LJ0I0 - ALFAJ0*ALFA1
            LK0I0 = LK0I0 - ALFAK0*ALFA1
            XM(PPQ(J0,I3),IELEM) = ALFAJ0*ALFA
            XM(PPQ(K0,I3),IELEM) = ALFAK0*ALFA
C
         ENDIF
C
         ALFA2 = MIN(LI0I3,SOM3)
         IF (ALFA2.GT.EPS) THEN
C
            BETA = 1.D0 / SOM3
            ALFAJ3 = LI3J3 * BETA
            ALFAK3 = LI3K3 * BETA
            ALFA = MAX(0.D0,ALFA2-SOM0)
C
            LI3J3 = LI3J3 - ALFAJ3*ALFA2
            LI3K3 = LI3K3 - ALFAK3*ALFA2
            XM(PPQ(I0,J3),IELEM) = ALFAJ3*ALFA
            XM(PPQ(I0,K3),IELEM) = ALFAK3*ALFA
C
         ENDIF
C
         ALFA = MIN(ALFA1,ALFA2)
         IF (ALFA.GT.0.D0) THEN
C
            ALFAJ3 = ALFAJ3 * ALFA
            ALFAK3 = ALFAK3 * ALFA
            ALFAJJ = ALFAJ3 * ALFAJ0
            ALFAKK = ALFAK3 * ALFAK0
            ALFAJK = ALFAJ3 * ALFAK0
            ALFAKJ = ALFAK3 * ALFAJ0
            ALFA = MIN(ALFAJK,ALFAKJ)
C
            XM(PPQ(J0,J3),IELEM) = ALFAJJ + ALFA
            XM(PPQ(K0,K3),IELEM) = ALFAKK + ALFA
            XM(PPQ(K0,J3),IELEM) = ALFAJK - ALFA
            XM(PPQ(J0,K3),IELEM) = ALFAKJ - ALFA
C
         ENDIF
C
         XM(PPQ(I0,I3),IELEM) = MAX(0.D0,LI0I3-MAX(SOM3,SOM0))
C
         LJ3I0 = 0.D0
         LK3I0 = 0.D0
         LI3J0 = 0.D0
         LK3J0 = 0.D0
         LI3K0 = 0.D0
         LJ3K0 = 0.D0
C
         SOM0 = LI0J0+LI0K0
         SOM3 = LJ3I3+LK3I3
C
         ALFA1 = MIN(LI3I0,SOM0)
         IF (ALFA1.GT.EPS) THEN
C
            BETA = 1.D0 / SOM0
            ALFAJ0 = LI0J0 * BETA
            ALFAK0 = LI0K0 * BETA
            ALFA = MAX(0.D0,ALFA1-SOM3)
C
            LI0J0 = LI0J0 - ALFAJ0*ALFA1
            LI0K0 = LI0K0 - ALFAK0*ALFA1
            LI3J0 = LI3J0 + ALFAJ0*ALFA
            LI3K0 = LI3K0 + ALFAK0*ALFA
C
         ENDIF
C
         XM(PPQ(I0,J0),IELEM) = LI0J0
         XM(PPQ(I0,K0),IELEM) = LI0K0
C
         ALFA2 = MIN(LI3I0,SOM3)
         IF (ALFA2.GT.EPS) THEN
C
            BETA = 1.D0 / SOM3
            ALFAJ3 = LJ3I3 * BETA
            ALFAK3 = LK3I3 * BETA
            ALFA = MAX(0.D0,ALFA2-SOM0)
C
            LJ3I3 = LJ3I3 - ALFAJ3*ALFA2
            LK3I3 = LK3I3 - ALFAK3*ALFA2
            LJ3I0 = LJ3I0 + ALFAJ3*ALFA
            LK3I0 = LK3I0 + ALFAK3*ALFA
C
         ENDIF
C
         XM(PPQ(J3,I3),IELEM) = LJ3I3
         XM(PPQ(K3,I3),IELEM) = LK3I3
C
         ALFA = MIN(ALFA1,ALFA2)
         IF (ALFA.GT.0.D0) THEN
C
            ALFAJ0 = ALFAJ0 * ALFA
            ALFAK0 = ALFAK0 * ALFA
            ALFAJJ = ALFAJ0 * ALFAJ3
            ALFAKK = ALFAK0 * ALFAK3
            ALFAJK = ALFAJ0 * ALFAK3
            ALFAKJ = ALFAK0 * ALFAJ3
            ALFA = MIN(ALFAJK,ALFAKJ)
C
            LJ3J0 = LJ3J0 + ALFAJJ + ALFA
            LK3K0 = LK3K0 + ALFAKK + ALFA
            LK3J0 = LK3J0 + ALFAJK - ALFA
            LJ3K0 = LJ3K0 + ALFAKJ - ALFA
C
         ENDIF
C
         LI3I0 = MAX(0.D0,LI3I0-MAX(SOM0,SOM3))
C
         SOM0 = LJ0I0+LJ0K0
         SOM3 = LI3J3+LK3J3
C
         ALFA1 = MIN(LJ3J0,SOM0)
         IF (ALFA1.GT.EPS) THEN
C
            BETA = 1.D0 / SOM0
            ALFAI0 = LJ0I0 * BETA
            ALFAK0 = LJ0K0 * BETA
            ALFA = MAX(0.D0,ALFA1-SOM3)
C
            LJ0I0 = LJ0I0 - ALFAI0*ALFA1
            LJ0K0 = LJ0K0 - ALFAK0*ALFA1
            LJ3I0 = LJ3I0 + ALFAI0*ALFA
            LJ3K0 = LJ3K0 + ALFAK0*ALFA
C
         ENDIF
C
         XM(PPQ(J0,I0),IELEM) = LJ0I0
         XM(PPQ(J0,K0),IELEM) = LJ0K0
C
         ALFA2 = MIN(LJ3J0,SOM3)
         IF (ALFA2.GT.EPS) THEN
C
            BETA = 1.D0 / SOM3
            ALFAI3 = LI3J3 * BETA
            ALFAK3 = LK3J3 * BETA
            ALFA = MAX(0.D0,ALFA2-SOM0)
C
            LI3J3 = LI3J3 - ALFAI3*ALFA2
            LK3J3 = LK3J3 - ALFAK3*ALFA2
            LI3J0 = LI3J0 + ALFAI3*ALFA
            LK3J0 = LK3J0 + ALFAK3*ALFA
C
         ENDIF
C
         XM(PPQ(I3,J3),IELEM) = LI3J3
         XM(PPQ(K3,J3),IELEM) = LK3J3
C
         ALFA = MIN(ALFA1,ALFA2)
         IF (ALFA.GT.0.D0) THEN
C
            ALFAI0 = ALFAI0 * ALFA
            ALFAK0 = ALFAK0 * ALFA
            ALFAII = ALFAI0 * ALFAI3
            ALFAKK = ALFAK0 * ALFAK3
            ALFAIK = ALFAI0 * ALFAK3
            ALFAKI = ALFAK0 * ALFAI3
            ALFA = MIN(ALFAIK,ALFAKI)
C
            LI3I0 = LI3I0 + ALFAII + ALFA
            LK3K0 = LK3K0 + ALFAKK + ALFA
            LK3I0 = LK3I0 + ALFAIK - ALFA
            LI3K0 = LI3K0 + ALFAKI - ALFA
C
         ENDIF
C
         LJ3J0 = MAX(0.D0,LJ3J0-MAX(SOM0,SOM3))
C
         SOM0 = LK0I0+LK0J0
         SOM3 = LI3K3+LJ3K3
C
         ALFA1 = MIN(LK3K0,SOM0)
         IF (ALFA1.GT.EPS) THEN
C
            BETA = 1.D0 / SOM0
            ALFAI0 = LK0I0 * BETA
            ALFAJ0 = LK0J0 * BETA
            ALFA = MAX(0.D0,ALFA1-SOM3)
C
            LK0I0 = LK0I0 - ALFAI0*ALFA1
            LK0J0 = LK0J0 - ALFAJ0*ALFA1
            LK3I0 = LK3I0 + ALFAI0*ALFA
            LK3J0 = LK3J0 + ALFAJ0*ALFA
C
         ENDIF
C
         XM(PPQ(K0,I0),IELEM) = LK0I0
         XM(PPQ(K0,J0),IELEM) = LK0J0
C
         ALFA2 = MIN(LK3K0,SOM3)
         IF (ALFA2.GT.EPS) THEN
C
            BETA = 1.D0 / SOM3
            ALFAI3 = LI3K3 * BETA
            ALFAJ3 = LJ3K3 * BETA
            ALFA = MAX(0.D0,ALFA2-SOM0)
C
            LI3K3 = LI3K3 - ALFAI3*ALFA2
            LJ3K3 = LJ3K3 - ALFAJ3*ALFA2
            LI3K0 = LI3K0 + ALFAI3*ALFA
            LJ3K0 = LJ3K0 + ALFAJ3*ALFA
C
         ENDIF
C
         XM(PPQ(I3,K3),IELEM) = LI3K3
         XM(PPQ(J3,K3),IELEM) = LJ3K3
C
         ALFA = MIN(ALFA1,ALFA2)
         IF (ALFA.GT.0.D0) THEN
C
            ALFAI0 = ALFAI0 * ALFA
            ALFAJ0 = ALFAJ0 * ALFA
            ALFAII = ALFAI0 * ALFAI3
            ALFAJJ = ALFAJ0 * ALFAJ3
            ALFAIJ = ALFAI0 * ALFAJ3
            ALFAJI = ALFAJ0 * ALFAI3
            ALFA = MIN(ALFAIJ,ALFAJI)
C
            LI3I0 = LI3I0 + ALFAII + ALFA
            LJ3J0 = LJ3J0 + ALFAJJ + ALFA
            LJ3I0 = LJ3I0 + ALFAIJ - ALFA
            LI3J0 = LI3J0 + ALFAJI - ALFA
C
         ENDIF
C
         LK3K0 = MAX(0.D0,LK3K0-MAX(SOM0,SOM3))
C
         XM(PPQ(I3,I0),IELEM) = LI3I0
         XM(PPQ(J3,I0),IELEM) = LJ3I0
         XM(PPQ(K3,I0),IELEM) = LK3I0
         XM(PPQ(I3,J0),IELEM) = LI3J0
         XM(PPQ(J3,J0),IELEM) = LJ3J0
         XM(PPQ(K3,J0),IELEM) = LK3J0
         XM(PPQ(I3,K0),IELEM) = LI3K0
         XM(PPQ(J3,K0),IELEM) = LJ3K0
         XM(PPQ(K3,K0),IELEM) = LK3K0
C
1     CONTINUE
C
C-----------------------------------------------------------------------
C
C     CALCUL DE LA SOMME DE CHAQUE LIGNE (AVEC UN SIGNE -)
C     LES TERMES DIAGONAUX SONT NULS
C
      IF(LEGO) THEN
C
        DO IELEM = 1,NELEM
C
            T(IELEM,1) = -XM(01,IELEM)-XM(02,IELEM)
     *                   -XM(03,IELEM)-XM(04,IELEM)-XM(05,IELEM)
            T(IELEM,2) = -XM(16,IELEM)-XM(06,IELEM)
     *                   -XM(07,IELEM)-XM(08,IELEM)-XM(09,IELEM)
            T(IELEM,3) = -XM(17,IELEM)-XM(21,IELEM)
     *                   -XM(10,IELEM)-XM(11,IELEM)-XM(12,IELEM)
            T(IELEM,4) = -XM(18,IELEM)-XM(22,IELEM)
     *                   -XM(25,IELEM)-XM(13,IELEM)-XM(14,IELEM)
            T(IELEM,5) = -XM(19,IELEM)-XM(23,IELEM)
     *                   -XM(26,IELEM)-XM(28,IELEM)-XM(15,IELEM)
            T(IELEM,6) = -XM(20,IELEM)-XM(24,IELEM)
     *                   -XM(27,IELEM)-XM(29,IELEM)-XM(30,IELEM)
C
        ENDDO
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
