C                       *******************
                        SUBROUTINE MT02AA_2
C                       *******************
C
     *( A11 , A12 , A13 ,
     *        A22 , A23 ,
     *              A33 ,
     *  XMUL,SU,SV,U,V,
     *  XEL,YEL,SURFAC,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.8           28/11/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : CONSTRUCTION DU TERME DE DIFFUSION POUR ESTEL2D
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     A11,A12... |<-- |  ELEMENTS DE LA MATRICE
C |     XMUL       | -->|  FACTEUR MULTIPLICATIF
C |     SF,SG,SH   | -->|  STRUCTURES DE F,G,H
C |     SU,SV,SW   | -->|  STRUCTURES DE U,V,W
C |     F,G,H      | -->|  FONCTIONS INTERVENANT DANS LE CALCUL DE LA
C |                |    |  MATRICE.
C |     U,V,W      | -->|  COMPOSANTES D'UN VECTEUR INTERVENANT DANS LE
C |                |    |  CALCUL DE LA MATRICE.
C |     XEL,YEL    | -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |     SURFAC     | -->|  SURFACE DES TRIANGLES.
C |     IKLE1,2,3  | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
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
      USE BIEF, EX_MT02AA_2 => MT02AA_2
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*)
      DOUBLE PRECISION, INTENT(INOUT) ::        A22(*),A23(*)
      DOUBLE PRECISION, INTENT(INOUT) ::               A33(*)
      DOUBLE PRECISION, INTENT(IN) :: XMUL
      DOUBLE PRECISION, INTENT(IN) :: U(*),V(*)
C     STRUCTURE DE U ET V
      TYPE(BIEF_OBJ)  , INTENT(IN) :: SU,SV
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,3),YEL(NELMAX,3)
      DOUBLE PRECISION, INTENT(IN) :: SURFAC(NELMAX)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DECLARATIONS SPECIFIQUES A CE SOUS-PROGRAMME
C
      INTEGER IELEM,IELMNU,IELMNV,ISOU,ISOV
C
      DOUBLE PRECISION X2,X3,Y2,Y3
      DOUBLE PRECISION KSAT1,KSAT2,KSAT3
      DOUBLE PRECISION SOM,XSUR12
C
C=======================================================================
C
C     EXTRACTION DU TYPE D'ELEMENT DE LA VISCOSITE
C
      IELMNU = SU%ELM
      ISOU   = SU%DIM2
C
      IELMNV = SV%ELM
      ISOV   = SV%DIM2
C
      XSUR12 = XMUL / 12.D0
C      
C-----------------------------------------------------------------------
C TEST SUR LES TYPES DE U ET V
C U (Kr) : P0 et de dim 3 (car P1 discontinu) - V (Ks) : P0 et de dim 3
C-----------------------------------------------------------------------
C
      IF(IELMNU.EQ.10.AND.ISOU.EQ.3.AND.SU%DIMDISC.EQ.11
     &   .AND.
     &   IELMNV.EQ.10.AND.ISOV.EQ.3) THEN
C
      DO 4 IELEM = 1 , NELEM
C
C SORTIE DES 3 TERMES DE LA MATRICE V (Ks EST SYMETRIQUE)
C
        KSAT1=SV%R(IELEM)
        KSAT2=SV%R(IELEM+NELEM)
        KSAT3=SV%R(IELEM+2*NELEM)
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
         X2  =  XEL(IELEM,2)
         X3  =  XEL(IELEM,3)
C
         Y2  =  YEL(IELEM,2)
         Y3  =  YEL(IELEM,3)
C
C   INITIALISATION DES INTERMEDIAIRES DE CALCULS
C
         SOM = ( SU%R(IELEM+2*NELEM)
     &      +   SU%R(IELEM+NELEM)
     &      +   SU%R(IELEM) ) * XSUR12 / SURFAC(IELEM) 
C
C  TERMES DIAGONAUX
C
      A11(IELEM) = (KSAT1*Y2**2-2*KSAT1*Y2*Y3+KSAT1*Y3**2+KSAT2*X2**2-
     &  2*KSAT2*X2*X3+KSAT2*X3**2-2*KSAT3*Y2*X2+2*KSAT3*Y2*X3+
     &  2*KSAT3*X2*Y3-2*KSAT3*Y3*X3)*SOM
C        
      A22(IELEM) = (KSAT1*Y3**2+KSAT2*X3**2-2*KSAT3*Y3*X3)*SOM
C
      A33(IELEM) = (KSAT1*Y2**2+KSAT2*X2**2-2*KSAT3*Y2*X2)*SOM

C
C  TERMES EXTRADIAGONAUX
C
      A12(IELEM) = -(-KSAT1*Y2*Y3+KSAT1*Y3**2-KSAT2*X2*X3+KSAT2*X3**2+
     &          KSAT3*X2*Y3-2*KSAT3*Y3*X3+KSAT3*Y2*X3)*SOM
C
      A13(IELEM) = -(KSAT1*Y2**2-KSAT1*Y2*Y3+KSAT2*X2**2-KSAT2*X2*X3-
     &          2*KSAT3*Y2*X2+KSAT3*Y2*X3+KSAT3*X2*Y3)*SOM
C
      A23(IELEM) = (-KSAT1*Y2*Y3-KSAT2*X2*X3+KSAT3*Y2*X3+KSAT3*X2*Y3)*
     &          SOM
C
C   FIN DE LA BOUCLE SUR LES ELEMENTS
C
4     CONTINUE
C-----------------------------------------------------------------------
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,10) 
        IF (LNG.EQ.2) WRITE(LU,11)
10      FORMAT(1X,'MT02AA_2 (BIEF) : TYPES NON PREVUS')
11      FORMAT(1X,
     *  'MT02AA_2 (BIEF) : TYPES NOT AVAILABLE')
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
