!                       *****************
                        SUBROUTINE INTERP
!                       *****************
!
     &(XINIT , YINIT , IKINIT , NPINIT , NEINIT ,
     & X , Y , NPOIN , NPMAX , SHP , ELT)
!
!***********************************************************************
! PROGICIEL : STBTEL  V5.2        24/04/91    J-C GALLAND  (LNH)
!                               09/11/94    P LANG / TRIGRID (LHF)
!***********************************************************************
!
! FONCTION : INTERPOLATION DES FONDS SUR LE MAILLAGE
!
!----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________
! |      NOM       |MODE|                   ROLE
! |________________|____|______________________________________________
! |    X,Y         | -->|  COORDONNEES DES POINTS DU MAILLAGE
! |    ZF          |<-- |  COTES DU FOND
! |    XRELV,YRELV | -->|  COORDONNEES DES POINTS DE BATHY
! |    ZRELV       | -->|  COTES DES POINTS DE BATHY
! |    NBAT        | -->|  NOMBRE DE POINTS DE BATHY
! |    NBOR        | -->|  NUMEROTATION DES ELEMENTS DE BORD
! |    NPTFR       | -->|  NOMBRE DE POINTS FRONTIERE
! |    NFOND       | -->|  CANAUX DES FICHIERS DES FONDS
! |    NBFOND      | -->|  NOMBRE DE FICHIERS FONDS DONNES PAR
! |                |    |  L'UTILISATEUR (5 MAXI)
! |    FOND        | -->|  NOM DES FICHIERS DES FONDS
! |    DM          | -->|  DISTANCE MINIMALE A LA FRONTIERE
! |                |    |  POUR L'INTERPOLATION DES FONDS
! |    FONTRI      | -->|  INDICATEUR DE LECTURE DES FONDS DANS TRIGRID
! |    CORTRI      | -->|  CORRECTION DES FONDS POUR TRIGRID
! |                |    |
! | COMMON:        |    |
! |  GEO:          |    |
! |    MESH        | -->| TYPE DES ELEMENTS DU MAILLAGE
! |    NDP         | -->| NOMBRE DE NOEUDS PAR ELEMENTS
! |    NPOIN       | -->| NOMBRE TOTAL DE NOEUDS DU MAILLAGE
! |    NELEM       | -->| NOMBRE TOTAL D'ELEMENTS DU MAILLAGE
! |    NPMAX       | -->| DIMENSION EFFECTIVE DES TABLEAUX X ET Y
! |                |    | (NPMAX = NPOIN + 0.1*NELEM)
! |    NELMAX      | -->| DIMENSION EFFECTIVE DES TABLEAUX CONCERNANT
! |                |    | LES ELEMENTS (NELMAX = NELEM + 0.2*NELEM)
! |                |    |
! |________________|____|______________________________________________
! MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
!----------------------------------------------------------------------
!
! APPELE PAR : STBTEL
! APPEL DE : LECFON, FASP
!
!**********************************************************************
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
      INTEGER IELEM , JELEM , NPINIT , NEINIT , IPOIN , NPOIN , NPMAX
      INTEGER IKINIT(NEINIT,3) , ELT(NPMAX)
!
      DOUBLE PRECISION XINIT(NPINIT) , YINIT(NPINIT)
      DOUBLE PRECISION X(NPMAX) , Y(NPMAX) , SHP(NPMAX,3)
      DOUBLE PRECISION XP,YP,A1,A2,A3,C1,C2,X1,X2,X3,Y1,Y2,Y3
!
!=======================================================================
!
      IF (LNG.EQ.1) WRITE(LU,2)
      IF (LNG.EQ.2) WRITE(LU,4)
!
2     FORMAT(//,1X,'INTERPOLATION DES DONNEES',/,
     &          1X,'-------------------------',/)
4     FORMAT(//,1X,'DATA INTERPOLATION',/,
     &          1X,'------------------',/)
!
      DO 10 IPOIN = 1,NPOIN
!
         XP = X(IPOIN)
         YP = Y(IPOIN)
         C1 = -999999.D0
!
         DO 20 IELEM = 1,NEINIT
            X1 = XINIT(IKINIT(IELEM,1))
            X2 = XINIT(IKINIT(IELEM,2))
            X3 = XINIT(IKINIT(IELEM,3))
            Y1 = YINIT(IKINIT(IELEM,1))
            Y2 = YINIT(IKINIT(IELEM,2))
            Y3 = YINIT(IKINIT(IELEM,3))
            A1 = (X3-X2)*(YP-Y2) - (Y3-Y2)*(XP-X2)
            A2 = (X1-X3)*(YP-Y3) - (Y1-Y3)*(XP-X3)
            A3 = (X2-X1)*(YP-Y1) - (Y2-Y1)*(XP-X1)
            IF (A1.GE.0.AND.A2.GE.0.AND.A3.GE.0) GOTO 30
            C2 = MIN(A1,A2,A3) / ((X3-X2)*(Y1-Y2)-(Y3-Y2)*(X1-X2))
            IF (C2.GT.C1) THEN
               C1 = C2
               JELEM = IELEM
            ENDIF
20       CONTINUE
!
         IF (LNG.EQ.1) WRITE(LU,*) 'EXTRAPOLATION NECESSAIRE POUR ',
     &                             'LE POINT :',IPOIN
         IF (LNG.EQ.2) WRITE(LU,*) 'EXTRAPOLATION REQUIRED FOR ',
     &                             'THE NODE :',IPOIN
         IELEM = JELEM
         X1 = XINIT(IKINIT(IELEM,1))
         X2 = XINIT(IKINIT(IELEM,2))
         X3 = XINIT(IKINIT(IELEM,3))
         Y1 = YINIT(IKINIT(IELEM,1))
         Y2 = YINIT(IKINIT(IELEM,2))
         Y3 = YINIT(IKINIT(IELEM,3))
         A1 = (X3-X2)*(YP-Y2) - (Y3-Y2)*(XP-X2)
         A2 = (X1-X3)*(YP-Y3) - (Y1-Y3)*(XP-X3)
         A3 = (X2-X1)*(YP-Y1) - (Y2-Y1)*(XP-X1)
!
30       CONTINUE
         C1 = (X3-X2)*(Y1-Y2)-(Y3-Y2)*(X1-X2)
         SHP(IPOIN,1) = A1/C1
         SHP(IPOIN,2) = A2/C1
         SHP(IPOIN,3) = A3/C1
         ELT(IPOIN) = IELEM
!
10    CONTINUE
!
!-----------------------------------------------------------------------
!
      RETURN
      END