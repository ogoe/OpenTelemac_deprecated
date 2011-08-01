C                       *****************
                        SUBROUTINE SEGBOR
C                       *****************
C
     *(NSEGBOR,IKLES,NELEM,NELMAX,NPOIN)
C
C***********************************************************************
C BIEF VERSION 5.9      19/06/2008    J-M HERVOUET (LNHE) 01 30 87 80 18
C COPYRIGHT 2009
C***********************************************************************
C
C    FONCTION : RECHERCHE DU NOMBRE DE SEGMENTS DE BORD DU MAILLAGE
C               EN INCLUANT EN PARALLELISME LES FRONTIERES INTERNES
C
C               ON REPREND ICI LE PRINCIPE DE VOISIN 
C               QUI SERA APPELE PLUS TARD
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    NSEGBOR     |<-- | NOMBRE DE SEGMENTS DE BORD|
C |    IKLE        | -->| NUMEROS GLOBAUX DES POINTS DE CHAQUE ELEMENT.
C |    X,Y         | -->| COORDONNEES DES ELEMENTS
C |    T1,2,3      | -->| TABLEAUX DE TRAVAIL ENTIERS.
C |    NELEM       | -->| NOMBRE D'ELEMENTS
C |    NELMAX      | -->| NOMBRE MAXIMUM D'ELEMENTS
C |    NPOIN       | -->| NOMBRE DE POINTS
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
      USE BIEF, EX_SEGBOR => SEGBOR
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)            :: NPOIN,NELMAX,NELEM
      INTEGER, INTENT(OUT)           :: NSEGBOR
      INTEGER, INTENT(IN)            :: IKLES(3,NELEM)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NFACE,NDP,KEL,IMAX,IFACE,IELEM,M1,M2,IV,IELEM2,IFACE2
      INTEGER I,ERR,I1,I2,IDIMAT
      INTEGER SOMFAC(2,4,2)
      DATA SOMFAC / 1,2 , 2,3 , 3,1 , 0,0   ,  1,2 , 2,3 , 3,4 , 4,1 /
C
C     TABLEAUX DE TRAVAIL ALLOUES DYNAMIQUEMENT
C 
      INTEGER, ALLOCATABLE :: IFABOR(:,:),MAT1(:),MAT2(:),MAT3(:)
      INTEGER, ALLOCATABLE :: NVOIS(:),IADR(:)
C
C-----------------------------------------------------------------------
C
      NFACE = 3
C     NOMBRE DE POINTS PAR ELEMENT
      NDP = 3
C     ADRESSE DANS SOMFAC
      KEL = 1
C
C     IDIMAT EST UNE MAJORATION DE LA SOMME DES NOMBRES DE VOISINS DE
C     TOUS LES POINTS (VOISIN = RELIE PAR UN SEGMENT)
C
      IDIMAT = NDP*2*NELEM
C
      ALLOCATE(MAT1(IDIMAT),STAT=ERR)
      ALLOCATE(MAT2(IDIMAT),STAT=ERR)
      ALLOCATE(MAT3(IDIMAT),STAT=ERR)
      ALLOCATE(IFABOR(NELEM,3),STAT=ERR)
      ALLOCATE(NVOIS(NPOIN),STAT=ERR)
      ALLOCATE(IADR(NPOIN),STAT=ERR)
C
      IF(ERR.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,1000) ERR
        IF(LNG.EQ.2) WRITE(LU,2000) ERR
1000    FORMAT(1X,'SEGBOR : ERREUR A L''ALLOCATION DE MEMOIRE : ',/,1X,
     *            'CODE D''ERREUR : ',1I6)
2000    FORMAT(1X,'SEGBOR: ERROR DURING ALLOCATION OF MEMORY: ',/,1X,
     *            'ERROR CODE: ',1I6)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C  CALCUL DU TABLEAU NVOIS POUR CHAQUE POINT
C  ATTENTION : NVOIS N'EST QU'UNE MAJORATION DU NOMBRE DE VOISINS
C              DONT LA SOMME VA FAIRE IDIMAT
C
      DO I=1,NPOIN
        NVOIS(I) = 0
      ENDDO
C
      DO IFACE = 1,NFACE
        DO IELEM=1,NELEM
          I1 = IKLES( SOMFAC(1,IFACE,KEL) , IELEM )
          I2 = IKLES( SOMFAC(2,IFACE,KEL) , IELEM )
          NVOIS(I1) = NVOIS(I1) + 1
          NVOIS(I2) = NVOIS(I2) + 1
        ENDDO
      ENDDO
C
C-----------------------------------------------------------------------
C
C  CALCUL DES ADRESSES DE CHAQUE POINT DANS UNE STRUCTURE DE TYPE
C  MATRICE COMPACTE
C
      IADR(1) = 1
      DO 50 I= 2,NPOIN
        IADR(I) = IADR(I-1) + NVOIS(I-1)
50    CONTINUE
C
      IMAX = IADR(NPOIN) + NVOIS(NPOIN) - 1
      IF(IMAX.GT.IDIMAT) THEN
        IF(LNG.EQ.1) WRITE(LU,51) IDIMAT,IMAX
        IF(LNG.EQ.2) WRITE(LU,52) IDIMAT,IMAX
51      FORMAT(1X,'SEGBOR: TAILLE DE MAT1,2,3 (',1I9,') INSUFFISANTE',/,
     *         1X,'IL FAUT AU MOINS : ',1I9)
52      FORMAT(1X,'SEGBOR: SIZE OF MAT1,2,3 (',1I9,') TOO SHORT',/,
     *         1X,'MINIMUM SIZE: ',1I9)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C  INITIALISATION A ZERO DE LA MATRICE COMPACTE
C
      DO I=1,IMAX
        MAT1(I) = 0
      ENDDO
C
C-----------------------------------------------------------------------
C
C  BOUCLE SUR LES FACES DE CHAQUE ELEMENT :
C
      DO 60 IFACE = 1 , NFACE
      DO 70 IELEM = 1 , NELEM
C
         IFABOR(IELEM,IFACE) = 0
C
C        NUMEROS GLOBAUX DES POINTS DE LA FACE :
C
         I1 = IKLES( SOMFAC(1,IFACE,KEL) , IELEM )
         I2 = IKLES( SOMFAC(2,IFACE,KEL) , IELEM )
C
C        NUMEROS GLOBAUX ORDONNES :
C
         M1 = MIN(I1,I2)
         M2 = MAX(I1,I2)
C
         DO 80 IV = 1,NVOIS(M1)
C
           IF(MAT1(IADR(M1)+IV-1).EQ.0) THEN
              MAT1(IADR(M1)+IV-1)=M2
              MAT2(IADR(M1)+IV-1)=IELEM
              MAT3(IADR(M1)+IV-1)=IFACE
              GO TO 81
           ELSEIF(MAT1(IADR(M1)+IV-1).EQ.M2) THEN
              IELEM2 = MAT2(IADR(M1)+IV-1)
              IFACE2 = MAT3(IADR(M1)+IV-1)
              IFABOR(IELEM,IFACE) = IELEM2
              IFABOR(IELEM2,IFACE2) = IELEM
              GO TO 81
           ENDIF
C
80       CONTINUE
C
         IF(LNG.EQ.1) WRITE(LU,82)
         IF(LNG.EQ.2) WRITE(LU,83)
82       FORMAT(1X,'SEGBOR : ERREUR DANS LE MAILLAGE       ',/,1X,
     *             '         PEUT-ETRE DES POINTS CONFONDUS')
83       FORMAT(1X,'SEGBOR : ERROR IN THE MESH             ',/,1X,
     *             '         MAYBE SUPERIMPOSED POINTS     ')
         CALL PLANTE(1)
         STOP
C
81       CONTINUE
C
70    CONTINUE
60    CONTINUE
C
      NSEGBOR = 0
      DO IFACE=1,NFACE
        DO IELEM=1,NELEM
          IF(IFABOR(IELEM,IFACE).EQ.0) NSEGBOR=NSEGBOR+1
        ENDDO
      ENDDO
C
      IF (LNG.EQ.1) WRITE(LU,500) NSEGBOR
      IF (LNG.EQ.2) WRITE(LU,501) NSEGBOR
500   FORMAT(1X,'SEGBOR (BIEF) : NOMBRE DE SEGMENTS DE BORD = ',1I6,/,
     *       1X,'EN COMPTANT CEUX DUS A LA DECOMPOSITION DE DOMAINE')
501   FORMAT(1X,'SEGBOR (BIEF) : NUMBER OF BOUNDARY SEGMENTS = ',1I6,/,
     *       1X,'INCLUDING THOSE DUE TO DOMAIN DECOMPOSITION')
C
C-----------------------------------------------------------------------
C
      DEALLOCATE(MAT1)
      DEALLOCATE(MAT2)
      DEALLOCATE(MAT3)
      DEALLOCATE(IFABOR)
      DEALLOCATE(NVOIS) 
      DEALLOCATE(IADR)  
C
C-----------------------------------------------------------------------
C
      RETURN
      END
