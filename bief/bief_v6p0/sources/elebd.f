C                       ****************
                        SUBROUTINE ELEBD
C                       ****************
C
     *(NELBOR,NULONE,KP1BOR,IFABOR,NBOR,IKLE,SIZIKL,IKLBOR,NELEM,NELMAX,
     * NPOIN,NPTFR,IELM,LIHBOR,KLOG,IFANUM,OPTASS,ISEG,T1,T2,T3)
C
C***********************************************************************
C BIEF VERSION 5.9      23/06/2008    J-M HERVOUET (LNHE) 01 30 87 80 18
C COPYRIGHT 1999
C***********************************************************************
C
C    PRISMES DECOUPES EN TETRAEDRES
C
C    FONCTION : 1) CONSTRUCTION DES TABLEAUX NELBOR ET NULONE
C               2) CONSTRUCTION DU TABLEAU KP1BOR
C               3) DISTINCTION DANS LE TABLEAU IFABOR ENTRE
C                  LES FACES DE BORD SOLIDES ET LES FACES LIQUIDES
C               4) CALCUL DE IKLBOR. 
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    NELBOR      |<-- | NUMERO DE L'ELEMENT ADJACENT AU KIEME SEGMENT|
C |    NULONE      |<-- | NUMERO LOCAL D'UN POINT DE BORD DANS         |
C |                |    | L'ELEMENT ADJACENT DONNE PAR NELBOR          |
C |    KP1BOR      |<-- | NUMERO DU POINT SUIVANT LE POINT DE BORD K.  |
C |    IFABOR      | -->| TABLEAU DES VOISINS DES FACES.
C |    NBOR        | -->| NUMERO GLOBAL DU POINT DE BORD K.
C |    IKLE        | -->| NUMEROS GLOBAUX DES POINTS DE CHAQUE ELEMENT.
C |    NELEM       | -->| NOMBRE TOTAL D'ELEMENTS DANS LE MAILLAGE.
C |    T1,2,3      | -->| TABLEAUX DE TRAVAIL ENTIERS.
C |    NPOIN       | -->| NOMBRE TOTAL DE POINTS DU DOMAINE.
C |    NPTFR       | -->| NOMBRE DE POINTS FRONTIERES.
C |    IELM        | -->| TYPE D'ELEMENT.
C |                |    | 11 : TRIANGLES.
C |                |    | 21 : QUADRILATERES.
C |    LIHBOR      | -->| TYPES DE CONDITIONS AUX LIMITES SUR H
C |    KLOG        | -->| CONVENTION POUR LA CONDITION LIMITE DE PAROI
C |    MXPTVS      | -->| NOMBRE MAXIMUM DE VOISINS D'UN POINT
C |    MXELVS      | -->| NOMBRE MAXIMUM D'ELEMENTS AUTOUR D'UN POINT
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
      USE BIEF, EX_ELEBD => ELEBD
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: KLOG,NELMAX,NELEM,SIZIKL
      INTEGER, INTENT(IN)    :: NPOIN,NPTFR,IELM,OPTASS
      INTEGER, INTENT(OUT)   :: NELBOR(NPTFR),NULONE(NPTFR,2)
      INTEGER, INTENT(OUT)   :: KP1BOR(NPTFR,2)
      INTEGER, INTENT(INOUT) :: NBOR(*)
      INTEGER, INTENT(INOUT) :: IFABOR(NELMAX,*)
      INTEGER, INTENT(IN)    :: IKLE(SIZIKL,*)
      INTEGER, INTENT(IN)    :: LIHBOR(NPTFR)
      INTEGER, INTENT(OUT)   :: IKLBOR(NPTFR,2)
      INTEGER, INTENT(INOUT) :: IFANUM(NELMAX,*)
      INTEGER, INTENT(IN)    :: ISEG(NPTFR)
      INTEGER, INTENT(OUT)   :: T1(NPOIN),T2(NPOIN),T3(NPOIN) 
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,NFACE,NPT,KEL,IPOIN
      INTEGER K,IFACE,I1,I2,N1,N2,IPT,IEL,I,K1,K2
C
      INTEGER SOMFAC(2,4,2)
C
      DATA SOMFAC / 1,2 , 2,3 , 3,1 , 0,0   ,  1,2 , 2,3 , 3,4 , 4,1 /
C
C-----------------------------------------------------------------------
C
      IF(IELM.EQ.11.OR.IELM.EQ.41.OR.IELM.EQ.51) THEN
C       TRIANGLES
        NFACE = 3
        NPT = 3
        KEL = 1
      ELSE
        IF(LNG.EQ.1) WRITE(LU,900) IELM
        IF(LNG.EQ.2) WRITE(LU,901) IELM
900     FORMAT(1X,'ELEBD : IELM=',1I6,' TYPE D''ELEMENT INCONNU')
901     FORMAT(1X,'ELEBD: IELM=',1I6,' UNKNOWN TYPE OF ELEMENT')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C  INITIALISATION DE T1,2,3 A ZERO
C
      DO IPOIN=1,NPOIN
        T1(IPOIN) = 0
        T2(IPOIN) = 0
        T3(IPOIN) = 0
      ENDDO
C
C  ON STOCKE K DANS TRAV(*,3) A L'ADRESSE NBOR(K)
C  CE QUI PERMET DE PASSER DE NUMERO GLOBAL A NUMERO DE BORD
C
      DO K = 1, NPTFR
        T3(NBOR(K)) = K
      ENDDO
C
C  BOUCLE SUR TOUTES LES FACES DE TOUS LES ELEMENTS :
C
      DO 20 IFACE = 1 , NFACE
      DO 10 IELEM = 1 , NELEM
C
      IF(IFABOR(IELEM,IFACE).EQ.-1) THEN
C
C      C'EST UNE VRAIE FACE DE BORD (LES FACES INTERNES EN PARALLELISME
C                                    SONT MARQUEES AVEC DES -2).
C      NUMEROS GLOBAUX DES POINTS DE LA FACE :
C
       I1 = IKLE( IELEM , SOMFAC(1,IFACE,KEL) )
       I2 = IKLE( IELEM , SOMFAC(2,IFACE,KEL) )
C
C      ON STOCKE DANS T1 ET T2 A L'ADRESSE I1 : I2 ET IELEM
C
       T1(I1) = I2
       T2(I1) = IELEM
C
C      UNE FACE LIQUIDE EST RECONNUE AVEC LA CONDITION LIMITE SUR H
C
       IF(NPTFR.GT.0) THEN
       IF(LIHBOR(T3(I1)).NE.KLOG.AND.LIHBOR(T3(I2)).NE.KLOG) THEN
C        FACE LIQUIDE : IFABOR=0  FACE SOLIDE : IFABOR=-1
         IFABOR(IELEM,IFACE)=0
       ENDIF
       ENDIF
C
      ENDIF
C
10    CONTINUE
20    CONTINUE
C
C     BOUCLE SUR TOUS LES POINTS:
C
      IF(NPTFR.GT.0) THEN
        DO I = 1 , NPOIN
          IF(T1(I).NE.0) THEN
C           POINT SUIVANT
            KP1BOR(T3(I),1)=T3(T1(I))
C           POINT PRECEDENT
            KP1BOR(T3(T1(I)),2)=T3(I)
            NELBOR(T3(I))=T2(I)
          ENDIF
        ENDDO
      ENDIF
C
C     VALEURS BIDON MISES POUR KP1BOR QUAND LE POINT SUIVANT EST DANS UN AUTRE
C     SOUS-DOMAINE. NELBOR ET NULONE MIS A ZERO
C
      IF(NCSIZE.GT.1) THEN
        DO 49 K1=1,NPTFR
          IF(ISEG(K1).GT.0) THEN
            KP1BOR(K1,1)=K1
            NELBOR(K1)=0
            NULONE(K1,1)=0
            NULONE(K1,2)=0
          ELSEIF(ISEG(K1).EQ.-9999) THEN
            KP1BOR(K1,1)=K1
            KP1BOR(K1,2)=K1
            NELBOR(K1)  =0
            NULONE(K1,1)=0
            NULONE(K1,2)=0
          ELSEIF(ISEG(K1).LT.0) THEN
            KP1BOR(K1,2)=K1
          ENDIF
49      CONTINUE
      ENDIF
C
C CALCUL DU TABLEAU NULONE
C
      DO 50 K1=1,NPTFR
C
      IF(NCSIZE.GT.1) THEN
        IF(ISEG(K1).GT.0.OR.ISEG(K1).EQ.-9999) GO TO 50
      ENDIF
C
      K2=KP1BOR(K1,1)
      IEL = NELBOR(K1)
      N1  = NBOR(K1)
      N2  = NBOR(K2)
C
      I1 = 0
      I2 = 0
C
      DO IPT=1,NPT
        IF(IKLE(IEL,IPT).EQ.N1) THEN
          NULONE(K1,1) = IPT
          I1 = 1
        ENDIF
        IF(IKLE(IEL,IPT).EQ.N2) THEN
          NULONE(K1,2) = IPT
          I2 = 1
        ENDIF
      ENDDO
C
      IF(I1.EQ.0.OR.I2.EQ.0) THEN
        IF(LNG.EQ.1) WRITE(LU,810) IEL
        IF(LNG.EQ.2) WRITE(LU,811) IEL
810     FORMAT(1X,'ELEBD: ERREUR DE NUMEROTATION DANS L''ELEMENT:',I6,/,
     *         1X,'       CAUSE POSSIBLE :                       '   ,/,
     *         1X,'       LE FICHIER DES CONDITIONS AUX LIMITES NE'  ,/,
     *         1X,'       CORRESPOND PAS AU FICHIER DE GEOMETRIE  ')
811     FORMAT(1X,'ELEBD: ERROR OF NUMBERING IN THE ELEMENT:',I6,
     *         1X,'       POSSIBLE REASON:                       '   ,/,
     *         1X,'       THE BOUNDARY CONDITION FILE IS NOT      '  ,/,
     *         1X,'       RELEVANT TO THE GEOMETRY FILE           ')
        CALL PLANTE(1)
        STOP
      ENDIF
C
50    CONTINUE
C
C  COMPUTING IKLBOR : LIKE IKLE FOR BOUNDARY POINTS, WITH BOUNDARY
C                     POINTS NUMBERING
C
      DO K=1,NPTFR
        IKLBOR(K,1) = K
        IKLBOR(K,2) = KP1BOR(K,1)
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
