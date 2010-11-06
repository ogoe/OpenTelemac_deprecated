C                       *****************
                        SUBROUTINE INIP4D
C                       *****************
C
     *( U , V , T , W , X , Y , SHP1 ,SHP2 , SHP3 , SHT , SHF , ELT ,
     * ETA , FRE , XCONV , YCONV , TCONV, FCONV , TETA , FREQ ,IKLE2 , 
     * NPOIN2 , NELEM2 , NPLAN  , IFF , NF  ,IFABOR,GOODELT)
C
C***********************************************************************
C  TOMAWAC  VERSION 1.0       1/02/93         F MARCOS (LNH) 30 87 72 66
C***********************************************************************
C
C      FONCTION:
C
C   - FIXE, POUR LES "HYPER PRISMES" DE TOMAWAC ET,
C     AVANT LA REMONTEE DES COURBES CARACTERISTIQUES,
C     LES COORDONNEES BARYCENTRIQUES DE TOUS LES NOEUDS DU
C     MAILLAGE DANS L'ELEMENT VERS OU POINTE CETTE COURBE.
C     (ROUTINE INSPIREE DE GSHP41 DE BIEF)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    U,V,T,W     ! -->! COMPOSANTES DU CHAMP CONVECTEUR              !
C !    X,Y         ! -->! COORDONNEES DES POINTS DU MAILLAGE.          !
C !    SHP1-2-3    !<-- ! COORDONNEES BARYCENTRIQUES DES NOEUDS DANS   !
C !                !    ! LEURS ELEMENTS 2D "ELT" ASSOCIES.            !
C !    SHT         !<-- ! COORDONNEES BARYCENTRIQUES SUIVANT Z DES     !
C !                !    ! NOEUDS DANS LEURS ETAGES "ETA" ASSOCIES.     !
C !    ELT         !<-- ! NUMEROS DES ELEMENTS 2D CHOISIS POUR CHAQUE  !
C !                !    ! NOEUD.                                       !
C !    ETA         !<-- ! NUMEROS DES DIREC. CHOISIS POUR CHAQUE NOEUD.!
C !    FRE         !<-- ! NUMEROS DES FREQ. CHOISIES POUR CHAQUE NOEUD.!
C !    XCONV       !<-- ! POSITION INITIALE DES DERIVANT EN X          !
C !    YCONV       !<-- ! POSITION INITIALE DES DERIVANT EN Y          !
C !    TCONV       !<-- ! POSITION INITIALE DES DERIVANT EN TETA       !
C !    FCONV       !<-- ! POSITION INITIALE DES DERIVANT EN F          !
C !    TETA        ! -->! DIRECTIONS DE PROPAGATION                    !
C !    FREQ        ! -->! FREQUENCES DE PROPAGATION                    !
C !    IKLE2       ! -->! TRANSITION ENTRE LES NUMEROTATIONS LOCALE    !
C !                ! -->! ET GLOBALE                                   !
C !    NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE 2D.             !
C !    NELEM2      ! -->! NOMBRE D'ELEMENTS DU MAILLAGE 2D.            !
C !    NPLAN       ! -->! NOMBRE DE DIRECTIONS                         !
C !    NF          ! -->! NOMBRE DE FREQUENCES                         !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : WAC
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER NPOIN2,NELEM2,NPLAN,NF
      INTEGER N1,N2,N3,IPOIN,IELEM,IPLAN,IP,IFF
C
      DOUBLE PRECISION U(NPOIN2,NPLAN),V(NPOIN2,NPLAN)
      DOUBLE PRECISION T(NPOIN2,NPLAN),W(NPOIN2,NPLAN)
      DOUBLE PRECISION TETA(*),FREQ(NF)
      DOUBLE PRECISION X(NPOIN2),Y(NPOIN2)
      DOUBLE PRECISION XCONV(NPOIN2,NPLAN),YCONV(NPOIN2,NPLAN)
      DOUBLE PRECISION TCONV(NPOIN2,NPLAN),FCONV(NPOIN2,NPLAN)
      DOUBLE PRECISION SHP1(NPOIN2,NPLAN),SHP2(NPOIN2,NPLAN)
      DOUBLE PRECISION SHP3(NPOIN2,NPLAN),SHT(NPOIN2,NPLAN)
      DOUBLE PRECISION SHF(NPOIN2,NPLAN)
      DOUBLE PRECISION DET1,DET2,EPS
      INTEGER IFABOR(NELEM2,7),GOODELT(NPOIN2,NPLAN)
C
      INTEGER IKLE2(NELEM2,3),ELT(NPOIN2,NPLAN),ETA(NPOIN2,NPLAN)
      INTEGER FRE(NPOIN2,NPLAN)
C
!      DATA EPS / 1.D-6 /
!      DATA EPS / 1.D-12 /
       DATA EPS /0.d0 /
C-----------------------------------------------------------------------
C
C  INITIALISATION DES POINTS A CONVECTER
C
      GOODELT = 0
      DO 60 IP=1,NPLAN
        DO 90 IPOIN=1,NPOIN2
          XCONV(IPOIN,IP)=X(IPOIN)
          YCONV(IPOIN,IP)=Y(IPOIN)
          TCONV(IPOIN,IP)=TETA(IP)
          FCONV(IPOIN,IP)=FREQ(IFF)
90        CONTINUE
60      CONTINUE
C
C-----------------------------------------------------------------------
        DO 10 IPLAN=1,NPLAN
C
C-----------------------------------------------------------------------
C  REMPLISSAGE INITIAL DES SHP ET DES ELT (POUR LES POINTS DE BORD
C  LATERAUX ON N'EST PAS SUR DE TROUVER UN ELEMENT VERS LEQUEL POINTE
C  -(U,V)).
C
         DO 20 IELEM = 1,NELEM2
C
            N1=IKLE2(IELEM,1)
               ELT(N1,IPLAN) = IELEM
               SHP1(N1,IPLAN) = 1.D0
               SHP2(N1,IPLAN) = 0.D0
               SHP3(N1,IPLAN) = 0.D0
            N1=IKLE2(IELEM,2)
               ELT(N1,IPLAN) = IELEM
               SHP1(N1,IPLAN) = 0.D0
               SHP2(N1,IPLAN) = 1.D0
               SHP3(N1,IPLAN) = 0.D0
            N1=IKLE2(IELEM,3)
               ELT(N1,IPLAN) = IELEM
               SHP1(N1,IPLAN) = 0.D0
               SHP2(N1,IPLAN) = 0.D0
               SHP3(N1,IPLAN) = 1.D0
20       CONTINUE
C
C-----------------------------------------------------------------------
C  REMPLISSAGE ELEMENT PAR ELEMENT DES SHP ET DES ELT DES POINTS DE
C  L'ELEMENT POUR LESQUELS -(U,V) POINTE VERS CET ELEMENT.
C
        DO 450 IELEM=1,NELEM2
C
          N1=IKLE2(IELEM,1)
          N2=IKLE2(IELEM,2)
          N3=IKLE2(IELEM,3)
C
C LOCALIZED ON THE BORDER OF EACH PROC
C ----------------------------------------------
C
          IF ((IFABOR(IELEM,1)==-2)) THEN
             ELT(N1,IPLAN) = IELEM
             SHP1(N1,IPLAN) = 1.D0
             SHP2(N1,IPLAN) = 0.D0
             SHP3(N1,IPLAN) = 0.D0 
             ELT(N2,IPLAN) = IELEM
             SHP1(N2,IPLAN) = 0.D0
             SHP2(N2,IPLAN) = 1.D0
             SHP3(N2,IPLAN) = 0.D0 
          ENDIF
C
          IF ((IFABOR(IELEM,2)==-2)) THEN
             ELT(N2,IPLAN) = IELEM
             SHP1(N2,IPLAN) = 0.D0
             SHP2(N2,IPLAN) = 1.D0
             SHP3(N2,IPLAN) = 0.D0 
             ELT(N3,IPLAN) = IELEM
             SHP1(N3,IPLAN) = 0.D0
             SHP2(N3,IPLAN) = 0.D0
             SHP3(N3,IPLAN) = 1.D0
          ENDIF
C
          IF ((IFABOR(IELEM,3)==-2)) THEN
             ELT(N3,IPLAN) = IELEM
             SHP1(N3,IPLAN) = 0.D0
             SHP2(N3,IPLAN) = 0.D0
             SHP3(N3,IPLAN) = 1.D0
             ELT(N1,IPLAN) = IELEM
             SHP1(N1,IPLAN) = 1.D0
             SHP2(N1,IPLAN) = 0.D0
             SHP3(N1,IPLAN) = 0.D0 
          ENDIF 
C
450       CONTINUE


        DO 50 IELEM=1,NELEM2
          N1=IKLE2(IELEM,1)
          N2=IKLE2(IELEM,2)
          N3=IKLE2(IELEM,3)
C
C DET1 = (NINI+1,UNILAG)  DET2 = (UNILAG,NINI-1)
C ----------------------------------------------
C
      DET1=(X(N2)-X(N1))*V(N1,IPLAN)-(Y(N2)-Y(N1))*U(N1,IPLAN)
      DET2=(Y(N3)-Y(N1))*U(N1,IPLAN)-(X(N3)-X(N1))*V(N1,IPLAN)
          IF (DET1.LE.EPS.AND.DET2.LE.EPS) THEN
             ELT(N1,IPLAN) = IELEM
             SHP1(N1,IPLAN) = 1.D0
             SHP2(N1,IPLAN) = 0.D0
             SHP3(N1,IPLAN) = 0.D0
             GOODELT(N1,IPLAN) = 1
          ENDIF
C
      DET1=(X(N3)-X(N2))*V(N2,IPLAN)-(Y(N3)-Y(N2))*U(N2,IPLAN)
      DET2=(Y(N1)-Y(N2))*U(N2,IPLAN)-(X(N1)-X(N2))*V(N2,IPLAN)
          IF (DET1.LE.EPS.AND.DET2.LE.EPS) THEN
             ELT(N2,IPLAN) = IELEM
             SHP1(N2,IPLAN) = 0.D0
             SHP2(N2,IPLAN) = 1.D0
             SHP3(N2,IPLAN) = 0.D0
             GOODELT(N2,IPLAN) = 1
          ENDIF
C
      DET1=(X(N1)-X(N3))*V(N3,IPLAN)-(Y(N1)-Y(N3))*U(N3,IPLAN)
      DET2=(Y(N2)-Y(N3))*U(N3,IPLAN)-(X(N2)-X(N3))*V(N3,IPLAN)
          IF (DET1.LE.EPS.AND.DET2.LE.EPS) THEN
             ELT(N3,IPLAN) = IELEM
             SHP1(N3,IPLAN) = 0.D0
             SHP2(N3,IPLAN) = 0.D0
             SHP3(N3,IPLAN) = 1.D0
             GOODELT(N3,IPLAN) = 1
          ENDIF 
C
50       CONTINUE
         DO 230 IELEM = 1,NELEM2
            N1=IKLE2(IELEM,1)
            N2=IKLE2(IELEM,2)
            N3=IKLE2(IELEM,3)
          IF (IFABOR(IELEM,1)==0) GOODELT(N1,IPLAN)=
     *                                        GOODELT(N1,IPLAN)+10
          IF (IFABOR(IELEM,1)==0) GOODELT(N2,IPLAN)=
     *                                        GOODELT(N2,IPLAN)+10
          IF (IFABOR(IELEM,2)==0) GOODELT(N2,IPLAN)=
     *                                        GOODELT(N2,IPLAN)+10
          IF (IFABOR(IELEM,2)==0) GOODELT(N3,IPLAN)=
     *                                        GOODELT(N3,IPLAN)+10
          IF (IFABOR(IELEM,3)==0) GOODELT(N3,IPLAN)=
     *                                        GOODELT(N3,IPLAN)+10
          IF (IFABOR(IELEM,3)==0) GOODELT(N1,IPLAN)=
     *                                        GOODELT(N1,IPLAN)+10
          IF (IFABOR(IELEM,1)==-1) GOODELT(N1,IPLAN)=
     *                                        GOODELT(N1,IPLAN)+100
          IF (IFABOR(IELEM,1)==-1) GOODELT(N2,IPLAN)=
     *                                        GOODELT(N2,IPLAN)+100
          IF (IFABOR(IELEM,2)==-1) GOODELT(N2,IPLAN)=
     *                                        GOODELT(N2,IPLAN)+100
          IF (IFABOR(IELEM,2)==-1) GOODELT(N3,IPLAN)=
     *                                        GOODELT(N3,IPLAN)+100
          IF (IFABOR(IELEM,3)==-1) GOODELT(N3,IPLAN)=
     *                                        GOODELT(N3,IPLAN)+100
          IF (IFABOR(IELEM,3)==-1) GOODELT(N1,IPLAN)=
     *                                        GOODELT(N1,IPLAN)+100
          IF (IFABOR(IELEM,1)==-2) GOODELT(N1,IPLAN)=
     *                                       GOODELT(N1,IPLAN)+1000
          IF (IFABOR(IELEM,1)==-2) GOODELT(N2,IPLAN)=
     *                                       GOODELT(N2,IPLAN)+1000
          IF (IFABOR(IELEM,2)==-2) GOODELT(N2,IPLAN)=
     *                                       GOODELT(N2,IPLAN)+1000
          IF (IFABOR(IELEM,2)==-2) GOODELT(N3,IPLAN)=
     *                                       GOODELT(N3,IPLAN)+1000
          IF (IFABOR(IELEM,3)==-2) GOODELT(N1,IPLAN)=
     *                                       GOODELT(N1,IPLAN)+1000
          IF (IFABOR(IELEM,3)==-2) GOODELT(N3,IPLAN)=
     *                                       GOODELT(N3,IPLAN)+1000
230      CONTINUE



C
C-----------------------------------------------------------------------
C  REMPLISSAGE POINT PAR POINT DES SHT ETA SHF ET FRE
C
         DO 70 IPOIN=1,NPOIN2
C
           IF (T(IPOIN,IPLAN).GT.0.D0) THEN
            IF (IPLAN.EQ.1) THEN
                ETA(IPOIN,1) = NPLAN
                SHT(IPOIN,1) = 1.D0
                TCONV(IPOIN,1)=TETA(NPLAN+1)
              ELSE
                ETA(IPOIN,IPLAN) = IPLAN-1
                SHT(IPOIN,IPLAN) = 1.D0
              ENDIF 
           ELSE
              ETA(IPOIN,IPLAN) = IPLAN
              SHT(IPOIN,IPLAN) = 0.D0
           ENDIF 
C
           IF (((W(IPOIN,IPLAN).GT.0.D0).AND.(IFF.NE.1)).OR.
     *             (IFF.EQ.NF)) THEN
            FRE(IPOIN,IPLAN) = IFF-1
            SHF(IPOIN,IPLAN) = 1.D0
           ELSE
            FRE(IPOIN,IPLAN) = IFF
            SHF(IPOIN,IPLAN) = 0.D0
           ENDIF
            
70       CONTINUE
C
C-----------------------------------------------------------------------
C
10      CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
