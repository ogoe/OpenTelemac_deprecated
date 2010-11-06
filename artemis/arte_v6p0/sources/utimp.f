C                       ****************
                        SUBROUTINE UTIMP
C                       ****************
C
     *(PHIR,PHII,C,CG,K,X,Y,ZF,H,
     * HHO,U0,V0,PHAS,S,TRA01,TRA02,TRA03,TRA04,INCI,
     * GRAV,PER,OMEGA,IKLE,NBOR,KP1BOR,
     * NELEM,NELMAX,IELM,IELMB,NPTFR,NPOIN,PRIVE)
C
C***********************************************************************
C
C  ARTEMIS VERSION 5.1   04/06/99    D. AELBRECHT (LNH) 01 30 87 74 12
C
C  LINKED TO BIEF VERS. 5.0          J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C      FONCTION:    PROGRAMME UTILISATEUR OU IL DISPOSE
C                   D'A PEU PRES TOUS LES ARGUMENTS DE SON CALCUL POUR
C                   FAIRE DES IMPRESSIONS, CALCULER DES SOLUTIONS
C                   ANALYTIQUES...(NE FAIT RIEN EN STANDARD)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    PHIR,PHII   | -->|  COMPOSANTES DU POTENTIEL                    |
C |    C,CG        | -->|  VITESSE DE PHASE ET DE GROUPE               |
C |      K         | -->|  NOMBRE D'ONDE                               |
C |      X,Y       | -->|  COORDONNEES DES POINTS DU MAILLAGE          |
C |      ZF        | -->|  COTE DU FOND                                |
C |      H         | -->|  HAUTEUR D'EAU AU REPOS                      |
C |      HHO       | -->|  HAUTEUR DE LA HOULE                         |
C |      U0,V0     | -->|  VITESSES EN SURFACE (A T=0)                 |
C |      PHAS      | -->|  PHASE DE LA HOULE                           |
C |      S         | -->|  COTE DE LA SURFACE LIBRE                    |
C |    TRA01,...,4 |<-->|  TABLEAUX DE TRAVAIL                         |
C |      INCI      | -->|  INCIDENCE DE LA HOULE                       |
C |     GRAV       | -->|  GRAVITE                                     |
C |     PER        | -->|  PERIODE DE LA HOULE                         |
C |     OMEGA      | -->|  PULSATION DE LA HOULE                       |
C |     IKLE       | -->|  TABLE DE CONNECTIVITE                       |
C |     NBOR       | -->|  ADRESSES DES POINTS DE BORD                 |
C |     KP1BOR     | -->|  NUMERO DE BORD DU POINT SUIVANT             |
C |     W1         |<-->|  TABLEAU DE TRAVAIL                          |
C |     NELEM      | -->|  NOMBRE D'ELEMENTS                           |
C |     IELM       | -->|  TYPE D'ELEMENT                              |
C |     IELMB      | -->|  TYPE D'ELEMENT DE BORD                      |
C |     NPTFR      | -->|  NOMBRE DE POINTS FRONTIERE                  |
C |     NPOIN      | -->|  NOMBRE DE POINTS                            |
C |     PRIVE      |<-->|  TABLEAUX RESERVE A L'UTILISATEUR            |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C**********************************************************************
C
      USE INTERFACE_ARTEMIS, EX_UTIMP=> UTIMP 
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER NELEM,NELMAX,IELM,IELMB,NPTFR,NPOIN
      INTEGER IKLE(NELMAX,3),NBOR(NPTFR),KP1BOR(NPTFR)
C
      DOUBLE PRECISION PHIR(NPOIN),PHII(NPOIN)
      DOUBLE PRECISION C(NPOIN),CG(NPOIN),K(NPOIN),X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION ZF(NPOIN),H(NPOIN),HHO(NPOIN),U0(NPOIN),V0(NPOIN)
      DOUBLE PRECISION INCI(NPOIN)
      DOUBLE PRECISION PHAS(NPOIN),S(NPOIN)
      DOUBLE PRECISION TRA01(NPOIN),TRA02(NPOIN)
      DOUBLE PRECISION TRA03(NPOIN),TRA04(NPOIN)
C
      TYPE(BIEF_OBJ) :: PRIVE
C
      DOUBLE PRECISION GRAV,PER,OMEGA
C
C------------------------------------------------------------------
C     EXEMPLE : ON MET LES ANCIENNES VARIABLES U1 ET V1
C               (VITESSES HORIZONTALES A TEMPS = T/4)
C               DANS PRIVE(I,1) ET PRIVE(I,2)
C------------------------------------------------------------------
C
C      CALL VECTOR(TRA02, '=' , 'GRADF          X' , IELM ,
C     *            1.D0 , PHII , BID , BID , BID , BID , BID ,
C     *            MESH , MESH , MSK , MASKEL )
C
C      CALL VECTOR(TRA03 , '=' , 'GRADF          Y' , IELM ,
C     *            1.D0 , PHII , BID , BID , BID , BID , BID ,
C     *            MESH , MESH , MSK , MASKEL )
C     *            MESH , XMESH ,
C
C      CALL VECTOR(TRA01 , '=' , 'MASBAS          ' , IELM ,
C     *            1.D0 , BID , BID , BID , BID , BID , BID ,
C     *            MESH , MESH , MSK , MASKEL )
C     *            MESH , XMESH ,
C
C      CALL OS( 'X=Y/Z   ' , TRA02 , TRA02 , TRA01 , BID )
C      CALL OS( 'X=Y/Z   ' , TRA03 , TRA03 , TRA01 , BID )
C
C      DO 25 I = 1,NPOIN
C         PRIVE%ADR(1)%P%R(1) = TRA02(I)
C         PRIVE%ADR(1)%P%R(2) = TRA03(I)
C 25   CONTINUE
C
      RETURN
      END
