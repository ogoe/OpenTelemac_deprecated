C                       *****************
                        SUBROUTINE MARAST
C                       *****************
C
     *(MARDAT,MARTIM,PHI0,NPOIN,AT,FU1,FV1,X,SINLAT,COSLAT,GRAV)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.2    01/03/94      E. DAVID    (LHF) 76 33 42 08
C                                         F LEPEINTRE (LNH) 30 87 78 54
C                                         J-M JANIN   (LNH) 30 87 72 84
C***********************************************************************
C
C      FONCTION:
C      =========
C
C      CALCULE LA FORCE DE LA MAREE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C | MARDAT         | -->| TABLEAU CONTENANT LES INFOS DE DATE          |
C | MARTIM         | -->| TABLEAU CONTENANT LES INFOS DE TEMPS (TU)    |
C | PHI0           | -->| LONGITUDE DU POINT ORIGINE                   |
C | NPOIN          | -->| NOMBRE DE POINTS DU MAILLAGE                 |
C | AT             | -->| TEMPS
C | FU1,FV1        |<-- | FORCES GENERATRICES CALCULEES                |
C | X              | -->| ABSCISSES DU MAILLAGE                        |
C | SINLAT         | -->| SIN DE LA LATITUDE EN COORD SPHERIQUE        |
C | COSLAT         | -->| COS DE LA LATITUDE EN COORD SPHERIQUE        |
C | GRAV           | -->| PESANTEUR                                    |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C SOUS-PROGRAMME APPELE PAR : PROSOU
C SOUS-PROGRAMMES APPELES   : ASTRO, TSLOC
C
C***********************************************************************
C
      USE INTERFACE_TELEMAC2D, EX_MARAST => MARAST
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: MARDAT(3),MARTIM(3),NPOIN
      DOUBLE PRECISION, INTENT(INOUT) :: FU1(NPOIN),FV1(NPOIN)
      DOUBLE PRECISION, INTENT(IN) :: COSLAT(NPOIN),SINLAT(NPOIN)
      DOUBLE PRECISION, INTENT(IN) :: X(NPOIN),GRAV,AT,PHI0
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER YEAR,MONTH,DAY,HOUR,MIN,SEC,I
C     
      DOUBLE PRECISION ARL,ARS,DL,DS,AL,AS
      DOUBLE PRECISION RT,LONG,AHL,AHS,MLT,MST,LONG0
      DOUBLE PRECISION F0L,F0S,FXL,FYL,FXS,FYS
      DOUBLE PRECISION TLOC,TSLOC,TLOC1
      DOUBLE PRECISION K2,H2
C
      INTRINSIC ACOS,COS,SIN
C
C-----------------------------------------------------------------------
C
C ECLATEMENT DES TABLEAUX MARDAT ET MARTIM
C
      YEAR  = MARDAT(1)
      MONTH = MARDAT(2)
      DAY   = MARDAT(3)
      HOUR  = MARTIM(1)
      MIN   = MARTIM(2)
      SEC   = MARTIM(3)
C
C RAPPEL : HOUR EN TEMPS UNIVERSEL |
C REMARQUE GENERALE : T=TERRE, L=LUNE , S=SOLEIL
C
C LONG0 : LONGITUDE DE REFERENCE EN RADIAN (0,2PI)
C
      LONG0=PHI0*ACOS(-1.D0)/180.D0
C
C APPEL FONCTION PRINCIPALE CALCULANT LES ANGLES LUNAIRES ET SOLAIRES
C
      CALL ASTRO(YEAR,MONTH,DAY,HOUR,MIN,SEC,AT,ARL,ARS,DL,DS,AL,AS)
C
C RT : RAYON DE LA TERRE EN M
C
      RT   = 6378000.D0
C
C RAPPORT MASSE LUNE SUR MASSE TERRE
C
      MLT  = 1.D0 / 81.53D0
C
C RAPPORT MASSE SOLEIL SUR MASSE TERRE
C
      MST  = 331954.D0
C
C AMPLITUDE DE LA FORCE INDUITE PAR LE POTENTIEL ASTRONOMIQUE POUR :
C
C     - LA LUNE
C
      F0L  = GRAV * MLT * ARL**2
C
C     - LE SOLEIL
C
      F0S  = GRAV * MST * ARS**2
C
C TEMPS SIDERAL
C
      TLOC1 = TSLOC(YEAR,MONTH,DAY,HOUR,MIN,SEC,AT)
C
      DO 10 I=1,NPOIN
C
C LONGITUDE DU NOEUD CONSIDERE
C
        LONG = X(I)/RT+LONG0
C
C TEMPS SIDERAL LOCAL
C
        TLOC = TLOC1 + LONG
C
C ANGLE HORAIRE DE LA LUNE
C
        AHL  = TLOC - AL
C
C ANGLE HORAIRE DU SOLEIL
C
        AHS  = TLOC - AS
C
C FORCE INDUITE PAR LE POTENTIEL ASTRONOMIQUE SEUL
C
C    FXL : FORCE SELON X DE LA LUNE
C     Y  : SELON Y
C     S  : IDEM POUR LE SOLEIL
C
        FXL  = F0L * COS(DL) * SIN(AHL) *
     *         ( ( 1.D0-2*ARL*(SINLAT(I)*SIN(DL)+
     *           COSLAT(I)*COS(DL)*COS(AHL))+ARL*ARL )**(-1.5D0) -1.D0 )
C
        FXS  = F0S * COS(DS) * SIN(AHS) *
     *         ( ( 1.D0-2*ARS*(SINLAT(I)*SIN(DS)+
     *           COSLAT(I)*COS(DS)*COS(AHS))+ARS*ARS )**(-1.5D0) -1.D0 )
C
        FYL  = F0L*(COSLAT(I)*SIN(DL)-SINLAT(I)*COS(DL)*COS(AHL))*
     *         ( ( 1.D0-2*ARL*(SINLAT(I)*SIN(DL)+
     *           COSLAT(I)*COS(DL)*COS(AHL))+ARL*ARL )**(-1.5D0) -1.D0 )
C
        FYS  = F0S*( COSLAT(I)*SIN(DS)-SINLAT(I)*COS(DS)*COS(AHS))*
     *         ( ( 1.D0-2*ARS*(SINLAT(I)*SIN(DS)+
     *           COSLAT(I)*COS(DS)*COS(AHS))+ARS*ARS )**(-1.5D0) -1.D0 )

C
C PRISE EN COMPTE DE :
C
C    - LA MAREE TERRESTRE (NOMBRE DE LOVE H2)
C
        H2=0.61D0
C
C    - LES PERTUBATION STATIQUES (NOMBRE DE LOVE K2)
C
        K2=0.30D0
C
C MANQUE :
C
C    - PERTUBATION DYNAMIQUE D'AUTO-ATTRACTION
C    - PERTUBATION DYNAMIQUE DES EFFETS DE CHARGES
C
C FORCE FINALE
C
        FU1(I)=FU1(I)+(1.D0+K2-H2)*(FXL+FXS)
        FV1(I)=FV1(I)+(1.D0+K2-H2)*(FYL+FYS)
C
10    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
