C                       *****************
                        SUBROUTINE DISMOY
C                       *****************
C
     *(NPOIN,NELEM,X,Y,IKLE,K,LISHHO)
C
C***********************************************************************
C
C  ARTEMIS VERSION 5.1   04/06/99    D. AELBRECHT (LNH) 01 30 87 74 12
C                                    P. THELLIER  (LNH) Merci Paul
C
C  LINKED TO BIEF VERS. 5.0          J-M HERVOUET (LNH) 01 30 87 80 18
C
C***********************************************************************
C
C     FONCTION  : CALCULE LE NOMBRE DE LISSAGES A PRIORI NECESSAIRE
C                 A OPERER SUR LA HAUTEUR DE HOULE (LISHHO)
C                 POUR FILTRER LES OSCILLATIONS PARASITES EN HOULE 
C                 REGULIERE, DETERMINE A PARTIR D'UNE DISTANCE
C                 MOYENNE ENTRE POINTS ET DU NOMBRE DE NOEUDS MOYEN
C                 PAR DEMI-LONGUEUR D'ONDE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     LISHHO     |<-- |  NOMBRE DE LISSAGES SUR HHO                  |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : ARTEMI
C
C  SOUS-PROGRAMME APPELE : OS, VECTOR
C
C***********************************************************************
C
      USE INTERFACE_ARTEMIS, EX_DISMOY => DISMOY
C                   
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C
C variables pour le calcul des distances moyennes
C
      INTEGER ng, isom, ielem
      DOUBLE PRECISION som, d1, d2, dmoy, pi, somd
C
C
      INTEGER NPOIN,NELEM,LISHHO
      INTEGER IKLE(NELEM,*)
C
      DOUBLE PRECISION X(NPOIN),Y(NPOIN),K(NPOIN)
C
C-----------------------------------------------------------------------
C
      pi = 3.1415926535897932384626433D0                    
c
c     calcul de la distance moyenne entre un point et ses voisins
c     -----------------------------------------------------------
c
      somd = 0.d0
      do 300 ng=1,npoin
         isom = 0
         som = 0.d0
c
         do 250 ielem=1,nelem
c
            if (ikle(ielem,1).eq.ng) then
c           --                       ---
c
c            pt 1 commun, calcul de la distance des pts 2 et 3
             d1 = sqrt ( (x(ng)-x(ikle(ielem,2))) **2.d0 +
     *                   (y(ng)-y(ikle(ielem,2))) **2.d0 )
             d2 = sqrt ( (x(ng)-x(ikle(ielem,3))) **2.d0 +
     *                   (y(ng)-y(ikle(ielem,3))) **2.d0 )
             som = som + d1 + d2
             isom = isom + 2
c                   
             elseif (ikle(ielem,2).eq.ng) then
c            ------                       ----
c
c            pt 2 commun, calcul de la distance des pts 1 et 3
             d1 = sqrt ( (x(ng)-x(ikle(ielem,1))) **2.d0 +
     *                   (y(ng)-y(ikle(ielem,1))) **2.d0 )
             d2 = sqrt ( (x(ng)-x(ikle(ielem,3))) **2.d0 +
     *                   (y(ng)-y(ikle(ielem,3))) **2.d0 )
             som = som + d1 + d2
             isom = isom + 2
c
             elseif (ikle(ielem,3).eq.ng) then
c            ------                       ----
c
c            pt 3 commun, calcul de la distance des pts 1 et 2
             d1 = sqrt ( (x(ng)-x(ikle(ielem,1))) **2.d0 +
     *                   (y(ng)-y(ikle(ielem,1))) **2.d0 )
             d2 = sqrt ( (x(ng)-x(ikle(ielem,2))) **2.d0 +
     *                   (y(ng)-y(ikle(ielem,2))) **2.d0 )
             som = som + d1 + d2
             isom = isom + 2
c
            endif
c           -----
c
 250     continue
c
         dmoy = som / float(isom)
         somd = somd + (pi / (K(ng)*dmoy))
c
 300  continue
c
c     calcul du nombre de lissage en fonction de la moyenne des distances
c     -------------------------------------------------------------------
c
      LISHHO = int((somd/float(npoin))) * 10
c
      RETURN
      END
