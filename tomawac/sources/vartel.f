C                       *****************
                        SUBROUTINE VARTEL
C                       *****************
C
     *( VAR, X, Y, DEPTH, UC, VC, ZREPOS, TRA01, F, NPLAN, NF, NPOIN2) 
C
C***********************************************************************
C  TOMAWAC VERSION 1.0    09/06/95       F. MARCOS (LNH) 30 87 72 66
C***********************************************************************
C
C     FONCTION  : PERMET L'UTILISATION D'UNE VARIABLE RECUPEREE
C                 DANS UN FICHIER TELEMAC 
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    VAR         ! -->! VARIABLE RECUPEREE DANS LE FICHIER TELEMAC   !
C !    X,Y         ! -->! COORDONNEES DES POINTS DU MAILLAGE 2D        !
C !    DEPTH       !<-->! HAUTEUR D'EAU                                !
C !    UC,VC       !<-->! CHAMPS DE COURANT                            !
C !    ZREPOS      !<-->! COTE INITIALE DU PLAN D'EAU AU REPOS         !
C !    TRA01       !<-->! TABLEAU DE TRAVAIL                           !
C !    F           !<-->! SPECTRE DE VARIANCE                          !
C !    NPLAN       ! -->! NOMBRE DE DIRECTIONS                         !
C !    NF          ! -->! NOMBRE DE FREQUENCES                         !
C !    NPOIN2      ! -->! NOMBRE DE POINTS 2D                          !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : WAC
C
C  SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER  NPOIN2,NPLAN,NF
C
      DOUBLE PRECISION F (NPOIN2,NPLAN,NF) , TRA01(NPOIN2)
      DOUBLE PRECISION X (NPOIN2) , Y (NPOIN2) , DEPTH(NPOIN2)
      DOUBLE PRECISION UC(NPOIN2) , VC(NPOIN2) , VAR(NPOIN2)
      DOUBLE PRECISION ZREPOS
C
C-----------------------------------------------------------------------
C
C     UTILISER ICI LA VARIABLE VAR POUR INITIALISER LES
C              TABLEAUX VOULUS
C     PAR EXEMPLE :
C
C     CALL OV( 'X=Y     ' , DEPTH , VAR , X , 0.D0 , NPOIN2)
C
      RETURN
      END
