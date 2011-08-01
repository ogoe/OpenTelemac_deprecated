C                       *****************
                        SUBROUTINE LECUTI
C                       *****************
C
     *(X,Y,NPOIN,NCOU,XRELV,YRELV,UR,VR,TRA01,NP,NPMAX)
C
C***********************************************************************
C  COWADIS VERSION 1.0    30/08/95        F.MARCOS     (LNH) 30 87 72 66
C***********************************************************************
C
C   FONCTION : CE SOUS-PROGRAMME PERMET LA LECTURE DES COURANTS
C               OU DES VENTS A UN FORMAT DEFINI PAR L'UTILISATEUR
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    X,Y         ! -->!  COORDONNEES DU MAILLAGE                     !
C !    NPOIN       ! -->!  NOMBRE DE POINTS DU MAILLAGE                !
C !    NCOU        ! -->!  NUMERO D'UNITE LOGIQUE DU FICHIER DES COURA.!
C !    XRELV       !<-->!  TABLEAU DES ABSCISSES DES POINTS RELEVES    !
C !    YRELV       !<-->!  TABLEAU DES ORDONNEES DES POINTS RELEVES    !
C !    UR,VR       !<-->!  TABLEAU DES COURANTS OU VENTS RELEVES       !
C !    TRA01       !<-->!  TABLEAU DE TRAVAIL DE DIMENSION NPMAX       !
C !    NP          !<-- !  NOMBRE DE POINTS RELEVES                    !
C !    NPMAX       ! -->!  NOMBRE DE POINTS RELEVES MAXIMUM            !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : LECDON
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/ LNG,LU
C
      INTEGER NPMAX,NP
C
      INTEGER NCOU,NPOIN
C
      DOUBLE PRECISION X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION XRELV(NPMAX),YRELV(NPMAX),UR(NPMAX),VR(NPMAX)
      DOUBLE PRECISION TRA01(NPMAX)
C
C-----------------------------------------------------------------------
C     PROGRAMMER ICI LA LECTURE DE VOTRE FICHIER DE COURANT OU DE VENT
C     LES DONNEES LUES SERONT STOCKEES DANS LES TABLEAUX
C     XRELV, YRELV, UR ET VR
C     L'INTERPOLATION SUR LE MAILLAGE DE CALCUL EST FAITE 
C     AUTOMATIQUEMENT A LA SORTIE DE LA PROCEDURE
C      (PENSEZ A RENSEIGNER NP POUR L'INTERPOLATION ET
C       A RETIRER LES LIGNES FORTRAN SUIVANTES)
C-----------------------------------------------------------------------
C     PROGRAM HERE THE READING OF YOUR CURRENT OR WIND VELOCITIES FILE
C     DATA READ WILL BE HOLD IN ARRAYS
C     XRELV, YRELV, UR AND VR
C     INTERPOLATION ON COMPUTATION MESH IS DONE 
C     AUTOMATICALLY AFTER THE END OF THE SUBROUTINE
C      (THINK TO INITIATE NP FOR THE INTERPOLATION AND
C       ERASE THE FOLLOWING FORTRAN LINES)
C-----------------------------------------------------------------------
C
      WRITE(LU,*)'***********************************************'
      IF (LNG.EQ.1) THEN
        WRITE(LU,*)'  VOUS FAITES APPEL A LA PROCEDURE LECUTI      '
        WRITE(LU,*)'  (FORMAT DU FICHIER DES COURANTS OU VENTS = 4)'
        WRITE(LU,*)'     MAIS VOUS NE L''AVEZ PAS MODIFIEE         '
      ELSE
        WRITE(LU,*)'     YOU CALL SUBROUTINE LECUTI                '
        WRITE(LU,*)'  (FORMAT OF CURRENT OF WIND FILE = 4)         '
        WRITE(LU,*)'     BUT YOU DID NOT CHANGE IT                 '
      ENDIF
      WRITE(LU,*)'***********************************************'
      CALL PLANTE(0)
      
      RETURN
      END
