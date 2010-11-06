C                        ********************
                         SUBROUTINE INIT_ZERO
C                        ********************
C 
C
C***********************************************************************
C SISYPHE VERSION 5.9                             C. VILLARET    
C                                                
C 
C 04/06/2008 : JMH INITIALISATION DE MASDEP
C                                                
C*********************************************************************** 
C
C     FONCTION  : INITIALISATIONS
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |                |--> | 
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C     - PROGRAMME APPELANT : SISYPHE  
C     - SOUS-PROGRAMMES APPELES : OS
C
C***********************************************************************
C
      USE DECLARATIONS_SISYPHE
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C-----------------------------------------------------------------------
C
      INTEGER I
C
C-----------------------------------------------------------------------
C
C========================================================================
C                         INITIALISATIONS 
C =======================================================================
C
C     DES TABLEAUX de TRAVAIL
C
      IF(NPRIV > 0) CALL OS ('X=0     ', X=PRIVE)
      CALL OS('X=0     ', X=T12   )   
      CALL OS('X=0     ', X=COEFPN)
C
C     INITIALISATION DES VARIABLES SEDIMENTAIRES :
C
      CALL OS('X=0     ', X=QS)
      CALL OS('X=0     ', X=QSX)
      CALL OS('X=0     ', X=QSY)
      CALL OS('X=0     ', X=QSCLXC )
      CALL OS('X=0     ', X=QSCLYC )
      CALL OS('X=0     ', X=QSCLXS )
      CALL OS('X=0     ', X=QSCLYS )
C
C 7 FOLLOWING LINES ADDED BY JMH 22/04/2005
C PROVISIONAL INITIALISATION FOR FIRST OUTPUT IN RESULTS FILE
C
      CALL OS('X=0     ', X=QSCL )  ! BLOCK OF SIZE NSICLA
      CALL OS('X=0     ', X=QS_S )
      CALL OS('X=0     ', X=QSXS )
      CALL OS('X=0     ', X=QSYS )
      CALL OS('X=0     ', X=QS_C )
      CALL OS('X=0     ', X=QSXC )
      CALL OS('X=0     ', X=QSYC )
C
C FOLLOWING LINE ADDED BY JMH 04/05/2005
C PROBABLY USEFUL ONLY IF(CHARR) AND WITH FINITE ELEMENTS
C
      CALL OS('X=0     ', X=ZFCL_C )
C
C     DEPOSITION MASSES FOR EVERY CLASS IN SUSPENSION
C
      DO I=1,NSICLA
        MASDEP(I)=0.D0
      ENDDO
C
      CALL OS('X=0     ', X=E    )
      CALL OS('X=0     ', X=ESOMT)
      CALL OS('X=0     ', X=CS   )
C
C---- INITIALISATION DES VARIABLES  HYDRODYNAMIQUES :
C
      CALL OS('X=0     ', X=QU )
      CALL OS('X=0     ', X=QV )
      CALL OS('X=0     ', X=U2D )
      CALL OS('X=0     ', X=V2D )
      CALL OS('X=0     ', X=HN )
      CALL OS('X=0     ', X=Q  )
      CALL OS('X=0     ', X=TOB)
C
C---- PARAMETRE DE HOULE SI BESOIN
C
C     TOUTES LES INITIALISATIONS DE LA HOULE SONT A SUPPRIMER
C     QUAND TOUS LES CONTROLES SERONT FAITS
C     VOIR BEDLOAD_BAILARD,DIBWAT,BIJKER ET SOULSBY
C
C
C     FW=0.3 CORRESPONDS TO NO WAVES, 0 WOULD DO THEN A LOG(0) 
      CALL OS('X=C     ', X=FW ,C=0.3D0   )   ! 
      CALL OS('X=0     ', X=HW    )   ! 
      CALL OS('X=0     ', X=TW    )   ! 
      CALL OS('X=C     ', X=THETAW, C=90.D0)  ! 
      CALL OS('X=0     ', X=UW    )   ! 
      CALL OS('X=0     ', X=TOBW)     ! 
C
C-----------------------------------------------------------------------
C
      RETURN
      END
