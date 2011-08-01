C                       ****************                               
                        SUBROUTINE VFCFL                              
C                       ****************                               
C                                                                       
     *(NUBO,VNOIN,NSEG,NPOIN,X,Y,G,H,EPS,QU,QV,DT,CFLWTD,DTVARI,LISTIN)                   
C                                                                       
C***********************************************************************
C TELEMAC 2D VERSION 5.2     18/03/98            N.GOUTAL           
C
C-----------------------------------------------------------------------
C                                                                       
C  FONCTION  : . OPTIMISATION DU PAS DE TEMPS                                        
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C ! . VNOIN        ! -->! NORMALE DU SEGMENT INTERNE                   !
C !                !    ! (2 PREMIERES COMPOSANTES) ET                 !
C !                !    ! LONGUEUR DE CE SEGMENT (3IEME COMPOSANTE)    !
C !    UI,VI       ! -->!  COMPOSANTES DE LA VITESSE DE LA             !         
C !                !    !  CELLULE i                                   !
C !    X,Y         ! -->!  COORDONNEES DES POINTS DU MAILLAGE          !
C !    H           ! -->!  HAUTEUR AU TEMPS N                          !
C !________________!____!______________________________________________!
C !   /COMMON/     !    !                                              !
C ! . NSEG         ! -->! NOMBRE TOTAL DE SEGMENTS DU MAILLAGE         !
C !________________!____!______________________________________________!
C***********************************************************************
C                                                                       
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C   
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NSEG,NPOIN
      INTEGER, INTENT(IN) :: NUBO(2,*)
      LOGICAL, INTENT(IN) :: DTVARI,LISTIN
      DOUBLE PRECISION, INTENT(IN) :: X(NPOIN),Y(NPOIN),H(NPOIN)
      DOUBLE PRECISION, INTENT(IN) :: QU(NPOIN),QV(NPOIN),VNOIN(3,*)
      DOUBLE PRECISION, INTENT(IN) :: EPS,G,CFLWTD
      DOUBLE PRECISION, INTENT(INOUT) :: DT  
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                       
      INTEGER ISEGIN,IEL1,IEL2                                                                            
C                                                                       
      DOUBLE PRECISION W1,DIJ,UI,VI,UN,XNC,UNI,UNJ,UJ,VJ,HI,HJ,HMAX                                                        
C                                                                       
      INTRINSIC SQRT,MAX,ABS                                            
C                                                                       
C-----------------------------------------------------------------------
C                                                                      
C------
C 1. INITIALISATION
C------
C
      XNC = 0.D0
C
C------                                                                 
C 2. BOUCLE SUR LES SEGMENTS INTERIEURS                                                    
C------                                                                 
C                                                
         
      DO 100 ISEGIN = 1 , NSEG                                        
C                                                                       
         IEL1 = NUBO(1,ISEGIN)                                          
         IEL2 = NUBO(2,ISEGIN)                                          
C                                                                       
C                                                                       
C   --->    CALCUL DE DIJ                             
C           -------------                           
C                                                                             
       DIJ = SQRT((X(IEL1)-X(IEL2))**2+(Y(IEL1)-Y(IEL2))**2)                      
C                                                                       
C                                                                       
C   --->    CALCUL DES COMPOSANTES DE LA VITESSE INTERNE                        
C           --------------------------------------------                        
C
       HI = H(IEL1)
       IF (H(IEL1).GE.EPS) THEN
          UI = QU(IEL1)/H(IEL1)
          VI = QV(IEL1)/H(IEL1)
       ELSE
          UI = 0.D0
          VI = 0.D0
       ENDIF                                                                    
C                           
C
       IF (H(IEL2).GE.EPS) THEN
          HJ = H(IEL2)
          UJ = QU(IEL2)/H(IEL2)
          VJ = QV(IEL2)/H(IEL2)
       ELSE
          HJ = 0.D0
          UJ = 0.D0
          VJ = 0.D0
       ENDIF                                                                    
C       
C
C   --->    PROJECTION DU VECTEUR VITESSE SUR LA NORMALE
C           --------------------------------------------
C
          UNI = ABS(UI*VNOIN(1,ISEGIN)+VI*VNOIN(2,ISEGIN))
          UNJ = ABS(UJ*VNOIN(1,ISEGIN)+VJ*VNOIN(2,ISEGIN))
          UN = (UNI+UNJ)/2.D0
          HMAX = MAX(HJ,HI)
C
C   --->    CALCUL DU NOMBRE DE COURANT
C           ---------------------------
C        
          W1 = (UN+SQRT(G*HMAX))/DIJ             
C   
          IF (W1.GT.XNC) XNC = W1
C
100    CONTINUE
C       
       IF(DTVARI) THEN 
         DT = CFLWTD/XNC
         IF(LISTIN.AND.LNG.EQ.1) WRITE(LU,*) 'PAS DE TEMPS : ',DT
         IF(LISTIN.AND.LNG.EQ.2) WRITE(LU,*) 'TIME-STEP: ',DT
       ENDIF         
C
       RETURN                                                            
       END        
