C                       ***********************                         
                        LOGICAL FUNCTION INPOLY                         
C                       ***********************                         
C                                                                       
     *( X , Y , XSOM , YSOM , NSOM )                                    
C                                                                       
C***********************************************************************
C BIEF VERSION 5.2     18/06/96         J.-M. HERVOUET (LNH) 30 87 83 81
C                        IDEE ET CODE INITIAL FOURNIS PAR E. DAVID (LHF)
C                                                           MERCI ERIC !
C                                        CORRECTION D'UN CAS PARTICULIER
C                     PAR JEAN-PHILIPPE RENAUD (CSN BRISTOL) LE 27/07/99 
C                                                  MERCI JEAN-PHILIPPE !
C***********************************************************************
C                                                                       
C    FONCTION  : INDIQUE SI UN POINT DE COORDONNEES X ET Y EST SITUE    
C                DANS UN POLYGONE DE SOMMETS DONNES.                    
C                                                                       
C                EXEMPLE : POUR UN TRIANGLE NSOM=3                      
C                                                                       
C    LES SOMMETS DU POLYGONE DOIVENT ETRE DISTINCTS.                    
C                                                                       
C    PRINCIPE  : ON PREND UNE DEMI-DROITE ISSUE DU POINT ET ON COMPTE LE
C                NOMBRE D'INTERSECTIONS AVEC LE POLYGONE.               
C                                                                       
C    CETTE FONCTION MARCHE AUSSI SI LE POLYGONE EST NON CONVEXE.        
C                                                                       
C    L'INTERSECTION EST CHERCHE A L'AIDE D'EQUATIONS PARAMETRIQUES      
C    DES DROITES :                                                      
C                                                                       
C    X + A * MU = XDEP + (XARR-XDEP) * LAMBDA                           
C    Y + B * MU = YDEP + (YARR-YDEP) * LAMBDA                           
C                                                                       
C    LA DEMI-DROITE EST CARACTERISEE PAR LE CHOIX DE A ET B, ET LE SIGNE
C    DE MU. IL Y A INTERSECTION SI MU >0 ET  0 < LAMBDA < 1             
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____._______________________________________________
C |      NOM       :MODE:                   ROLE                        
C |________________:____:_______________________________________________
C |   X,Y          :<-->: COORDONNEES DU POINT.                         
C |   XSOM         : -->: TABLEAU DES ABSCISSES DES SOMMETS DU POLYGONE 
C |   YSOM         : -->: TABLEAU DES ORDONNEES DES SOMMETS DU POLYGONE 
C |   NSOM         : -->: NOMBRE DE SOMMETS                             
C |________________:____:_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE) 
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C  APPELE PAR :                                                         
C                                                                       
C  SOUS-PROGRAMME APPELE :                                              
C                                                                       
C***********************************************************************
C                                                                       
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NSOM 
      DOUBLE PRECISION, INTENT(IN) :: X , Y                           
      DOUBLE PRECISION, INTENT(IN) :: XSOM(NSOM) , YSOM(NSOM)   
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                        
      INTEGER N,NSECT                                              
C                                                                       
      DOUBLE PRECISION A,B,ANGLE,XDEP,YDEP,XARR,YARR,DET,MU,LAMBDA,EPS            
C                                                                       
      INTRINSIC COS,SIN,ABS,MOD                                         
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      EPS = 1.D-9                                                       
      ANGLE = -1.D0                                                     
C                                                                       
C BOUCLE DE CHOIX DE A ET B POUR EVITER LES CAS PARTICULIERS            
C                                                                       
1000  CONTINUE                                                          
      ANGLE = ANGLE + 1.D0                                              
      IF(ANGLE.GT.360.D0) THEN
C       CAS PARTICULIER D'UN POINT SUR LE CONTOUR                                          
        INPOLY=.TRUE.
        RETURN                                               
      ENDIF                                                             
      A = COS(ANGLE*3.141592653D0/180.D0)                                 
      B = SIN(ANGLE*3.141592653D0/180.D0)                                 
      NSECT=0                                                           
C                                                                       
C BOUCLE SUR TOUS LES SEGMENTS DU POLYGONE                              
C                                                                       
      DO 10 N=1,NSOM                                                    
C                                                                       
C     DEP : PREMIER POINT SU SEGMENT    ARR : SECOND POINT              
C                                                                       
      XDEP=XSOM(N)                                                      
      YDEP=YSOM(N)                                                      
      IF(N.LT.NSOM) THEN                                                
        XARR=XSOM(N+1)                                                  
        YARR=YSOM(N+1)                                                  
      ELSE                                                              
        XARR=XSOM(1)                                                    
        YARR=YSOM(1)                                                    
      ENDIF                                                             
C                                                                       
C CAS DE DEUX POINTS SUCCESSIFS CONFONDUS                               
C                                                                       
      IF(ABS(XDEP-XARR)+ABS(YDEP-YARR).LT.EPS) THEN                     
        WRITE(LU,*) ' '                                                 
        WRITE(LU,*) ' '                                                 
        IF(LNG.EQ.1) THEN                                               
          WRITE(LU,*) 'INPOLY : POINTS CONFONDUS DANS LE POLYGONE'      
          WRITE(LU,*) 'AU POINT DE COORDONNEES : ',XDEP,'  ET  ',YDEP   
          WRITE(LU,*) 'DE NUMERO ',N                                    
          IF(N.EQ.NSOM) THEN                                            
            WRITE(LU,*) 'LE DERNIER POINT NE DOIT PAS ETRE EGAL AU'     
            WRITE(LU,*) 'PREMIER (POUR UN TRIANGLE, PAR EXEMPLE,'       
            WRITE(LU,*) 'IL FAUT DONNER TROIS POINTS ET NON QUATRE)'    
          ENDIF                                                         
          WRITE(LU,*) 'INPOLY EST PROBABLEMENT APPELE PAR FILPOL'       
        ENDIF                                                           
        IF(LNG.EQ.2) THEN                                               
          WRITE(LU,*) 'INPOLY: SUPERIMPOSED POINTS IN THE POLYGON'      
          WRITE(LU,*) 'AT POINT: ',XDEP,'  AND  ',YDEP,' WITH NUMBER ',N
          IF(N.EQ.NSOM) THEN                                            
            WRITE(LU,*) 'THE LAST POINT MUST NOT BE EQUAL TO THE FIRST' 
            WRITE(LU,*) 'FOR EXAMPLE, GIVE 3 POINTS FOR A TRIANGLE'     
          ENDIF                                                         
          WRITE(LU,*) 'INPOLY IS PROBABLY CALLED BY FILPOL'             
        ENDIF                                                           
        WRITE(LU,*) ' '                                                 
        WRITE(LU,*) ' '                                                 
        STOP                                                            
      ENDIF                                                             
C                                                                       
C     CAS OU LE POINT EST UN SOMMET (L'ALGORITHME GENERAL LE DONNERAIT  
C                                    EXTERIEUR AVEC 2 INTERSECTIONS)    
C                                                                       
      IF(ABS(X-XDEP).LE.EPS.AND.ABS(Y-YDEP).LE.EPS) THEN                
        NSECT=1                                                         
        GO TO 2000                                                      
      ENDIF                                                             
C                                                                       
C     DETERMINANT DU SYSTEME DE KRAMER                                  
C                                                                       
      DET = A*(YDEP-YARR)-B*(XDEP-XARR)                                 
      IF(ABS(DET).LT.EPS) GO TO 1000                                    
C                                                                       
      MU     = ( (XDEP-X)*(YDEP-YARR)-(YDEP-Y)*(XDEP-XARR) ) / DET      
      LAMBDA = (    A    *(YDEP-Y   )-    B   *(XDEP-X   ) ) / DET
C
C-------------------------------------------------------
C CORRECTION JP RENAUD (CSN BRISTOL) POUR EVITER QUE LE 
C POINT D'INTERSECTION NE SOIT UN DES SOMMETS
C
C SI LE POINT D'INTERSECTION EST UN SOMMET, ON INCREMENTE L'ANGLE
C CAR LE POINT SERAIT COMPTE DEUX FOIS AU LIEU D'UNE SEULE.
C
      IF ((ABS(X+A*MU-XDEP).LE.EPS.AND.ABS(Y+B*MU-YDEP).LE.EPS)
     *    .OR.
     *    (ABS(X+A*MU-XARR).LE.EPS.AND.ABS(Y+B*MU-YARR).LE.EPS))
     *    GOTO 1000
C
C FIN DE CORRECTION JP RENAUD
C-------------------------------------------------------
C      
      IF(MU.GE.-EPS.AND.LAMBDA.GE.-EPS.AND.LAMBDA.LE.1.D0+EPS) THEN     
        NSECT=NSECT+1                                                   
      ENDIF                                                             
C                                                                       
10    CONTINUE                                                          
C                                                                       
2000  CONTINUE                                                          
C                                                                       
      INPOLY=(MOD(NSECT,2).EQ.1)                                        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END                                                               
 
