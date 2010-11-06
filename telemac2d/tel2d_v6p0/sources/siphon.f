C                       *****************                               
                        SUBROUTINE SIPHON                               
C                       *****************                               
C                                                                       
     *(RELAXS,NSIPH,ENTSIP,SORSIP,GRAV,                                 
     * H,ZF,ISCE,DSCE,SECSCE,ALTSCE,CSSCE,CESCE,DELSCE,ANGSCE,LSCE,     
     * NTRAC,T,TSCE,USCE,VSCE,U,V,ENTET,MAXSCE)                                       
C                                                                       
C***********************************************************************
C  TELEMAC 2D VERSION 5.9     19/04/96     V. GUINOT   (LHF)            
C              MODIFIE LE     03/10/96  J.-M. HERVOUET (LNH) 30 87 80 18
C                             03/2000    MODIF E. DAVID (SOGREAH)
C
C 16/02/2009 JMH : CORRECTION DE LA CORRECTION DU 03/2000 (EN //)
C
C***********************************************************************
C                                                                       
C      FONCTIONS: TRAITEMENT DES SIPHONS.                               
C      ==========                                                       
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   RELAXS       | -->| COEFFICIENT DE RELAXATION.                    
C |   NPSIPH       | -->| NOMBRE DE SIPHONS.                            
C |   ENTSIP       | -->| NUMERO DE L'ENTREE D'UNE BUSE DANS LA         
C |                |    | NUMEROTATION DES SOURCES.                     
C |   SORSIP       | -->| NUMERO DE LA SORTIE D'UNE BUSE DANS LA        
C |                |    | NUMEROTATION DES SOURCES.                     
C |   GRAV         | -->| PESANTEUR.                                    
C |   H            | -- | HAUTEUR D'EAU.                                
C |   ZF           | -- | COTES DU FOND.                                
C |   ISCE         | -->| NUMERO GLOBAL DES POINTS SOURCES.             
C |   DSCE         |<-- | DEBIT DES SOURCES.                            
C |   SECSCE       | -->| SECTION DES BUSES.                            
C |   ALTSCE       | -->| COTE DES BUSES.                               
C |   CSSCE        |<-- | COEFFICIENTS DE PERTE DE CHARGE               
C |                |    | LORS D'UN FONCTIONNEMENT EN SORTIE.           
C |   CESCE        |<-- | COEFFICIENTS DE PERTE DE CHARGE               
C |                |    | LORS D'UN FONCTIONNEMENT EN ENTREE.           
C |   DELSCE       | -->| ANGLE DES BUSES AVEC LA VERTICALE             
C |   ANGSCE       | -->| ANGLE DES BUSES AVEC L'AXE OX.                
C |   LSCE         |<-- | PERTE DE CHARGE LINEAIRE DE LA CONDUITE.      
C |   TRAC         | -->| LOGIQUE INDIQUANT S'IL Y A UN TRACEUR.        
C |   T            | -->| TRACEUR.                                      
C |   TSCE         | -->| VALEUR DU TRACEUR AUX SOURCES.                
C |   USCE         | -->| VITESSE U DU COURANT AUX SOURCES.             
C |   VSCE         | -->| VITESSE V DU COURANT AUX SOURCES.             
C |   U            | -->| VITESSE U.                                    
C |   V            | -->| VITESSE V.                                    
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE) 
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C APPELE PAR : TELMAC                                                   
C                                                                       
C SOUS-PROGRAMMES APPELES : NEANT                                       
C                                                                       
C********************************************************************** 
C 
      USE BIEF
C                                                                      
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C                                                                                                            
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NSIPH,NTRAC,MAXSCE
      INTEGER, INTENT(IN)             :: ENTSIP(*),SORSIP(*),ISCE(*)
      LOGICAL, INTENT(IN)             :: ENTET 
      DOUBLE PRECISION, INTENT(IN)    :: RELAXS,GRAV
      DOUBLE PRECISION, INTENT(INOUT) :: USCE(*),VSCE(*),DSCE(*)
      DOUBLE PRECISION, INTENT(INOUT) :: TSCE(MAXSCE,NTRAC)
      DOUBLE PRECISION, INTENT(IN)    :: ANGSCE(*),LSCE(*),CESCE(*)                        
      DOUBLE PRECISION, INTENT(IN)    :: CSSCE(*),DELSCE(*)
      DOUBLE PRECISION, INTENT(IN)    :: SECSCE(*),ALTSCE(*)
      DOUBLE PRECISION, INTENT(IN)    :: H(*),ZF(*),U(*),V(*)
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: T
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                       
      INTEGER N,I1,I2,IR1,IR2,ITRAC                                                                                                                             
C                                                                                        
      DOUBLE PRECISION SEC,L                                   
      DOUBLE PRECISION D1,D2,S1,S2,CE1,CE2,CS1,CS2,Q,HIR1,HIR2        
C                                                                       
      INTRINSIC SQRT,COS,SIN                                        
C                                                                       
      DOUBLE PRECISION P_DMAX                                          
      EXTERNAL         P_DMAX                                                   
C                                                                       
C-----------------------------------------------------------------------
C                                                                     
C BOUCLE SUR LES SIPHONS                                                
C                                                                       
      DO 10 N=1,NSIPH                                                   
C                                                                       
C     ON REPERE LES NOEUDS D'ENTREE, DE SORTIE                          
C 
C     NUMEROS DES SOURCES CORRESPONDANTES                                                                      
      I1=ENTSIP(N)                                                      
      I2=SORSIP(N) 
C     NUMEROS DES POINTS DE CES SOURCES                                                     
      IR1=ISCE(I1)                                                      
      IR2=ISCE(I2)                                                      
C                                                                       
C     CHARGES, ASSIMILEES A LA COTE DE SURFACE LIBRE                    
C                                                                       
      IF(IR1.GT.0) THEN                                                 
        S1=H(IR1)+ZF(IR1)
        HIR1=H(IR1)                                               
      ELSE                                                              
        S1=-1.D10
        HIR1=-1.D10                                                       
      ENDIF                                                             
      IF(IR2.GT.0) THEN                                                 
        S2=H(IR2)+ZF(IR2)
        HIR2=H(IR2)                                               
      ELSE                                                              
        S2=-1.D10 
        HIR2=-1.D10                                                      
      ENDIF                                                             
C     CAS OU UNE DES EXTREMITES N'EST PAS DANS LE SOUS-DOMAINE          
      IF(NCSIZE.GT.1) THEN                                              
        S1=P_DMAX(S1)                                                   
        S2=P_DMAX(S2)
        HIR1=P_DMAX(HIR1)
        HIR2=P_DMAX(HIR2)                                                   
      ENDIF                                                             
C                                                                       
C     COEFFICIENTS POUR LES CALCULS DE PERTE DE CHARGE                  
C                                                                       
      D1=DELSCE(I1)                                                     
      D2=DELSCE(I2)                                                     
      CE1=CESCE(I1)                                                     
      CE2=CESCE(I2)                                                     
      CS1=CSSCE(I1)                                                     
      CS2=CSSCE(I2)                                                     
      SEC=SECSCE(I1)                                                    
      L  =LSCE(I1)                                                      
C                                                                       
C     CALCUL DU DEBIT EN FONCTION DU DELTAH                             
C     SI LA PERTE DE CHARGE LINEAIRE EST NEGLIGEABLE, ON PEUT SE        
C     PERMETTRE D'AVOIR DES SECTIONS DIFFERENTES EN ENTREE ET EN SORTIE 
C                                                                       
      IF(S1.GE.S2) THEN                                                 
C EDd + CCt 03/2000 (CORRIGE PAR JMH 16/02/2009 H(IR1) ET H(IR2)
C                    NON GARANTIS EN PARALLELISME, D'OU HIR1 ET HIR2) 
C        IF(S1.GT.ALTSCE(I1).AND.S1.GT.ALTSCE(I2)) THEN    
        IF(S1.GT.ALTSCE(I1).AND.S1.GT.ALTSCE(I2).AND.
     *     HIR1.GT.0.02D0) THEN                                
C
          Q = SEC * SQRT( 2.D0*GRAV*(S1-S2)/(CE1+L+CS2) )               
        ELSE                                                            
          Q=0.D0                                                        
        ENDIF                                                           
      ELSE                                                              
C EDd + CCt 03/2000 
C        IF(S2.GT.ALTSCE(I1).AND.S2.GT.ALTSCE(I2)) THEN                  
        IF(S2.GT.ALTSCE(I1).AND.S2.GT.ALTSCE(I2).AND.
     *     HIR2.GT.0.02D0) THEN                  
C
          Q = - SEC * SQRT( 2.D0*GRAV*(S2-S1)/(CS1+L+CE2) )             
        ELSE                                                            
          Q=0.D0                                                        
        ENDIF                                                           
      ENDIF                                                             
C                                                                       
C     SI LES CHARGES DES DEUX COTES SONT INFERIEURES AUX                
C     COTES D'ORIFICES DU SIPHON, IL NE PASSE RIEN                      
C                                                                       
      IF(S1.LT.ALTSCE(I1).AND.S2.LT.ALTSCE(I2)) Q=0.D0                  
C                                                                       
C  ON REMPLIT LE TABLEAU DSCE EN EFFECTUANT UNE RELAXATION.             
C                                                                       
      DSCE(I2)= RELAXS * Q + (1.D0-RELAXS) * DSCE(I2)                   
      DSCE(I1)=-DSCE(I2)                                                
C
      IF(ENTET) THEN
        WRITE(LU,*) ' '
        IF(LNG.EQ.1) THEN
          WRITE(LU,*) 'SIPHON ',N,' DEBIT DE ',DSCE(I2),' M3/S'
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*) 'CULVERT ',N,' DISCHARGE OF ',DSCE(I2),' M3/S'
        ENDIF
        WRITE(LU,*) ' '
      ENDIF 
C                                                                      
C  TRAITEMENT DES VITESSES DES SOURCES.                                 
C  ON FAIT POUR LES VITESSES COMME POUR LE TRACEUR                      
C                                                                       
      IF(DSCE(I1).GT.0.D0) THEN                                         
        USCE(I1) = ( COS(D1)*DSCE(I1)/SECSCE(I1) ) * COS(ANGSCE(I1))    
        VSCE(I1) = ( COS(D1)*DSCE(I1)/SECSCE(I1) ) * SIN(ANGSCE(I1))    
      ELSE                                                              
        IF(IR1.GT.0) THEN                                               
          USCE(I1) = U(IR1)                                             
          VSCE(I1) = V(IR1)                                             
        ENDIF                                                           
      ENDIF                                                             
      IF(DSCE(I2).GT.0.D0) THEN                                         
        USCE(I2) = ( COS(D2)*DSCE(I2)/SECSCE(I2) ) * COS(ANGSCE(I2))    
        VSCE(I2) = ( COS(D2)*DSCE(I2)/SECSCE(I2) ) * SIN(ANGSCE(I2))    
      ELSE                                                              
        IF(IR2.GT.0) THEN                                               
          USCE(I2) = U(IR2)                                             
          VSCE(I2) = V(IR2)                                             
        ENDIF                                                           
      ENDIF                                                             
C                                                                       
C  TRAITEMENT DU TRACEUR :                                              
C                                                                       
      IF(NTRAC.GT.0) THEN
        DO ITRAC=1,NTRAC                                                
          IF(DSCE(I1).GT.0.D0) THEN                                       
            IF(IR2.GT.0) THEN                                             
              TSCE(I1,ITRAC)=T%ADR(ITRAC)%P%R(IR2)                                             
            ELSE                                                          
              TSCE(I1,ITRAC)=-1.D10                                             
            ENDIF                                                         
            IF(NCSIZE.GT.1) TSCE(I1,ITRAC)=P_DMAX(TSCE(I1,ITRAC))                   
          ENDIF                                                           
          IF(DSCE(I2).GT.0.D0) THEN                                       
            IF(IR1.GT.0) THEN                                             
              TSCE(I2,ITRAC)=T%ADR(ITRAC)%P%R(IR1)                                             
            ELSE                                                          
              TSCE(I2,ITRAC)=-1.D10                                        
            ENDIF                                                         
            IF(NCSIZE.GT.1) TSCE(I2,ITRAC)=P_DMAX(TSCE(I2,ITRAC))                   
          ENDIF
        ENDDO                                                           
      ENDIF                                                             
C                                                                       
C  FIN DE LA BOUCLE SUR LES SIPHONS.                                    
C                                                                       
10    CONTINUE                                                          
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END
