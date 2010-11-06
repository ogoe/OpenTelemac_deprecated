C                       **************
                        SUBROUTINE MAJ
C                       **************
C
     *(NS,NSEG,NPTFR,G,DT,AIRS,
     * H,QU,QV,UA,CE,NBOR,LIMPRO,XNEBOR,YNEBOR,KNEU,SMH,KFROT,CF)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.4                                           INRIA
C
C***********************************************************************
C
C     FONCTION  : CALCUL DES SOLUTIONS HYDRO A L'INSTANT N+1 
C                      PAR SCHEMA EXPLICITE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |  NS            | -->|  NOMBRE DE POINTS DU MAILLAGE                |
C |  NSEG          | -->|  NOMBRE D'ARETES DU MAILLAGE                 |
C |  NPTFR         | -->|  NOMBRE DE POINTS FRONTIERE                  |
C |  G             | -->|  CONSTANTE DE GRAVITE                        |
C |  DT            | -->|  PAS DE TEMPS                                |
C |  AIRS          | -->|  AIRES DES CELLULES                          |
C |  H             | -->|  HAUTEURS D'EAU AU TEMPS N                   |
C |  QU,QV         | -->|  COMPOSANTES DU DEBIT AU TEMPS N             |
C |  UA            |<-- |  H, QU, QV  AU TEMPS N+1                     |
C |  CE            | -->|  FLUX                                        |
C |  NBOR          | -->|  NUMEROS GLOBAUX DES POINTS DE BORD          |
C |  LIMPRO        | -->|  TYPES DE CONDITIONS AUX LIMITES             |
C |  XNEBOR,YNEBOR | -->|  NORMALE AUX POINTS FRONTIERE                |
C |  KNEU          | -->|  CONVENTION POUR LES POINTS NEUMANN          |
C |  SMH           | -->|  TERMES SOURCES DE L'EQUATION DE CONTINUITE  |
C |  KFROT         | -->|  LOI DE FROTTEMENT SUR LE FOND               |
C |  CF            | -->|  COEFFICIENT DE FROTTEMENT                   !
C !________________!____!______________________________________________!
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)   
C        -- (TABLEAU DE TRAVAIL)                                        
C-----------------------------------------------------------------------
C     - SOUS PROGRAMME(S) APPELANT : RESOLU                             
C     - SOUS PROGRAMME(S) APPELES  : CDLPROJ, FRICTION
C 
C***********************************************************************
C
      USE INTERFACE_TELEMAC2D, EX_MAJ => MAJ
C
      IMPLICIT NONE
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NS,NSEG,NPTFR,KNEU,KFROT
      INTEGER, INTENT(IN) :: NBOR(NPTFR),LIMPRO(NPTFR,6)
      DOUBLE PRECISION, INTENT(IN) :: XNEBOR(2*NPTFR),YNEBOR(2*NPTFR)
      DOUBLE PRECISION, INTENT(IN) :: H(NS),QU(NS),QV(NS),SMH(NS)
      DOUBLE PRECISION, INTENT(IN) :: CE(3,NS),CF(NS),AIRS(NS)
      DOUBLE PRECISION, INTENT(IN) :: G,DT  
      DOUBLE PRECISION, INTENT(INOUT) :: UA(3,NS) 
C     
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IS
C                  
      DOUBLE PRECISION USAIS
C
      INTRINSIC ABS
C
C-----------------------------------------------------------------------
C      
C EXPLICIT RESOLUTION
C	
C     UPDATING THE PHYSICAL SOLUTION 
C
      DO IS =1,NS
C
        USAIS = DT/AIRS(IS)
        UA(1,IS)  = H(IS) + USAIS*(CE(1,IS)+SMH(IS))
        UA(2,IS)  = QU(IS) + USAIS*CE(2,IS)
        UA(3,IS)  = QV(IS) + USAIS*CE(3,IS)
C
      ENDDO
C
C   PROJECTION ON THE SLIPPING BOUNDARY CONDITIONS
C
      CALL CDLPROJ(NS,NPTFR,NBOR,LIMPRO,XNEBOR,YNEBOR,KNEU,UA)
C
      DO IS =1,NS
        IF(UA(1,IS).LE.1.D-12) UA(1,IS)=0.D0
        IF(ABS(UA(2,IS)).LE.1.D-12) UA(2,IS)=0.D0
        IF(ABS(UA(3,IS)).LE.1.D-12) UA(3,IS)=0.D0
      ENDDO
C
      IF(KFROT.NE.0) CALL FRICTION(NS,G,DT,UA,H,QU,QV,CF)
C
C-----------------------------------------------------------------------
C 
      RETURN
      END
