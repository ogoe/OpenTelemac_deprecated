C                       *****************                               
                        SUBROUTINE PORO11                               
C                       *****************                               
C                                                                       
     *(TETA,ZF,HN,IKLE,NELEM,NELMAX)                                    
C                                                                       
C***********************************************************************
C TELEMAC 2D VERSION 5.2     01/08/97    J-M HERVOUET (LNH) 30 87 80 18 
C                                  PAUL BATES (BRISTOL) 44 117 928 9108
C***********************************************************************
C                                                                       
C FONCTION : MARQUAGE DES BANCS DECOUVRANTS   
C MODIFIED TO IMPLEMENT WETTING/DRYING ALGORITHM
C OF DELFINA ET AL                          
C                                                                       
C            PARTIALLY WET ELEMENT : TETA = 0 < NU < 1                           
C            WET ELEMENT           : TETA = NU = 1.0
C            DRY ELEMENT           : TETA = NU = 1.0                             
C                                                                       
C            LE CRITERE DE DECOUVREMENT EST CELUI DE J.-M. JANIN :      
C            LE FOND D'UN POINT DE L'ELEMENT EST PLUS HAUT QUE LA       
C            SURFACE LIBRE D'UN AUTRE.                                  
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      TETA      |<-- |  NU (PAR ELEMENT)                            |
C |      SL,ZF     | -->|  SURFACE LIBRE ET FOND                       |
C |      MESH      | -->|  STRUCTURE DE MAILLAGE                       |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE) 
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C  APPELE PAR : SETNU                                                  
C                                                                       
C  SOUS-PROGRAMME APPELE : OS                                           
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
      INTEGER, INTENT(IN)             :: NELEM,NELMAX
      INTEGER, INTENT(IN)             :: IKLE(NELMAX,*)
      DOUBLE PRECISION, INTENT(IN)    :: ZF(*),HN(*)
      DOUBLE PRECISION, INTENT(INOUT) :: TETA(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                       
      INTEGER IELEM,OPTNU                                                                                          
C                                            
      DOUBLE PRECISION Y(1),Z(1)                                
      DOUBLE PRECISION SL1,SL2,SL3,ZF1,ZF2,ZF3,H1,H2,H3,EPS1,EPS2                      
C
      LOGICAL TEST1, TEST2
C                                                                       
      INTRINSIC MAX,MIN                                                 
C                                                                       
C-----------------------------------------------------------------------
C                                                                 
      CALL OV( 'X=C     ' , TETA , Y , Z , 1.D0 , NELEM )                   
C                                                                       
C--------------------------------------------------------------
C
C  OPTION FOR CALCULATION OF NU
C
C  1 = EXACT CALCULATION BASED ON PROJECTED INUNDATED AREA
C  2 = ROUGH APPROXIMATION TO INUNDATED AREA
C  3 = CONSTANT VALUE
C  4 = USER SPECIFIED
C   
      EPS1=1.D-6
      EPS2=0.3D0
C
C-----------------------------------------------------------------------
C
      DO 4 IELEM = 1 , NELEM                                         
C                                                                            
        ZF1 = ZF(IKLE(IELEM,1))                                      
        ZF2 = ZF(IKLE(IELEM,2))                                      
        ZF3 = ZF(IKLE(IELEM,3))                                      
C 
        H1 = HN(IKLE(IELEM,1))
        H2 = HN(IKLE(IELEM,2))
        H3 = HN(IKLE(IELEM,3))
C
        SL1 = H1 + ZF1                                      
        SL2 = H2 + ZF2                                      
        SL3 = H3 + ZF3                                      
C     
        TEST1 = .FALSE.
        TEST2 = .FALSE.
C
C       IF TIDAL FLAT SUSPECTED AND SOME DEPTH
C       (POROSITY LEFT TO 1 IN FULLY DRY ELEMENTS)
C             
        IF(     MAX(ZF1,ZF2,ZF3).GT.MIN(SL1,SL2,SL3)
     *     .AND.MAX(H1,H2,H3).GT.EPS1                ) TEST1 = .TRUE.
C
C       
C
        IF(MAX(H1,H2,H3).EQ.H1.AND.SL1.GT.MIN(ZF1,ZF2,ZF3)
     *                        .AND.SL1.LT.MAX(ZF1,ZF2,ZF3).OR.
     *     MAX(H1,H2,H3).EQ.H2.AND.SL2.GT.MIN(ZF1,ZF2,ZF3)
     *                        .AND.SL2.LT.MAX(ZF1,ZF2,ZF3).OR.
     *     MAX(H1,H2,H3).EQ.H3.AND.SL3.GT.MIN(ZF1,ZF2,ZF3)
     *                    .AND.SL3.LT.MAX(ZF1,ZF2,ZF3)) TEST2=.TRUE.
C
C---------------------------------------------------------------------
C
        IF(TEST1.AND.TEST2) THEN
C         ROUGH CALCULATION OF POROSITY BY PERCENTAGE OF ELEMENT INUNDATED
          TETA(IELEM)=MAX(H1,H2,H3)/(MAX(ZF1,ZF2,ZF3)-MIN(ZF1,ZF2,ZF3))
          TETA(IELEM) = MAX(TETA(IELEM),EPS2) 
        ENDIF
C
C---------------------------------------------------------------------      
C
4     CONTINUE                                                       
C                                                                       
C---------------------------------------------------------------------
C                                                              
      RETURN                                                            
      END
