C                       *******************                       
                        SUBROUTINE MAXSLOPE
C                       *******************                       
C                                                                       
     *(SLOPE,ZF,ZR,XEL,YEL,NELEM,NELMAX,NPOIN,IKLE,EVOL,UNSV2D,MESH)                                                   
C                                                                       
C***********************************************************************
C SISYPHE VERSION 5.8      16/11/07    J-M HERVOUET (LNH) 01 30 87 80 18
C
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT 
C***********************************************************************
C                                                                       
C   FONCTION : COLLAPSE OF SAND WITH A SLOPE GREATER THAN A STABILITY
C              CRITERION
C                                                                       
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   SLOPE        | -->| MAXIMUM SLOPE IN DEGREES
C |   ZF           |<-->| BOTTOM THAT WILL BE MODIFIED
C |   ZR           | -->| NON ERODABLE BED
C |   XEL,YEL      | -->| MESH COORDINATES PER ELEMENT
C |   NELEM        | -->| NUMBER OF ELEMENTS
C |   NELMAX       | -->| MAXIMUM NUMBER OF ELEMENTS
C |   NPOIN        | -->| NUMBER OF POINTS IN THE MESH
C |   IKLE         | -->| CONNECTIVITY TABLE
C |   EVOL         |<-->| WORK ARRAY, THEN EVOLUTION DUE TO SLIDE
C |   UNSV2D       | -->| INVERSE OF INTEGRAL OF BASES
C |   MESH         | -->| MESH STRUCTURE
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE) 
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C APPELE PAR : SUBIEF                                                   
C                                                                       
C SOUS-PROGRAMME APPELE : OV                                            
C                                                                       
C***********************************************************************
C
      USE BIEF
C
      USE INTERFACE_SISYPHE, EX_MAXSLOPE => MAXSLOPE
C                                                                       
      IMPLICIT NONE                                                     
      INTEGER LNG,LU                                                    
      COMMON/INFO/LNG,LU                                                
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX,NPOIN
      INTEGER, INTENT(IN) :: IKLE(NELMAX,3)
C                                                                       
      DOUBLE PRECISION, INTENT(IN   ) :: SLOPE
      DOUBLE PRECISION, INTENT(INOUT) :: ZF(NPOIN)   
      DOUBLE PRECISION, INTENT(IN)    :: ZR(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: XEL(NELMAX,3),YEL(NELMAX,3)
C
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: EVOL
      TYPE(BIEF_OBJ), INTENT(IN)      :: UNSV2D
      TYPE(BIEF_MESH) :: MESH 
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C                                                                      
      INTEGER IELEM,I1,I2,I3,I
      DOUBLE PRECISION X2,X3,Y2,Y3,Z2,Z3,A,B,L,ZC,DEUXSURF,TANSL                                         
C  
      INTRINSIC SQRT,MIN,MAX,TAN
C                                                                     
C-----------------------------------------------------------------------
C
      TANSL=TAN(3.141592653589D0*SLOPE/180.D0)
C
C     INITIALIZING THE RIGHT-HAND SIDE EVOL TO ZERO
C
      CALL CPSTVC(UNSV2D,EVOL)
      CALL OS('X=0     ',X=EVOL)
C
C     LOOP ON ELEMENTS
C
      DO IELEM=1,NELEM
C
        I1=IKLE(IELEM,1)
        I2=IKLE(IELEM,2)
        I3=IKLE(IELEM,3)
C
        X2=XEL(IELEM,2)
        X3=XEL(IELEM,3)
        Y2=YEL(IELEM,2)
        Y3=YEL(IELEM,3)
        Z2=ZF(I2)-ZF(I1)
        Z3=ZF(I3)-ZF(I1)
C
C       TWICE THE TRIANGLE AREA
C
        DEUXSURF=X2*Y3-X3*Y2
C
C       AVERAGE BOTTOM IN THE ELEMENT
C
        ZC=(ZF(I1)+ZF(I2)+ZF(I3))/3.D0
C
C       COMPONENTS OF BOTTOM GRADIENT
C
        A=(Z2*Y3-Z3*Y2)/DEUXSURF
        B=(Z3*X2-Z2*X3)/DEUXSURF
C
C       CORRECTING FACTOR ON SLOPE 
C  
        L=MIN(1.D0,TANSL/MAX(SQRT(A**2+B**2),1.D-8))
C
C       L LIMITED DUE TO NON ERODABLE BEDS : ZF MUST NOT GO BELOW ZR
C
        IF(ZF(I1).GT.ZC) L=MAX(L,(ZR(I1)-ZC)/MAX(ZF(I1)-ZC,1.D-8))
        IF(ZF(I2).GT.ZC) L=MAX(L,(ZR(I2)-ZC)/MAX(ZF(I2)-ZC,1.D-8))
        IF(ZF(I3).GT.ZC) L=MAX(L,(ZR(I3)-ZC)/MAX(ZF(I3)-ZC,1.D-8))
C
C       BUILDING THE RIGHT-HAND SIDE       
C
        EVOL%R(I1)=EVOL%R(I1)+(1.D0-L)*(ZC-ZF(I1))*DEUXSURF/6.D0
        EVOL%R(I2)=EVOL%R(I2)+(1.D0-L)*(ZC-ZF(I2))*DEUXSURF/6.D0
        EVOL%R(I3)=EVOL%R(I3)+(1.D0-L)*(ZC-ZF(I3))*DEUXSURF/6.D0
C
      ENDDO
C             
C-----------------------------------------------------------------------
C
C     FINAL RESOLUTION
C
      IF(NCSIZE.GT.1) THEN
        CALL PARCOM(EVOL,2,MESH)
      ENDIF
C
      DO I=1,NPOIN
        EVOL%R(I)=EVOL%R(I)*UNSV2D%R(I)
      ENDDO           
C             
C-----------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END
