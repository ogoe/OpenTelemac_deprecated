      ! *************************** !
        SUBROUTINE BEDLOAD_SOLVS_FE !
      ! *************************** !

     &(MESH,S,EBOR,MASKEL,MASK,  
     & QSX,QSY,IELMT,NPOIN,NPTFR,KENT,KDIR,LIMTEC,DT,
     & MSK, ENTET, T1, T4, T8, 
     & ZFCL,HZ,HZN,GLOSEG,DIMGLO,FLODEL,FLULIM,NSEG,UNSV2D)
C
C***********************************************************************
C SISYPHE VERSION 6.0  05/09/2009  J.-M. HERVOUET (NEW METHOD) 
C SISYPHE VERSION 5.8  29/10/2007  J.-M. HERVOUET 
C SISYPHE VERSION 5.5  14/09/2004  F.    HUVELIN                    
C SISYPHE VERSION 5.3  --/--/2002  B.    MINH DUC                     
C SISYPHE VERSION 5.1  11/09/1995  E.    PELTIER                      
C SISYPHE VERSION 5.1  11/09/1995  C.    LENORMANT                    
C SISYPHE VERSION 5.1  11/09/1995  J.-M. HERVOUET                     
C***********************************************************************
C
C  FONCTION: 
C
C                  ! =============================== !
C                  !               D(HZ)             !
C                  ! RESOLUTION DE ---- + DIV(T) = 0 !
C                  !                DT               !
C                  ! =============================== !
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   ZFCL         |<-- | ZFCL=HZ-HZN
C |   HZ           |<-- | NEW AVAILABLE LAYER OF SEDIMENT
C |   HZN          | -->| OLD AVAILABLE LAYER OF SEDIMENT
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : BEDLOAD_EVOL
C
C SOUS-PROGRAMMES APPELES : Cpstvc
C                           positive_depths
C
C***********************************************************************
C
      USE BIEF
      USE INTERFACE_SISYPHE, EX_BEDLOAD_SOLVS_FE => BEDLOAD_SOLVS_FE
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
      TYPE(BIEF_OBJ),   INTENT(IN)    :: S,LIMTEC,MASKEL,MASK,QSX,QSY
      INTEGER,          INTENT(IN)    :: IELMT,NPOIN,NPTFR,KENT,KDIR
      INTEGER,          INTENT(IN)    :: DIMGLO,NSEG
      INTEGER,          INTENT(IN)    :: GLOSEG(DIMGLO,2)
      DOUBLE PRECISION, INTENT(IN)    :: DT
      DOUBLE PRECISION, INTENT(INOUT) :: FLULIM(NSEG)
      LOGICAL,          INTENT(IN)    :: MSK,ENTET
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: FLODEL,T1,T4,T8
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: HZ,EBOR
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: ZFCL
      TYPE(BIEF_OBJ),   INTENT(IN)    :: HZN,UNSV2D
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER K
C
C-----------------------------------------------------------------------
C  
C     BOUNDARY FLUXES
C    
      CALL CPSTVC(QSX,T4)
      CALL OS('X=C     ',X=T4,C=1.D0)
      CALL VECTOR(T8, '=', 'FLUBDF          ',IELBOR(IELMT,1),1.D0,
     &            T4,S,S,QSX,QSY,S,MESH,.TRUE.,MASK)
C
C     HERE THE VARIABLE WILL BE THE LAYER DEPTH OF THE SEDIMENT CLASS
C     NOT THE EVOLUTION
C    
      DO K=1,NPTFR
        IF(LIMTEC%I(K).EQ.KDIR) THEN
          EBOR%R(K)=EBOR%R(K)+HZN%R(MESH%NBOR%I(K))
        ENDIF
      ENDDO      
C 
!     CALL VECTOR(T1,'=','VGRADP          ',QSX%ELM,-1.D0,
!    *            S,S,S,QSX,QSY,S,MESH,MSK,MASKEL)
C                 T1 AS HUGRADP IS NOT USED AS AN ASSEMBLED VECTOR
C                 BUT TO GET THE NON ASSEMBLED FORM MESH%W
C     JUST LIKE CALL VECTOR BUT WITHOUT ASSEMBLING T1 BECAUSE LEGO IS SET
C     TO FALSE (ONLY NON ASSEMBLED MESH%W%R IS USED AFTER)
      CALL VECTOS(T1%R,'=','VGRADP          ',-1.D0,
     *            S%R,S%R,S%R,QSX%R,QSY%R,S%R,
     *            S,S,S,QSX,QSY,S,
C                           LEGO
     *            MESH%W%R,.FALSE.,
     *            MESH%XEL%R  , MESH%YEL%R  , MESH%ZEL%R  ,
     *            MESH%SURFAC%R,MESH%IKLE%I,MESH%NBOR%I,
     *            MESH%XSGBOR%R, MESH%YSGBOR%R, MESH%ZSGBOR%R,
     *            NBPTS(QSX%ELM),MESH%NELEM,MESH%NELMAX,
     *            QSX%ELM,MESH%LV,MSK,MASKEL%R,MESH)
C
      CALL POSITIVE_DEPTHS(T1,T4,HZ,HZN,MESH,
     *                     FLODEL,.TRUE.,T8,DT,UNSV2D,NPOIN,
     *                     GLOSEG(1:DIMGLO,1),GLOSEG(1:DIMGLO,2),
     *                     MESH%NBOR%I,NPTFR,.FALSE.,T8,.FALSE.,
     *                     1,FLULIM,
     *                     LIMTEC%I,EBOR%R,KDIR,ENTET,MESH%W%R)
C 
      CALL OS('X=Y-Z   ' ,X=ZFCL,Y=HZ,Z=HZN) 
C
C     DIRICHLET CONDITIONS       
C
      DO K=1,NPTFR
        IF(LIMTEC%I(K).EQ.KDIR) THEN
          EBOR%R(K)=EBOR%R(K)-HZN%R(MESH%NBOR%I(K))
          ZFCL%R(MESH%NBOR%I(K)) = EBOR%R(K)
        ENDIF
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
