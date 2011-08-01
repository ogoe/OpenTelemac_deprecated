C                       **************************
                        SUBROUTINE SHARE_3D_FLUXES
C                       **************************
C
     *(FLUX,XMUL,NPLAN,MESH2,MESH3,OPT)
C
C***********************************************************************
C BIEF VERSION 6.0        14/04/2010  J-M HERVOUET (LNHE) 01 30 87 80 18
C COPYRIGHT EDF 2009              
C***********************************************************************
C
C   FONCTION : SHARING ASSEMBLED FLUXES BETWEEN SUB-DOMAINS
C              THIS IS IN SOME SENSE THE CONTRARY OF PARCOM2_SEG, BUT
C              THE FLUXES ON THE HORIZONTAL SEGMENTS WILL BE DIVIDED
C              BY XMUL*2 AND THE FLUXES OF VERTICAL INTERFACE SEGMENTS
C              WILL BE DIVIDED BY XMUL*MESH%FAC%R
C 
C   BEWARE: MUST NOT BE CALLED WHEN NCSIZE = 0            
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C |      FLUX      |<-->| FLUXES TO BE SHARED
C |      XMUL      | -->| MULTIPLICATING FACTOR
C |      NPLAN     | -->| NUMBER OF PLANES
C |      MESH2     | -->| 2D MESH
C |      MESH3     | -->| 3D MESH
C |      OPT       | -->| 1 : HORIZONTAL AND VERTICAL SEGMENTS ONLY
C |                |    | 2 : ALL SEGMENTS
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   :  PLANTE
C
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NPLAN,OPT
C
      TYPE(BIEF_MESH) , INTENT(INOUT) :: MESH2,MESH3
C
      DOUBLE PRECISION, INTENT(INOUT) :: FLUX(*)
      DOUBLE PRECISION, INTENT(INOUT) :: XMUL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IPLAN,NSEG,NSEGH,NSEGV,NPOIN2,I,I3D,I2D,IAD
C
C-----------------------------------------------------------------------
C
      NSEG=MESH2%NSEG
      NPOIN2=MESH2%NPOIN
      NSEGH=NSEG*NPLAN
      NSEGV=NPOIN2*(NPLAN-1)
C
C     HORIZONTAL FLUXES
C
      DO IPLAN=1,NPLAN
        CALL MULT_INTERFACE_SEG(FLUX(1+(IPLAN-1)*NSEG:IPLAN*NSEG),
     *                          MESH2%NH_COM_SEG%I,
     *                          MESH2%NH_COM_SEG%DIM1,
     *                          MESH2%NB_NEIGHB_SEG,
     *                          MESH2%NB_NEIGHB_PT_SEG%I,
     *                          0.5D0*XMUL,NSEG)
      ENDDO
C
C     VERTICAL FLUXES (SAME NUMBERING AS POINTS, SO FAC%R(I))
C
      IAD=1
      DO I=1,NPTIR
C       I2D=NACHB(1,I) WITH NACHB OF SIZE NACHB(NBMAXNSHARE,NPTIR)
C       IAD IS (I-1)*NBMAXNSHARE+1
        I2D=MESH2%NACHB%I(IAD)
        DO IPLAN=1,NPLAN-1
          I3D=(IPLAN-1)*NPOIN2+I2D
          FLUX(NSEGH+I3D)=FLUX(NSEGH+I3D)*MESH3%FAC%R(I3D)*XMUL
        ENDDO
        IAD=IAD+NBMAXNSHARE
      ENDDO
C
C     ALTERNATE WAY, SIMPLER BUT LESS EFFICIENT FOR VERTICAL FLUXES:
C
!     DO I3D=1,NSEGV
!       FLUX(NSEGH+I3D)=FLUX(NSEGH+I3D)*MESH3%FAC%R(I3D)*XMUL
!     ENDDO
C
C     CROSSED FLUXES (SEE STOSEG41 FOR STORAGE). THERE ARE 2*NESG
C     PER LAYER AND NPLAN-1 LAYER. HERE ORISEG=1 AND ORISEG=2 SEGMENTS
C     ARE MULTIPLIED BY THE SAME NUMBER, SO GIVEN HOW THE NUMBERING
C     IS BUILT IT IS LIKE IF WE HAVE 2*NPLAN-2 LAYERS OF HORIZONTAL SEGMENTS
C      
      IF(OPT.EQ.2) THEN
        DO IPLAN=1,2*(NPLAN-1)
          CALL MULT_INTERFACE_SEG(FLUX(1+NSEGH+NSEGV+(IPLAN-1)*NSEG:
     *                                   NSEGH+NSEGV+ IPLAN   *NSEG),
     *                            MESH2%NH_COM_SEG%I,
     *                            MESH2%NH_COM_SEG%DIM1,
     *                            MESH2%NB_NEIGHB_SEG,
     *                            MESH2%NB_NEIGHB_PT_SEG%I,
     *                            0.5D0*XMUL,NSEG)
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
