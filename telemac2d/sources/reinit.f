C                       *****************                            
                        SUBROUTINE REINIT
C                       *****************                             
C                                                                       
     *(NS,NSEG,NPTFR,H,SMTR,HSTOK,HC,HCSTOK,FLUXT,FLUHBOR,DTT,NTRAC)
C                                                                       
C***********************************************************************
C  TELEMAC 2D VERSION 5.8                                         INRIA
C***********************************************************************
C
C   INITIALISATION D'UN PAS DE TEMPS TRACEUR
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C |  NS            | -->|  NOMBRE DE POINTS DU MAILLAGE                |
C |  NSEG          | -->|  NOMBRE D'ARETES DU MAILLAGE                 |
C |  NPTFR         | -->|  NOMBRE DE POINTS FRONTIERE                  |
C |  H             | -->|  HAUTEURS D'EAU                              |
C |  SMTR          |<-- !  TERMES SOURCES DU TRACEUR                   !
C |  HSTOK         |<-- |  HAUTEURS D'EAU  STOCKEES                    |
C !  HC            ! -->!  H RECONSTRUIT ORDRE 2   CORRIGE             !
C !  HCSTOK        !<-- !  H RECONSTRUIT ORDRE 2   CORRIGE  STOCKE     !
C |  FLUXT         |<-- |  FLUX  TRACEUR  REINITIALISE                 |
C |  FLUHBOR       |<-- |  FLUX  TRACEUR FRONTIERE REINITIALISE        |
C !  DTT           !<-- !  PAS DE TEMPS TRACEUR                        !
C !________________!____!______________________________________________!
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)   
C        -- (TABLEAU DE TRAVAIL)                                        
C-----------------------------------------------------------------------
C
C     - SOUS PROGRAMME(S) APPELANT : RESOLU                             
C 
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NS,NSEG,NPTFR,NTRAC
      DOUBLE PRECISION, INTENT(INOUT) :: DTT
      DOUBLE PRECISION, INTENT(INOUT) :: HSTOK(*),HCSTOK(2,*)
      DOUBLE PRECISION, INTENT(IN)    :: H(*),HC(2,*)
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: SMTR,FLUXT,FLUHBOR
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IS,NSG,ITRAC
C
      DTT = 0.D0
C     
      DO IS=1,NS
        HSTOK(IS)=H(IS)
      ENDDO
C
      DO ITRAC=1,NTRAC
        DO IS=1,NS
          SMTR%ADR(ITRAC)%P%R(IS)=0.D0
        ENDDO
      ENDDO
C
      DO NSG=1,NSEG
        HCSTOK(1,NSG) = HC(1,NSG)
        HCSTOK(2,NSG) = HC(2,NSG)
      ENDDO
C
      DO ITRAC=1,NTRAC
        DO NSG=1,NSEG
          FLUXT%ADR(ITRAC)%P%R(NSG)=0.D0
        ENDDO
      ENDDO
C
      DO ITRAC=1,NTRAC
        DO IS=1,NPTFR
          FLUHBOR%ADR(ITRAC)%P%R(IS)=0.D0
        ENDDO
      ENDDO
C
C-----------------------------------------------------------------------
C 
      RETURN
      END
