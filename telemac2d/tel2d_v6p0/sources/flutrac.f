C                       ******************                            
                        SUBROUTINE FLUTRAC
C                       ******************                             
C                                                                       
     *(NSEG,NPTFR,DT,FLUXT,FLUHBOR,FLUXTEMP,FLUHBTEMP,DTT)
C                                                                       
C***********************************************************************
C  TELEMAC 2D VERSION 5.4                                         INRIA
C-----------------------------------------------------------------------
C
C      INCREMENTE LES FLUX TRACEUR DES FLUX D'UN PAS DE TEMPS HYDRO
C-----------------------------------------------------------------------
C                             ARGUMENTS                                 
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C |  NSEG          | -->|  NOMBRE D'ARETES DU MAILLAGE                 |
C |  NPTFR         | -->|  NOMBRE DE POINTS FRONTIERE                  |
C |  DT            | -->|  PAS DE TEMPS HYDRO                          |
C |  FLUXT         |<-->|  FLUX  TRACEUR INCREMENTE                    |
C |  FLUHBOR       |<-->|  FLUX  TRACEUR FRONTIERE INCREMENTE          |
C |  FLUXTEMP      | -->|  FLUX D'UN PAS DE TEMPS HYDRO                |
C |  FLUHBTEMP     | -->|  FLUX FRONTIERE D'UN PAS DE TEMPS HYDRO      |
C !  DTT           !<-->!  PAS DE TEMPS TRACEUR                        !
C !________________!____!______________________________________________!
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)   
C        -- (TABLEAU DE TRAVAIL)                                        
C-----------------------------------------------------------------------
C     - SOUS PROGRAMME(S) APPELANT : RESOLU                             
C 
C***********************************************************************
C
      IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NSEG, NPTFR
      DOUBLE PRECISION, INTENT(IN)    :: DT,FLUXTEMP(*),FLUHBTEMP(*)
      DOUBLE PRECISION, INTENT(INOUT) :: FLUHBOR(*),FLUXT(*),DTT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IS,NSG
C
C-----------------------------------------------------------------------
C
      DTT=DTT+DT
C
      DO NSG=1,NSEG 
        FLUXT(NSG) = FLUXT(NSG) + DT * FLUXTEMP(NSG)
      ENDDO
C
      DO IS=1,NPTFR
        FLUHBOR(IS) = FLUHBOR(IS) + DT * FLUHBTEMP(IS)
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
