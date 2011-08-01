C                       *****************
                        SUBROUTINE CONLIT
C                       *****************
C 
     *(NBOR)
C
C
C***********************************************************************
C SISYPHE VERSION 5.9                             E. PELTIER    11/09/95
C                                                 C. LENORMANT
C                                                 J.-M. HERVOUET
C                                                 C. MACHET     07/06/02
C 
C 19/06/2008 : CV : PRISE EN COMPTE DE CBOR_VASE ET CBOR_SABLE
C                                               
C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT   
C***********************************************************************
C
C     FONCTION  : TRAITE LES POINTS D'ENTREE ET DE DIRICHLET.
C
C                 PERMET DE MODIFIER LES TABLEAUX DE CONDITIONS   
C                 AUX LIMITES DANS LE CAS OU ELLES SONT VARIABLES 
C                 EN TEMPS. 
C                     
C                 PERMET D'IMPOSER LE DEBIT SOLIDE EN CERTAINS POINTS
C                 FRONTIERES (TABLEAUX QBOR ET LIQBOR). ATTENTION, ON
C                 DOIT ALORS SPECIFIER EN CES POINTS LIEBOR = KSORT !
C 
C      FUNCTION : IMPOSED BOUNDARY CONDITIONS 
C                 THIS SUBROUTINE ALLOWS TO IMPOSE 
C                             TIME VARYING BOUNDARY CONDITIONS
C                 (CONSTANT VALUES CAN BE DIRECTLY IMPOSED IN CONDIM INPUT FILE)
C          
C                  IT IS ALSO POSSIBLE TO IMPOSE A SAND TRANSPORT RATE IN
C                  A NUMBER OF BOUNDARY POINTS - IT IS THEN NECESSARY TO  
C                  ALSO IMPOSE LIEBOR = KSORT !
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |     EBOR       |<-->| IMPOSED BED EVOLUTION AT THE BOUNDARY
C |     QBOR       |<-->| IMPOSED SOLID TRANSPORT AT THE BOUNDARY
C |     CBOR       |<-->| IMPOSED SUSPENDED SAND CONC AT THE BOUNDARY
C
C |     LIEBOR     |<-->| TYPE OF BOUNDARY CONDITIONS ON BED EVOLUTION 
C |     LIQBOR     |<-->| TYPE OF BOUNDARY CONDITIONS ON SAND TRANSPORT RATE 
C |     LICBOR     |<-->| TYPE OF BOUNDARY CONDITIONS ON SUSPENDED SAND CONC 
C 
C |     NBOR       | -->| GLOBAL NUMBER OF BOUNDARY POINT
C |                |    | 
C |     NPOIN      | -->| NUMBER OF 2D POINTS 
C |     NPTFR      | -->| NUMBER OF BOUNDARY POINTS
C 
C |   KENT,KSORT   | -->| TYPES OF
C |   KADH,KLOG    | -->| BOUNDARY 
C | KNEU,KDIR,KDDL | -->| CONDITIONS
C
C |   
C |________________|____|______________________________________________
C MODE : -->(INPUT), <--(OUTPUT), <-->(MODIFIED DATA)
C-----------------------------------------------------------------------
C CALLED BY  : SISYPHE
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_SISYPHE
      USE DECLARATIONS_TELEMAC
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN):: NBOR(NPTFR)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I,K,IFRLIQ,IRANK
C
C-----------------------------------------------------------------------
C
      DO  K=1,NPTFR
C
        I = NBOR(K)
C
C
C  DIRICHLET CONDITIONS
C  +++++++++++++++++++++
C
        IF(LIEBOR%I(K).EQ.KADH) THEN
          LIEBOR%I(K)= KLOG
        ENDIF
C
C IMPOSED SOLID DISCHARGE - FREE BED EVOLUTION
C ++++++++++++++++++++++++++++++++++++++++++++
C QBOR%ADR(J)%P%R(K) IS THE SOLID DISCHARGE IMPOSED AT THE FRONTIER 
C                   NODE K , CLASS OF SEDIMENT J
C  
C               LIEBOR%I(K)=KSORT
c               LIQBOR%I(K)=KENT
C
c               QBOR%ADR(1)%P%R(K)=1.D-4
C               QBOR%ADR(2)%P%R(K)=1.D-4 .....
c 
c  IMPOSED BED EVOLUTON 
C +++++++++++++++++++++
c          IF (LIEBOR%I(K).EQ.KENT) THEN
c               EBOR%ADR(1)%P%R(K)=1.D-4
c               EBOR%ADR(2)%P%R(K)=1.D-4.....
c         ENDIF
C
       ENDDO
C
C-----------------------------------------------------------------------
C     LICBOR : BOUNDARY CONDITION FOR SEDIMENT CONCENTRATION
C-----------------------------------------------------------------------
     
      IF(SUSP) THEN 
C
        DO K=1,NPTFR
C
C         SO FAR LICBOR=LIEBOR (WITH KADH CHANGED INTO KLOG, SEE ABOVE,
C                               BUT CAN BE CHANGED)
C
          LICBOR%I(K) = LIEBOR%I(K)
C
C         ENTRANCE : IMPOSED CONCENTRATION
C         -------------------------------
C
C         NOTE JMH: KSORT MUST BE TREATED ALSO BECAUSE SUBROUTINE DIFFIN
C                   MAY CHANGE A KSORT INTO KENT, DEPENDING OF FLOW
C
          IFRLIQ=NUMLIQ%I(K)        
          IF(LIEBOR%I(K).EQ.KENT.OR.LIEBOR%I(K).EQ.KSORT) THEN 
            DO I=1,NSICLA
               IRANK=I+(IFRLIQ-1)*NSICLA
               CBOR%ADR(I)%P%R(K) = CBOR_CLASSE(IRANK)             
            ENDDO
          ENDIF
C
        ENDDO  
C                 
      ENDIF        
C
C-----------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE CONLIT
