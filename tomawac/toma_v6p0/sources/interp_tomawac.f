C                       *************************
                        SUBROUTINE INTERP_TOMAWAC
C                       *************************
C
     * ( F , B , SHP1 , SHP2 , SHP3 , SHZ , ELT , ETA , IKLE2,
     *   ETAP1, NPOIN2 , NELEM2 , NPLAN , TRA01)
C
C***********************************************************************
C  TOMAWAC  VERSION 1.0       01/02/95        F MARCOS (LNH) 30 87 72 66
C***********************************************************************
C
C  FONCTION :
C
C   INTERPOLE AU PIED DES CARACTERISTIQUES
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !    F           !<-->! DENSITE SPECTRALE D'ACTION D'ONDE            !
C !    B           ! -->! FACTEUR DE PROPORTIONNALITE                  !
C !    SHP1-2-3    ! -->! COORDONNEES BARYCENTRIQUES DES NOEUDS DANS   !
C !                !    ! LEURS ELEMENTS 2D "ELT" ASSOCIES.            !
C !    SHZ         ! -->! COORDONNEES BARYCENTRIQUES SUIVANT Z DES     !
C !                !    ! NOEUDS DANS LEURS ETAGES "ETA" ASSOCIES.     !
C !    ELT         ! -->! NUMEROS DES ELEMENTS 2D CHOISIS POUR CHAQUE  !
C !                !    ! NOEUD.                                       !
C !    ETA         ! -->! NUMEROS DES ETAGES CHOISIS POUR CHAQUE NOEUD.!
C !    IKLE2       ! -->! TRANSITION ENTRE LES NUMEROTATIONS LOCALE    !
C !                ! -->! ET GLOBALE                                   !
C !    ETAP1       ! -->! TABLEAU DES ETAGES SUPERIEURS                !
C !    NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE 2D.             !
C !    NELEM2      ! -->! NOMBRE D'ELEMENTS DU MAILLAGE 2D.            !
C !    NPLAN       ! -->! NOMBRE DE DIRECTIONS                         !
C !    TRA01       !<-->! TABLEAU DE TRAVAIL                           !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C APPELE PAR : WAC
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      USE DECLARATIONS_TOMAWAC ,ONLY : MESH
      USE TOMAWAC_MPI
      USE BIEF
C
      IMPLICIT NONE
C
      INTEGER ETAG,ETAGP1,IPOIN1,IPOIN2,IPOIN3,IPLAN
      INTEGER NELEM2,NPOIN2,NPLAN,IP
      INTEGER I,I3D,TYPE
      INTEGER IPOIN
      DOUBLE PRECISION F(NPOIN2,NPLAN),TRA01(NPOIN2,NPLAN)
      DOUBLE PRECISION SHP3(NPOIN2,NPLAN),B(NPOIN2)
      DOUBLE PRECISION SHP1(NPOIN2,NPLAN),SHP2(NPOIN2,NPLAN)
      DOUBLE PRECISION SHZ(NPOIN2,NPLAN),UMSHZ
C
      INTEGER IKLE2(NELEM2,3),ELT(NPOIN2,NPLAN),ETA(NPOIN2,NPLAN)
      INTEGER ETAP1(NPLAN)
      INTEGER ELT2,ETA2
      INTEGER IP2,IP3,IP4,NNSEND,NNRECV
C
C-----------------------------------------------------------------------
C
C
      DO 2 IPLAN=1,NPLAN
        DO 3 IP=1,NPOIN2
      	TRA01(IP,IPLAN)=F(IP,IPLAN)*B(IP)
3       CONTINUE
2     CONTINUE
C
       IF (NCSIZE.GT.1) THEN
C
       IF (.NOT.ALLOCATED(F_SEND)) ALLOCATE(F_SEND(SUM(
     *                                      RECVCOUNTS(:,IFREQ))))
       IF (.NOT.ALLOCATED(F_RECV)) ALLOCATE(F_RECV(SUM(
     *                                      SENDCOUNTS(:,IFREQ))))  
        NNSEND = SUM(SENDCOUNTS(:,IFREQ))
        NNRECV = SUM(RECVCOUNTS(:,IFREQ))
C
        DO IP2 = 1,SUM(RECVCOUNTS(:,IFREQ))
           ELT2 = SH_LOC(IFREQ)%ELT(IP2)
           ETA2 = SH_LOC(IFREQ)%ETA(IP2)
            F_SEND(IP2)%F(1) = TRA01(IKLE2(ELT2,1),ETA2)
            F_SEND(IP2)%F(2) = TRA01(IKLE2(ELT2,2),ETA2)
            F_SEND(IP2)%F(3) = TRA01(IKLE2(ELT2,3),ETA2)
            F_SEND(IP2)%F(4) = TRA01(IKLE2(ELT2,1),ETAP1(ETA2))
            F_SEND(IP2)%F(5) = TRA01(IKLE2(ELT2,2),ETAP1(ETA2))
            F_SEND(IP2)%F(6) = TRA01(IKLE2(ELT2,3),ETAP1(ETA2))
            F_SEND(IP2)%IOR   = RECVCHAR(IP2,IFREQ)%IOR
            F_SEND(IP2)%MYPID = RECVCHAR(IP2,IFREQ)%MYPID
            F_SEND(IP2)%NEPID = RECVCHAR(IP2,IFREQ)%NEPID
            F_SEND(IP2)%SHP1  = SH_LOC(IFREQ)%SHP1(IP2)
            F_SEND(IP2)%SHP2  = SH_LOC(IFREQ)%SHP2(IP2)
            F_SEND(IP2)%SHP3  = SH_LOC(IFREQ)%SHP3(IP2)
            F_SEND(IP2)%SHZ   = SH_LOC(IFREQ)%SHZ(IP2)
            F_SEND(IP2)%BP    = (F_SEND(IP2)%F(1)*F_SEND(IP2)%SHP1 + 
     *                  F_SEND(IP2)%F(2)*F_SEND(IP2)%SHP2 +
     *    F_SEND(IP2)%F(3)*F_SEND(IP2)%SHP3)*(1.d0-F_SEND(IP2)%SHZ) +
     *                  (F_SEND(IP2)%F(4)*F_SEND(IP2)%SHP1 + 
     *                  F_SEND(IP2)%F(5)*F_SEND(IP2)%SHP2 +
     *    F_SEND(IP2)%F(6)*F_SEND(IP2)%SHP3)*F_SEND(IP2)%SHZ
       ENDDO

       CALL GLOB_FONCTION_COMM ()

       ENDIF
C    
      DO 10 IPLAN = 1 , NPLAN
         DO 20 IP = 1 , NPOIN2
         ETAG=ETA(IP,IPLAN)
          ETAGP1=ETAP1(ETAG)
         IPOIN1=IKLE2(ELT(IP,IPLAN),1)
         IPOIN2=IKLE2(ELT(IP,IPLAN),2)
         IPOIN3=IKLE2(ELT(IP,IPLAN),3)
         UMSHZ=1.D0-SHZ(IP,IPLAN)
           F(IP,IPLAN)=
     *      ((TRA01(IPOIN1,ETAG) * SHP1(IP,IPLAN)
     *      + TRA01(IPOIN2,ETAG) * SHP2(IP,IPLAN) 
     *      + TRA01(IPOIN3,ETAG) * SHP3(IP,IPLAN)) * UMSHZ
     *     +( TRA01(IPOIN1,ETAGP1) * SHP1(IP,IPLAN)
     *      + TRA01(IPOIN2,ETAGP1) * SHP2(IP,IPLAN)
     *      + TRA01(IPOIN3,ETAGP1) * SHP3(IP,IPLAN)) 
     *       * SHZ(IP,IPLAN) ) /B(IP)
         IF (F(IP,IPLAN).LT.0.D0) THEN
            F(IP,IPLAN)=0.D0 
         ENDIF    

20      CONTINUE
10    CONTINUE

      IF (NCSIZE.GT.1) THEN

        DO 110 IP = 1 ,NNSEND

           UMSHZ=1.D0-F_RECV(IP)%SHZ
           IPLAN = F_RECV(IP)%IOR/NPOIN2+1
           IP2 = F_RECV(IP)%IOR-(IPLAN-1)*NPOIN2            
           IF (IP2==0) IPLAN = IPLAN-1
           IF (IP2==0) IP2 = NPOIN2
           F(IP2,IPLAN) = 
     *         F_RECV(IP)%BP/B(IP2)
           IF (F(IP2,IPLAN).LT.0.d0) THEN
              F(IP2,IPLAN) = 0.d0
           ENDIF
C
110     CONTINUE
       DEALLOCATE(F_SEND,F_RECV)

         DO IPLAN=1,NPLAN 
           CALL PARCOM2
     *     ( F(:,IPLAN) , 
     *       F(:,IPLAN) , 
     *       F(:,IPLAN) ,
     *       NPOIN2 , 1 , 3 , 1 , MESH )
         ENDDO

       ENDIF
C
      RETURN
      END
