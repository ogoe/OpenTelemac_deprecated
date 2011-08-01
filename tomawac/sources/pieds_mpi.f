      SUBROUTINE PIED_MPI(CX,CY,CT,DT,NRK,X,Y,TETA,IKLE2,IFABOR,ETAP1,
     *                    TRA01,SHP1,SHP2,SHP3,SHZ,JF,ELT,ETA,ITR01,
     *                    GOODELT,NPLAN,NPOIN2,NPOIN3,NF,NELEM2,MESH)
      USE BIEF
!BD_INCKA modif pour //
      USE TOMAWAC_MPI_TOOLS
      USE TOMAWAC_MPI, ONLY : SH_AGAIN,RECVAGAIN,SH_LOC,RECVCHAR,
     *                        NARRV,NCHARA,NLOSTCHAR,NSEND,TEST,
     *                        NCHDIM,NFREQ,IFREQ,ISPDONE,INIT_TOMAWAC,
     *                        PIEDS_TOMAWAC,PIEDS_TOMAWAC_MPI,
     *                        WIPE_HEAPED_CHAR,PREP_INITIAL_SEND,
     *                        GLOB_CHAR_COMM
!BD_INCKA fin modif
      IMPLICIT NONE

      DOUBLE PRECISION CX(NPOIN3,NF) , CY(NPOIN3,NF)
      DOUBLE PRECISION CT(NPOIN3,NF) 
      DOUBLE PRECISION SHP1(NPOIN3,NF) , SHP2(NPOIN3,NF)
      DOUBLE PRECISION SHP3(NPOIN3,NF) , SHZ(NPOIN3,NF)
      DOUBLE PRECISION SHF(NPOIN3,NF) 
      DOUBLE PRECISION X(NPOIN2),Y(NPOIN2)
      DOUBLE PRECISION TETA(NPLAN)
      DOUBLE PRECISION SURDET(NELEM2)
      DOUBLE PRECISION DT,TRA01(NPOIN3,8),PROMIN
      INTEGER ELT(NPOIN3,NF),ETA(NPOIN3,NF)
      INTEGER IKLE2(NELEM2,3),IFABOR(NELEM2,7),ETAP1(NPLAN)
      INTEGER ITR01(NPOIN3,3),JF
      INTEGER NPLAN,NPOIN2,NPOIN3,NF,NELEM2,NRK
      INTEGER LAST_NOMB,NLOSTAGAIN,NUMBER,IER,NRECV,NUMBERLOST
      INTEGER ITE,IP,IPLAN,NBB,IPOIN,GOODELT(NPOIN2,NPLAN)
      INTEGER NARRSUM
      INTEGER P_ISUM,P_IMAX
      EXTERNAL P_ISUM,P_IMAX
      DOUBLE PRECISION :: TEST2(NPOIN3,NF)
!      DOUBLE PRECISION :: TES(NPOIN2,NPLAN)
      TYPE(BIEF_MESH)  ::  MESH


         CALL CORRECT_GOODELT(GOODELT,NPOIN2,NPLAN,MESH)
C
         IF (.not.allocated(NCHARA)) allocate(NCHARA(NF),NLOSTCHAR(NF),
     *                                        NSEND(NF))
         CALL INIT_TOMAWAC(NCHARA(JF),NCHDIM,1,
     *                                       NPOIN3,LAST_NOMB)

C
         IF(.not.allocated(TEST)) allocate(TEST(NPOIN3,NF))
         IFREQ=JF

           CALL PIEDS_TOMAWAC
     *(CX,CY,CT,DT,NRK,X,Y,TETA,IKLE2,IFABOR,
     *  ETAP1,TRA01,TRA01(1,2),
     *  TRA01(1,3),TRA01(1,4),TRA01(1,5),TRA01(1,6),SHP1(1,JF),
     *  SHP2(1,JF),SHP3(1,JF),SHZ(1,JF),ELT(1,JF),ETA(1,JF),
     *  ITR01(1,1),NPOIN3,NPOIN2,
     *  NELEM2,NPLAN,JF,SURDET,-1,ITR01(1,2),MESH%IFAPAR%I,TEST(1,JF), 
     *  NCHDIM,NCHARA(JF),MESH,GOODELT)
!On regarde si une carcateristique en bord de domaine sort et celle 
!de l'autre ne sort pas, alors on ne prendra que la contribution 
!maximal des deux et ainsi on ne traitera pas la caracteristique sortie
!          DO IP = 1,NPOIN2
!              DO IPLAN = 1,NPLAN
!                 TES(IP,IPLAN)  =TEST(IP+NPOIN2*(IPLAN-1),JF)
!              ENDDO
!          ENDDO
         where (TEST(:,JF).LT.0.5d0)
             SHP1(:,JF)=0.d0
             SHP2(:,JF)=0.d0
             SHP3(:,JF)=0.d0
             SHZ(:,JF) = 0.d0
         end where
!          DO IPLAN = 1,NPLAN
!          CALL PARCOM2
!      * ( TES(1,IPLAN) , 
!      *   TES(1,IPLAN) , 
!      *   TES(1,IPLAN) ,
!      *   NPOIN2 , 1 , 2 , 1 , MESH )
!          ENDDO
!          DO IP = 1,NPOIN2
!             DO IPLAN = 1,NPLAN
!                TEST(IP+NPOIN2*(IPLAN-1),JF)=TES(IP,IPLAN)
!             ENDDO
!          ENDDO  
!          where (TEST(:,JF).GT.1.5d0)
!             SHP1(:,JF)=SHP1(:,JF)/TEST(:,JF)
!             SHP2(:,JF)=SHP2(:,JF)/TEST(:,JF)
!             SHP3(:,JF)=SHP3(:,JF)/TEST(:,JF)
!          end where
!Ici on ressort heapchar(nchara,nfreq) et heapcount(ncsize,nfreq)
!heapcount=> nombre de caracteristiques a envoyer sur chaques procs
         CALL WIPE_HEAPED_CHAR(TEST(1,JF),NPOIN3,.TRUE.,NSEND(JF),
     *                        NLOSTCHAR(JF),NCHDIM,
     &                        NCHARA(JF)) 
         
!Pas forcement utile, regarde si Test==1 alors on retire de la liste
! des caracteristiques à envoyer en affectant heapcahr%nepid==-1
!        DO WHILE(P_IMAX(NLOSTCHAR(JF))>0)! THERE ARE -REALLY- LOST TRACEBACKS SOMEWHERE
          CALL PREP_INITIAL_SEND(NSEND,NLOSTCHAR,NCHARA)

!Creation du tabelaux SDISP et organisation des donnees sous forme croissante de SENDCHAR
          CALL GLOB_CHAR_COMM ()
!Envoie de SENDCHAR et ecritrue dans RECVCHAR


!
         IF(.not.allocated(ISPDONE)) allocate(ISPDONE(NPOIN3,NF))
         IF(.not.allocated(NARRV)) allocate(NARRV(NF))
         CALL ALLOC_LOCAL(NARRV(IFREQ),IFREQ,NF,NLOSTAGAIN,
     *                      NUMBERLOST,NARRSUM)

       TEST2(:,JF) = 1.d0
         IF (NUMBERLOST>0) THEN
       CALL PIEDS_TOMAWAC_MPI
     *(CX,CY,CT,DT,NRK,X,Y,TETA,IKLE2,IFABOR,ETAP1,
     *  TRA01,TRA01(1,2),
     *  TRA01(1,3),TRA01(1,4),TRA01(1,5),TRA01(1,6),SH_LOC(JF)%SHP1,
     *  SH_LOC(JF)%SHP2,SH_LOC(JF)%SHP3,SH_LOC(JF)%SHZ,
     *  SH_LOC(JF)%ELT,SH_LOC(JF)%ETA,
     *  NARRV(JF),NPOIN2,
     *  NELEM2,NPLAN,JF,SURDET,-1,MESH%IFAPAR%I, 
     *  2,NCHARA(JF),RECVCHAR(1,JF))

        CALL ALLOC_AGAIN(NARRV(IFREQ),IFREQ,NLOSTAGAIN,NUMBERLOST,
     *                   NUMBER)
        CALL ORGANIZE_SENDAGAIN()

        CALL SUPP_ENVOI_AGAIN(IFREQ,NUMBER)

!
           ITE = 0
          DO WHILE((NUMBERLOST>0).AND.(ITE.LE.20)) 
           ITE= ITE + 1
          CALL ORGANIZE_SENDAGAIN()
          CALL ENVOI_AGAIN(NRECV)
          TEST2(:,JF)=1.d0
          CALL PIEDS_TOMAWAC_MPI
     *(CX,CY,CT,DT,NRK,X,Y,TETA,IKLE2,IFABOR,
     *  ETAP1,TRA01,TRA01(1,2),
     *  TRA01(1,3),TRA01(1,4),TRA01(1,5),TRA01(1,6),SH_AGAIN%SHP1,
     *  SH_AGAIN%SHP2,SH_AGAIN%SHP3,SH_AGAIN%SHZ,
     *  SH_AGAIN%ELT,SH_AGAIN%ETA,
     *  NRECV,NPOIN2,
     *  NELEM2,NPLAN,JF,SURDET,-1,MESH%IFAPAR%I, 
     *  2,NCHARA(JF),RECVAGAIN)
        CALL INCREM_ENVOI_RECV(IFREQ,NUMBER,NLOSTAGAIN,NUMBERLOST,
     *                         NRECV)
        ENDDO ! fin sur la boucle dowhile
         CALL FINAL_ORGA_RECV(NARRV(IFREQ),IFREQ)
          ELSE
           CALL RESET_COUNT(IFREQ)
          ENDIF
      END SUBROUTINE PIED_MPI

C                       *****************
      SUBROUTINE PRELEO_MPI
C                       *****************
     *(XLEO,YLEO,NLEO,X,Y,IKLE,SURDET,NPOIN2,NELEM2,NOLEO,ISLEO)
C
C***********************************************************************
C  COWADIS  VERSION 1.0     1/02/95    F. MARCOS (LNH) 30 87 72 66
C***********************************************************************
C
C     FONCTION  : CHOISI LES POINTS DE CALCUL LES PLUS PROCHES
C                 DES POINTS DE SORTIE DEMANDES
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C !      NOM       !MODE!                   ROLE                       !
C !________________!____!______________________________________________!
C !   XLEO         ! -->! TABLEAU DES ABSCISSES DES POINTS DE SORTIE   !
C !   YLEO         ! -->! TABLEAU DES ORDONNEES DES POINTS DE SORTIE   !
C !   NLEO         ! -->! NOMBRE DE POINTS DE SORTIE                   !
C !   X            ! -->! ABSCISSES DES POINTS                         !
C !   Y            ! -->! ORDONNEES DES POINTS                         !
C !   IKLE         ! -->! CONNECTIVITE NOEUD ELEMENTS                  !
C !   SURDET       ! -->! 1/SUPERFICIE ELEMENTS                        !
C !   NPOIN2       ! -->! NOMBRE DE POINTS 2D                          !
C !   NELEM2       ! -->! NOMBRE D'ELEMENTS 2D                         !
C !   NCSIZE       ! -->! NOMBRE D'ELEMENTS 2D                         !
C !   NOLEO        !<-->! TABLEAU DES NUMERO DES POINTS CHOISIS        !
C !   NOPID        ! <--! TABLEAU NUMERO POINT / PROCS                 !
C !________________!____!______________________________________________!
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : COW
C
C***********************************************************************
C
      USE BIEF
      IMPLICIT NONE

      COMMON/ECRSPE_MPI/SPE_SEND
      INTEGER SPE_SEND
C
      INTEGER I,ILEO,NLEO,NPOIN2,NELEM2,IELEM,NOELEM
C
      DOUBLE PRECISION X(NPOIN2)  , Y(NPOIN2)
      DOUBLE PRECISION XLEO(NLEO)  , YLEO(NLEO)
      DOUBLE PRECISION SURDET(NELEM2)
      DOUBLE PRECISION DIST,DIST2,SHP1,SHP2,SHP3
C
      INTEGER NOLEO(NLEO),N1G,N2G,N3G
      INTEGER IKLE(NELEM2,3)
      LOGICAL ISLEO(NLEO)
C
C-----------------------------------------------------------------------
C
!       DO 10 ILEO=1,NLEO
!         DIST=1.D99
!         DO 20 I=1,NPOIN2
!          DIST2=(XLEO(ILEO)-X(I))**2+(YLEO(ILEO)-Y(I))**2
!          IF (DIST2.LT.DIST) THEN
!              DIST=DIST2
!              NOLEO(ILEO)=I
!          ENDIF
! 20      CONTINUE
! 10    CONTINUE
       SPE_SEND = 0
       ISLEO = .FALSE.
       NOLEO = 1
       DO ILEO = 1,NLEO
          NOELEM = 0
          DO 20 IELEM = 1,NELEM2
             N1G=IKLE(IELEM,1)
             N2G=IKLE(IELEM,2)
             N3G=IKLE(IELEM,3)
               SHP1 = ((X(N3G)-X(N2G))*(YLEO(ILEO)-Y(N2G))
     *               -(Y(N3G)-Y(N2G))*(XLEO(ILEO)-X(N2G)))*SURDET(IELEM)
               SHP2 = ((X(N1G)-X(N3G))*(YLEO(ILEO)-Y(N3G))
     *               -(Y(N1G)-Y(N3G))*(XLEO(ILEO)-X(N3G)))*SURDET(IELEM)
               SHP3 = ((X(N2G)-X(N1G))*(YLEO(ILEO)-Y(N1G))
     *               -(Y(N2G)-Y(N1G))*(XLEO(ILEO)-X(N1G)))*SURDET(IELEM)
             IF ((SHP1.GE.0.d0).AND.(SHP2.GE.0.d0)
     *                                        .AND.(SHP3.GE.0.d0)) THEN
               ISLEO(ILEO) = .TRUE.
               NOELEM = IELEM
               IF (SHP2>SHP1) THEN
                  NOLEO(ILEO) = N2G
                  IF (SHP3>SHP2) NOLEO(ILEO) = N3G
               ELSE
                  NOLEO(ILEO) = N1G
                  IF (SHP3>SHP1) NOLEO(ILEO) = N3G
               ENDIF
               SPE_SEND=SPE_SEND+1 
               GOTO 30
              ENDIF
20         CONTINUE
30       CONTINUE            
       ENDDO

C
C-----------------------------------------------------------------------
C
      RETURN
      END


