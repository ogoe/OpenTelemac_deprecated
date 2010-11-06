!                       ******************
                        SUBROUTINE ELEB3DT                              
!                       ******************
!
     *(IKLE3 , NBOR   , KP1BOR, NELBOR, IKLBOR, NULONE,                              
     * NELEM2, NELMAX2, NPOIN2, NPLAN , NETAGE, NPTFR  )
!
!***********************************************************************
! BIEF VERSION 5.3         23/08/99      J-M HERVOUET(LNH) 30 87 80 18
!                                                    
!***********************************************************************
!
!    FONCTION: CAS DES PRISMES DECOUPES EN TETRAEDRES
!    =========
!
!    CONSTRUCTION DU MAILLAGE 3D : A L'ENTREE, LES TABLEAUX DU MAILLAGE
!    3D REMPLIS PAR UN APPEL PREALABLE DE ELEBD. A LA SORTIE, LES
!    TABLEAUX COMPLETES EN 3D.
!
!----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________.
! !  NOM           !MODE!                  ROLE                        !
! !________________!____!______________________________________________!
! !  IKLE3         !<-- ! CORRESPONDANCE LOCALE - GLOBALE EN 3D        !
! !  SURFAC        !<-->! SURFACE DES TRIANGLES ETENDUE EN 3D          !
! !  NBOR3         !<-- ! CORRESPONDANCE ENTRE LA NUMEROTATION DE BORD !
! !                !    ! ET LA NUMEROTATION GLOBALE (3D)              !
! !  NBOR          ! -->! CORRESPONDANCE ENTRE LA NUMEROTATION DE BORD !
! !                !    ! ET LA NUMEROTATION GLOBALE (2D)              !
! !  KP1BOR        ! -->! PT FRONT. SUIVANT LE PT FRONT. CONSIDERE     !
! !  NELBOR        ! -->! NUMERO GLOBAUX DES ELEMENTS DE BORD          !
! !  IKLBOR        !<-- ! TABLE DE CONNECTIVITE ELEMENTS DE BORD       !
! !  NELBO3        !<-- ! ASSOCIE A CHAQUE FACE DE BORD L'ELEMENT 3D   !
! !                !    ! AUQUEL ELLE APPARTIENT                       !
! !  NULONE        !<-- ! ASSOCIE LA NUMEROTATION LOCALE DE BORD A LA  !
! !                !    ! NUMEROTATION LOCALE 3D                       !
! !  IKLE2         ! -->! CORRESPONDANCE LOCALE - GLOGALE EN 2D        !
! !  NELEM2        ! -->! NOMBRE D'ELEMENTS EN 2D                      !
! !  NPOIN2        ! -->! NOMBRE DE POINTS 2D                          !
! !  NPLAN         ! -->! NOMBRE DE PLANS HORIZONTAUX                  !
! !  NETAGE        ! -->! NPLAN - 1                                    !
! !  NPTFR         ! -->! NOMBRE DE POINTS DE BORD 2D                  !
! !________________!____!______________________________________________!
! MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
!
!-----------------------------------------------------------------------
!
! SOUS-PROGRAMME APPELE PAR : MITRID
! SOUS-PROGRAMMES APPELES : OV , PLANTE
!
!***********************************************************************
!
      USE BIEF, EX_ELEB3DT => ELEB3DT
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C                              
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NELEM2,NPOIN2,NPLAN,NETAGE,NPTFR,NELMAX2
      INTEGER, INTENT(INOUT) :: IKLE3(NELEM2,3,NETAGE,4)
      INTEGER, INTENT(INOUT) :: IKLBOR(NPTFR,2,NETAGE,3)
      INTEGER, INTENT(INOUT) :: NULONE(NPTFR,2,NETAGE,3)
      INTEGER, INTENT(INOUT) :: NELBOR(NPTFR,2,NETAGE),NBOR(NPTFR*NPLAN)
      INTEGER, INTENT(INOUT) :: KP1BOR(NPTFR)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      LOGICAL OK(2)
!
      INTEGER IELEM,IPOIN,T(3)
      INTEGER IETAGE,IPTFR,IL1,IL2,IL3,IL4,IG(2,2,3),IL(2,2,3)
      INTEGER IG1,IG2,IG3,IG4
      INTEGER NUM1(12),NUM2(12),NUM3(12),K,L,M,N
!
      DATA NUM1 / 1 , 2 , 4 , 1 , 3 , 2 , 2 , 3 , 4 , 3 , 1 , 4 /
      DATA NUM2 / 2 , 4 , 1 , 3 , 2 , 1 , 3 , 4 , 2 , 1 , 4 , 3 /
      DATA NUM3 / 4 , 1 , 2 , 2 , 1 , 3 , 4 , 2 , 3 , 4 , 3 , 1 /
!
!***********************************************************************
!
! TABLES DE CONNECTIVITE POUR LES FACES DE BORDS --> IKLBOR , NBOR3 ,
! CORRESPONDANCE NUMEROTATION LOCALE BORD NUMEROTATION LOCALE 3D
! --> NULONE
! ET CALCUL DE NELBO3
!
! BORDS LATERAUX
! POUR CHAQUE FACE RECTANGULAIRE DECOUPEE EN DEUX TRIANGLES
! LE TRIANGLE BAS EST LE NUMERO 1, LE HAUT NUMERO 2.
!
      DO IPTFR = 1,NPTFR
!
!        NUMERO DU TRIANGLE TOUCHANT LE COTE EN 2D
         IELEM = NELBOR(IPTFR,1,1)
!
!        NUMERO GLOBAL DU POINT BAS-GAUCHE DE LA FACE RECTANGULAIRE
!        DU PREMIER ETAGE (LE MOTIF SE REPETE ENSUITE)
         IPOIN = NBOR(IPTFR)
!
         DO IETAGE = 1,NETAGE
!
!           NUMEROTATION DE BORD 3D DES 4 POINTS DE LA FACE RECTANGULAIRE
!           
            IL1 =       IPTFR   + (IETAGE-1)*NPTFR
            IL2 = KP1BOR(IPTFR) + (IETAGE-1)*NPTFR
            IL3 = IL2 + NPTFR
            IL4 = IL1 + NPTFR
!
!           NUMEROTATION GLOBALE 3D DES 4 POINTS DE LA FACE RECTANGULAIRE
!           
            IG1 =        NBOR(IPTFR)  + (IETAGE-1)*NPOIN2
            IG2 = NBOR(KP1BOR(IPTFR)) + (IETAGE-1)*NPOIN2
            IG3 = IG2 + NPOIN2
            IG4 = IG1 + NPOIN2
!
! NUMEROS DES 3 TETRAEDRES POUVANT TOUCHER LA FACE
!
            T(1) = (IETAGE-1)*3*NELEM2+IELEM
            T(2) = T(1) + NELEM2
            T(3) = T(2) + NELEM2
!
! RECHERCHE DU TRIANGLE BAS (PEUT ETRE 1-2-4 OU 1-2-3)
!
!           2 FORMES POSSIBLES DU TRIANGLE BAS (EN GLOBAL ET DE BORD)
            IG(1,1,1)=IG1
            IG(1,1,2)=IG2
            IG(1,1,3)=IG4
            IG(1,2,1)=IG1
            IG(1,2,2)=IG2
            IG(1,2,3)=IG3
            IL(1,1,1)=IL1
            IL(1,1,2)=IL2
            IL(1,1,3)=IL4
            IL(1,2,1)=IL1
            IL(1,2,2)=IL2
            IL(1,2,3)=IL3
!           2 FORMES POSSIBLES DU TRIANGLE HAUT (EN GLOBAL ET DE BORD)
            IG(2,1,1)=IG1
            IG(2,1,2)=IG3
            IG(2,1,3)=IG4
            IG(2,2,1)=IG2
            IG(2,2,2)=IG3
            IG(2,2,3)=IG4
            IL(2,1,1)=IL1
            IL(2,1,2)=IL3
            IL(2,1,3)=IL4
            IL(2,2,1)=IL2
            IL(2,2,2)=IL3
            IL(2,2,3)=IL4
!
            OK(1)=.FALSE.
            OK(2)=.FALSE.
!
!           K=1 TRIANGLE BAS   K=2 TRIANGLE HAUT
            DO K=1,2
!           2 DECOUPAGES POSSIBLES
            DO L=1,2
!           12 FACONS POUR UN TETRAEDRE DE PRESENTER SES FACES
            DO M=1,12
!           3 TETRAEDRES POSSIBLES
            DO N=1,3
              IF(IG(K,L,1).EQ.IKLE3(IELEM,N,IETAGE,NUM1(M)).AND.
     *           IG(K,L,2).EQ.IKLE3(IELEM,N,IETAGE,NUM2(M)).AND.
     *           IG(K,L,3).EQ.IKLE3(IELEM,N,IETAGE,NUM3(M))) THEN
!
                  IKLBOR(IPTFR,K,IETAGE,1) = IL(K,L,1)
                  IKLBOR(IPTFR,K,IETAGE,2) = IL(K,L,2)
                  IKLBOR(IPTFR,K,IETAGE,3) = IL(K,L,3)
                  NELBOR(IPTFR,K,IETAGE)   = T(N)
                  NULONE(IPTFR,K,IETAGE,1) = NUM1(M) 
                  NULONE(IPTFR,K,IETAGE,2) = NUM2(M) 
                  NULONE(IPTFR,K,IETAGE,3) = NUM3(M)
!
                  OK(K) = .TRUE.
!
              ENDIF
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            IF(.NOT.OK(1).OR..NOT.OK(2)) THEN
            WRITE(LU,*) 'PB IN ELEB3DT IELEM=',IELEM,' IPTFR=',IPTFR
            CALL PLANTE(1)
            STOP
            ENDIF
!
!           NUMEROS GLOBAUX DES POINTS DE BORDS LATERAUX
!           TOUS PLANS SAUF SURFACE           
            NBOR(IPTFR +(IETAGE-1)*NPTFR)=IPOIN+(IETAGE-1)*NPOIN2
!
         ENDDO
      ENDDO
!
! NUMEROS GLOBAUX DES POINTS DE BORDS LATERAUX : EN SURFACE
!
      DO IPTFR = 1,NPTFR
         NBOR(IPTFR+(NPLAN-1)*NPTFR)=NBOR(IPTFR)+(NPLAN-1)*NPOIN2
      ENDDO
!
!-----------------------------------------------------------------------
!
      RETURN
      END
