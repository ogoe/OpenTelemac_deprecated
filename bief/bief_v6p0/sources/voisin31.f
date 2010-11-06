!                       *******************
                        SUBROUTINE VOISIN31
!                       *******************
!
     *(IFABOR,NELEM,NELMAX,IELM,IKLE,SIZIKL,
     * NPOIN,NACHB,NBOR,NPTFR,LIHBOR,KLOG,IKLESTR,NELEMTOTAL,NELEB2)
!
!***********************************************************************
! BIEF VERSION 5.8      22/01/08    REGINA NEBAUER (LNHE) 01 30 87 83 93
!
!***********************************************************************
!
!    FONCTION : CONSTRUCTION DU TABLEAU IFABOR, OU IFABOR(IELEM,IFACE)
!               EST LE NUMERO GLOBAL DU VOISIN DE LA FACE IFACE DE
!               L'ELEMENT IELEM SI CE VOISIN EXISTE ET 0 SI LA FACE EST
!               SUR LA FRONTIERE DU DOMAINE.
!
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________.
! |      NOM       |MODE|                   ROLE                       |
! |________________|____|______________________________________________|
! |    IFABOR      |<-- | TABLEAU DES VOISINS DES FACES.
! |                |    | (CAS DES MAILLAGES ADAPTATIFS)
! |    IELM        | -->| 31: TETRAEDRES NON STRUCTURES
! |    NELEM       | -->| NOMBRE D'ELEMENTS DANS LE MAILLAGE.
! |    NELMAX      | -->| NOMBRE MAXIMUM D'ELEMENTS DANS LE MAILLAGE.
! |    NPOIN       | -->| NOMBRE TOTAL DE POINTS DU DOMAINE
! |    NPTFR       | -->| NOMBRE DE POINTS DE BORD
! |    IKLE        | -->| Table de connectivite domaine
! |    SIZIKLE     | -->| ??
! |    NBOR        | -->| Correspondance no noeud de bord/no global
! |    LIHBOR      | -->| Type de CL par noeud
! |    KLOG        | -->| ????
! |    NACHB       | -->| TABLEAU DE VOISINAGE POUR PARALLELISME
! !  IKLETR,NELEB2 | -->/ CONNECTIVITE DES TRIA DE BORD POUR ESTEL3D
! |________________|____|_______________________________________________
!  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
!***********************************************************************
!
      USE BIEF !, EX_VOISIN31 => VOISIN31
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN   ) :: IELM 
      INTEGER, INTENT(IN   ) :: NPTFR
      INTEGER, INTENT(IN   ) :: NELEM
      INTEGER, INTENT(IN   ) :: NELMAX
      INTEGER, INTENT(IN   ) :: NPOIN
      INTEGER, INTENT(IN   ) :: SIZIKL
      INTEGER, INTENT(IN   ) :: NBOR(NPTFR)
      INTEGER, INTENT(IN   ) :: NACHB(NBMAXNSHARE,NPTIR)
      ! NOTE : on donne explicitement la deuxieme dimension de IFABOR et
      ! IKLE, car il s'agit ici toujours de tetraedres!
      INTEGER, INTENT(INOUT) :: IFABOR(NELMAX,4)
      INTEGER, INTENT(IN   ) :: IKLE(SIZIKL,4)
      INTEGER, INTENT(IN   ) :: LIHBOR(NPTFR)
      INTEGER, INTENT(IN   ) :: KLOG
      INTEGER, INTENT(IN   ) :: NELEMTOTAL
      INTEGER, INTENT(IN   ) :: IKLESTR(NELEMTOTAL,3)
      INTEGER, INTENT(IN   ) :: NELEB2
!
! VARIABLES LOCALES 
!-----------------------------------------------------------------------

      ! Le tableau qui est l'inverse de NBOR (ca donne pour chaque noeud
      ! du domaine le numero de noeud e bord, ou zero si le noeud est a
      ! l'interieur du domaine.
      INTEGER, DIMENSION(:  ), ALLOCATABLE :: NBOR_INV
      ! Le tableau definissant le nombre d'element (tetraedres) voisins
      ! d'un noeud.
      INTEGER, DIMENSION(:  ), ALLOCATABLE :: NVOIS
      ! Le tableau definissant les identifiants des elements voisins de
      ! chaque noeud. 
      INTEGER, DIMENSION(:  ), ALLOCATABLE :: NEIGH

      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IKLE_TRI

      INTEGER, DIMENSION(:,:), ALLOCATABLE :: VOIS_TRI

      ! Un tableau definissant les adresses des differents entrees dans
      ! le tableau NEIGH
      INTEGER, DIMENSION(NPOIN)            :: IADR
      ! La valeur d'une entree dans ce tableau.
      INTEGER                              :: ADR
!
      ! Le nombre maximal d'elements voisin d'un noeud.
      INTEGER :: NMXVOISIN
      INTEGER :: IMAX       ! Dimensionnement tableau IADR
!
      INTEGER :: NFACE      ! Nombre de faces de l'element (tetra : 4)
      INTEGER :: NBTRI      ! Le nombre de triangles definis
!
      INTEGER :: IELEM      ! Compteur elements
      INTEGER :: IELEM2     ! Compteur elements
      INTEGER :: IPOIN      ! Compteur noeuds domaine
      INTEGER :: INOEUD     ! Compteur noeuds tetraedres/triangles
      INTEGER :: IFACE      ! Compteur face
      INTEGER :: IFACE2     ! Compteur face
      INTEGER :: ITRI       ! Compteur trianlges
      INTEGER :: IVOIS      ! Compteur voisins
      INTEGER :: NV         ! Nombre de voisins
!
      INTEGER :: ERR        ! Code d'erreur allocation memoire
!
      LOGICAL :: found      ! trouve ou pas ...
!
      INTEGER :: I1, I2, I3 ! Les trois noeuds d'un triangle   
      INTEGER :: M1, M2, M3 ! La meme chose en ordonne.
!
      INTEGER :: I,J,K      ! ca sert ...
!
      INTEGER :: IR1,IR2,IR3,IR4,IR5,IR6,COMPT
      LOGICAL :: BORD
!
!   ~~~~~~~~~~~~~~~~~~~~~~~   
!     Definition des quatre triangles du tetraedre : la premiere
!     dimension du tableau est le numero du triangle, la deuxieme donne
!     les numeros des noeuds de tetraedres qui le definissent.
      INTEGER SOMFAC(3,4)
      DATA SOMFAC /  1,2,3 , 4,1,2 , 2,3,4 , 3,4,1   /
!-----------------------------------------------------------------------
! Debut du code 
!-----------------------------------------------------------------------
!
!  
!
! D'abord on verifie qu'on est bien dans le cas des tetraedres. Sinon,
! goodbye.
!
      IF(IELM.EQ.31) THEN
       NFACE = 4
      ELSE
       IF(LNG.EQ.1) WRITE(LU,98) IELM
       IF(LNG.EQ.2) WRITE(LU,99) IELM
98     FORMAT(1X,'VOISIN31: IELM=',1I6,' TYPE D''ELEMENT NON PREVU')
99     FORMAT(1X,'VOISIN31: IELM=',1I6,' TYPE OF ELEMENT NOT AVAILABLE')
       CALL PLANTE(1)
       STOP
      ENDIF
!
! Allocation of the temporary arrays
      allocate(nbor_inv(npoin),stat=err)
      if(err.ne.0) then
        if(lng.eq.1) then
          write(lu,*) 'voisin31 : allocation de nbor_inv defectueuse'
        endif
        if(lng.eq.2) then
          write(lu,*) 'voisin31 : wrong allocation of nbor_inv'
        endif
        stop
      endif
!
! Allocation of the temporary arrays
      allocate(nvois(npoin),stat=err)
      if(err.ne.0) then
        if(lng.eq.1) then
          write(lu,*) 'voisin31 : allocation de nvois defectueuse'
        endif
        if(lng.eq.2) then
          write(lu,*) 'voisin31 : wrong allocation of nvois'
        endif
        stop
      endif
!
!-----------------------------------------------------------------------
! ETAPE 1 : Comptage du nombre d'elements voisins d'un noeud.
!-----------------------------------------------------------------------
! Calcul du nombre d'elements voisins pour chaque noeud du maillage.
! RESULTAT : NVOIS(INOEUD) donne le nombre d'elements voisins pour le
! noeud INOEUD

      ! On commence avec l'initialisation a 0 du compteur des elements
      ! voisins.
      DO I = 1, NPOIN
        NVOIS(I) = 0
      END DO   
      ! Puis on compte les elements voisins.
      ! En parcourant la table de connectivite, on incremente le
      ! compteur a chaque fois qu'un element reference le noeud IPOIN

      ! Boucle sur les 4 noeuds de l'element
      DO INOEUD = 1, 4   
        ! Boucle sur les elements
        DO IELEM = 1,NELEM
          ! L'id du noeud I de l'element IELEM
          IPOIN        = IKLE( IELEM , INOEUD )
          ! Incrementer le compteur.
          NVOIS(IPOIN) = NVOIS(IPOIN) + 1
        END DO
      END DO

!-----------------------------------------------------------------------
! ETAPE 2 : Determination de la taille du tableau NEIGH() et de la
! table auxiliaire pour indexer NEIGH. allocation de NEIGH
!-----------------------------------------------------------------------
! On va dans la suite creer un tableau qui va contenir les identifiant
! des elements voisins de chaque noeud. Comme le nombre de voisins est
! a priori different pour chaque noeud, et comme on ne veut pas maximise
! le tableau (trop gros), on a besoin d'un tableau auxiliaire qui va
! nous donner l'adresse des entrees pour un noeud donnee. Ce tableau a
! autant d'entree que de noeud.
! Et a l'occasion, on va aussi calculer le nombre maximum de voisins.

      ! La premiere entree dans le tableau des id des voisins est 1.
      ADR       = 1
      IADR(1)   = ADR
      ! Le nombre max d'elements voisins 
      NV        = NVOIS(1)
      NMXVOISIN = NV

      DO IPOIN = 2,NPOIN
          ! Calcul de l'adresse des autres entrees:
          ADR         = ADR + NV
          IADR(IPOIN) = ADR
          NV          = NVOIS(IPOIN)
          ! Reperage du nombre de voisins max.
          NMXVOISIN   = MAX(NMXVOISIN,NV)
      END DO

      ! Le nombre total d'elements voisins pour tous les noeuds donne la
      ! taille du tableau des voisins :

      IMAX = IADR(NPOIN) + NVOIS(NPOIN)

      ! Allocation de la table contenant les identifiants des elements
      ! voisins pour chaque noeud.
      ALLOCATE(NEIGH(IMAX),STAT=ERR)
      IF(ERR.NE.0) GOTO 999
!
!-----------------------------------------------------------------------
! ETAPE 3 : initialisation de NEIGH
!-----------------------------------------------------------------------
! Apres allocation du tableu NEIGH, il faut le remplir : donc on
! recommence la boucle sur les quatre noeuds de chaque element et cette
! fois ci, on note en meme temps l'identifiant dans le tableau NEIGH.
!
!
      ! Reinitialisation du compteur des elements voisins a 0, pour
      ! savoir ou on en est.
      NVOIS(:) = 0

      ! Pour chaque noeud des elements, on recupere son identifiant.
      DO inoeud = 1, 4  ! Boucle sur les noeuds de l'element
        DO IELEM=1,NELEM ! Boucle sur les elements
          IPOIN     = IKLE( IELEM , INOEUD )
          ! On a un voisin en plus.
          NV           = NVOIS(IPOIN) + 1
          NVOIS(IPOIN) = NV
          ! On note l'identifiant de l'element voisin dans le tableau
          NEIGH(IADR(IPOIN)+NV) = IELEM
        END DO ! fin boucle elements
      END DO  ! fin boucle noeuds
!
!-----------------------------------------------------------------------
! ETAPE 4 : Reperer les faces communes des tetraedres et remplir le
! tableau IFABOR. 
!-----------------------------------------------------------------------
! Pour reperer les faces communes aux elements :
! Parmis les elements qui partagent un noeud, il y a au moins
! deux qui partagent une face. (s'il ne s'agit pas d'un noeud de
! bord). 
! Le principe de l'algorithme : en reperant les tetraedres
! partagent le noeud IPOIN, on reconstruit les triangles de leur
! faces. 
! Si on trouve deux triangles qui partagent les memes noeuds, ca
! veut dire que les tetraedres qui les definissent sont voisins.
! Si on ne trouve pas de voisin, cela veut dire que le triangles
! est une face de bord.
! On part du principe qu'un triangles ne peut etre definit par
! plus de deux tetraedres. si c'etait le cas, il y a un probleme
! de maillage, et ca, on ne s'en occupe pas ...
!
! Avantages : on economise pas mal de memoire, en memorisant que les
! triangles autour d'un noeud. 
! Desavantages : on risque de faire des calculs en trop 
! (pour arriver a l'etape ou on definit la table de connectivite des
! triangles)
! on peut peut-etre sauter cette etape en regardant si IFABOR contient
! deja qqchose ou pas ...
!
! On va definir une table de connectivite pour les triangles.
! Cette table de connectivite n'a pas pour but de recenser la
! totalite des triangles, mais uniquement ceux autour d'un noeud.
! On connait le nombre de (tetraedres) voisins maximal d'un
! noeuds. Dans le pire des cas, il s'agit d'un noeud de bord 
! On va maximiser (beaucoup) en disant que le nombre maximal de
! triangels autour d'un noeud peut etre le nombre de tetraedres
! voisins.
! On cree aussi un tableau VOIS_TRI, qui contient l'id de l'element
! tetraedre qui l'a definit en premier (et qui sera le voisin du
! tetraedre qui va trouver qu'il y a deja un autre qui le defini)
! Ce tableau a deux entrees : l'id de l'element et l'id de la face.
!
      NBTRI = NMXVOISIN * 3
!
      allocate(IKLE_TRI(NBTRI,3),STAT=ERR)
      IF(ERR.NE.0) GOTO 999
      allocate(VOIS_TRI(NBTRI,2),STAT=ERR)
      IF(ERR.NE.0) GOTO 999
!
      IFABOR(:,:) = 0
!
      ! Boucle sur tous les noeuds du maillage.
      DO IPOIN = 1, NPOIN
          ! Pour chaque noeud, on regarde les elements tetraedres
          ! voisins (plus precisement : les faces triangles qu'il
          ! constitue)
          ! On reinitialise la table de connectivite des triangles des
          ! ces tetraedres a 0 ainsi que le nombre de triangles qu'on a
          ! trouve :
          IKLE_TRI(:,:) = 0
          ! La meme chose pour le tableau qui dit quel element a deja
          ! defini le triangle :
          VOIS_TRI(:,:) = 0
          ! On recommence a compter les triangles :
          NBTRI         = 0
          NV            = NVOIS(IPOIN)
          ADR           = IADR(IPOIN)
          DO IVOIS = 1, NV
              ! L'identifiant de l'element voisin no ivois du noeud
              ! ipoin :
              IELEM = NEIGH(ADR+IVOIS)
              ! On boucle sur les quatre faces de cet element.
              DO IFACE = 1 , NFACE
                  ! Si on a deja un voisin pour cette face, on va pas
                  ! plus loin et prend la prochaine.
                  ! Si on n'a pas encore de voisin, on le cherche...
                  IF ( IFABOR(IELEM,IFACE) .EQ. 0 ) THEN 
                  ! Chaque face definit un triangle. Le triangle est
                  ! donne par trois noeuds.
                  I1 = IKLE(IELEM,SOMFAC(1,IFACE))
                  I2 = IKLE(IELEM,SOMFAC(2,IFACE))
                  I3 = IKLE(IELEM,SOMFAC(3,IFACE))
                  ! On ordonne ces trois noeuds, M1 est le noeud avec 
                  ! l'identifiant le plus petit, M3 celui avec
                  ! l'identifiant le plus grand et M2 est au milieu :
                  M1 = MAX(I1,(MAX(I2,I3)))
                  M3 = MIN(I1,(MIN(I2,I3)))
                  M2 = I1+I2+I3-M1-M3

                  ! On parcourt le tableau des triangles deja definis
                  ! pour voir s'il en a un qui commence deja par M1.
                  ! Si c'est le cas, on verifie s'il a aussi M2 et M3
                  ! comme noeud. Si oui, c'est gagne, on a trouve un
                  ! voisin. Si la recherche echoue, on cree un nouveau
                  ! triangle.

                  found = .FALSE.
                  DO ITRI = 1, NBTRI
                      IF ( IKLE_TRI(ITRI,1) .EQ. M1 ) THEN
                          IF ( IKLE_TRI(ITRI,2) .EQ. M2 .AND.
     &                         IKLE_TRI(ITRI,3) .EQ. M3 ) THEN
                               ! on a trouve ! c'est tout bon.
                               ! On recupere l'info qu'on a dans
                               ! vois_tri. (cad l'element qui a deja
                               ! defini le triangle et la face)
                               IELEM2 = VOIS_TRI(ITRI,1)
                               IFACE2 = VOIS_TRI(ITRI,2)
                               IF ( IELEM2 .EQ. IELEM ) THEN
                                  IF(LNG.EQ.1) WRITE(LU,908) IELM
                                  IF(LNG.EQ.2) WRITE(LU,909) IELM
908                               FORMAT(1X,'VOISIN: IELM=',1I6,', 
     &                            PROBLEME DE VOISIN')
909                               FORMAT(1X,'VOISIN: IELM=',1I6,',
     &                            NEIGHBOUR PROBLEM')
                                  CALL PLANTE(1)
                                  STOP
                               END IF
                               ! pour etre sur :
                               IF ( IELEM2 .EQ. 0 .OR.
     &                              IFACE2 .EQ. 0 ) THEN
                                IF(LNG.EQ.1) WRITE(LU,918) IELEM2,IFACE2
                                IF(LNG.EQ.2) WRITE(LU,919) IELEM2,IFACE2
918                            FORMAT(1X,'VOISIN31:TRIANGLE NON DEFINI,
     &                         IELEM=',1I6,'IFACE=',1I6)
919                            FORMAT(1X,'VOISIN31:UNDEFINED TRIANGLE,
     &                         IELEM=',1I6,'IFACE=',1I6)
                                CALL PLANTE(1)
                                STOP
                               END IF
                               ! l'element et son voisin : on note la 
                               ! correspondance dans IFABOR.
                               IFABOR(IELEM ,IFACE ) = IELEM2
                               IFABOR(IELEM2,IFACE2) = IELEM
                               found = .TRUE.
                          END IF
                      END IF
                  END DO
                  ! Et non, on n'a pas encore trouve ce triangle, alors
                  ! on en cree un nouveau.
                  IF ( .NOT. FOUND) THEN
                      NBTRI             = NBTRI + 1
                      IKLE_TRI(NBTRI,1) = M1
                      IKLE_TRI(NBTRI,2) = M2
                      IKLE_TRI(NBTRI,3) = M3
                      VOIS_TRI(NBTRI,1) = IELEM
                      VOIS_TRI(NBTRI,2) = IFACE
                  END IF
              END IF ! IFABOR zero 
              END DO ! Fin boucle faces des elements voisins
!
          END DO ! fin boucle elements voisins du noeud
      END DO ! fin boucle noeuds
!
      DEALLOCATE(NEIGH)
      DEALLOCATE(IKLE_TRI)
      DEALLOCATE(VOIS_TRI)                      
!
!-----------------------------------------------------------------------
! ETAPE 5 : faces entre differents domaine de calcul (decomposition de
! domaine) : on va laisser ca aux specialistes du domaine !
!-----------------------------------------------------------------------
!
!  ON POURRAIT ESSAYER AVEC UN ALGORITHME PLUS LEGER.
!  PAR EXEMPLE EN UTILISANT INDPU
!
      IF (NCSIZE.GT.1) THEN
!

        DO 61 IFACE=1,NFACE
          DO 71 IELEM=1,NELEM
!
!  CERTAINES FACES DE BORD SONT EN FAIT DES INTERFACES ENTRE
!  SOUS-DOMAINES : ON LEUR MET UNE VALEUR -2 AU LIEU DE 0
!
!      IF(IFABOR(IELEM,IFACE).EQ.-1) THEN
            IF (IFABOR(IELEM,IFACE).EQ.0) THEN
!
         I1 = IKLE( IELEM , SOMFAC(1,IFACE) )
         I2 = IKLE( IELEM , SOMFAC(2,IFACE) )
         I3 = IKLE( IELEM , SOMFAC(3,IFACE) )
!
         IR1=0
         IR2=0
         IR3=0
!
         DO J=1,NPTIR
           IF(I1.EQ.NACHB(1,J)) IR1=1
           IF(I2.EQ.NACHB(1,J)) IR2=1
           IF(I3.EQ.NACHB(1,J)) IR3=1
         ENDDO
! 
! Affectation eventuelle à -2 fr IFABOR --> maille d'interface
! sinon, sa valeur reste à 0 --> maille de bord 
              IF (IR1.EQ.1.AND.IR2.EQ.1.AND.IR3.EQ.1) THEN
!               print*,'ce sont des noeuds d''interface'

! Ces trois noeuds sont des noeuds d'interface; Correspondent ils
! a un triangle d'interface (virtuel) ou a un triangle de bord ?
                BORD=.FALSE.
                IR4=0
                IR5=0
                IR6=0
                DO 55 J=1,NPTFR
                  IF (I1.EQ.NBOR(J)) IR5=1
                  IF (I2.EQ.NBOR(J)) IR4=1
                  IF (I3.EQ.NBOR(J)) IR6=1
55              CONTINUE
! Ce sont aussi des noeuds de bord
                IF (IR5.EQ.1.AND.IR4.EQ.1.AND.IR6.EQ.1) THEN
!
                  DO J=1,NELEB2
! C'est un triangle de bord
                    COMPT=0
                    DO I=1,3
!                     IF (IKLETR(J,I)==I1) COMPT=COMPT+1
!                     IF (IKLETR(J,I)==I2) COMPT=COMPT+10
!                     IF (IKLETR(J,I)==I3) COMPT=COMPT+100
! Bon scénario mais pb de compil avec bief
                      IF (IKLESTR(J,I)==I1) COMPT=COMPT+1
                      IF (IKLESTR(J,I)==I2) COMPT=COMPT+10
                      IF (IKLESTR(J,I)==I3) COMPT=COMPT+100
                    ENDDO
! Ces trois noeuds appartiennent bien au meme triangle de bord
                    IF (COMPT==111) THEN
                      BORD=.TRUE.
!                     print*,'SOMMETS D''UN TRIANGLE DE BORD'
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF
!                IF (IR5.EQ.0.OR.IR4.EQ.0.OR.IR6.EQ.0) THEN
                IF (.NOT.BORD) THEN
! Ces trois noeuds appartiennent donc a une maille d'interface
!                 print*,'NOEUDS D''INTERFACE'
                  IFABOR(IELEM,IFACE)=-2
                ENDIF
              ENDIF
!
            ENDIF
!
71        CONTINUE
61      CONTINUE
!     
      ENDIF
!
!-----------------------------------------------------------------------
!  
!  DISTINCTION DANS IFABOR ENTRE LES FACES DE BORD ET LES FACES LIQUIDES
!
!  INITIALISATION DE NBOR_INV A ZERO
!
      DO IPOIN=1,NPOIN
        NBOR_INV(IPOIN) = 0
      ENDDO      
!
!  PERMET DE PASSER DU NUMERO GLOBAL AU NUMERO DE BORD
!
      DO K = 1, NPTFR
        NBOR_INV(NBOR(K)) = K
      ENDDO           
!
!  BOUCLE SUR TOUTES LES FACES DE TOUS LES ELEMENTS :
!
      DO 90 IFACE = 1 , NFACE
      DO 100 IELEM = 1 , NELEM
!
      IF(IFABOR(IELEM,IFACE).EQ.-1) THEN
!
!      C'EST UNE VRAIE FACE DE BORD (LES FACES INTERNES EN PARALLELISME
!                                    SONT MARQUEES AVEC DES -2).
!      NUMEROS GLOBAUX DES POINTS DE LA FACE :
!
       I1 = IKLE( IELEM , SOMFAC(1,IFACE) )
       I2 = IKLE( IELEM , SOMFAC(2,IFACE) )
       I3 = IKLE( IELEM , SOMFAC(3,IFACE) )
!
!      UNE FACE LIQUIDE EST RECONNUE AVEC LA CONDITION LIMITE SUR H
!       
       IF(NPTFR.GT.0) THEN
       IF(LIHBOR(NBOR_INV(I1)).NE.KLOG.AND.LIHBOR(NBOR_INV(I2)).NE.KLOG
     *     .AND.LIHBOR(NBOR_INV(I3)).NE.KLOG  ) THEN
!        FACE LIQUIDE : IFABOR=0  FACE SOLIDE : IFABOR=-1
         IFABOR(IELEM,IFACE)=0       
       ENDIF
       ENDIF
!
      ENDIF
!
100    CONTINUE
90    CONTINUE    
!   
!-----------------------------------------------------------------------
!  
      RETURN
!
!-----------------------------------------------------------------------
!  
999   IF(LNG.EQ.1) WRITE(LU,1000) ERR
      IF(LNG.EQ.2) WRITE(LU,2000) ERR
1000  FORMAT(1X,'VOISIN31 : ERREUR A L''ALLOCATION DE MEMOIRE :',/,1X,
     *            'CODE D''ERREUR : ',1I6)
2000  FORMAT(1X,'VOISIN31: ERROR DURING ALLOCATION OF MEMORY: ',/,1X,
     *            'ERROR CODE: ',1I6)
      CALL PLANTE(1)
      STOP
!
!-----------------------------------------------------------------------
!  
      END
