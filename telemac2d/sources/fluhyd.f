C                       *****************
                        SUBROUTINE FLUHYD
C                       *****************
C
     *(NS,NT,NSEG,NPTFR,NUBO,G,DT,X,Y,AIRS,NU,AIRE,
     * UA,ZF,VNOIN,CE,NBOR,LIMPRO,XNEBOR,YNEBOR,KDIR,KNEU,
     * KDDL,HBOR,UBOR,VBOR,FLUENT,FLUSORT,NORDRE,CMI,JMI,
     * DJX,DJY,DX,DY,DTHAUT,CFLWTD,
     * DPX,DPY,IVIS,CVIS,FLUHBTEMP,BETA,DSZ,AIRST,HC,FLUXTEMP,NTRAC)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.8                                           INRIA
C
C***********************************************************************
C
C     FONCTION  : CALCUL DES FLUX A L'INSTANT N 
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |  NS            | -->|  NOMBRE DE POINTS DU MAILLAGE                |
C |  NT            | -->|  NOMBRE D'ELEMENTS DU MAILLAGE               |
C |  NSEG          | -->|  NOMBRE D'ARETES DU MAILLAGE                 |
C |  NPTFR         | -->|  NOMBRE DE POINTS FRONTIERE                  |
C !  NUBO          ! -->!  NUMEROS GLOBAUX DES EXTREMITES DES ARETES   !
C |  G             | -->|  CONSTANTE DE GRAVITE                        |
C |  DT            |<-->|  PAS DE TEMPS                                |
C |  X,Y           | -->|  COORDONNEES DES NOEUDS DU MAILLAGE          |
C |  AIRS          | -->|  AIRES DES CELLULES                          |
C |  NU            | -->|  NUMEROS DES NOEUDS PAR TRIANGLE             |
C |  AIRE          | -->|  AIRES DES TRIANGLES                         |
C |  UA            | -->|  UA(1,IS) = H,  UA(2,IS)=U  ,UA(3,IS)=V      |
C |  ZF            | -->|  COTES DU FOND                               |
C !  VNOIN         ! -->!  NORMALE A l'INTERFACE                       |
C !                !    !   (2 PREMIERES COMPOSANTES) ET               |
C !                !    !   LONGUEUR DE CE SEGMENT (3IEME COMPOSANTE)  |
C |  CE            |<-- |  FLUX   +  TERMES DIFFUSION                  |
C |  NBOR          | -->|  NUMEROS GLOBAUX DES POINTS DE BORD          |
C |  LIMPRO        | -->|  TYPES DE CONDITIONS AUX LIMITES             |
C |  XNEBOR,YNEBOR | -->|  NORMALE AUX POINTS FRONTIERE                |
C |  KDIR          | -->|  CONVENTION POUR LES POINTS DIRICHLET        |
C |  KNEU          | -->|  CONVENTION POUR LES POINTS NEUMANN          |
C |  KDDL          | -->|  CONVENTION POUR LES POINTS LIBRES           |
C |  HBOR          | -->|  VALEURS IMPOSEES DE H                       |
C |  UBOR          | -->|  VALEURS IMPOSEES DE U                       |
C |  VBOR          | -->|  VALEURS IMPOSEES DE V                       |
C !  FLUENT,FLUSORT|<-- |  FLUX MASSE ENTREE ET SORTIE DE TN A TN+1    |
C |  NORDRE        | -->|  ORDRE DU SCHEMA                             |
C !  CMI           ! -->!  COORDONNEES DES POINTS MILIEUX D'INTERFACE  !
C !  JMI           ! -->!  NUMERO DU TRIANGLE AUQUEL APPARTIENT LE     !
C |                !    !  POINT MILIEU DE L'INTERFACE                 !
C |  DJX,DJY       | -- |  GRADIENTS PAR TRIANGLES                     |
C |  DX,DY         | -- |  GRADIENTS PAR NOEUDS                        |
C |  DTHAUT        ! -->!  UTILISE POUR CONDITION CFL                  !
C !  CFLWTD        ! -->!  NOMBRE DE CFL                               !
C |  DPX, DPY      ! -->!  GRADIENTS DES FONCTIONS DE BASE             !
C !  IVIS          ! -->!  OPTION DIFFUSION DES VITESSES               !
C !  CVIS          ! -->!  COEFFICIENT DE DIFFUSION DES VITESSES       !
C |  FLUHBTEMP     !<-- !  FLUX BORD POUR TRACEUR                      !
C |  BETA          ! -- !  COEFFICIENT EXTRAPOLATION POUR ORDRE 2      !
C |  DSZ           ! -->!  VARIATIONS DE Z POUR ORDRE 2                !
C |  AIRST         ! -->!  AIRES DES SOUS-TRIANGLES DANS CELLULES      !
C !  HC            !<-- !  H RECONSTRUIT ORDRE 2   CORRIGE             !
C |  FLUXTEMP      !<-- !  FLUX POUR TRACEUR                           !
C |  TRAC          | -->|  LOGIQUE INDIQUANT LA PRESENCE D'UN TRACEUR  |
C !________________!____!______________________________________________!
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)   
C        -- (TABLEAU DE TRAVAIL)                                        
C-----------------------------------------------------------------------
C     - SOUS PROGRAMME(S) APPELANT : RESOLU                             
C     - SOUS PROGRAMME(S) APPELES  : GRADNOD, CDL, FLUCIN
C 
C***********************************************************************
C
      USE BIEF
      USE INTERFACE_TELEMAC2D, EX_FLUHYD => FLUHYD
C
      IMPLICIT NONE
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NS,NT,NSEG,NPTFR,KDIR,KNEU,KDDL,NORDRE
      INTEGER, INTENT(IN) :: NBOR(NPTFR),LIMPRO(NPTFR,6),NU(NT,3)
      INTEGER, INTENT(IN) :: NUBO(2,NSEG),JMI(*),IVIS,NTRAC
      DOUBLE PRECISION, INTENT(IN) :: XNEBOR(2*NPTFR),YNEBOR(2*NPTFR)
      DOUBLE PRECISION, INTENT(IN) :: HBOR(NPTFR),G,CFLWTD,DTHAUT(*)
      DOUBLE PRECISION, INTENT(IN) :: UBOR(NPTFR),VBOR(NPTFR),CMI(2,*)
      DOUBLE PRECISION, INTENT(IN) :: AIRST(2,*),CVIS
      DOUBLE PRECISION, INTENT(IN) :: X(NS),Y(NS),AIRS(NS),AIRE(NT)
      DOUBLE PRECISION, INTENT(INOUT) :: BETA,DT,HC(2,*)
      DOUBLE PRECISION, INTENT(INOUT) :: CE(3,NS),FLUENT,FLUSORT
      DOUBLE PRECISION, INTENT(IN) :: UA(3,NS),ZF(NS),VNOIN(3,NSEG)
      DOUBLE PRECISION, INTENT(IN) :: DSZ(2,*),DPX(3,NT),DPY(3,NT)
      DOUBLE PRECISION, INTENT(INOUT) :: DJX(3,*),DJY(3,*)
      DOUBLE PRECISION, INTENT(INOUT) :: DX(3,*),DY(3,*)
      TYPE(BIEF_OBJ), INTENT(INOUT) :: FLUXTEMP,FLUHBTEMP 
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C  
      INTEGER IS,IVAR
C
C-----------------------------------------------------------------------
C       
C     EXPLICIT RESOLUTION
C	
      DO IS=1,NS
        DO IVAR=1,3
          CE(IVAR,IS) = 0.D0
        ENDDO
      ENDDO
C
C   CALCUL DES GRADIENTS PAR NOEUD ET TERMES DE DIFFUSION
C
      IF(NORDRE.EQ.2.OR.IVIS.EQ.1) CALL GRADNOD(NS,NT,NU,AIRE,AIRS,
     *                        UA,DPX,DPY,DJX,DJY,DX,DY,IVIS,CVIS,CE,ZF)
C
      CALL FLUCIN(NS,NSEG,NUBO,G,X,Y,CFLWTD,DT,UA,ZF,VNOIN,CE,NORDRE,
     *            CMI,JMI,DJX,DJY,DX,DY,BETA,DSZ,AIRS,
     *            AIRST,HC,FLUXTEMP,NPTFR,NBOR,XNEBOR,YNEBOR,NTRAC)
C
C        BOUNDARY CONDITIONS TREATMENT
C
      CALL CDL(NS,NPTFR,NBOR,LIMPRO,XNEBOR,YNEBOR,KDIR,KNEU,
     *         G,HBOR,UBOR,VBOR,UA,CE,FLUENT,FLUSORT,DTHAUT,DT,CFLWTD,
     *         FLUHBTEMP,NTRAC)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
