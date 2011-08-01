C                       ******************
                        SUBROUTINE VC11TT0
C                       ******************
C
     *( XMUL,SF,SG,F,G,X,Y,Z,IKLE1,IKLE2,IKLE3,IKLE4,NELEM,NPOIN,
     *  W,ICOORD)
C
C***********************************************************************
C BIEF VERSION 5.5         06/04   J-M HERVOUET (LNH) 01 30 87 80 18
C                                        
C***********************************************************************
C
C FONCTION : CALCUL DU VECTEUR SUIVANT EN ELEMENTS FINIS :
C
C            (EXEMPLE DE LA COMPOSANTE X QUI CORRESPOND A ICOORD = 1)
C
C
C                       /            DF
C    VEC(I)  =  XMUL   /  ( G  P  *( --  )) D(OMEGA)
C                     /OMEGA    I    DX
C
C
C    P   EST UNE BASE LINEAIRE
C     I
C
C    F EST UN VECTEUR DE DISCRETISATION P1
C    W est un vecteur de discretisation P0.
C    G est aussi un vecteur de discretisation P0.
C    De ce fait, un assemblage du vecteur resultat n'est pas necessaire.
C    Le vecteur W contient directement le resultat.
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      XMUL      | -->|  COEFFICIENT MULTIPLICATEUR.
C |      SF,SG,SH  | -->|  STRUCTURES DES FONCTIONS F,G ET H
C |      SU,SV,SW  | -->|  STRUCTURES DES FONCTIONS U,V ET W
C |      F,G,H     | -->|  FONCTIONS INTERVENANT DANS LA FORMULE.
C |      U,V,W     | -->|  COMPOSANTES D'UN VECTEUR
C |                |    |  INTERVENANT DANS LA FORMULE.
C |      X, Y, Z   | -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |      SURFAC    | -->|  SURFACE DES ELEMENTS.
C |      IKLE1,... | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE.
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE.
C |      NPOIN     | -->|  NOMBRE de points dans le maillage 
C |      W         |<-- |  VECTEUR RESULTAT
C |      ICOORD    | -->|  COORDONNEE SUIVANT LAQUELLE ON DERIVE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C  PROGRAMMES APPELES : ASSVEC
C
C-----------------------------------------------------------------------
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C-----------------------------------------------------------------------
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM          ! number of elements
      INTEGER, INTENT(IN) :: NPOIN          ! number of points
      INTEGER, INTENT(IN) :: ICOORD         ! Direction of grad :
                                            ! 1 - X
                                            ! 2 - Y
                                            ! 3 - Z
C
      INTEGER, INTENT(IN) :: IKLE1(NELEM)  ! Node number 1 of elements
      INTEGER, INTENT(IN) :: IKLE2(NELEM)  ! Node number 2 of elements
      INTEGER, INTENT(IN) :: IKLE3(NELEM)  ! Node number 3 of elements
      INTEGER, INTENT(IN) :: IKLE4(NELEM)  ! Node number 4 of elements
C
      DOUBLE PRECISION, DIMENSION(NPOIN), TARGET, INTENT(IN) :: X,Y,Z 
      DOUBLE PRECISION, INTENT(IN)         :: XMUL        ! constant factor
      DOUBLE PRECISION, INTENT(OUT)        :: W(NELEM)   ! Resultat
C
C     STRUCTURES DE F,G,H,U,V,W ET DONNEES REELLES
C
      TYPE(BIEF_OBJ)  , INTENT(IN)   :: SF,SG ! 
      DOUBLE PRECISION, INTENT(IN)   :: G(*)  ! Vector to multiply by.
      DOUBLE PRECISION, INTENT(IN)   :: F(*)  ! Vector we calculate the
                                              ! grad of.
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C      
C VARIABLES LOCALES
C
      INTEGER          :: IELEM,IELMF,IELMG  ! Element types
      DOUBLE PRECISION :: F1,F2,F3,F4        ! The 4 values of F at the
                                             ! nodes of an element
      DOUBLE PRECISION :: X2,X3,X4,Y2,Y3,Y4  ! Delta_X, Delta_Y
      DOUBLE PRECISION :: X1, Y1             ! Coord of the first node
      INTEGER          :: I1,I2,I3,I4        ! The numbers of the nodes
                                             ! of an element 
C
      DOUBLE PRECISION, DIMENSION(:), POINTER :: PX, PY ! Pointer to
                                                        ! the coord
      DOUBLE PRECISION          :: XSUR24
      DOUBLE PRECISION          :: F2MF1,F3MF1,F4MF1
      DOUBLE PRECISION          :: DET               ! Jacobien
C
C-----------------------------------------------------------------------
C INITIALISATIONS
C
      XSUR24 = XMUL/24.D0
C
      IELMF = SF%ELM
      IELMG = SG%ELM      
C
C-----------------------------------------------------------------------
C Test sur la composante a deriver : 
C 1 pour X, 2 pour Y 3 pour Z. le reste n'est pas permis.
C On met le pointeur sur les tableaux de coordonnees qui seront
C utilise. 
         
      SELECT CASE (ICOORD )

       CASE ( 1 )   

         PX => Y
         PY => Z
         
       CASE ( 2 )
          
         PX => Z
         PY => X
     
       CASE ( 3 )
       
         PX => X
         PY => Y

       CASE DEFAULT 

         IF (LNG.EQ.1) WRITE(LU,202) ICOORD
         IF (LNG.EQ.2) WRITE(LU,203) ICOORD
 202     FORMAT(1X,'VC11TT0 (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *        1I6,' VERIFIER ICOORD')
 203     FORMAT(1X,'VC11TT0 (BIEF) : IMPOSSIBLE COMPONENT ',
     *        1I6,' CHECK ICOORD')
         CALL PLANTE(1)

      END SELECT


      IF(IELMF.EQ.31.AND.IELMG.EQ.30) THEN
      
       
C Loop over elements.
      DO  IELEM = 1 , NELEM

C Get the id of the four nodes of the element.

         I1 = IKLE1(IELEM)
         I2 = IKLE2(IELEM)
         I3 = IKLE3(IELEM)
         I4 = IKLE4(IELEM)

! Get the four nodal values of the vector we want to derive.
 
         F1 = F(I1)
         F2 = F(I2)
         F3 = F(I3)
         F4 = F(I4)

! The differences of the nodal values of F 

         F2MF1 = F2-F1
         F3MF1 = F3-F1
         F4MF1 = F4-F1

!
!  COORDONNEES REELLES DES POINTS DE L'ELEMENT ( ORIGINE EN 1 )
!
         X1  =  PX(I1)
         X2  =  PX(I2) - X1
         X3  =  PX(I3) - X1
         X4  =  PX(I4) - X1
         Y1  =  PY(I1)
         Y2  =  PY(I2) - Y1
         Y3  =  PY(I3) - Y1
         Y4  =  PY(I4) - Y1
       
         DET =  (X3*Y4-X4*Y3)*F2MF1 + (Y2*X4-X2*Y4)*F3MF1+
     &           (X2*Y3-Y2*X3)*F4MF1 
      
! The result

         W(IELEM) = DET* G(IELEM) * XSUR24

      ENDDO 
    
 
C
C=======================================================================
C     ERREUR sur les types d'elements.      
      
      
      ELSE
C-----------------------------------------------------------------------
C
         IF (LNG.EQ.1) WRITE(LU,1100) IELMF,SF%NAME
         IF (LNG.EQ.1) WRITE(LU,1200) IELMG,SG%NAME
         IF (LNG.EQ.1) WRITE(LU,1300)
         IF (LNG.EQ.2) WRITE(LU,1101) IELMF,SF%NAME
         IF (LNG.EQ.2) WRITE(LU,1201) IELMG,SG%NAME
         IF (LNG.EQ.2) WRITE(LU,1301)
         CALL PLANTE(1)
         STOP
 1100    FORMAT(1X,'VC11TT0 (BIEF) :',/,
     *          1X,'DISCRETISATION DE F : ',1I6,
     *          1X,'NOM REEL : ',A6)
 1200    FORMAT(1X,'DISCRETISATION DE G : ',1I6,
     *          1X,'NOM REEL : ',A6)
 1300    FORMAT(1X,'CAS NON PREVU')
 1101    FORMAT(1X,'VC11TT0 (BIEF) :',/,
     *          1X,'DISCRETIZATION OF F:',1I6,
     *          1X,'REAL NAME: ',A6)
 1201    FORMAT(1X,'DISCRETIZATION OF G:',1I6,
     *          1X,'REAL NAME: ',A6)
 1301    FORMAT(1X,'CASE NOT IMPLEMENTED')
 
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE VC11TT0
