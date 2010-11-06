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
          NOLEO(ILEO) = 1
          ISLEO(ILEO) = .FALSE.
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


