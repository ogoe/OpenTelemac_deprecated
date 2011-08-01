!                    **************************
                     MODULE INTERFACE_TELEMAC2D
!                    **************************
!
!
!***********************************************************************
! TELEMAC2D 6.1
!***********************************************************************
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF_DEF
!
!-----------------------------------------------------------------------
!
!     DEFINITION OF INTERFACES
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE AKEPIN
     &(AK,EP,U,V,H,NPOIN,KFROT,CMU,C2,ESTAR,SCHMIT,KMIN,EMIN,CF)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPOIN,KFROT
      DOUBLE PRECISION, INTENT(INOUT) :: AK(NPOIN),EP(NPOIN)
      DOUBLE PRECISION, INTENT(IN) :: KMIN,EMIN,CMU,C2,ESTAR,SCHMIT
      DOUBLE PRECISION, INTENT(IN) :: U(NPOIN),V(NPOIN),H(NPOIN),CF(*)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE ASSIGNSTR(CHESTR,SETSTR,PZONE,NZONE,NPOIN)
      USE BIEF_DEF
      IMPLICIT NONE
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: CHESTR
      TYPE(BIEF_OBJ), INTENT(IN)      :: SETSTR
      INTEGER, INTENT(IN)             :: PZONE(*)
      INTEGER, INTENT(IN)             :: NZONE
      INTEGER, INTENT(IN)             :: NPOIN
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE ASTRO
     &(YEAR,MONTH,DAY,HOUR,MIN,SEC,AT,ARL,ARS,DL,DS,AL,AS)
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: YEAR,MONTH,DAY,HOUR,MIN,SEC
      DOUBLE PRECISION, INTENT(IN)    :: AT
      DOUBLE PRECISION, INTENT(INOUT) :: ARL,ARS,DL,DS,AL,AS
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE BILANT1
     &(H,UCONV,VCONV,HPROP,WORK1,WORK2,WORK3,WORK4,WORK5,DT,LT,NIT,INFO,
     & MASKTR,T,TN,TETAT,MASSOU,MSK,MASKEL,MESH,FLUSOR,FLUENT,EQUA,LTT,
     & ITRAC)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)            :: LT,NIT,LTT,ITRAC
      DOUBLE PRECISION, INTENT(IN)   :: DT,TETAT,MASSOU,FLUSOR,FLUENT
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: WORK1,WORK2,WORK3,WORK4,WORK5
      TYPE(BIEF_OBJ), INTENT(IN)     :: HPROP,UCONV,VCONV,H,T,TN,MASKEL
      TYPE(BIEF_OBJ), INTENT(IN)     :: MASKTR
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      LOGICAL, INTENT(IN)            :: MSK,INFO
      CHARACTER(LEN=20), INTENT(IN)  :: EQUA
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE BILAN(MESH,H,WORK,MASK,
     & AT,DT,LT,NIT,INFO,MASSES,MSK,MASKEL,EQUA,POROS,
     & OPTBAN,NPTFR,FLBOR,FLUX_BOUNDARIES,NUMLIQ,NFRLIQ)
        USE BIEF_DEF
        IMPLICIT NONE
        INTEGER, INTENT(IN)            :: LT,NIT,OPTBAN,NPTFR,NFRLIQ
        INTEGER, INTENT(IN)            :: NUMLIQ(*)
        CHARACTER(LEN=20), INTENT(IN)  :: EQUA
        LOGICAL, INTENT(IN)            :: INFO,MSK
        TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
        TYPE(BIEF_OBJ), INTENT(INOUT)  :: WORK,FLBOR
        TYPE(BIEF_OBJ), INTENT(IN)     :: H,MASKEL,POROS,MASK
        DOUBLE PRECISION, INTENT(IN)   :: AT,DT
        DOUBLE PRECISION, INTENT(INOUT):: MASSES,FLUX_BOUNDARIES(*)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE BORD
     &(HBOR,UBOR,VBOR,TBOR,U,V,H,
     & ZF,NBOR,TRA05,TRA06,LIHBOR,LIUBOR,LITBOR,
     & XNEBOR,YNEBOR,NPOIN,NPTFR,NPTFR2,TEMPS,
     & NDEBIT,NCOTE,NVITES,
     & NTRAC,NTRACE,NFRLIQ,NUMLIQ,KENT,KENTU,PROVEL,MASK,MESH,EQUA,
     & NOMIMP)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPOIN,NPTFR,NDEBIT,NCOTE,NVITES,NTRACE
      INTEGER, INTENT(IN) :: KENT,KENTU,NFRLIQ,NTRAC,NPTFR2
      INTEGER, INTENT(IN) :: PROVEL(*)
      INTEGER, INTENT(IN) :: LIHBOR(NPTFR),LIUBOR(NPTFR2)
      INTEGER, INTENT(IN) :: NUMLIQ(NPTFR),NBOR(NPTFR2)
      DOUBLE PRECISION, INTENT(IN) :: TEMPS
      DOUBLE PRECISION, INTENT(IN) :: ZF(NPOIN)
      DOUBLE PRECISION, INTENT(IN) :: XNEBOR(NPTFR),YNEBOR(NPTFR)
      CHARACTER(LEN=20), INTENT(IN) :: EQUA
      CHARACTER(LEN=144), INTENT(IN) :: NOMIMP
      DOUBLE PRECISION, INTENT(INOUT) :: UBOR(NPTFR2,2),VBOR(NPTFR2,2)
      DOUBLE PRECISION, INTENT(INOUT) :: HBOR(NPTFR)
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: H,U,V,TRA05,TRA06,TBOR
      TYPE(BIEF_OBJ), INTENT(IN)  :: MASK,LITBOR
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE BORD_TIDAL_BC
     &(NBOR,LIHBOR,LIUBOR,NPTFR,
     & KENT,KENTU,MESH,GEOSYST,NUMZONE,TIDALTYPE,BOUNDARY_COLOUR,MAXFRO,
     & NFO2,NBI2,NRFO,XSHIFT,YSHIFT)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)            :: NPTFR,NFO2,NBI2,NRFO
      INTEGER, INTENT(IN)            :: KENT,KENTU,MAXFRO
      INTEGER, INTENT(IN)            :: GEOSYST,NUMZONE,TIDALTYPE
      INTEGER, INTENT(IN)            :: LIHBOR(NPTFR),LIUBOR(NPTFR)
      INTEGER, INTENT(IN)            :: NBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)   :: XSHIFT,YSHIFT
      TYPE(BIEF_OBJ), INTENT(IN)     :: BOUNDARY_COLOUR
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE BORD_TIDE
     &(ZF,NBOR,LIHBOR,LIUBOR,NPOIN,NPTFR,TEMPS,NCOTE,NVITES,
     & NUMLIQ,KENT,KENTU,NOMIMP,TIDALTYPE,CTIDE,MSL,NODALCORR,NFOT,
     & BOUNDARY_COLOUR,HBTIDE,UBTIDE,VBTIDE,NUMTIDE,ICALHW,
     & MARDAT,MARTIM)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NPOIN,NPTFR,NCOTE,NVITES,NFOT
      INTEGER, INTENT(IN)             :: KENT,KENTU,NODALCORR
      INTEGER, INTENT(IN)             :: LIHBOR(NPTFR),LIUBOR(NPTFR)
      INTEGER, INTENT(IN)             :: NUMLIQ(NPTFR),NBOR(NPTFR)
      INTEGER, INTENT(IN)             :: TIDALTYPE,MARDAT(3),MARTIM(3)
      INTEGER, INTENT(INOUT)          :: ICALHW
      DOUBLE PRECISION, INTENT(IN)    :: TEMPS,CTIDE,MSL
      DOUBLE PRECISION, INTENT(IN)    :: ZF(NPOIN)
      TYPE(BIEF_OBJ), INTENT(IN)      :: BOUNDARY_COLOUR
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: NUMTIDE,UBTIDE,VBTIDE,HBTIDE
      CHARACTER(LEN=144), INTENT(IN)  :: NOMIMP
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CALDT(NS,G,H,U,V,DTHAUT,DT,CFL,ICIN,DTVARI,LISTIN)
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NS,ICIN
      DOUBLE PRECISION, INTENT(INOUT) :: DT
      DOUBLE PRECISION, INTENT(IN)    :: H(NS),U(NS),V(NS),DTHAUT(NS)
      DOUBLE PRECISION, INTENT(IN)    :: G,CFL
      LOGICAL, INTENT(IN) :: DTVARI,LISTIN
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CARAFR
     & ( U,V,H,T,UCONV,VCONV,X , Y , SHP ,
     &   SURDET , DT , IKLE , IFABOR , ELT ,
     &   NBOR , NELBOR , NULONE , IELM , NELEM , NELMAX ,
     &   NPOIN , NDP , NPTFR ,
     &   MSK , MASKEL , MASKPT , NPT , LISPFR, NTRAC ,
     &   HBTIL , UBTIL , VBTIL , TBTIL , ZBTIL , ZF, T5  )
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM,NELMAX,NPOIN,NDP,NPTFR,IELM,NPT,NTRAC
      INTEGER, INTENT(IN) :: LISPFR(NPTFR)
      INTEGER, INTENT(IN) :: IKLE(NELMAX,NDP),IFABOR(NELMAX,*)
      INTEGER, INTENT(IN) :: NBOR(NPTFR),NELBOR(NPTFR),NULONE(NPTFR)
      INTEGER, INTENT(INOUT)          :: ELT(NPOIN)
      LOGICAL, INTENT(IN)             :: MSK
      DOUBLE PRECISION, INTENT(INOUT) :: HBTIL(NPTFR),UBTIL(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: VBTIL(NPTFR),T5(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: ZBTIL(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: U(NPOIN),V(NPOIN),H(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: UCONV(NPOIN),VCONV(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: SHP(NDP,NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: ZF(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: MASKEL(NELMAX),MASKPT(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: SURDET(NELEM)
      DOUBLE PRECISION, INTENT(IN)    :: DT
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: TBTIL
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: T
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CDLPROJ(NS,NPTFR,NBOR,LIMPRO,XNEBOR,YNEBOR,KNEU,UA)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NS,NPTFR,KNEU
      INTEGER, INTENT(IN) :: NBOR(NPTFR),LIMPRO(NPTFR,6)
      DOUBLE PRECISION, INTENT(IN) :: XNEBOR(2*NPTFR),YNEBOR(2*NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: UA(3,NS)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CDL
     &(NS,NPTFR,NBOR,LIMPRO,XNEBOR,YNEBOR,KDIR,KNEU,G,HBOR,
     & UBOR,VBOR,UA,CE,FLUENT,FLUSORT,FLBOR,
     & DTHAUT,DT,CFL,FLUHBTEMP,NTRAC)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NS,NPTFR,KDIR,KNEU,NTRAC
      INTEGER, INTENT(IN)             :: NBOR(NPTFR),LIMPRO(NPTFR,6)
      DOUBLE PRECISION, INTENT(IN)    :: XNEBOR(2*NPTFR),YNEBOR(2*NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: HBOR(NPTFR),UA(3,NS),DTHAUT(*)
      DOUBLE PRECISION, INTENT(IN)    :: UBOR(NPTFR),VBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: G,CFL
      DOUBLE PRECISION, INTENT(INOUT) :: DT
      DOUBLE PRECISION, INTENT(INOUT) :: CE(3,NS),FLUENT,FLUSORT
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: FLUHBTEMP,FLBOR
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CDL_TCH
     &(NS,NPTFR,NBOR,LIMPRO,XNEBOR,YNEBOR,KDIR,KNEU,KDDL,G,
     & HBOR,UBOR,VBOR,W,CE,FLUENT,FLUSORT,
     & FLBOR,DTHAUT,DT,CFL,EPS,ZF,WINF)
      USE BIEF
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)             :: NS,NPTFR,KDIR,KNEU,KDDL
      INTEGER, INTENT(IN)             :: NBOR(NPTFR),LIMPRO(NPTFR,6)
      DOUBLE PRECISION, INTENT(IN)    :: XNEBOR(2*NPTFR),YNEBOR(2*NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: HBOR(NPTFR),W(3,NS),DTHAUT(*)
      DOUBLE PRECISION, INTENT(IN)    :: UBOR(NPTFR),VBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: G,CFL,EPS,ZF(NS)
      DOUBLE PRECISION, INTENT(IN)    :: DT,WINF(3,NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: CE(3,NS),FLUENT,FLUSORT
      TYPE(BIEF_OBJ) , INTENT(INOUT)  :: FLBOR
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CDLZZ
     &(NS,NPTFR,NBOR,LIMPRO,XNEBOR,YNEBOR,KDIR,KNEU,KDDL,G,
     & HBOR,UBOR,VBOR,W,CE,FLUENT,FLUSORT,
     & FLBOR,DTHAUT,DT,CFL,EPS,ZF,WINF)
      USE BIEF
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)             :: NS,NPTFR,KDIR,KNEU,KDDL
      INTEGER, INTENT(IN)             :: NBOR(NPTFR),LIMPRO(NPTFR,6)
      DOUBLE PRECISION, INTENT(IN)    :: XNEBOR(2*NPTFR),YNEBOR(2*NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: HBOR(NPTFR),W(3,NS),DTHAUT(*)
      DOUBLE PRECISION, INTENT(IN)    :: UBOR(NPTFR),VBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: G,CFL,EPS,ZF(NS)
      DOUBLE PRECISION, INTENT(IN)    :: DT,WINF(3,NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: CE(3,NS),FLUENT,FLUSORT
      TYPE(BIEF_OBJ) , INTENT(INOUT)  :: FLBOR
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CHPCON(UCONV,VCONV,U,V,UN,VN,TETAU)
      USE BIEF_DEF
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: TETAU
      TYPE(BIEF_OBJ), INTENT(IN)    :: U,UN,V,VN
      TYPE(BIEF_OBJ), INTENT(INOUT) :: UCONV,VCONV
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CLHUVT
     &(NWEIRS,NPSING,NPSMAX,NUMDIG,ZDIG,
     & X,Y,ZF,IOPTAN,UNORM,CHESTR,
     & NKFROT,KARMAN,T,NTRAC,H,UBOR,VBOR,TBOR,NBOR,
     & LIHBOR,LIUBOR,LIVBOR,LITBOR)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: NWEIRS,NPSMAX,IOPTAN,NTRAC
      INTEGER, INTENT(IN)    :: NPSING(NWEIRS),NUMDIG(2,NWEIRS,NPSMAX)
      INTEGER, INTENT(IN)    :: NBOR(*),NKFROT(*)
      INTEGER, INTENT(INOUT) :: LIUBOR(*),LIHBOR(*),LIVBOR(*)
      DOUBLE PRECISION, INTENT(IN)    :: ZDIG(NWEIRS,NPSMAX)
      DOUBLE PRECISION, INTENT(IN)    :: X(*),Y(*),ZF(*),CHESTR(*),H(*)
      DOUBLE PRECISION, INTENT(IN)    :: UNORM(*)
      DOUBLE PRECISION, INTENT(INOUT) :: UBOR(*),VBOR(*)
      DOUBLE PRECISION, INTENT(IN)    :: KARMAN
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: LITBOR,TBOR
      TYPE(BIEF_OBJ), INTENT(IN)      :: T
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CLSING
     &(NWEIRS,NPSING,NPSMAX,NUMDIG,X,Y,ZF,CHESTR,NKFROT,KARMAN,
     & ZDIG,PHIDIG,NBOR,H,T,NTRAC,IOPTAN,UNORM,
     & UBOR,VBOR,TBOR,LIHBOR,LIUBOR,LIVBOR,LITBOR,GRAV)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NWEIRS,NPSMAX,IOPTAN,NTRAC
      INTEGER, INTENT(IN) :: NKFROT(*),NBOR(*)
      INTEGER, INTENT(IN) :: NPSING(NWEIRS),NUMDIG(2,NWEIRS,NPSMAX)
      INTEGER, INTENT(INOUT) :: LIUBOR(*),LIVBOR(*),LIHBOR(*)
      DOUBLE PRECISION, INTENT(IN)    :: PHIDIG(NWEIRS,NPSMAX)
      DOUBLE PRECISION, INTENT(IN)    :: ZDIG(NWEIRS,NPSMAX),H(*)
      DOUBLE PRECISION, INTENT(IN)    :: X(*),Y(*),ZF(*),CHESTR(*)
      DOUBLE PRECISION, INTENT(IN)    :: KARMAN,GRAV
      DOUBLE PRECISION, INTENT(INOUT) :: UBOR(*),VBOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: UNORM(*)
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: TBOR,LITBOR
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: T
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CLTRAC
     &(NWEIRS,NPSING,NPSMAX,NUMDIG,ZF,ZDIG,H,T,NBOR,LITBOR,TBOR,NTRAC)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NWEIRS,NPSMAX,NTRAC
      INTEGER, INTENT(IN) :: NPSING(NWEIRS),NUMDIG(2,NWEIRS,NPSMAX)
      INTEGER, INTENT(IN) :: NBOR(*)
      DOUBLE PRECISION, INTENT(IN)  :: ZDIG(NWEIRS,NPSMAX)
      DOUBLE PRECISION, INTENT(IN)  :: ZF(*),H(*)
      TYPE(BIEF_OBJ), INTENT(INOUT) :: TBOR,LITBOR
      TYPE(BIEF_OBJ), INTENT(IN)    :: T
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE COEFMAT(PERIAF,DT,M,AM,NPERIAF)
      IMPLICIT NONE
      INTEGER,          INTENT(IN   ) :: NPERIAF,M
      DOUBLE PRECISION, INTENT(IN   ) :: DT,PERIAF(NPERIAF)
      DOUBLE PRECISION, INTENT(INOUT) :: AM(2*NPERIAF,2*NPERIAF)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE COEFRO(CF,H,U,V,KARMAN,KFROT,CHESTR,GRAV,MESH,T1)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)            :: KFROT
      DOUBLE PRECISION, INTENT(IN)   :: GRAV,KARMAN
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: CF,T1
      TYPE(BIEF_OBJ), INTENT(IN)     :: CHESTR,H,U,V
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CONDIN_ADJ(ALIRE,NRES,TROUVE)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: ALIRE(*),NRES
      INTEGER, INTENT(INOUT) :: TROUVE(*)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CORNOR
     &(XNEBOR,YNEBOR,XSGBOR,YSGBOR,KP1BOR,NPTFR,KLOG,
     & LIHBOR,T1,T2,MESH)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NPTFR,KLOG
      INTEGER, INTENT(IN)             :: LIHBOR(NPTFR)  ,KP1BOR(NPTFR,2)
      DOUBLE PRECISION, INTENT(IN)    :: XSGBOR(NPTFR,4),YSGBOR(NPTFR,4)
      DOUBLE PRECISION, INTENT(INOUT) :: XNEBOR(NPTFR,2),YNEBOR(NPTFR,2)
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: T1,T2
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CORPOR(POROS)
      USE BIEF_DEF
      IMPLICIT NONE
      TYPE(BIEF_OBJ), INTENT(INOUT) :: POROS
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CORRECTION_DEPTH_2D(GLOSEG,DIMGLO,YAFLODEL,YASMH,
     &                                 YAFLULIM)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: DIMGLO
      INTEGER, INTENT(IN)    :: GLOSEG(DIMGLO,2)
      LOGICAL, INTENT(IN)    :: YAFLODEL,YASMH
      LOGICAL, INTENT(INOUT) :: YAFLULIM
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE COSAKE
     &(KARMAN,CMU,C1,C2,SIGMAK,SIGMAE,ESTAR,SCHMIT,KMIN,KMAX,EMIN,EMAX)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT) :: KMIN,KMAX,EMIN,EMAX
      DOUBLE PRECISION, INTENT(OUT) :: KARMAN,CMU,C1,C2
      DOUBLE PRECISION, INTENT(OUT) :: SIGMAK,SIGMAE,ESTAR,SCHMIT
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE COST_FUNCTION(JCOUT,OPTION,MODE)
      IMPLICIT NONE
      DOUBLE PRECISION , INTENT(INOUT) :: JCOUT
      INTEGER , INTENT(IN)             :: OPTION
      CHARACTER(LEN=3) , INTENT(IN)    :: MODE
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CUBEEQUATION(ACOF, BCOF, CCOF, DCOF, REALS, X)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: ACOF, BCOF, CCOF, DCOF
      INTEGER,          INTENT(OUT) :: REALS
      DOUBLE PRECISION, INTENT(OUT) :: X(3)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE DEBIMP
     &(Q,UBOR,VBOR,U,V,H,NUMLIQ,IFRLIQ,WORK1,WORK2,NPTFR,MASK,MESH,
     & KP1BOR,EQUA)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NPTFR,IFRLIQ
      INTEGER, INTENT(IN)             :: NUMLIQ(NPTFR),KP1BOR(NPTFR,2)
      CHARACTER(LEN=20), INTENT(IN)   :: EQUA
      DOUBLE PRECISION, INTENT(INOUT) :: UBOR(NPTFR),VBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: MASK(*),Q
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
      TYPE(BIEF_OBJ), INTENT(IN)      :: H,U,V
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: WORK1,WORK2
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        DOUBLE PRECISION FUNCTION DEBSCE( TIME , I , DISCE )
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: TIME,DISCE(*)
      INTEGER         , INTENT(IN) :: I
        END FUNCTION
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE DIFSOU
     &(TEXP,TIMP,YASMI,TSCEXP,HPROP,TN,TETAT,NREJTR,ISCE,DSCE,TSCE,
     & MAXSCE,MAXTRA,AT,DT,MASSOU,NTRAC,FAC)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: ISCE(*),NREJTR,NTRAC
      INTEGER, INTENT(IN)             :: MAXSCE,MAXTRA
      LOGICAL, INTENT(INOUT)          :: YASMI(*)
      DOUBLE PRECISION, INTENT(IN)    :: AT,DT,TETAT,DSCE(*),FAC(*)
      DOUBLE PRECISION, INTENT(IN)    :: TSCE(MAXSCE,MAXTRA)
      DOUBLE PRECISION, INTENT(INOUT) :: MASSOU(*)
      TYPE(BIEF_OBJ), INTENT(IN)      :: TN,HPROP
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: TSCEXP,TEXP,TIMP
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE DISPER( VISC , U , V , H , CF , ELDER , PROPNU )
      USE BIEF_DEF
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: ELDER(2),PROPNU
      DOUBLE PRECISION, INTENT(IN)  :: H(*),CF(*),U(*),V(*)
      TYPE(BIEF_OBJ), INTENT(INOUT) :: VISC
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE DRAGCOEFF(V,D,VK,CW)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: V, D, VK
      DOUBLE PRECISION, INTENT(OUT) :: CW
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE DRAGFO(FUDRAG,FVDRAG)
      USE BIEF_DEF
      IMPLICIT NONE
      TYPE(BIEF_OBJ), INTENT(INOUT) :: FUDRAG,FVDRAG
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE ENTETE(IETAPE,AT,LT)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: AT
      INTEGER, INTENT(IN)          :: LT,IETAPE
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        DOUBLE PRECISION FUNCTION EXLIM(ILIM,BETA,GRI,GRIJ)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ILIM
      DOUBLE PRECISION, INTENT(IN) :: GRI,GRIJ,BETA
        END FUNCTION
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FILTER_H
     &(VEC,T1,MESH,MSK,MASKEL,N,FLODEL,YAFLODEL,DT,W1,UNSV2D)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)           :: N
      DOUBLE PRECISION, INTENT(IN)  :: DT
      LOGICAL, INTENT(IN)           :: MSK,YAFLODEL
      TYPE(BIEF_MESH), INTENT(INOUT):: MESH
      TYPE(BIEF_OBJ), INTENT(INOUT) :: VEC,T1,FLODEL,W1
      TYPE(BIEF_OBJ), INTENT(IN)    :: MASKEL,UNSV2D
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FLOT
     &(XFLOT,YFLOT,NFLOT,NITFLO,FLOPRD,X,Y,NPOIN,DEBFLO,FINFLO,NIT)
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NPOIN,NIT,NFLOT,NITFLO,FLOPRD
      INTEGER, INTENT(INOUT)          :: DEBFLO(NFLOT),FINFLO(NFLOT)
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: XFLOT(NITFLO,NFLOT)
      DOUBLE PRECISION, INTENT(INOUT) :: YFLOT(NITFLO,NFLOT)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FLUCINT
     &(NS,NSEG,DIMT,NUBO,G,X,Y,UA,TN,ZF,VNOCL,CE,
     & NORDRE,CMI,JMI,DJX,DJY,DX,DY,DJXT,DJYT,DXT,DYT,EPSWL)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NS,NSEG,DIMT,NORDRE
      INTEGER, INTENT(IN) :: NUBO(2,NSEG)
      DOUBLE PRECISION, INTENT(IN) :: G,X(NS),Y(NS),VNOCL(3,NSEG),ZF(NS)
      DOUBLE PRECISION, INTENT(INOUT) :: CE(NS,3)
      DOUBLE PRECISION, INTENT(IN)    :: UA(3,NS),TN(DIMT),CMI(2,*)
      DOUBLE PRECISION, INTENT(IN)    :: JMI(*),EPSWL
      DOUBLE PRECISION, INTENT(IN)    :: DJX(3,*),DJY(3,*)
      DOUBLE PRECISION, INTENT(IN)    :: DX(3,*),DY(3,*)
      DOUBLE PRECISION, INTENT(IN)    :: DJXT(*),DJYT(*),DXT(*),DYT(*)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FLUCIN
     &(NS,NSEG,NUBO,G,X,Y,CFL,DT,UA,ZF,VNOCL,CE,NORDRE,CMI,JMI,
     & DJX,DJY,DX,DY,BETA,DSZ0,AIRS,AIRST,HC,FLUXTEMP,NPTFR,NBOR,
     & XNEBOR,YNEBOR,NTRAC)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NS,NSEG,NPTFR,NORDRE,NTRAC
      INTEGER, INTENT(IN) :: NBOR(*),NUBO(2,NSEG),JMI(*)
      DOUBLE PRECISION, INTENT(IN)    :: XNEBOR(*),YNEBOR(*),X(NS),Y(NS)
      DOUBLE PRECISION, INTENT(IN)    :: ZF(NS),VNOCL(3,NSEG),AIRS(*)
      DOUBLE PRECISION, INTENT(IN)    :: G,CFL,UA(3,NS),AIRST(2,*)
      DOUBLE PRECISION, INTENT(IN)    :: DSZ0(2,*),CMI(2,*)
      DOUBLE PRECISION, INTENT(INOUT) :: BETA,DT,CE(NS,3),HC(2,*)
      DOUBLE PRECISION, INTENT(IN)    :: DJX(3,*),DJY(3,*)
      DOUBLE PRECISION, INTENT(IN)    :: DX(3,*),DY(3,*)
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: FLUXTEMP
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FLUHYD
     &(NS,NT,NSEG,NPTFR,NUBO,G,DT,X,Y,AIRS,NU,AIRE,
     & UA,ZF,VNOIN,CE,NBOR,LIMPRO,XNEBOR,YNEBOR,KDIR,KNEU,
     & KDDL,HBOR,UBOR,VBOR,FLUENT,FLUSORT,NORDRE,CMI,JMI,
     & DJX,DJY,DX,DY,DTHAUT,CFLWTD,FLBOR,
     & DPX,DPY,IVIS,CVIS,FLUHBTEMP,BETA,DSZ,AIRST,HC,FLUXTEMP,NTRAC)
      USE BIEF_DEF
      IMPLICIT NONE
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
      TYPE(BIEF_OBJ), INTENT(INOUT) :: FLUXTEMP,FLUHBTEMP,FLBOR
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FLUROE
     &(W,FLUSCE,NUBO,VNOIN,WINF,FLUX,FLUSORT,FLUENT,
     & NELEM,NSEG,NPTFR,NPOIN,X,Y,AIRS,ZF,EPS,DDMIN,G,
     & XNEBOR,YNEBOR,LIMPRO,NBOR,KDIR,KNEU,KDDL,FLBOR)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPOIN,NELEM,NSEG,NPTFR,KDIR,KNEU,KDDL
      INTEGER, INTENT(IN) :: NUBO(2,*),LIMPRO(NPTFR,6),NBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN) :: X(NPOIN),Y(NPOIN),W(3,NPOIN)
      DOUBLE PRECISION, INTENT(IN) :: AIRS(NPOIN),ZF(NPOIN),VNOIN(3,*)
      DOUBLE PRECISION, INTENT(IN) :: XNEBOR(2*NPTFR),YNEBOR(2*NPTFR)
      DOUBLE PRECISION, INTENT(IN) :: WINF(3,NPTFR),G,EPS
      DOUBLE PRECISION, INTENT(INOUT) :: FLUX(NPOIN,3),DDMIN
      DOUBLE PRECISION, INTENT(INOUT) :: FLUSCE(3,NPOIN),FLUSORT,FLUENT
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: FLBOR
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FLUSEW
     &(AMINF,UBOR,VBOR,NPOIN,EPS,G,W,
     & XNEBOR,YNEBOR,XSGBOR,YSGBOR,
     & NPTFR,LIMPRO,NBOR,KDIR,KNEU,KDDL)
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NPOIN,NPTFR,KDIR,KNEU,KDDL
      INTEGER, INTENT(IN)             :: LIMPRO(NPTFR,6),NBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: XNEBOR(NPTFR),YNEBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: XSGBOR(NPTFR,4),YSGBOR(NPTFR,4)
      DOUBLE PRECISION, INTENT(IN)    :: W(3,NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: UBOR(NPTFR),VBOR(NPTFR),EPS,G
      DOUBLE PRECISION, INTENT(INOUT) :: AMINF(3,NPTFR) 
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FLUSRC
     &(IEL1,IEL2,ISEGIN,VNOIN,W,FLUSCE,X,Y,AIRS,NPOIN,NSEG,ZF,EPS,G)
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NPOIN,NSEG,ISEGIN,IEL1,IEL2
      DOUBLE PRECISION, INTENT(IN)    :: G,EPS,VNOIN(3,NSEG),ZF(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: AIRS(NPOIN),X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: W(3,NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: FLUSCE(3,NPOIN)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FLUTRAC
     &(NSEG,NPTFR,DT,FLUXT,FLUHBOR,FLUXTEMP,FLUHBTEMP,DTT)
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NSEG, NPTFR
      DOUBLE PRECISION, INTENT(IN)    :: DT,FLUXTEMP(*),FLUHBTEMP(*)
      DOUBLE PRECISION, INTENT(INOUT) :: FLUHBOR(*),FLUXT(*),DTT
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FLUXE
     &(HJ,UJ,VJ,HI,UI,VI,XN,YN,RNORM,G,FLULOC)
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(INOUT) :: FLULOC(3)
         DOUBLE PRECISION, INTENT(IN) :: G,HI,HJ,UI,UJ,VI,VJ,RNORM
         DOUBLE PRECISION, INTENT(IN) :: XN,YN
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FLUXZZ
     &(NS,NSEG,NUBO,G,X,Y,W,ZF,VNOCL,CE,AIRS)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NS,NSEG
      INTEGER, INTENT(IN) :: NUBO(2,NSEG)
      DOUBLE PRECISION, INTENT(IN)    :: X(NS),Y(NS)
      DOUBLE PRECISION, INTENT(IN)    :: ZF(NS),VNOCL(3,NSEG),AIRS(*)
      DOUBLE PRECISION, INTENT(IN)    :: G,W(3,NS)
      DOUBLE PRECISION, INTENT(INOUT) :: CE(3,NS)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FLU_TCH
     &(NS,NSEG,NUBO,G,X,Y,W,ZF,VNOCL,CE,AIRS)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NS,NSEG
      INTEGER, INTENT(IN) :: NUBO(2,NSEG)
      DOUBLE PRECISION, INTENT(IN)    :: X(NS),Y(NS)
      DOUBLE PRECISION, INTENT(IN)    :: ZF(NS),VNOCL(3,NSEG),AIRS(*)
      DOUBLE PRECISION, INTENT(IN)    :: G,W(3,NS)
      DOUBLE PRECISION, INTENT(INOUT) :: CE(3,NS)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FLU_TCHAMEN
     &(H1,H2,ETA1,ETA2,U1,U2,V1,V2,
     & XNN,YNN,FLXI,FLXJ,G,EPS)
      USE BIEF_DEF
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)    :: G,H1,H2,ETA1,ETA2,U1,U2
      DOUBLE PRECISION, INTENT(IN)    :: V1,V2,XNN,YNN,EPS
      DOUBLE PRECISION, INTENT(INOUT) :: FLXI(3),FLXJ(3)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FLU_ZOKAGOA
     *(H1,H2,ETA1,ETA2,U1,U2,V1,V2,
     * XNN,YNN,FLXI,FLXJ,G)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)    :: G,H1,H2,ETA1,ETA2,U1,U2
      DOUBLE PRECISION, INTENT(IN)    :: V1,V2,XNN,YNN
      DOUBLE PRECISION, INTENT(INOUT) :: FLXI(3),FLXJ(3)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FRICTION(NS,G,DT,UA,H,QU,QV,CF)
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NS
      DOUBLE PRECISION, INTENT(IN)    :: G,DT
      DOUBLE PRECISION, INTENT(IN)    :: CF(NS)
      DOUBLE PRECISION, INTENT(IN)    :: H(NS),QU(NS),QV(NS)
      DOUBLE PRECISION, INTENT(INOUT) :: UA(3,NS)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FRICTION_BUBBLE
     & (IKLE, NPOIN, NELEM, NELMAX, LINDNER, NKFROT, CHESTR, NDEFMA,
     &  LINDDP, LINDSP)
      USE BIEF_DEF
      IMPLICIT NONE
      TYPE(BIEF_OBJ), INTENT(IN)    :: IKLE
      INTEGER,        INTENT(IN)    :: NPOIN, NELEM, NELMAX
      LOGICAL,        INTENT(IN)    :: LINDNER
      TYPE(BIEF_OBJ), INTENT(INOUT) :: NKFROT, CHESTR, NDEFMA
      TYPE(BIEF_OBJ), INTENT(INOUT) :: LINDDP, LINDSP
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FRICTION_CALC
     &(N_START, N_END, KFROT, NDEF, VK, GRAV,
     & KARMAN, CHESTR, DW_MESH, HC, VRES, CF)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER,          INTENT(IN)    :: N_START, N_END, KFROT
      DOUBLE PRECISION, INTENT(IN)    :: NDEF, VK, GRAV, KARMAN
      TYPE(BIEF_OBJ),   INTENT(IN)    :: CHESTR, DW_MESH, HC, VRES
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: CF
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FRICTION_CHOICE(FRICTION_PASS, KARMAN)
      IMPLICIT NONE
      INTEGER,          INTENT(IN) :: FRICTION_PASS
      DOUBLE PRECISION, INTENT(IN) :: KARMAN
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FRICTION_LINDNER(VA,HA,CF,VK,G,DP,SP,CP)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: VA,HA,CF,VK,G,DP,SP
      DOUBLE PRECISION, INTENT(OUT) :: CP
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FRICTION_READ
     & (NCOF, NZONMX, ITURB, LISRUG, LINDNER, NOMCOF, NZONES, FRTAB)
      USE FRICTION_DEF
      IMPLICIT NONE
      INTEGER,            INTENT(IN)    :: NCOF, NZONMX
      INTEGER,            INTENT(IN)    :: ITURB, LISRUG
      LOGICAL,            INTENT(IN)    :: LINDNER
      CHARACTER(LEN=144), INTENT(IN)    :: NOMCOF
      INTEGER,            INTENT(OUT)   :: NZONES
      TYPE(FRICTION_OBJ), INTENT(INOUT) :: FRTAB
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FRICTION_SCAN(NCOF,NOMCOF,TYP,LINE)
      IMPLICIT NONE
      INTEGER,       INTENT(IN)    :: NCOF
      CHARACTER(LEN=144), INTENT(IN)    :: NOMCOF
      INTEGER,       INTENT(INOUT) :: LINE
      INTEGER,       INTENT(OUT)   :: TYP
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FRICTION_UNIF
     & (MESH,H,U,V,CHESTR,S,KFROT,KFROTL,ITURB,LISRUG,LINDNER,
     &  SB,NDEF,DP,SP,VK,KARMAN,GRAV,T1,T2,CHBORD,CF,CFBOR)
      USE BIEF_DEF
      IMPLICIT NONE
      TYPE(BIEF_MESH),  INTENT(IN)      :: MESH
      TYPE(BIEF_OBJ),   INTENT(IN)      :: H,U,V,CHESTR,CHBORD,S
      INTEGER,          INTENT(IN)      :: KFROT,KFROTL,ITURB,LISRUG
      LOGICAL,          INTENT(IN)      :: LINDNER
      DOUBLE PRECISION, INTENT(IN)      :: NDEF,DP,SP
      DOUBLE PRECISION, INTENT(IN)      :: VK,KARMAN,GRAV
      DOUBLE PRECISION, INTENT(INOUT)   :: SB
      TYPE(BIEF_OBJ),   INTENT(INOUT)   :: T1,T2
      TYPE(BIEF_OBJ),   INTENT(INOUT)   :: CF,CFBOR
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FRICTION_ZONES
     & (MESH, H, U, V, S, CHESTR, CHBORD, NKFROT, NDEFMA, LINDDP,
     &  LINDSP, KFRO_B, NDEF_B, ITURB, LISRUG, LINDNER, VK,
     &  KARMAN, GRAV, T1, T2, CF, CFBOR)
      USE BIEF_DEF
      IMPLICIT NONE
      TYPE(BIEF_MESH),    INTENT(IN)    :: MESH
      TYPE(BIEF_OBJ),     INTENT(IN)    :: H, U, V, S
      TYPE(BIEF_OBJ),     INTENT(IN)    :: CHESTR
      TYPE(BIEF_OBJ),     INTENT(IN)    :: CHBORD
      TYPE(BIEF_OBJ),     INTENT(IN)    :: NKFROT
      TYPE(BIEF_OBJ),     INTENT(IN)    :: NDEFMA, LINDDP, LINDSP
      TYPE(BIEF_OBJ),     INTENT(IN)    :: KFRO_B, NDEF_B
      INTEGER,            INTENT(IN)    :: ITURB, LISRUG
      LOGICAL,            INTENT(IN)    :: LINDNER
      DOUBLE PRECISION,   INTENT(IN)    :: VK, KARMAN, GRAV
      TYPE(BIEF_OBJ),     INTENT(INOUT) :: CF, CFBOR
      TYPE(BIEF_OBJ),     INTENT(INOUT) :: T1, T2
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FRICTI
     &(FU_IMP,FV_IMP,FUDRAG,FVDRAG,UN,VN,HN,CF,MESH,T1,T2,VERTIC,
     & UNSV2D,MSK,MASKEL,HFROT)
      USE BIEF_DEF
      IMPLICIT NONE
      LOGICAL, INTENT(IN)                 :: VERTIC,MSK
      INTEGER, INTENT(IN)                 :: HFROT
      TYPE(BIEF_OBJ),  INTENT(IN)         :: UN,VN,CF,UNSV2D,MASKEL
      TYPE(BIEF_OBJ),  INTENT(IN), TARGET :: HN
      TYPE(BIEF_OBJ),  INTENT(INOUT)      :: FU_IMP,FV_IMP,FUDRAG,FVDRAG
      TYPE(BIEF_OBJ),  INTENT(INOUT), TARGET :: T1,T2
      TYPE(BIEF_MESH), INTENT(INOUT)      :: MESH
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE GESTIO
     &(U,V,C,T,AK,EP,UTILD,VTILD,CTILD,TTILD,AKTILD,EPTILD,
     & TRAC,PROPA,CONVV,ITURB,IETAPE)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)           :: ITURB,IETAPE
      LOGICAL, INTENT(IN)           :: TRAC,CONVV(4),PROPA
      TYPE(BIEF_OBJ), INTENT(IN)    :: T,AK,EP
      TYPE(BIEF_OBJ), INTENT(INOUT) :: U,V,C
      TYPE(BIEF_OBJ), INTENT(INOUT) :: UTILD,VTILD,CTILD,TTILD
      TYPE(BIEF_OBJ), INTENT(INOUT) :: AKTILD,EPTILD
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE GRADNODT(NS,NT,NU,AIRT,AIRS,H,T,DPX,DPY,DJX,DJY,
     &                      DX,DY,DIFT,CVIST,CE,DTT)
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NS,NT
      INTEGER, INTENT(IN)             :: NU(NT,3)
      LOGICAL, INTENT(IN)             :: DIFT
      DOUBLE PRECISION, INTENT(IN)    :: DPX(3,NT),DPY(3,NT)
      DOUBLE PRECISION, INTENT(IN)    :: AIRT(NT),AIRS(NS),H(NS),T(NS)
      DOUBLE PRECISION, INTENT(INOUT) :: DJX(NT),DJY(NT),DX(NS),DY(NS)
      DOUBLE PRECISION, INTENT(INOUT) :: CE(NS)
      DOUBLE PRECISION, INTENT(IN)    :: DTT,CVIST
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE GRADNOD
     &(NS,NT,NU,AIRT,AIRS,UA,DPX,DPY,DJX,DJY,DX,DY,IVIS,CVIS,CE,ZF)
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NS,NT,IVIS
      INTEGER, INTENT(IN)             :: NU(NT,3)
      DOUBLE PRECISION, INTENT(IN)    :: AIRT(NT),AIRS(NS),CVIS
      DOUBLE PRECISION, INTENT(INOUT) :: DJX(3,NT),DJY(3,NT)
      DOUBLE PRECISION, INTENT(INOUT) :: DX(3,NS),DY(3,NS)
      DOUBLE PRECISION, INTENT(INOUT) :: CE(NS,3)
      DOUBLE PRECISION, INTENT(IN)    :: UA(3,NS),ZF(NS)
      DOUBLE PRECISION, INTENT(IN)    :: DPX(3,NT),DPY(3,NT)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE GRADZ
     &(NS,NT,NSEG,NU,NUBO,X,Y,AIRT,AIRS,CMI,JV,
     & ZF,DPX,DPY,DSZ,BETA,AIRST,DXIZ,DYIZ,DSP,DSM,CORR)
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NS,NT,NSEG
      INTEGER, INTENT(IN)             :: NU(NT,3),NUBO(2,NSEG),JV(*)
      DOUBLE PRECISION, INTENT(INOUT) :: DSZ(2,*)
      DOUBLE PRECISION, INTENT(IN)    :: X(NS),Y(NS),AIRT(NT),AIRS(NS)
      DOUBLE PRECISION, INTENT(INOUT) :: DXIZ(NS),DYIZ(NS)
      DOUBLE PRECISION, INTENT(INOUT) :: DSP(NS),DSM(NS),CORR(NS),BETA
      DOUBLE PRECISION, INTENT(IN)    :: DPX(3,NT),DPY(3,NT)
      DOUBLE PRECISION, INTENT(IN)    :: CMI(2,*),AIRST(2,*),ZF(NS)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE HPROPA(HPROP,HN,H,PROLIN,HAULIN,TETA,NSOUSI)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)           :: NSOUSI
      LOGICAL, INTENT(IN)           :: PROLIN
      DOUBLE PRECISION, INTENT(IN)  :: TETA,HAULIN
      TYPE(BIEF_OBJ), INTENT(IN)    :: HN,H
      TYPE(BIEF_OBJ), INTENT(INOUT) :: HPROP
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE INCIDE
     &(COTOND,H,C0,PATMOS,ATMOS,ZF,MESH,LT,AT,GRAV,ROEAU,PRIVE)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)            :: LT
      LOGICAL, INTENT(IN)            :: ATMOS
      DOUBLE PRECISION, INTENT(IN)   :: AT,GRAV,ROEAU
      TYPE(BIEF_OBJ), INTENT(IN)     :: PATMOS,H,C0,ZF,PRIVE
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: COTOND
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE INIT_KIN  
     * (W,FLUSCE,NUBO,VNOIN,AT,DT,LT,NIT,
     *  NELEM,NSEG,NPTFR,FLUX,AIRS,AIRE,
     *  X,Y,IKLE,ZF,CF,NPOIN,HN,H,U,V,QU,QV,G,LISTIN,
     *  XNEBOR,YNEBOR,XSGBOR,YSGBOR,
     *  LIMPRO,NBOR,KDIR,KNEU,KDDL, 
     *  HBOR,UBOR,VBOR,FLUSORT,FLUENT,CFLWTD,DTVARI,NELMAX,KFROT,  
     *  NREJET,ISCE,TSCE2,MAXSCE,MAXTRA,YASMH,SMH,MASSES,
     *  NTRAC,DIMT,T,HTN,TN,DLIMT,LIMTRA,
     *  TBOR,MASSOU,FLUTENT,FLUTSOR,DTHAUT,DPX,DPY,DJX,DJY,CMI,JMI,
     *  SMTR,DXT,DYT,DJXT,DJYT,
     *  DIFVIT,ITURB,PROPNU,DIFT,DIFNU,
     *  DX,DY,OPTVF,FLUSORTN,FLUENTN,
     *  DSZ,AIRST,HSTOK,HCSTOK,FLUXT,FLUHBOR,FLBOR,
     *  LOGFR,LTT,DTN,FLUXTEMP,FLUHBTEMP,
     *  HC,TMAX,DTT,T1,T2,T3,T4,T5,ICIN,EPS,NORDRE,IVIS)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM,NPOIN,NSEG,NPTFR,LT,NIT,NREJET,DIMT
      INTEGER, INTENT(IN) :: MAXSCE,MAXTRA,ICIN,NORDRE,IVIS
      INTEGER, INTENT(IN) :: DLIMT,OPTVF,JMI(*)
      INTEGER, INTENT(IN) :: KDIR,KNEU,KDDL,ITURB,NELMAX,KFROT,NTRAC
      INTEGER, INTENT(IN) :: NUBO(2,*),LIMPRO(NPTFR,6),NBOR(NPTFR)
      INTEGER, INTENT(IN) :: IKLE(NELMAX,3),ISCE(NREJET),LIMTRA(DLIMT)
      INTEGER, INTENT(INOUT) :: LTT,LOGFR(*)
      LOGICAL, INTENT(IN) :: LISTIN,DTVARI,YASMH,DIFVIT,DIFT
      DOUBLE PRECISION, INTENT(INOUT) :: T1(*),T2(*),T3(*),T4(*),T5(*)
      DOUBLE PRECISION, INTENT(IN)    :: XNEBOR(2*NPTFR),YNEBOR(2*NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: XSGBOR(NPTFR,4),YSGBOR(NPTFR,4)	
      DOUBLE PRECISION, INTENT(INOUT) :: DT
      DOUBLE PRECISION, INTENT(IN)    :: EPS,AT,VNOIN(3,*)
      DOUBLE PRECISION, INTENT(IN)    :: TSCE2(MAXSCE,MAXTRA)   
      DOUBLE PRECISION, INTENT(INOUT) :: W(3,NPOIN),FLUSORTN,FLUENTN   
      DOUBLE PRECISION, INTENT(IN)    :: AIRE(NPOIN),DTHAUT(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: HBOR(NPTFR),UBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: VBOR(NPTFR),HN(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: SMH(NPOIN),ZF(NPOIN),CF(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: U(NPOIN),V(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: H(NPOIN),QU(NPOIN),QV(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: DPX(3,NELMAX),DPY(3,NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN),Y(NPOIN),AIRS(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: FLUSCE(3,NPOIN),FLUX(NPOIN,3)
      DOUBLE PRECISION, INTENT(INOUT) :: FLUSORT,FLUENT,MASSES
      DOUBLE PRECISION, INTENT(INOUT) :: FLUTENT(*),FLUTSOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: MASSOU(*)
      DOUBLE PRECISION, INTENT(IN)    :: G,CFLWTD,AIRST(2,NSEG)
      DOUBLE PRECISION, INTENT(INOUT) :: HSTOK(*),HCSTOK(2,*),DTT
      DOUBLE PRECISION, INTENT(INOUT) :: CMI(2,NSEG)
      DOUBLE PRECISION, INTENT(IN)    :: PROPNU,DIFNU,TMAX
      DOUBLE PRECISION, INTENT(INOUT) :: DJX(3,NELMAX),DJY(3,NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: DX(3,NPOIN),DY(3,NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: DJXT(NELMAX),DJYT(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: DXT(NPOIN),DYT(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: DSZ(2,NSEG)
      DOUBLE PRECISION, INTENT(INOUT) :: HC(2,NSEG),DTN 
      TYPE(BIEF_OBJ) , INTENT(IN)     :: TBOR,TN
      TYPE(BIEF_OBJ) , INTENT(INOUT)  :: T,HTN,SMTR,FLUHBOR,FLUHBTEMP
      TYPE(BIEF_OBJ) , INTENT(INOUT)  :: FLUXTEMP,FLUXT,FLBOR         
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE INITSTR(CHESTR,SETSTR,PZONE,NZONE,NPOIN,T1)
      USE BIEF_DEF
      IMPLICIT NONE
      TYPE(BIEF_OBJ), INTENT(IN)      :: CHESTR
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: SETSTR,T1
      INTEGER, INTENT(IN)             :: PZONE(*)
      INTEGER, INTENT(IN)             :: NZONE
      INTEGER, INTENT(IN)             :: NPOIN
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE INTEMP
     &(W,FLUX,FLUX_OLD,AIRS,DT,NPOIN,ZF,CF,EPS,KFROT,SMH,
     & HN,QU,QV,LT,GAMMA) 
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NPOIN,KFROT,LT
      DOUBLE PRECISION, INTENT(INOUT) :: W(3,NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: FLUX(NPOIN,3),DT,EPS
      DOUBLE PRECISION, INTENT(IN)    :: AIRS(NPOIN),ZF(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: FLUX_OLD(NPOIN,3)
      DOUBLE PRECISION, INTENT(IN)    :: CF(NPOIN),SMH(NPOIN),GAMMA 
      DOUBLE PRECISION, INTENT(IN)    :: HN(NPOIN),QU(NPOIN),QV(NPOIN)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE INTERPOL(RO,R02,R03,JCOUT1,JCOUT2,JCOUT3)
      IMPLICIT NONE
      DOUBLE PRECISION , INTENT(IN)    :: R02,R03,JCOUT1,JCOUT2,JCOUT3
      DOUBLE PRECISION , INTENT(INOUT) :: RO
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE ISITOK
     &(H,NPH,U,NPU,V,NPV,NTRAC,T,NPT,X,Y,BORNES,ARRET)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)          :: NPH,NPU,NPV,NPT,NTRAC
      LOGICAL, INTENT(INOUT)       :: ARRET
      DOUBLE PRECISION, INTENT(IN) :: H(NPH),U(NPU),V(NPV)
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*),BORNES(8)
      TYPE(BIEF_OBJ)  , INTENT(IN) :: T
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE KEPSCL
     &(KBOR,EBOR,AUBOR,CF,CFBOR,DISBOR,
     & UN,VN,HN,LIMKEP,LIUBOR,LIMPRO,NBOR,NPTFR,
     & KARMAN,CMU,C2,ESTAR,SCHMIT,LISRUG,PROPNU,KMIN,EMIN,
     & KNEU,KDIR,KENT,KENTU,KADH,KLOG,UETUTA)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPTFR,LISRUG
      INTEGER, INTENT(IN) :: KNEU,KDIR,KENT,KADH,KLOG,KENTU
      INTEGER, INTENT(IN) :: LIMPRO(NPTFR,6),NBOR(NPTFR)
      INTEGER, INTENT(IN) :: LIMKEP(NPTFR,2),LIUBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: KMIN,EMIN
      DOUBLE PRECISION, INTENT(IN)    :: CF(*),CFBOR(*),UETUTA(*)
      DOUBLE PRECISION, INTENT(IN)    :: UN(*),VN(*),HN(*),DISBOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: KBOR(*),EBOR(*),AUBOR(*)
      DOUBLE PRECISION, INTENT(IN) :: KARMAN,CMU,C2,ESTAR,SCHMIT,PROPNU
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE KEPSIL
     &(AK,EP,AKTILD,EPTILD,AKN,EPN,VISC,CF,U,V,HN,UCONV,VCONV,
     & KBOR,EBOR,LIMKEP,IELMK,IELME,
     & SMK,SME,TM1,MAK,MAE,CM2,TE1,TE2,NPTFR,DT,
     & MESH,T1,T2,T3,TB,CMU,C1,C2,SIGMAK,SIGMAE,ESTAR,SCHMIT,
     & KMIN,KMAX,EMIN,EMAX,
     & INFOKE,KDIR,MSK,MASKEL,MASKPT,S,SLVK,SLVEP,ICONV,OPTSUP)
      USE BIEF_DEF
      IMPLICIT NONE
      TYPE(SLVCFG), INTENT(INOUT)  :: SLVK,SLVEP
      INTEGER, INTENT(IN)          :: ICONV,NPTFR,KDIR,LIMKEP(NPTFR,2)
      INTEGER, INTENT(IN)          :: OPTSUP,IELMK,IELME
      LOGICAL, INTENT(IN)          :: INFOKE,MSK
      DOUBLE PRECISION, INTENT(IN) :: KMIN,KMAX,EMIN,EMAX,SCHMIT
      DOUBLE PRECISION, INTENT(IN) :: CMU,C1,C2,SIGMAK,SIGMAE,ESTAR
      DOUBLE PRECISION, INTENT(IN) :: DT
      TYPE(BIEF_OBJ), INTENT(INOUT) :: TM1,MAK,MAE,CM2
      TYPE(BIEF_OBJ), INTENT(IN)    :: UCONV,VCONV,AKN,EPN,AKTILD,EPTILD
      TYPE(BIEF_OBJ), INTENT(IN)    :: HN,VISC,U,V,MASKEL,S,MASKPT,CF
      TYPE(BIEF_OBJ), INTENT(IN)    :: KBOR,EBOR
      TYPE(BIEF_OBJ), INTENT(INOUT) :: T1,T2,T3,AK,EP,SMK,SME,TE1,TE2
      TYPE(BIEF_MESH) :: MESH
      TYPE(BIEF_OBJ) :: TB
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE KEPSIN
     &(LIMKEP,LIUBOR,NPTFR,
     & KENT,KENTU,KSORT,KADH,KLOG,KINC,KNEU,KDIR)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: NPTFR,KENT,KSORT,KADH,KLOG
      INTEGER, INTENT(IN)    :: KINC,KNEU,KDIR,KENTU
      INTEGER, INTENT(INOUT) :: LIMKEP(NPTFR,2)
      INTEGER, INTENT(IN)    :: LIUBOR(NPTFR)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE LAGRAN(NLAG,DEBLAG,FINLAG)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: NLAG
      INTEGER, INTENT(INOUT) :: DEBLAG(NLAG) , FINLAG(NLAG)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE LECDON_TELEMAC2D(MOTCAR,FILE_DESC,PATH,NCAR)
      IMPLICIT NONE
      INTEGER, INTENT(IN)               :: NCAR
      CHARACTER(LEN=250), INTENT(IN)    :: PATH
      CHARACTER(LEN=144), INTENT(INOUT) :: FILE_DESC(4,300)
      CHARACTER(LEN=144), INTENT(INOUT) :: MOTCAR(300)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE LECSIP
     & (RELAXS,NSIPH,ENTSIP,SORSIP,SECSCE,
     &  ALTSCE,CSSCE,CESCE,DELSCE,ANGSCE,LSCE,MAXSCE,IFIC)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MAXSCE,IFIC
      INTEGER, INTENT(INOUT) :: ENTSIP(*),SORSIP(*),NSIPH
      DOUBLE PRECISION, INTENT(INOUT) :: RELAXS
      DOUBLE PRECISION, INTENT(INOUT) :: SECSCE(*),ALTSCE(*),CSSCE(*)
      DOUBLE PRECISION, INTENT(INOUT) :: DELSCE(*),ANGSCE(*)
      DOUBLE PRECISION, INTENT(INOUT) :: CESCE(*),LSCE(*)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE LECSNG
     &(NWEIRS,NWRMAX,NPSING,NUMDIG,ZDIG,PHIDIG,IOPTAN,NPSMAX,NPOIN,IFIC)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: NWRMAX,NPOIN,IFIC,NWEIRS
      INTEGER, INTENT(INOUT) :: NPSMAX,IOPTAN
      INTEGER, INTENT(INOUT) :: NPSING(NWEIRS),NUMDIG(2,NWEIRS,*)
      DOUBLE PRECISION, INTENT(INOUT) :: ZDIG(NWEIRS,*     )
      DOUBLE PRECISION, INTENT(INOUT) :: PHIDIG(NWEIRS,*     )
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE LOIDEN(YAM,YS,PHI,DEB,G)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(INOUT) :: DEB
      DOUBLE PRECISION, INTENT(IN)    :: G,YAM,PHI,YS
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE LOINOY(YAM,YAV,YS,PHI,DEB,G)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(INOUT) :: DEB
      DOUBLE PRECISION, INTENT(IN)    :: G,YAM,YAV,PHI,YS
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE MAJTRAC
     & (NS,NT,DIMT,DLIMT,NSEG,NPTFR,NUBO,
     & X,Y,AIRS,NU,AIRE,HT,HTN,TN,ZF,NBOR,
     & TBOR,FLUTENT,FLUTSOR, SMTR,NORDRE,CMI,JMI,
     & DJXT,DJYT,DXT,DYT,
     & DPX,DPY,DIFT,CVIST,BETA,DSZ,AIRST,HSTOK,
     & HCSTOK,FLUXT,FLUHBOR,MASSOU,DTT)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: DIFT
      INTEGER, INTENT(IN) :: NSEG,NPTFR,NORDRE,DIMT,DLIMT,NS,NT
      INTEGER, INTENT(IN) :: NUBO(2,NSEG),NU(NT,3)
      INTEGER, INTENT(IN) :: NBOR(NPTFR),JMI(*)
      DOUBLE PRECISION, INTENT(INOUT) :: HT(DIMT),FLUTENT,FLUTSOR
      DOUBLE PRECISION, INTENT(INOUT) :: MASSOU
      DOUBLE PRECISION, INTENT(IN)    :: TBOR(DLIMT),DSZ(2,*)
      DOUBLE PRECISION, INTENT(IN)    :: X(NS),Y(NS),AIRS(NS),AIRE(NT)
      DOUBLE PRECISION, INTENT(IN)    :: HTN(DIMT),TN(DIMT),ZF(*)
      DOUBLE PRECISION, INTENT(IN)    :: SMTR(DIMT),DPX(3,NT),DPY(3,NT)
      DOUBLE PRECISION, INTENT(IN)    :: CMI(2,*),AIRST(2,*),CVIST
      DOUBLE PRECISION, INTENT(INOUT) :: DJXT(*),DJYT(*),DXT(*),DYT(*)
      DOUBLE PRECISION, INTENT(INOUT) :: BETA
      DOUBLE PRECISION, INTENT(IN)    :: HSTOK(*)
      DOUBLE PRECISION, INTENT(IN)    :: HCSTOK(2,*),FLUXT(*)
      DOUBLE PRECISION, INTENT(IN)    :: FLUHBOR(*),DTT
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE MAJ
     &(NS,NSEG,NPTFR,G,DT,AIRS,
     & H,QU,QV,UA,CE,NBOR,LIMPRO,XNEBOR,YNEBOR,KNEU,SMH,KFROT,CF)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NS,NSEG,NPTFR,KNEU,KFROT
      INTEGER, INTENT(IN) :: NBOR(NPTFR),LIMPRO(NPTFR,6)
      DOUBLE PRECISION, INTENT(IN) :: XNEBOR(2*NPTFR),YNEBOR(2*NPTFR)
      DOUBLE PRECISION, INTENT(IN) :: H(NS),QU(NS),QV(NS),SMH(NS)
      DOUBLE PRECISION, INTENT(IN) :: CE(3,NS),CF(NS),AIRS(NS)
      DOUBLE PRECISION, INTENT(IN) :: G,DT
      DOUBLE PRECISION, INTENT(INOUT) :: UA(3,NS)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE MAJZZ    
     &(W,FLUX,FLUX_OLD,AIRS,DT,NPOIN,ZF,CF,EPS,KFROT,SMH,
     & HN,QU,QV,LT,GAMMA,NPTFR,NBOR,LIMPRO,XNEBOR,YNEBOR,KNEU)
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NPOIN,KFROT,LT,NPTFR,KNEU
      INTEGER, INTENT(IN)             :: NBOR(NPTFR),LIMPRO(NPTFR,6)
      DOUBLE PRECISION, INTENT(IN)    :: XNEBOR(2*NPTFR),YNEBOR(2*NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: W(3,NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: FLUX(NPOIN,3),DT,EPS
      DOUBLE PRECISION, INTENT(IN)    :: FLUX_OLD(NPOIN,3),GAMMA
      DOUBLE PRECISION, INTENT(IN)    :: AIRS(NPOIN),ZF(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: CF(NPOIN),SMH(NPOIN) 
      DOUBLE PRECISION, INTENT(IN)    :: HN(NPOIN),QU(NPOIN),QV(NPOIN)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE MARAST
     &(MARDAT,MARTIM,PHI0,NPOIN,AT,FU1,FV1,X,SINLAT,COSLAT,GRAV)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MARDAT(3),MARTIM(3),NPOIN
      DOUBLE PRECISION, INTENT(INOUT) :: FU1(NPOIN),FV1(NPOIN)
      DOUBLE PRECISION, INTENT(IN) :: COSLAT(NPOIN),SINLAT(NPOIN)
      DOUBLE PRECISION, INTENT(IN) :: X(NPOIN),GRAV,AT,PHI0
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE MASBAS2D
     &(VOLU2D,V2DPAR,UNSV2D,IELM,MESH,MSK,MASKEL,T1,S)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)            :: IELM
      LOGICAL, INTENT(IN)            :: MSK
      TYPE(BIEF_OBJ) , INTENT(INOUT) :: VOLU2D,V2DPAR,UNSV2D,T1
      TYPE(BIEF_OBJ) , INTENT(IN)    :: MASKEL,S
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE MASKOB(MASKEL,X,Y,IKLE,NELEM,NELMAX,NPOIN,AT,LT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM,NPOIN,NELMAX,LT
      INTEGER, INTENT(IN) :: IKLE(NELMAX,*)
      DOUBLE PRECISION, INTENT(INOUT) :: MASKEL(NELEM)
      DOUBLE PRECISION, INTENT(IN) :: X(NPOIN),Y(NPOIN),AT
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE MATBOU
     &(MESH,M1,M2,A11,A12,A21,A22,SMU,SMV,VR,VS,H0,MSK,MASKEL,S)
      USE BIEF_DEF
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: MSK
      TYPE(BIEF_OBJ), INTENT(INOUT) :: A11,A12,A21,A22,M1,M2
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ), INTENT(INOUT) :: SMU,SMV
      TYPE(BIEF_OBJ), INTENT(IN)    :: MASKEL,H0,S,VR,VS
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE MESURES(ITER,TT)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
      INTEGER, INTENT(IN) :: ITER
      DOUBLE PRECISION, INTENT(IN) :: TT
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE METEO
     &(PATMOS,WINDX,WINDY,FUAIR,FVAIR,X,Y,AT,LT,NPOIN,VENT,ATMOS,
     & HN,TRA01,GRAV,ROEAU,NORD,PRIVE)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LT,NPOIN
      LOGICAL, INTENT(IN) :: ATMOS,VENT
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN),Y(NPOIN),HN(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: WINDX(NPOIN),WINDY(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: PATMOS(NPOIN),TRA01(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: FUAIR,FVAIR,AT,GRAV,ROEAU,NORD
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: PRIVE
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE METGRA
     &(RO,ESTIME,GRADJ,GRADJN,JCOUT1,DESC,NPARAM,OPTID,RSTART,R02,R03)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER , INTENT(IN)             :: NPARAM,OPTID
      CHARACTER(LEN=72)                :: ESTIME
      DOUBLE PRECISION , INTENT(IN)    :: JCOUT1
      LOGICAL , INTENT(IN)             :: RSTART
      TYPE(BIEF_OBJ) , INTENT(IN)      :: GRADJ,GRADJN
      TYPE(BIEF_OBJ) , INTENT(INOUT)   :: DESC
      DOUBLE PRECISION , INTENT(INOUT) :: R02,R03,RO
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE NEWSTR
     &(SETSTR,SETSTR2,DESC,RO,RSTART,NPARAM,ESTIME,KFROT)
      USE BIEF_DEF
      IMPLICIT NONE
      DOUBLE PRECISION , INTENT(IN)    :: RO
      TYPE (BIEF_OBJ)  , INTENT(IN)    :: DESC
      TYPE (BIEF_OBJ)  , INTENT(IN)    :: SETSTR2
      TYPE (BIEF_OBJ)  , INTENT(INOUT) :: SETSTR
      LOGICAL          , INTENT(INOUT) :: RSTART
      INTEGER          , INTENT(IN)    :: NPARAM,KFROT
      CHARACTER(LEN=72), INTENT(IN)    :: ESTIME
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE NODALF_PUGH
     &(FFMN2,FFM4,NODALCORR,TEMPS,DEJA,MARDAT,MARTIM)
      IMPLICIT NONE
      INTEGER, INTENT(IN)           :: NODALCORR,MARDAT(3),MARTIM(3)
      DOUBLE PRECISION, INTENT(OUT) :: FFMN2,FFM4
      DOUBLE PRECISION, INTENT(IN)  :: TEMPS
      LOGICAL, INTENT(IN)           :: DEJA
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE NODALUPV_PUGH(UPVM2,UPVN2,MARDAT,MARTIM)
      IMPLICIT NONE
      INTEGER, INTENT(IN)           :: MARDAT(3),MARTIM(3)
      DOUBLE PRECISION, INTENT(OUT) :: UPVM2,UPVN2
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE NOMVAR_TELEMAC2D(TEXTE,TEXTPR,MNEMO,NPERIAF,NTRAC,
     &   NAMETRAC)
      IMPLICIT NONE
      CHARACTER(LEN=32), INTENT(INOUT) :: TEXTE(*),TEXTPR(*)
      CHARACTER(LEN=8),  INTENT(INOUT) :: MNEMO(*)
      INTEGER, INTENT(IN)              :: NPERIAF,NTRAC
      CHARACTER(LEN=32), INTENT(IN)    :: NAMETRAC(32)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE OUTPUT_TELEMAC2D(TIME)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: TIME
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE PORO11(TETA,ZF,HN,IKLE,NELEM,NELMAX)
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NELEM,NELMAX
      INTEGER, INTENT(IN)             :: IKLE(NELMAX,*)
      DOUBLE PRECISION, INTENT(IN)    :: ZF(*),HN(*)
      DOUBLE PRECISION, INTENT(INOUT) :: TETA(*)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE POROS(TETA,ZF,HN,MESH)
      USE BIEF_DEF
      IMPLICIT NONE
      TYPE(BIEF_OBJ), INTENT(IN)    :: ZF,HN
      TYPE(BIEF_OBJ), INTENT(INOUT) :: TETA
      TYPE(BIEF_MESH), INTENT(IN)   :: MESH
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE PREBOR
     &(HBOR,UBOR,VBOR,TBOR,U,V,H,HN,T,NBOR,NPOIN,NPTFR,NTRAC,NFRLIQ,
     & FRTYPE,NUMLIQ)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NPOIN,NPTFR,NFRLIQ,NTRAC
      INTEGER, INTENT(IN)             :: NBOR(NPTFR)
      INTEGER, INTENT(IN)             :: NUMLIQ(NPTFR),FRTYPE(NFRLIQ)
      DOUBLE PRECISION, INTENT(IN)    :: U(NPOIN),V(NPOIN),H(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: HBOR(NPTFR),UBOR(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: VBOR(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: HN(NPOIN)
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: T
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: TBOR
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE PROSOU(FU,FV,SMH,    UN,VN,HN,GRAV,NORD,
     &                    FAIR,WINDX,WINDY,VENT,HWIND,CORIOL,FCOR,
     &                    SPHERI,YASMH,COSLAT,SINLAT,AT,LT,
     &                    NREJET,NREJEU,DSCE,ISCE,T1,MESH,MSK,MASKEL,
     &                    MAREE,MARDAT,MARTIM,PHI0,OPTSOU,
     &                    COUROU,NPTH,VARCL,NVARCL,VARCLA,UNSV2D,
     &                    FXWAVE,FYWAVE)
      USE BIEF_DEF
      IMPLICIT NONE
      TYPE(BIEF_OBJ), INTENT(INOUT) :: T1
      TYPE(BIEF_OBJ), INTENT(INOUT) :: FU,FV,SMH,FXWAVE,FYWAVE
      TYPE(BIEF_OBJ), INTENT(IN)    :: MASKEL,UN,VN,HN,UNSV2D
      TYPE(BIEF_OBJ), INTENT(IN)    :: WINDX,WINDY,COSLAT,SINLAT
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      INTEGER, INTENT(IN)           :: NVARCL,LT,NREJET,NREJEU,OPTSOU
      INTEGER, INTENT(IN)           :: NPTH
      INTEGER, INTENT(IN)           :: MARDAT(3),MARTIM(3),ISCE(NREJET)
      DOUBLE PRECISION, INTENT(IN)  :: HWIND,AT,FAIR,FCOR,DSCE(NREJET)
      DOUBLE PRECISION, INTENT(IN)  :: GRAV,NORD,PHI0
      CHARACTER(LEN=32), INTENT(IN) :: VARCLA(NVARCL)
      LOGICAL, INTENT(IN)           :: VENT,MAREE,CORIOL,SPHERI,MSK
      LOGICAL, INTENT(IN)           :: COUROU
      LOGICAL, INTENT(INOUT)        :: YASMH
      TYPE(BIEF_OBJ), INTENT(INOUT) :: VARCL
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE PROPAG
     &(U,V,H,UCONV,VCONV,CONVV,H0,C0,COTOND,PATMOS,ATMOS,
     & HPROP,UN,VN,HN,UTILD,VTILD,HTILD,DH,DU,DV,DHN,VISC,VISC_S,FU,FV,
     & SMH,MESH,ZF,AM1,AM2,AM3,BM1,BM2,CM1,CM2,TM1,A23,A32,MBOR,
     & CV1,CV2,CV3,W1,UBOR,VBOR,AUBOR,HBOR,DIRBOR,
     & TE1,TE2,TE3,TE4,TE5,T1,T2,T3,T4,T5,T6,T7,T8,
     & LIMPRO,MASK,GRAV,ROEAU,CF,DIFVIT,IORDRH,IORDRU,LT,AT,DT,
     & TETAH,TETAHC,TETAU,TETAD,
     & AGGLOC,AGGLOU,KDIR,INFOGR,KFROT,ICONVF,
     & PRIVE,ISOUSI,BILMAS,MASSES,YASMH,OPTBAN,CORCON,
     & OPTSUP,MSK,MASKEL,MASKPT,RO,ROVAR,
     & MAT,RHS,UNK,TB,S,BD,PRECCU,SOLSYS,CFLMAX,OPDVIT,OPTSOU,
     & NFRLIQ,SLVPRO,EQUA,VERTIC,ADJO,ZFLATS,TETAZCOMP,
     & UDEL,VDEL,DM1,ZCONV,COUPLING,FLBOR,BM1S,BM2S,CV1S,
     & VOLU2D,V2DPAR,UNSV2D,NUMDIG,NWEIRS,NPSING,HFROT,FLULIM,YAFLULIM)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LT,OPTSUP(4),KDIR,KFROT,ICONVF(4),NWEIRS
      INTEGER, INTENT(IN) :: IORDRH,IORDRU,ISOUSI,OPTBAN,OPTSOU,SOLSYS
      INTEGER, INTENT(IN)             :: OPDVIT,NFRLIQ,HFROT
      DOUBLE PRECISION, INTENT(IN)    :: TETAU,TETAD,TETAH,AGGLOC,AGGLOU
      DOUBLE PRECISION, INTENT(IN)    :: TETAHC,AT,DT,GRAV,ROEAU
      DOUBLE PRECISION, INTENT(IN)    :: TETAZCOMP
      DOUBLE PRECISION, INTENT(INOUT) :: CFLMAX,MASSES
!                                                           NPSMAX
      INTEGER, INTENT(IN) :: NPSING(NWEIRS),NUMDIG(2,NWEIRS,*     )
      LOGICAL, INTENT(IN) :: BILMAS,ATMOS,DIFVIT,INFOGR,CONVV(4),MSK
      LOGICAL, INTENT(IN) :: YASMH,ROVAR,PRECCU,VERTIC,ADJO,CORCON
      LOGICAL, INTENT(IN) :: YAFLULIM
      TYPE(SLVCFG), INTENT(INOUT)    :: SLVPRO
      CHARACTER(LEN=20), INTENT(IN)  :: EQUA
      CHARACTER(LEN=*) , INTENT(IN)  :: COUPLING
      TYPE(BIEF_OBJ), INTENT(IN)     :: UCONV,VCONV,SMH,UN,VN,HN
      TYPE(BIEF_OBJ), INTENT(IN)     :: VOLU2D,V2DPAR,UNSV2D
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: RO,FLBOR
      TYPE(BIEF_OBJ), INTENT(IN)     :: UTILD,VTILD,PATMOS,CF,FLULIM
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: U,V,H,CV1,CV2,CV3,PRIVE,DH,DHN
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: CV1S
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: DU,DV,FU,FV,VISC,VISC_S,HTILD
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: UBOR,VBOR,HBOR,AUBOR,COTOND
      TYPE(BIEF_OBJ), INTENT(IN)     :: MASKEL,MASKPT,ZF
      TYPE(BIEF_OBJ), INTENT(IN)     :: HPROP,H0,C0,LIMPRO
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: TE1,TE2,TE3,TE4,TE5,ZFLATS
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: T1,T2,T3,T4,T5,T6,T7,T8
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: W1,UDEL,VDEL,DM1,ZCONV
      TYPE(BIEF_OBJ), INTENT(IN)     :: S
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: AM1,AM2,AM3,BM1,BM2,CM1,CM2,TM1
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: A23,A32,MBOR,BM1S,BM2S
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: MASK,MAT,RHS,UNK,TB,BD,DIRBOR
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
        END SUBROUTINE
      END INTERFACE
!
!=======================================================================
!
      INTERFACE
        SUBROUTINE PROPAG_ADJ
     &(UCONV,VCONV,CONVV,H0,C0,COTOND,PATMOS,ATMOS,
     & HPROP,UN,VN,HN,UTILD,VTILD,HTILD,DH,DU,DV,DHN,VISC,VISC_S,FU,FV,
     & SMH,MESH,ZF,AM1,AM2,AM3,BM1,BM2,CM1,CM2,TM1,A23,A32,MBOR,
     & CV1,CV2,CV3,W1,UBOR,VBOR,AUBOR,HBOR,DIRBOR,
     & TE1,TE2,TE3,TE4,TE5,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,
     & LIMPRO,MASK,GRAV,ROEAU,CF,DIFVIT,IORDRH,IORDRU,LT,AT,DT,
     & TETAH,TETAHC,TETAU,TETAD,
     & AGGLOC,AGGLOU,KDIR,INFOGR,KFROT,ICONVF,
     & PRIVE,ISOUSI,BILMAS,MASSES,YASMH,OPTBAN,CORCON,
     & OPTSUP,MSK,MASKEL,MASKPT,RO,ROVAR,
     & MAT,RHS,UNK,TB,S,BD,PRECCU,SOLSYS,CFLMAX,OPDVIT,OPTSOU,
     & NFRLIQ,SLVPRO,EQUA,VERTIC,
     & ADJO,UD,VD,HD,U,V,H,UU,VV,HH,UIT1,VIT1,HIT1,PP,QQ,RR,
     & TAM1,TAM2,TAM3,TBM1,TBM2,TCM1,
     & TCM2,MATADJ,UNKADJ,ALPHA1,ALPHA2,ALPHA3,ADJDIR,ESTIME,OPTCOST,
     & NIT,NVARRES,VARSOR,
     & NRES,NREF,ALIRE,TROUVE,MAXVAR,VARCL,VARCLA,
     & TEXTE,TEXREF,TEXRES,W,OUTINI,CHESTR,KARMAN,NDEF,
     & ITURB,LISRUG,LINDNER,SB,DP,SP,CHBORD,CFBOR,HFROT,UNSV2D)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LT,OPTSUP(4),KDIR,KFROT,ICONVF(4)
      INTEGER, INTENT(IN) :: IORDRH,IORDRU,ISOUSI,OPTBAN,OPTSOU,SOLSYS
      INTEGER, INTENT(IN) :: OPDVIT,NFRLIQ,LISRUG,ITURB,OPTCOST,HFROT
      INTEGER, INTENT(IN)    :: NIT,NRES,NREF,MAXVAR
      INTEGER, INTENT(INOUT) :: NVARRES,TROUVE(*),ALIRE(*)
      LOGICAL, INTENT(IN)    :: BILMAS,ATMOS,DIFVIT,INFOGR,CONVV(4),MSK
      LOGICAL, INTENT(IN)    :: YASMH,ROVAR,PRECCU,VERTIC,ADJO,CORCON
      LOGICAL, INTENT(IN)    :: OUTINI,LINDNER
      DOUBLE PRECISION, INTENT(IN)    :: TETAU,TETAD,TETAH,AGGLOC,AGGLOU
      DOUBLE PRECISION, INTENT(IN)    :: TETAHC,AT,DT,GRAV,ROEAU,CFLMAX
      DOUBLE PRECISION, INTENT(IN)    :: KARMAN,NDEF,DP,SP
      DOUBLE PRECISION, INTENT(INOUT) :: MASSES,SB
      TYPE(SLVCFG), INTENT(INOUT)     :: SLVPRO
      CHARACTER(LEN=20), INTENT(IN)   :: EQUA
      TYPE(BIEF_OBJ), INTENT(IN)      :: UCONV,VCONV,SMH,UN,VN,HN,UNSV2D
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: RO
      TYPE(BIEF_OBJ), INTENT(IN)      :: UTILD,VTILD,PATMOS,CF
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: U,V,H,CV1,CV2,CV3,PRIVE,DH,DHN
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: DU,DV,FU,FV,VISC,VISC_S,HTILD
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: UBOR,VBOR,HBOR,AUBOR
      TYPE(BIEF_OBJ), INTENT(IN)      :: MASKEL,MASKPT,ZF
      TYPE(BIEF_OBJ), INTENT(IN)      :: HPROP,H0,C0,COTOND,LIMPRO
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: T1,T2,T3,T4,T5,T6,T7,T8,T9,T10
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: T11
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: TE1,TE2,TE3,TE4,TE5
!     STRUCTURES OF MATRICES
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: TAM1,TAM2,TAM3,TBM1
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: TBM2,TCM1,TCM2
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: AM1,AM2,AM3,BM1,BM2,CM1,CM2,TM1
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: A23,A32,MBOR
!
      TYPE(BIEF_OBJ),  INTENT(INOUT)  :: MASK,MAT,RHS,UNK,TB,BD,DIRBOR
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
      TYPE(BIEF_OBJ),  INTENT(INOUT)  :: CHESTR
      TYPE (BIEF_OBJ), INTENT(INOUT)  :: HD,UD,VD,ALPHA1,ALPHA2,ALPHA3
      TYPE (BIEF_OBJ), INTENT(INOUT)  :: HH,UU,VV,UIT1,VIT1,HIT1
      TYPE (BIEF_OBJ), INTENT(INOUT)  :: PP,QQ,RR,CHBORD,CFBOR
      TYPE(BIEF_OBJ),  INTENT(INOUT)  :: W1
      TYPE(BIEF_OBJ), INTENT(IN)      :: S
      REAL,  INTENT(INOUT)            :: W(*)
      TYPE(BIEF_OBJ),  INTENT(INOUT)  :: VARSOR
      TYPE(BIEF_OBJ),  INTENT(INOUT)  :: MATADJ,UNKADJ,ADJDIR,VARCL
      CHARACTER(LEN=72), INTENT(IN)   :: ESTIME
      CHARACTER(LEN=32), INTENT(INOUT):: VARCLA(10),TEXTE(*)
      CHARACTER(LEN=32), INTENT(INOUT):: TEXREF(*),TEXRES(*)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE PROPIN_TELEMAC2D
     &(LIMPRO,LIMDIM,MASK,LIUBOR,LIVBOR,LIHBOR,KP1BOR,NBOR,NPTFR,
     & KENT,KENTU,KSORT,KADH,KLOG,KINC,KNEU,KDIR,KDDL,KOND,
     & CLH,CLU,CLV,IELMU,U,V,GRAV,H,LT,NPOIN,NELBOR,NELMAX,MSK,MASKEL,
     & NFRLIQ,THOMFR,NUMLIQ,FRTYPE,XNEBOR,YNEBOR,ENTET)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPOIN,NELMAX,NPTFR,KOND,KENTU,LT,NFRLIQ
      INTEGER, INTENT(IN) :: KENT,KSORT,KADH,KLOG,KINC,KNEU,KDIR,KDDL
      INTEGER, INTENT(IN) :: LIMDIM,IELMU
      INTEGER, INTENT(IN) :: NELBOR(NPTFR),LIVBOR(NPTFR),LIHBOR(NPTFR)
      INTEGER, INTENT(IN) :: LIUBOR(NPTFR),FRTYPE(*)
      INTEGER, INTENT(INOUT) :: LIMPRO(LIMDIM,6)
      INTEGER, INTENT(IN) :: KP1BOR(NPTFR),NBOR(NPTFR),NUMLIQ(NFRLIQ)
      INTEGER, INTENT(INOUT) :: CLH(NPTFR),CLU(NPTFR),CLV(NPTFR)
      LOGICAL, INTENT(IN) :: MSK,THOMFR,ENTET
      DOUBLE PRECISION, INTENT(IN)   :: XNEBOR(NPTFR),YNEBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)   :: U(NPOIN),V(NPOIN),H(NPOIN)
      DOUBLE PRECISION, INTENT(IN)   :: GRAV,MASKEL(NELMAX)
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: MASK
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        DOUBLE PRECISION FUNCTION Q(I)
      IMPLICIT NONE
      INTEGER         , INTENT(IN) :: I
        END FUNCTION
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE READ_FIC_CURVES
     &                           (NFIC,NFRLIQ,STA_DIS_CURVES,PTS_CURVES)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: NFIC,NFRLIQ
      INTEGER, INTENT(IN)    :: STA_DIS_CURVES(NFRLIQ)
      INTEGER, INTENT(INOUT) :: PTS_CURVES(NFRLIQ)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE READ_FIC_FRLIQ(Q, WHAT , AT , NFIC , LISTIN , STAT )
      IMPLICIT NONE
      CHARACTER*8     , INTENT(IN)       :: WHAT
      DOUBLE PRECISION, INTENT(IN)       :: AT
      DOUBLE PRECISION, INTENT(INOUT)    :: Q
      INTEGER         , INTENT(IN)       :: NFIC
      LOGICAL         , INTENT(IN)       :: LISTIN
      LOGICAL         , INTENT(OUT)      :: STAT
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE REINIT
     &(NS,NSEG,NPTFR,H,SMTR,HSTOK,HC,HCSTOK,FLUXT,FLUHBOR,DTT,NTRAC)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NS,NSEG,NPTFR,NTRAC
      DOUBLE PRECISION, INTENT(INOUT) :: DTT
      DOUBLE PRECISION, INTENT(INOUT) :: HSTOK(*),HCSTOK(2,*)
      DOUBLE PRECISION, INTENT(IN)    :: H(*),HC(2,*)
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: SMTR,FLUXT,FLUHBOR
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE RESCUE
     &(U,V,H,S,ZF,T,TRAC0,NTRAC,ITURB,NPOIN,AKEP,TROUVE)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: TROUVE(36),ITURB,NPOIN,NTRAC
      LOGICAL, INTENT(INOUT)          :: AKEP
      DOUBLE PRECISION, INTENT(INOUT) :: U(NPOIN),V(NPOIN),H(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: S(NPOIN),ZF(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: TRAC0(NTRAC)
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: T
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE RESOLU
     &(W,FLUSCE,NUBO,VNOIN,WINF,AT,DT,LT,NIT,
     & NELEM,NSEG,NPTFR,FLUX,AIRS,AIRE,
     & X,Y,IKLE,ZF,CF,NPOIN,HN,H,U,V,QU,QV,G,LISTIN,
     & XNEBOR,YNEBOR,XSGBOR,YSGBOR,LIMPRO,NBOR,KDIR,KNEU,KDDL,
     & HBOR,UBOR,VBOR,FLUSORT,FLUENT,CFLWTD,DTVARI,NELMAX,KFROT,  
     & NREJET,ISCE,TSCE2,MAXSCE,MAXTRA,YASMH,SMH,MASSES,
     & NTRAC,DIMT,T,HTN,TN,DLIMT,LIMTRA,
     & TBOR,MASSOU,FLUTENT,FLUTSOR,DTHAUT,DPX,DPY,DJX,DJY,CMI,JMI,
     & SMTR,DXT,DYT,DJXT,DJYT,DIFVIT,ITURB,PROPNU,DIFT,DIFNU,
     & DX,DY,OPTVF,FLUSORTN,FLUENTN,
     & DSZ,AIRST,HSTOK,HCSTOK,FLUXT,FLUHBOR,FLBOR,
     & LOGFR,LTT,DTN,FLUXTEMP,FLUHBTEMP,
     & HC,TMAX,DTT,T1,T2,T3,T4,T5,GAMMA,FLUX_OLD)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM,NPOIN,NSEG,NPTFR,LT,NIT,NREJET,DIMT
      INTEGER, INTENT(IN) :: MAXSCE,MAXTRA
      INTEGER, INTENT(IN) :: DLIMT,OPTVF,JMI(*)
      INTEGER, INTENT(IN) :: KDIR,KNEU,KDDL,ITURB,NELMAX,KFROT,NTRAC
      INTEGER, INTENT(IN) :: NUBO(2,*),LIMPRO(NPTFR,6),NBOR(NPTFR)
      INTEGER, INTENT(IN) :: IKLE(NELMAX,3),ISCE(NREJET),LIMTRA(DLIMT)
      INTEGER, INTENT(INOUT) :: LTT,LOGFR(*)
      LOGICAL, INTENT(IN) :: LISTIN,DTVARI,YASMH,DIFVIT,DIFT
      DOUBLE PRECISION, INTENT(INOUT) :: T1(*),T2(*),T3(*),T4(*),T5(*)
      DOUBLE PRECISION, INTENT(IN)    :: XNEBOR(2*NPTFR),YNEBOR(2*NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: XSGBOR(NPTFR,4),YSGBOR(NPTFR,4)	
      DOUBLE PRECISION, INTENT(INOUT) :: DT
      DOUBLE PRECISION, INTENT(IN)    :: AT,VNOIN(3,*),GAMMA
      DOUBLE PRECISION, INTENT(IN)    :: TSCE2(MAXSCE,MAXTRA)   
      DOUBLE PRECISION, INTENT(INOUT) :: W(3,NPOIN),FLUSORTN,FLUENTN   
      DOUBLE PRECISION, INTENT(IN)    :: AIRE(NPOIN),DTHAUT(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: HBOR(NPTFR),UBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: VBOR(NPTFR),HN(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: SMH(NPOIN),ZF(NPOIN),CF(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: U(NPOIN),V(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: H(NPOIN),QU(NPOIN),QV(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: DPX(3,NELMAX),DPY(3,NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: WINF(3,*)
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN),Y(NPOIN),AIRS(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: FLUSCE(3,NPOIN),FLUX(NPOIN,3)
      DOUBLE PRECISION, INTENT(INOUT) :: FLUSORT,FLUENT,MASSES
      DOUBLE PRECISION, INTENT(INOUT) :: FLUTENT(*),FLUTSOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: MASSOU(*)
      DOUBLE PRECISION, INTENT(IN)    :: G,CFLWTD,AIRST(2,NSEG)
      DOUBLE PRECISION, INTENT(INOUT) :: HSTOK(*),HCSTOK(2,*),DTT
      DOUBLE PRECISION, INTENT(INOUT) :: CMI(2,NSEG)
      DOUBLE PRECISION, INTENT(IN)    :: PROPNU,DIFNU,TMAX
      DOUBLE PRECISION, INTENT(INOUT) :: DJX(3,NELMAX),DJY(3,NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: DX(3,NPOIN),DY(3,NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: DJXT(NELMAX),DJYT(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: DXT(NPOIN),DYT(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: DSZ(2,NSEG),FLUX_OLD(NPOIN,3)
      DOUBLE PRECISION, INTENT(INOUT) :: HC(2,NSEG),DTN 
      TYPE(BIEF_OBJ) , INTENT(IN)     :: TBOR,TN
      TYPE(BIEF_OBJ) , INTENT(INOUT)  :: T,HTN,SMTR,FLUHBOR,FLUHBTEMP
      TYPE(BIEF_OBJ) , INTENT(INOUT)  :: FLUXTEMP,FLUXT,FLBOR            
!
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE ROTNE0
     &(MESH,M1,A11,A12,A21,A22,SMU,SMV,UN,VN,H0,MSK,MASKEL,S,DT)
      USE BIEF_DEF
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: MSK
      DOUBLE PRECISION, INTENT(IN)   :: DT
      TYPE(BIEF_OBJ), INTENT(IN)     :: MASKEL,H0,S,UN,VN
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: SMU,SMV
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: A11,A12,A21,A22,M1
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE SIPHON
     &(RELAXS,NSIPH,ENTSIP,SORSIP,GRAV,
     & H,ZF,ISCE,DSCE,SECSCE,ALTSCE,CSSCE,CESCE,DELSCE,ANGSCE,LSCE,
     & NTRAC,T,TSCE,USCE,VSCE,U,V,ENTET,MAXSCE)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NSIPH,NTRAC,MAXSCE
      INTEGER, INTENT(IN)             :: ENTSIP(*),SORSIP(*),ISCE(*)
      LOGICAL, INTENT(IN)             :: ENTET
      DOUBLE PRECISION, INTENT(IN)    :: RELAXS,GRAV
      DOUBLE PRECISION, INTENT(INOUT) :: USCE(*),VSCE(*),DSCE(*)
      DOUBLE PRECISION, INTENT(INOUT) :: TSCE(MAXSCE,NTRAC)
      DOUBLE PRECISION, INTENT(IN)    :: ANGSCE(*),LSCE(*),CESCE(*)
      DOUBLE PRECISION, INTENT(IN)    :: CSSCE(*),DELSCE(*)
      DOUBLE PRECISION, INTENT(IN)    :: SECSCE(*),ALTSCE(*)
      DOUBLE PRECISION, INTENT(IN)    :: H(*),ZF(*),U(*),V(*)
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: T
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE SMOOTHING_FLUX
     &(XMUL,SF,F,SURFAC,IKLE1,IKLE2,IKLE3,NELEM,NELMAX,W)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX),IKLE3(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: W(NELMAX,3)
      DOUBLE PRECISION, INTENT(IN)    :: SURFAC(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL
      TYPE(BIEF_OBJ), INTENT(IN)      :: SF
      DOUBLE PRECISION, INTENT(IN)    :: F(*)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        DOUBLE PRECISION FUNCTION SL( I , N )
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I,N
        END FUNCTION
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE SMAGOR
     &(VISC,CF,U,V,MESH,T1,T2,T3,T4,MSK,MASKEL,PROPNU)
      USE BIEF_DEF
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: MSK
      DOUBLE PRECISION, INTENT(IN)   :: PROPNU
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: VISC,T1,T2,T3,T4
      TYPE(BIEF_OBJ), INTENT(IN)     :: MASKEL,CF,U,V
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE SMTRAC
     &(NPOIN,DIMT,AT,DT,SMTR,SMH,NREJET,ISCE,TSCE2,MAXSCE,MAXTRA,ITRAC)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPOIN,NREJET,ISCE(*),DIMT,ITRAC
      INTEGER, INTENT(IN) :: MAXSCE,MAXTRA
      DOUBLE PRECISION, INTENT(IN)    :: AT,DT,SMH(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: TSCE2(MAXSCE,MAXTRA)
      DOUBLE PRECISION, INTENT(INOUT) :: SMTR(DIMT)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE SORFLO
     &(XFLOT,YFLOT,IKLFLO,DEBFLO,FINFLO,NFLOT,NITFLO,FLOPRD,
     & NRBI,TITCAS,BINRES,NOMRBI,NIT,MAXVAR,DATE,TIME,MESH,
     & I_ORIG,J_ORIG)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: I_ORIG,J_ORIG,FLOPRD,NITFLO
      INTEGER, INTENT(IN)             :: NFLOT,NRBI,NIT,MAXVAR
      INTEGER, INTENT(IN)             :: DATE(3),TIME(3)
      INTEGER, INTENT(INOUT)          :: IKLFLO(3*NITFLO*NFLOT)
      INTEGER, INTENT(IN)             :: DEBFLO(NFLOT),FINFLO(NFLOT)
      DOUBLE PRECISION, INTENT(INOUT) :: XFLOT(NITFLO*NFLOT)
      DOUBLE PRECISION, INTENT(INOUT) :: YFLOT(NITFLO*NFLOT)
      CHARACTER(LEN=72), INTENT(IN)   :: TITCAS,NOMRBI
      CHARACTER(LEN=3), INTENT(IN)    :: BINRES
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        DOUBLE PRECISION FUNCTION STA_DIS_CUR
     &(IFRLIQ,FLUX,PTS,QZ,NFRLIQ,ZN)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IFRLIQ,NFRLIQ,PTS
      DOUBLE PRECISION, INTENT(IN) :: ZN,FLUX,QZ(2,NFRLIQ,PTS)
        END FUNCTION
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE STEADY
     &(H1,H2,NPH,U1,U2,NPU,V1,V2,NPV,NTRAC,T1,T2,NPT,CRIPER,ARRET)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)          :: NPH,NPU,NPV,NPT,NTRAC
      LOGICAL, INTENT(INOUT)       :: ARRET
      DOUBLE PRECISION, INTENT(IN) :: H1(NPH),H2(NPH),U1(NPU),U2(NPU)
      DOUBLE PRECISION, INTENT(IN) :: V1(NPV),V2(NPV)
      DOUBLE PRECISION, INTENT(IN) :: CRIPER(3)
      TYPE(BIEF_OBJ)  , INTENT(IN) :: T1,T2
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE TCHAM_SMALL(HI,HJ,ETAI,ETAJ,UI,UJ,VI,VJ,G,FLX)
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
      DOUBLE PRECISION, INTENT(IN)    :: G,HI,HJ,ETAI,ETAJ,UI,UJ
      DOUBLE PRECISION, INTENT(IN)    :: VI,VJ
      DOUBLE PRECISION, INTENT(INOUT) :: FLX(3)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE TELEMAC2D(PASS,ATDEP,NITER,CODE,
     &                                            DTDEP,NEWTIME,DOPRINT)
        IMPLICIT NONE
        INTEGER,          INTENT(IN) :: PASS,NITER
        DOUBLE PRECISION, INTENT(IN) :: ATDEP
        CHARACTER(LEN=*), INTENT(IN) :: CODE
        DOUBLE PRECISION, INTENT(IN), OPTIONAL :: DTDEP
        LOGICAL,          INTENT(IN), OPTIONAL :: NEWTIME,DOPRINT
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE TEL4DEL
     &(NPOIN,NPOIN2,NELEM,NSEG,IKLE,ELTSEG,GLOSEG,ORISEG,MAXSEG,
     & X,Y,NPTFR,LIHBOR,NBOR,NOLAY,AAT,DDT,LLT,NNIT,HNEW,HPROP,ZNEW,
     & U,V,SALI,TEMP,VISC,TITRE,NOMGEO,NOMLIM,NSTEPA,
     & NSOU,NOMSOU,NMAB,NOMMAB,NCOU,NOMCOU,NINI,NOMINI,NVEB,NOMVEB,
     & NMAF,NOMMAF,NCOB,NOMCOB,NSAL,NOMSAL,NTEM,NOMTEM,NVEL,NOMVEL,
     & NVIS,NOMVIS,INFOGR,NELEM2,SALI_DEL,TEMP_DEL,VELO_DEL,DIFF_DEL,
     & MARDAT,MARTIM,FLOW,INIFLOW,W,YAFLULIM,FLULIM,V2DPAR,KNOLG,
     & MESH2D,MESH3D)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)            :: NPOIN,NPOIN2,NELEM,NSEG,NPTFR
      INTEGER, INTENT(IN)            :: NOLAY,LLT,NNIT,NELEM2
      INTEGER, INTENT(IN)            :: MAXSEG,MARDAT(3),MARTIM(3)
      INTEGER, INTENT(INOUT)         :: NSTEPA
      INTEGER, INTENT(IN)            :: IKLE(NELEM2,3),LIHBOR(*)
      INTEGER, INTENT(IN)            :: ELTSEG(NELEM2,3)
      INTEGER, INTENT(IN)            :: ORISEG(NELEM2,3)
      INTEGER, INTENT(IN)            :: GLOSEG(MAXSEG,2)
      INTEGER, INTENT(IN)            :: NBOR(*),KNOLG(NPOIN2)
      INTEGER, INTENT(IN)            :: NSOU,NMAB,NCOU,NINI,NVEL,NVIS
      INTEGER, INTENT(IN)            :: NVEB,NMAF,NCOB,NSAL,NTEM
      DOUBLE PRECISION  , INTENT(IN) :: X(NPOIN2),Y(NPOIN2),ZNEW(NPOIN)
      DOUBLE PRECISION  , INTENT(IN) :: HPROP(NPOIN2),HNEW(NPOIN2)
      DOUBLE PRECISION  , INTENT(IN) :: AAT,DDT,V2DPAR(NPOIN2)
      DOUBLE PRECISION  , INTENT(IN) :: U(NPOIN),V(NPOIN)
      DOUBLE PRECISION  , INTENT(IN) :: SALI(NPOIN),TEMP(NPOIN)
      DOUBLE PRECISION  , INTENT(IN) :: VISC(NPOIN),FLULIM(*)
      DOUBLE PRECISION  , INTENT(INOUT) :: FLOW(*),W(NELEM,*)
      CHARACTER(LEN=72) , INTENT(IN) :: TITRE
      CHARACTER(LEN=144), INTENT(IN) :: NOMSOU,NOMMAB,NOMCOU,NOMINI
      CHARACTER(LEN=144), INTENT(IN) :: NOMVEB,NOMMAF,NOMCOB,NOMSAL
      CHARACTER(LEN=144), INTENT(IN) :: NOMTEM,NOMGEO,NOMLIM,NOMVEL
      CHARACTER(LEN=144), INTENT(IN) :: NOMVIS
      LOGICAL           , INTENT(IN) :: INFOGR,SALI_DEL,TEMP_DEL,INIFLOW
      LOGICAL           , INTENT(IN) :: VELO_DEL,DIFF_DEL,YAFLULIM
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH2D,MESH3D
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE TESTEUR
     &(NS,NSEG,NPTFR,NUBO,DT,NBOR,NORDRE,AIRS,AIRST,HSTOK,
     & HCSTOK,FLUXT,FLUXTEMP,FLUHBOR,FLUHBTEMP,LOGFR,TEST,NTRAC)
      USE BIEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: NS,NSEG,NPTFR,NORDRE,NTRAC
      INTEGER, INTENT(IN)             :: NUBO(2,NSEG),NBOR(NPTFR)
      INTEGER, INTENT(IN)             :: LOGFR(*)
      DOUBLE PRECISION, INTENT(IN)    :: AIRS(NS),AIRST(2,*)
      DOUBLE PRECISION, INTENT(IN)    :: DT
      DOUBLE PRECISION, INTENT(INOUT) :: TEST
      DOUBLE PRECISION, INTENT(IN)    :: HSTOK(*)
      DOUBLE PRECISION, INTENT(IN)    :: HCSTOK(2,*)
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: FLUHBOR,FLUHBTEMP,FLUXT
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: FLUXTEMP
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE THOMPS
     &(HBOR,UBOR,VBOR,TBOR,U,V,H,T,ZF,X,Y,NBOR,FRTYPE,UNA,C,
     & UCONV,VCONV,T6,T7,XCONV,YCONV,FU,FV,LIHBOR,LIUBOR,LIVBOR,
     & LITBOR,LISPFR,T8,ITRAV2,W1R,W2R,W3R,
     & HBTIL,UBTIL,VBTIL,TBTIL,ZBTIL,SURDET,IKLE,CF,SMH,IFABOR,NELEM,
     & MESH,XNEBOR,YNEBOR,NPOIN,NPTFR,LT,TEMPS,DT,GRAV,
     & NTRAC,NFRLIQ,KSORT,KINC,KENT,KENTU,LV,MSK,MASKEL,
     & NELMAX,IELM,NORD,FAIR,WINDX,WINDY,
     & VENT,HWIND,CORIOL,FCOR,SPHERI,
     & MAREE,MARDAT,MARTIM,PHI0,OPTSOU,ISCE,DSCE,SHPP,
     & COUROU,NPTH,VARCL,NVARCL,VARCLA,NUMLIQ,SHP,UNSV2D,HFROT,
     & FXWAVE,FYWAVE,DX_T,DY_T,DZ_T,ELT_T,IT3,IT4,HFIELD,UFIELD,VFIELD,
     & ZS,GZSX,GZSY)
      USE BIEF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPTFR,LT,NPOIN,NELEM,NELMAX,NFRLIQ,LV
      INTEGER, INTENT(IN) :: NVARCL,NPTH,IELM,NTRAC,HFROT
      INTEGER, INTENT(IN) :: KSORT,KINC,KENT,KENTU
      INTEGER, INTENT(IN) :: MARDAT(3),MARTIM(3),OPTSOU,ISCE(*)
      INTEGER, INTENT(IN) :: NBOR(NPTFR),IKLE(*),IFABOR(*)
      INTEGER, INTENT(IN) :: LIHBOR(NPTFR),LIUBOR(NPTFR),LIVBOR(NPTFR)
      INTEGER, INTENT(IN) :: FRTYPE(NFRLIQ),NUMLIQ(NPTFR)
      INTEGER, INTENT(INOUT) :: LISPFR(NPTFR),ELT_T(NPTFR)
!     ITRAV2 : OF DIMENSION NPOIN
      INTEGER, INTENT(INOUT) :: ITRAV2(*)
      LOGICAL, INTENT(IN) :: VENT,MAREE,CORIOL,SPHERI,MSK,COUROU
      DOUBLE PRECISION, INTENT(INOUT) :: HBOR(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: UBOR(NPTFR),VBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: XNEBOR(NPTFR),YNEBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: SURDET(*),DSCE(*)
      DOUBLE PRECISION, INTENT(IN)    :: TEMPS,GRAV,DT,FAIR,FCOR,NORD
      DOUBLE PRECISION, INTENT(IN)    :: HWIND,PHI0
      DOUBLE PRECISION, INTENT(INOUT) :: W1R(NPTFR),W2R(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: W3R(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: SHPP(3,NPTFR),SHP(*)
      DOUBLE PRECISION, INTENT(INOUT) :: DX_T(NPTFR),DY_T(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: DZ_T(NPTFR)
      TYPE(BIEF_OBJ), INTENT(IN)      :: WINDX,WINDY,MASKEL
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: HBTIL,UBTIL,VBTIL,ZBTIL,T7
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: VARCL,FXWAVE,FYWAVE,XCONV,YCONV
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: FU,FV,T8,UNA,UCONV,VCONV,C,U,V
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: H,T,SMH,TBOR,TBTIL,T6,IT3,IT4
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: HFIELD,UFIELD,VFIELD,ZS
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: GZSX,GZSY
      TYPE(BIEF_OBJ), INTENT(IN)      :: ZF,CF,LITBOR,UNSV2D
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
      CHARACTER(LEN=32), INTENT(IN)   :: VARCLA(NVARCL)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        DOUBLE PRECISION FUNCTION TRSCE( TIME , I , ITRAC )
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: TIME
      INTEGER         , INTENT(IN) :: I,ITRAC
        END FUNCTION
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        DOUBLE PRECISION FUNCTION TR( I , ITRAC , N , IERR )
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: I,N,ITRAC
      INTEGER, INTENT(INOUT) :: IERR
        END FUNCTION
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE UTIMP_TELEMAC2D
     &                         (LTL,ATL,GRADEBL,GRAPRDL,LISDEBL,LISPRDL)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: ATL
      INTEGER, INTENT(IN) :: LTL,GRADEBL,GRAPRDL,LISDEBL,LISPRDL
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE VALRO(RO,S,ROEAU)
      USE BIEF_DEF
      IMPLICIT NONE
      TYPE(BIEF_OBJ), INTENT(IN)    :: S
      TYPE(BIEF_OBJ), INTENT(INOUT) :: RO
      DOUBLE PRECISION, INTENT(IN)  :: ROEAU
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE VISTUR(VISC,AK,EP,NPOIN,CMU,PROPNU)
      USE BIEF_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)           :: NPOIN
      DOUBLE PRECISION, INTENT(IN)  :: CMU,PROPNU
      TYPE(BIEF_OBJ), INTENT(IN)    :: AK,EP
      TYPE(BIEF_OBJ), INTENT(INOUT) :: VISC
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        DOUBLE PRECISION FUNCTION VIT( I , N )
          IMPLICIT NONE
          INTEGER, INTENT(IN) :: I,N
        END FUNCTION
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE VOLFIN 
     & (W1,AT,DT,LT,NIT,NELEM,NPTFR,
     &  TB,ZF,CF,NPOIN,HN,H,U,V,QU,QV,G,LISTIN,
     &  S,MSK,MASKEL,MESH,LIMPRO,NBOR,KDIR,KNEU,KDDL, 
     &  HBOR,UBOR,VBOR,MASSES,FLUENT,FLUSOR,CFLWTD,DTVARI,KFROT,
     &  NREJET,ISCE,TSCE2,MAXSCE,MAXTRA,YASMH,SMH,
     &  NTRAC,DIMT,T,HT,TN,DLIMT,LIMTRA,
     &  TBOR,MASSOU,FLUTENT,FLUTSOR,DTHAUT,DPX,DPY,DJX,DJY,CMI,JMI,
     &  DJXT,DJYT,DIFVIT,ITURB,PROPNU,DIFT,DIFNU,DX,DY,OPTVF,
     &  HSTOK,HCSTOK,LOGFR,DSZ,FLUXT,FLUHBOR,FLBOR,DTN,FLUSORTN,
     &  FLUENTN,LTT,
     &  FLUXTEMP,FLUHBTEMP,HC,SMTR,AIRST,TMAX,DTT,GAMMA,FLUX_OLD)
      USE BIEF_DEF
      IMPLICIT NONE 
      INTEGER, INTENT(IN)    :: NPTFR,KDIR,KNEU,KDDL,DIMT,KFROT,OPTVF
      INTEGER, INTENT(IN)    :: NELEM,NPOIN,LT,NIT,NREJET,ITURB,DLIMT
      INTEGER, INTENT(IN)    :: NTRAC,MAXSCE,MAXTRA
      INTEGER, INTENT(INOUT) :: LTT
      INTEGER, INTENT(IN)    :: LIMPRO(NPTFR,6),NBOR(NPTFR)
      INTEGER, INTENT(IN)    :: LIMTRA(DLIMT)
      INTEGER, INTENT(IN)    :: ISCE(NREJET)
      INTEGER, INTENT(INOUT) :: JMI(*),LOGFR(NPOIN)
      LOGICAL, INTENT(IN)    :: DIFVIT,DIFT,LISTIN,MSK,DTVARI,YASMH
      DOUBLE PRECISION, INTENT(IN) :: PROPNU,DIFNU,GAMMA
      DOUBLE PRECISION, INTENT(INOUT) :: AT,DT,MASSES,DTT
      DOUBLE PRECISION, INTENT(INOUT) :: H(NPOIN),QU(NPOIN),QV(NPOIN)           
      DOUBLE PRECISION, INTENT(INOUT) :: W1(*)
C                                              NSEG    NSEG               
      DOUBLE PRECISION, INTENT(INOUT) :: DSZ(2,*),HC(2,*)
      DOUBLE PRECISION, INTENT(INOUT) :: U(NPOIN),V(NPOIN) 
      DOUBLE PRECISION, INTENT(IN)    :: HN(NPOIN),SMH(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: CF(NPOIN),ZF(NPOIN),G 
      DOUBLE PRECISION, INTENT(INOUT) :: HSTOK(DIMT),HCSTOK(2,*)
      DOUBLE PRECISION, INTENT(IN)    :: HBOR(NPTFR),UBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: VBOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: TSCE2(MAXSCE,MAXTRA)
      DOUBLE PRECISION, INTENT(INOUT) :: DX(3,*),DY(3,*)
      DOUBLE PRECISION, INTENT(IN)    :: AIRST(2,*)
      DOUBLE PRECISION, INTENT(INOUT) :: DPX(3,*),DPY(3,*)
      DOUBLE PRECISION, INTENT(INOUT) :: CMI(2,*),DJX(3,*),DJY(3,*)
      DOUBLE PRECISION, INTENT(IN)    :: CFLWTD,DTHAUT(NPOIN),TMAX 
      DOUBLE PRECISION, INTENT(INOUT) :: FLUSOR,FLUENT,DTN,MASSOU(*)
      DOUBLE PRECISION, INTENT(INOUT) :: FLUSORTN,FLUENTN
      DOUBLE PRECISION, INTENT(INOUT) :: DJXT(*),DJYT(*)
      DOUBLE PRECISION, INTENT(INOUT) :: FLUTENT(*),FLUTSOR(*)    
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: TB
      TYPE(BIEF_OBJ), INTENT(IN)      :: S,MASKEL 
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
      TYPE(BIEF_OBJ) , INTENT(IN)     :: TBOR,TN
      TYPE(BIEF_OBJ) , INTENT(INOUT)  :: T,HT,SMTR,FLUHBOR,FLUHBTEMP
      TYPE(BIEF_OBJ) , INTENT(INOUT)  :: FLUXTEMP,FLUXT,FLBOR,FLUX_OLD   
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        DOUBLE PRECISION FUNCTION VUSCE( TIME , I )
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: TIME
      INTEGER         , INTENT(IN) :: I
        END FUNCTION
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE WETDRY
     &(ETA1,Z1,H1,U1,V1,ETA2,Z2,H2,U2,V2,EPS)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)    :: EPS
      DOUBLE PRECISION, INTENT(INOUT)   :: ETA1,Z1,H1,U1,V1
      DOUBLE PRECISION, INTENT(INOUT)   :: ETA2,Z2,H2,U2,V2
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        DOUBLE PRECISION FUNCTION VVSCE( TIME , I )
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: TIME
      INTEGER         , INTENT(IN) :: I
        END FUNCTION
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE WALL_FRICTION
     &(UETUTA,AUBOR,CFBOR,DISBOR,UN,VN,LIMPRO,NBOR,NPTFR,
     & KARMAN,PROPNU,LISRUG,KNEU,KDIR,KENT,KENTU,KADH,KLOG,IELMU,KP1BOR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPTFR,LISRUG,KNEU,KDIR,KENT,KADH,KLOG,KENTU
      INTEGER, INTENT(IN) :: IELMU
      INTEGER, INTENT(IN) :: LIMPRO(NPTFR,6),NBOR(NPTFR),KP1BOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: CFBOR(*),UN(*),VN(*),DISBOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: AUBOR(*),UETUTA(*)
      DOUBLE PRECISION, INTENT(IN)    :: KARMAN,PROPNU
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE WRIHYD
     &(TITRE , ITSTRT , ITSTOP , ITSTEP , NPOIN2 , MBND   ,
     & NSEG  , NOLAY  , NOMGEO , NOMLIM ,
     & F     , NSTEPA , NOMSOU , NOSUIS , NOMCOU ,
     & NOMINI, NOMVEB , NORSED , NOMSAL , NOMTEM , NOMVEL , NOMVIS ,
     & NHYD,
     & SALI_DEL,TEMP_DEL,VELO_DEL,DIFF_DEL,MARDAT,MARTIM)
      IMPLICIT NONE
      INTEGER,          INTENT(IN) :: NHYD,ITSTRT,ITSTOP,ITSTEP,NPOIN2
      INTEGER,          INTENT(IN) :: NSEG,NOLAY,NSTEPA,MBND
      INTEGER,          INTENT(IN) :: MARDAT(3),MARTIM(3)
      CHARACTER(*),     INTENT(IN) :: TITRE,NOMGEO,NOMLIM
      CHARACTER(*),     INTENT(IN) :: NOMSOU,NOSUIS,NOMCOU,NOMSAL,NOMTEM
      CHARACTER(*),     INTENT(IN) :: NOMINI,NOMVEB,NORSED,NOMVEL,NOMVIS
      DOUBLE PRECISION, INTENT(IN) :: F(NPOIN2,NOLAY)
      LOGICAL,          INTENT(IN) :: SALI_DEL,TEMP_DEL
      LOGICAL,          INTENT(IN) :: VELO_DEL,DIFF_DEL
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE ZEROPHI(X0,X,NIT,CA1)
      IMPLICIT NONE
      INTEGER, INTENT(INOUT)          :: NIT
      DOUBLE PRECISION, INTENT(IN)    :: X0,CA1
      DOUBLE PRECISION, INTENT(INOUT) :: X
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE ZEROPSI(X0,X,NIT,CA1,A2)
      IMPLICIT NONE
      INTEGER, INTENT(INOUT)          :: NIT
      DOUBLE PRECISION, INTENT(IN)    :: X0,A2,CA1
      DOUBLE PRECISION, INTENT(INOUT) :: X
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE ZOKA_SMALL(HI,HJ,ETAI,ETAJ,UI,UJ,VI,VJ,G,FLX)
      USE BIEF_DEF
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)    :: G,HI,HJ,ETAI,ETAJ,UI,UJ
      DOUBLE PRECISION, INTENT(IN)    :: VI,VJ
      DOUBLE PRECISION, INTENT(INOUT) :: FLX(3)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      END MODULE INTERFACE_TELEMAC2D
