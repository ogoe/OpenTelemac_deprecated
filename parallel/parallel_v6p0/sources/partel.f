!
! 29/01/2009:
! New version written by Christophe DENIS (EDF SINETICS) to decrease
! the computing time. This version requires more RAM but a parallel
! partel version is being designed to deacrease the amount of memory required 
!
! //// halo elements neighbourhood description added 
!      for parallel characteristics, follow four slashes 
!      jaj pinxit Thu Jul  3 09:55:31 CEST 2008 
!
! 08/08/2007
!
! Warning by JMH: there is a CALL EXIT(ICODE) which is a Fortran extension
!                 it will not work with some compilers, like Nag
!
!
!
!                       **************
                        program PARTEL
!                       **************
!
!======================================================================
! Telemac System V5P1-V2P2       Rebekka Kopmann rebekka.kopmann@baw.de
! (C) 2000-2002 BAW           Jacek A. Jankowski jacek.jankowski@baw.de
!                              Jean-Michel Hervouet j-m.hervouet@edf.fr
!======================================================================
!
!  Partitioning program for Telemac's base mesh of triangles
!
!  PARTEL (C) 2000 Copyright Bundesanstalt fuer Wasserbau, Karlsruhe
!  METIS 4.0.1 Copyright 1998, Regents of the University of Minnesota
!  BIEF (C) 2000 Electricite de France
!
!  based partially on program HANSEL, dated 12th July 1995
!  by Reinhard Hinkelmann et al., (C) University of Hannover
!
!  first  version January-March 2000 by RK (sel_metis & metis_sel)
!  second version jaj pinxit Tue Dec 12 10:48:41 MET 2000
!     (partitioning of geometry and 2D result files possible)
!  third  version jaj pinxit Fri Feb 22 15:46:23 MET 2002
!     (errors in BC values in decomposed BC files removed)
!     (erroneous treatment of islands debugged)
!  
!  fourth version Wed Apr 17 15:51:44 MDT 2002
!     (partitioning for 3D result files done by JMH)
!     (including both partitioning methods and beautifying by jaj)  
!
!  fifth version delivered Tue Jan 21 17:36:25 MDT 2002 by JMH
!     (corrected a wrong dimension of the array cut, an error
!      occuring by a larger number of processors)
!
!  sixth version delivered Mon Jan 27 11:47:07 MET 2003
!      by jaj with Matthieu Gonzales de Linares
!     (corrected a wrong dimension of the array allvar)
!  version corrected Wed Feb 19 14:17:58 MET 2003
!     by jaj for 1000 processors. 
!     BEWARE: no check of the subdomain topology
!             log must be carefully studied for Metis messages
!
!  seventh version delivered Wed Mar 12 11:47:07 MET 2003
!      by j-m hervouet
!      algorithm changed : a segment is in a subdomain if it belongs
!                          to an element in the subdomain
!                          not if the 2 points of the segment belong to the subdomain.
!      specific ELEBD included, all reference to MPI or BIEf removed
!
!  eighth version delivered 
!      by j-m hervouet in September 2003
!      ubor and vbor inverted line 613 when reading the cli file.
!
!  05/12/2006: modification by Charles Moulinec to avoid the pathological case
!              of a boundary point without following and preceding point in
!              the same subdomain. 
!              Also dimension of NBOR changed to NPTFRMAX*2
!              Look for "MOULINEC".
!          
!  19/02/2008: modifications by Olivier Boiteau (SINETICS) :
!              parameterization of array dimensions like NACHB
!              see NBSDOMVOIS AND NBMAXNSHARE, the latter must
!              be equal to its value in bief.f
!
!  15/05/2008: modification by Pascal Vezolle (IBM) :
!              Increasing of the  number of maximal partitions
!              from 1000 to 100000
!
!  16/06/2008: modification by Jean-Michel Hervouet (LNHE) :
!              Adapting VOISIN_PARTEL for meshes which are not really
!              finite element meshes (this is the case of sub-domains)
!
!  12/08/2008: received by JMH from JAJ:
!
!

! //// jaj pinxit Thu Jul  3 09:55:31 CEST 2008 
!
! :: appending partitioned BC files with elemental interface description
! :: required for parallel characteristics  
! :: explanations in the code, follow the //// marker (four slashes) 
! :: jaj notices: this program has gone too ugly, needs rewriting!
!
!
!  Fortran-90; requires linking a C subroutine
!  uses subroutine METIS_PartMeshDual from Metis library
!  available from http://www-users.cs.umn.edu/~karypis/metis/
!
!  reads: 
!        * geometry file in selafin format "<geo>"
!        * boundary conditions file "<cli>"
!        * number of required partitions (or subdomains, processors)
!          1: dual graph (each element becomes a vertex of a graph)
!          2: nodal graph (each node becomes a vertex of a graph)
!  writes:
!        * decomposed files mentioned above named with Telemac
!          convention for decomposed files in form <cli>00n-00i
!          <geo>00n-00i where n = number of partitions - 1
!           (i.e. nparts-1) and i = partition number in [0..n] 
! 
!        * metis output files; element and node partitioning:
!          named <geo>.epart.<nparts> and geo.npart.<nparts>
!          and <geo>.met (connectivity table)
!          [this files can be used by standalone Metis software]
!
! Notice:  Instead of the geometry file "<geo>" the results 
!          files of Telemac2D and Telemac3D (Selafin format) 
!          can be partitioned (as previous computation or reference
!          files. Only the last time-step results are kept when 
!          there are outputs for several time-steps.
!
! Compiler-check: SGI MIPSpro Fortran 90 Compiler Version 7.3.1.2m
!                 NAGWare Fortran 95 compiler Release 4.1(345)
!
! important assumption: no more than maxnproc processors 
!
!----------------------------------------------------------------------
!  calls Front2 (modified from BIEF, in this file), 
!        Extens (in this file), 
!        METIS_PartMeshDual, METIS_PartMeshNodal (from METIS library),
!        Voisin (in this file),
!        and Elebd (external, from BIEF)
!  [in Telemac System, BIEF library must be linked]
!----------------------------------------------------------------------
!
!
      implicit none
!
!     MAXIMUM GEOMETRICAL MULTIPLICITY OF A NODE (variable aussi
!     presente dans la BIEF, ne pas changer l'une sans l'autre)
      INTEGER, PARAMETER :: NBMAXNSHARE =  10
!     maximum number of HALO, in the parallel version the number of halo will be directly computed
      INTEGER, PARAMETER :: NBMAXHALO=100000
!
      integer, parameter :: maxnproc = 100000 ! max partition number [00000..99999]
      integer, parameter :: maxlensoft = 144 ! soft max file name length
      integer, parameter :: maxlenhard = 250 ! hard max file name length
      integer, parameter :: maxaddch = 10 ! max added suffix length
      integer, parameter :: maxvar = 100  ! max number of variables
      integer, parameter :: maxallvarlength = 3200 ! maxvar*32 for allvar
!
      integer pmethod
      integer nvar, nrec, nplan, nptfr, nptir, nptfrmax
      integer nelem, npoin, ndp, nelem2, npoin2, ndum
      integer ib(10)
!
      integer, allocatable :: ikles(:), ikles_p(:)
      integer, allocatable :: ikles3d(:),ikles3d_p(:,:,:)
      integer, allocatable :: irand(:), irand_p(:)
      integer, allocatable :: lihbor(:), liubor(:), livbor(:)
      integer, allocatable :: litbor(:)
      integer, allocatable :: npoin_p(:), nelem_p(:), nptfr_p(:)
      integer, allocatable :: nbor(:), nbor_p(:), nptir_p(:)
      integer, allocatable :: numliq(:), numsol(:)
      integer, allocatable :: knolg(:,:), knogl(:,:),check(:)
      integer, allocatable :: elelg(:,:), elegl(:)
      integer, allocatable :: cut(:), cut_p(:,:), sort(:)
      integer, allocatable :: part_p(:,:), part(:)
!
      real, allocatable    :: f(:,:), f_p(:,:,:)
      real, allocatable    :: hbor(:) 
      real, allocatable    :: ubor(:), vbor(:), aubor(:)
      real, allocatable    :: tbor(:), atbor(:), btbor(:)
!
      real times, timed
!
      integer :: ninp=10, ncli=11, nmet=12,ninpformat=52
      integer :: nepart=15, nnpart=16, nout=17, nclm=18
      integer time(3), date(3)
!
      character(len=80)  :: title
      character(len=32)  :: vari, variable(maxvar)
      character(len=maxallvarlength) :: allvar 
      character(len=maxlenhard)  :: nameinp, namecli, nameout, nameclm
      character(len=maxlenhard)  :: namemet,nameepart,namenpart,
     c     nameninpformat,nameoutforma
      character(len=5)   :: chch  
      character(len=12)  :: fmt4
!
      integer max_nelem_p, min_nelem_p
      integer  max_npoin_p,max_n_neigh
      integer i, j, k, l , m, n, p, err, iso, idum
      integer istop, istart, iseg, ii, iloop
      integer i_len, i_s, i_sp, i_lencli, i_leninp
      integer ielem_p, ipoin_p, iptfr_p
!
      real xseg, yseg, bal, rdum
      double precision area, x1, x2, x3, y1, y2, y3
      logical is, wrt, timecount
!
! Metisology
!
      integer nparts, etype, numflag, edgecut
      integer, allocatable :: epart(:), npart(:)
      character(len=10) fmt1, fmt2, fmt3
!
! for calling front2
!
      integer, parameter :: maxfro = 300   ! max number of boundaries
      integer nfrliq, nfrsol, debliq(maxfro), finliq(maxfro)
      integer debsol(maxfro), finsol(maxfro)
      integer, allocatable :: dejavu(:), kp1bor(:,:)
!
! for calling BIEF mesh subroutines (to be optimised soon):
      integer, allocatable :: ifabor(:,:), ifanum(:,:), nelbor(:)
      integer, allocatable :: nulone(:,:)
      integer, allocatable :: ikle(:,:), iklbor(:,:), isegf(:)
      integer, allocatable :: it1(:), it2(:), it3(:)

      integer npoin_tot

      integer lng,lu,li
      common /info/ lng,lu
!
! time measuring 
!
      integer  tdeb, tfin, tdebp, tfinp, temps, parsec
      integer  tdeb_glob, tfin_glob
      integer  time_in_seconds
      external time_in_seconds
!
! extens function
!
      character(len=11) :: extens
      external extens   
!
!----------------------------------------------------------------------
!
!jaj new for parallel characteristics ////
! halo elements: these adjacent to the interface edges having 
! neighbours behind a boundary 
!
      ! the elemental global->local numbering translation table 
      ! this is elegl saved from all partitions for further use
      INTEGER, ALLOCATABLE :: gelegl(:,:)
!
      ! the halo elements neighbourhood description for a halo cell 
      INTEGER, ALLOCATABLE :: ifapar(:,:,:)
!
      ! the number of halo cells pro partition 
      INTEGER, ALLOCATABLE :: nhalo(:) 
!
      ! work variables 
      INTEGER ifaloc(3)
      LOGICAL found
      INTEGER ndp_2d,ndp_3d
      INTEGER ef,pos
      INTEGER, ALLOCATABLE :: nbre_ef(:),nbre_ef_loc(:),ef_i(:),
     c     tab_tmp(:),ef_ii(:)
      LOGICAL trouve,halo
      INTEGER noeud,nbre_noeud_interne
      INTEGER nbre_noeud_interf
      INTEGER frontiere,nbre_ef_i
      LOGICAL interface      

! #### for sections 

      TYPE chain_type
        INTEGER :: npair(2)
        DOUBLE PRECISION :: xybeg(2), xyend(2)
        CHARACTER(LEN=24) :: descr
        INTEGER :: nseg
        INTEGER, POINTER :: liste(:,:) 
      END TYPE 
      TYPE (chain_type), ALLOCATABLE :: chain(:)
      INTEGER, PARAMETER :: nsemax=500 ! max number of segments in a section 
      INTEGER, ALLOCATABLE :: liste(:,:), anpbeg(:),anpend(:) 
      INTEGER :: nsec, ihowsec, isec, ielem, im(1), in(1), npbeg, npend 
      INTEGER :: ncp, pt, i1,i2,i3, arr,dep, ilprec,ilbest,elbest,igbest 
      DOUBLE PRECISION :: xa, ya, distb, diste, dminb, dmine
      DOUBLE PRECISION :: dist1, dist2, dist3, dist
      CHARACTER(len=maxlenhard) :: namesec
      LOGICAL :: with_sections=.FALSE.

!
!----------------------------------------------------------------------
!
      ndp_2d=3
      ndp_3d=6

      call system_clock (count=temps, count_rate=parsec)
      timecount = .true.
      if (parsec==0) timecount = .false.  ! count_rate == 0 : no clock
      if (timecount) tdeb = temps
!
      lng=2 ! je ne parle francais, je suis barbarien
      lu=6  ! fortran standard ouput channel
      li=5  ! fortran standard input channel
       

!----------------------------------------------------------------------
! names of the input file to eventually guide to PARES3D
! if parallel computation with ESTEL3D 
!
!
!=>Fabs
!      do 
!        read(li,'(a)')nameinp
!        if (nameinp /= ' ') exitabout:
!      enddo
!      if (nameinp(1:3)=='ES3') then
! PARTEL adapted to estel3d code
!        call pares3d(nameinp,li)
! back to the end of partelabout:
!        goto 299
!      else
! continue with telemac codes
!        rewind li
!      endif 
!<=Fabs
!
!----------------------------------------------------------------------
! introduce yourself
!
      write(lu,*) ' '
      write(lu,*) '+-------------------------------------------------+'
      write(lu,*) '  PARTEL: Telemac Selafin Metisologic Partitioner'
      write(lu,*) '                                                   '          
      write(lu,*) '  Rebekka Kopmann & Jacek A. Jankowski (BAW)'
      write(lu,*) '                 Jean-Michel Hervouet (LNHE)'
      write(lu,*) '                 Christophe Denis     (SINETICS) '
      write(lu,*) '  PARTEL (C) Copyright 2000-2002 '
      write(lu,*) '  Bundesanstalt fuer Wasserbau, Karlsruhe'
      write(lu,*) ' '
      write(lu,*) '  METIS 4.0.1 (C) Copyright 1998 '
      write(lu,*) '  Regents of the University of Minnesota '
      write(lu,*) ' '
      write(lu,*) '  BIEF 5.9 (C) Copyright 2008 EDF'
      write(lu,*) '+-------------------------------------------------+'
      write(lu,*) ' '
!jaj ////
      write(lu,*) '  => This is a preliminary development version '
      write(lu,*) '     Dated:  Tue Jan 27 11:11:20 CET 2009'
      write(lu,*) ' '
      write(lu,*) '  Maximum number of partitions: ',maxnproc
      write(lu,*) ' '
      write(lu,*) '+--------------------------------------------------+'
      write(lu,*) ' '
!
!----------------------------------------------------------------------
! names of the input files:
!
      do 
        write(lu, advance='no', fmt=
     &         '(/,'' Selafin input name <input_name>: '')')
        read(li,'(a)') nameinp
        if (nameinp.eq.' ') then
          write (lu,'('' no filename'')') 
        else
!=>Fabs
          if (nameinp(1:3)=='ES3') then
! PARTEL adapted to estel3d code
            call pares3d(nameinp,li)
            goto 299
          else
!<=Fabs
! continue with telemac codes
            write(lu,*) 'input: ',nameinp
            exit 
!=>Fabs
          endif
!<=Fabs
        end if  
      end do

      inquire (file=nameinp,exist=is)
      if (.not.is) then 
        write (lu,'('' file does not exist: '',a30)') nameinp
        call plante2(-1)
        stop
      end if  
!
      do
        write(lu, advance='no', fmt=
     &           '(/,'' Boundary conditions file name : '')')
        read(li,'(a)') namecli
        if (namecli.eq.' ') then
          write (lu,'('' no filename'')') 
        else
          write(lu,*) 'input: ',namecli
          exit
        end if
      end do
!  
      inquire (file=namecli,exist=is)
      if (.not.is) then 
        write (lu,'('' file does not exist: '',a30)') namecli
        call plante2(-1)
        stop
      end if  
!
      do 
        write(lu, advance='no',fmt=
     &    '(/,'' Number of partitions <nparts> [2 -'',i6,'']: '')') 
     &        maxnproc
        read(li,*) nparts
        if ( (nparts > maxnproc) .or. (nparts < 2) ) then
          write(lu,
     &    '('' Number of partitions must be in [2 -'',i6,'']'')') 
     &      maxnproc
        else
          write(lu,'('' input: '',i4)') nparts
          exit
        end if 
      end do
!
      write(lu,fmt='(/,'' Partitioning options: '')')
!      write(lu,*) '  1: DUAL  graph', 
!     & ' (each element of the mesh becomes a vertex of the graph)'
!      write(lu,*) '  2: NODAL graph', 
!     & ' (each node of the mesh becomes a vertex of the graph)'

      do 
        write(lu, advance='no',fmt=
     &    '(/,'' Partitioning method <pmethod> [1 or 2]: '')') 
        read(li,*) pmethod
        if ( (pmethod > 2) .or. (pmethod < 1) ) then
          write(lu,
     &    '('' Partitioning method must be 1 or 2'')') 
        else
          write(lu,'('' input: '',i3)') pmethod
          exit
        end if 
      end do
!
! #### the sections file name 

      DO
        WRITE(lu, ADVANCE='no',FMT=
     &    '(/,'' With sections? [1:YES 0:NO]: '')') 
        READ(li,*) i
        IF ( i<0 .OR. i>1 ) THEN
          WRITE(lu,
     &    '('' Please answer 1:YES or 0:NO '')') 
        ELSE
          WRITE(lu,'('' input: '',i4)') i
          EXIT
        END IF 
      END DO
      IF (i==1) with_sections=.TRUE.


      IF (with_sections) THEN 
        DO
          WRITE(lu, ADVANCE='no', FMT=
     &      '(/,'' Control sections file name (or RETURN) : '')')
          READ(li,'(a)') namesec
          IF (namesec.EQ.' ') THEN
            WRITE (lu,'('' no filename '')') 
          ELSE
            WRITE(lu,*) 'input: ',namesec
            EXIT
          ENDIF
        END DO
!  
        INQUIRE (FILE=namesec,EXIST=is)
        IF (.NOT.is) THEN
          WRITE (lu,'('' file does not exist: '',a30)') namesec
          CALL plante2(-1)
          STOP
        ENDIF  
      ENDIF
!
! find the input file core name length
!
      i_s  = len(nameinp)
      i_sp = i_s + 1
      do i=1,i_s
        if (nameinp(i_sp-i:i_sp-i) .NE. ' ') exit
      enddo
      i_len=i_sp - i
      i_leninp = i_len
!
      if (i_leninp > maxlensoft) then
        write(lu,*) ' '
        write(lu,*) 'ATTENTION:'
        write(lu,*) 'The name of the input file:'
        write(lu,*) nameinp
        write(lu,*) 'is longer than ',maxlensoft,' characters' 
        write(lu,*) 'which is the longest applicable name for Telemac '
        write(lu,*) 'input and output files. Stopped. '
        call plante2(-1)
        stop
      endif
!
      namemet = nameinp(1:i_leninp)//'.met'
!
      open(ninp,file=nameinp,status='OLD',form='UNFORMATTED')
      rewind ninp
!
!----------------------------------------------------------------------
! start reading the geometry or result file
!
!
      read (ninp) title
      read (ninp) i, j
      nvar = i + j 
    
      allvar(1:41) = 'X-COORDINATE----M---,Y-COORDINATE----M---'
      istart = 42
!
      write (lu,*) 'variables are: '
      do i=1,nvar
        read(ninp) vari
       
        variable(i) = vari
     
        do j=1,32
          if(vari(j:j).eq.' ') vari(j:j) = '-'
        end do
        istop = istart+20
        if (istop.gt.maxallvarlength) then
          write(lu,*) 'variable names too long for string allvar'
          write(lu,*) 'stopped.'
          call plante2(-1)
          stop
        endif
        allvar(istart:istart) = ','
        allvar(istart+1:istop) = vari
        istart=istop+1
      enddo 
!
! read the rest of the selafin file
! 10 integers, the first is the number of records (timesteps)
!
      read (ninp) (ib(i), i=1,10)
      if (ib(8).ne.0.or.ib(9).ne.0) then
        write(lu,*) 'this is a partial output file'
        write(lu,*) 'maybe meet Gretel before...'
      endif 
      nrec  = ib(1)
      nplan = ib(7) 
      if (ib(10).eq.1) then 
        read(ninp) date(1), date(2), date(3), time(1), time(2), time(3)
        
      endif 
!
      read (ninp) nelem,npoin,ndp,ndum
      npoin_tot=npoin
      if (nplan.gt.1) then 
        write(lu,*) ' '
        write(lu,*) '3D mesh detected.' 
        npoin2 = npoin/nplan
        nelem2 = nelem/(nplan-1)
        write(lu,*) 'ndp nodes per element:             ',ndp
        write(lu,*) 'nplan number of mesh levels:       ',nplan
        write(lu,*) 'npoin2 number of 2D mesh nodes:    ',npoin2
        write(lu,*) 'npoin number of 3D mesh nodes:     ',npoin
        write(lu,*) 'nelem2 number of 2D mesh elements: ',nelem2
        write(lu,*) 'nelem number of 3D mesh elements:  ',nelem
        if (mod(npoin,nplan).ne.0) then 
          write (lu,*) 'But npoin2 /= npoin3/nplan!'
          call plante2(-1)
          stop   
        endif
        if (mod(nelem,(nplan-1)).ne.0) then 
          write (lu,*) 'But nelem2 /= nelem3/nplan!'
          call plante2(-1)
          stop
        endif
        write(lu,*) ' '
      else
        write(lu,*) ' '
        write(lu,*) 'one-level mesh.'
        write(lu,*) 'ndp nodes per element:         ',ndp
        write(lu,*) 'npoin number of mesh nodes:    ',npoin
        write(lu,*) 'nelem number of mesh elements: ',nelem
        write(lu,*) ' '
        npoin2 = npoin
        nelem2 = nelem
      endif
!
      if (ndp.eq.3) then  
        write(lu,*) 'The input file assumed to be 2D Selafin'
      elseif (ndp.eq.6) then
        write(lu,*) 'The input file assumed to be 3D Selafin'
      else   
        write(lu,*) 'the elements are neither triangles nor prisms!'
        write(lu,*) 'ndp = ',ndp
        call plante2(-1)
        stop
      endif
!
! now let us allocate 
!
      allocate (ikles(nelem2*3),stat=err)
      if (err.ne.0) call ALLOER (lu, 'ikles')
      if(nplan.gt.1) then
        allocate (ikles3d(nelem*ndp),stat=err)
        if (err.ne.0) call ALLOER (lu, 'ikles3d')
      endif 
      allocate (irand(npoin),stat=err)
      if (err.ne.0) call ALLOER (lu, 'irand')
!     nvar+2 : first two functions are x and y
!     npoin is 3D here in 3D
      allocate (f(npoin,nvar+2),stat=err)
      if (err.ne.0) call ALLOER (lu, 'f')
!
! connectivity table:
!
      if(nplan.le.1) then
        read(ninp) ((ikles((k-1)*ndp+j),j=1,ndp),k=1,nelem)
        
      else
        read(ninp) ((ikles3d((k-1)*ndp+j),j=1,ndp),k=1,nelem)
!       building ikles
        do j=1,3
          do k=1,nelem2
            ikles((k-1)*3+j)=ikles3d((k-1)*6+j)
          enddo
        enddo
      endif
!
! boundary nodes indications
!
      read(ninp) (irand(j),j=1,npoin)
! computation of nptfr done later with the boundary conditions file
! (modification by j-m hervouet on 10/04/02)
! irand is not always correct and may lead to errors
!
! number of boundary points in 2D mesh
!      nptfr = 0
!      do j=1,npoin2
!        if(irand(j).ne.0) nptfr = nptfr+1
!      end do 
!      write (lu,*) ' '
!      write (lu,*) 'nptfr number of boundary nodes in 2D mesh',nptfr
!      write (lu,*) ' '
!
! x-, y-coordinates
!
      read(ninp) (f(j,1),j=1,npoin)
      read(ninp) (f(j,2),j=1,npoin)
!     
! now the loop over all records (timesteps) - for an initial 
! conditions file automatically the last time step values are 
! taken (!)
!
      iloop = 0
      do 
!
! read the time step
!
        read(ninp, end=111, err=300) times
        iloop = iloop + 1
!
        timed = times/3600
        write(lu,*) 'timestep: ',times,'s = ',timed,'h'
!
! read the time variables; no 1 and 2 are x,y
!
        do k=3,nvar+2
!          write(lu,*) 'now reading variable',k-2
          read(ninp, end=300, err=300) (f(j,k), j=1,npoin)
!          write(lu,*) 'reading variable',k-2,' successful'
        end do
      end do
 111  close (ninp)
 !     write(lu,*) ' '
 !     write(lu,*) 'There has been ',iloop,' time-dependent recordings'
 !     write(lu,*) 'Only the last one taken into consideration'
 !     write(lu,*) ' '
!
!-----------------------------------------------------------------------
! ...check if the area of the elements are negative...
! ... area = 0.5*abs(x1*y2 - y1*x2 + y1*x3 - x1*y3 + x2*y3 - y2*x3)
! notice: area and x1, y1, x2, y2, x3, y3 must be double precision
!
!        do j=1,nelem
!          x1 = f(ikles((j-1)*3+1),1)
!          y1 = f(ikles((j-1)*3+1),2)
!          x2 = f(ikles((j-1)*3+2),1)
!          y2 = f(ikles((j-1)*3+2),2)
!          x3 = f(ikles((j-1)*3+3),1)
!          y3 = f(ikles((j-1)*3+3),2)
!          area = x1*y2-y1*x2+y1*x3-x1*y3+x2*y3-y2*x3
!          if ( area < 0.0 ) then
!            write(lu,*) 'Global domain'
!            write(lu,*) 'Determinant of element',j,' is negative'
!            write(lu,*) '(local node orientation is clockwise!)'
!            write(lu,*) 'Det-value: ',area
!            write(lu,*) 'node nr 1, x1,y1: ',ikles((j-1)*3+1),x1,y1
!            write(lu,*) 'node nr 2, x2,y2: ',ikles((j-1)*3+2),x2,y2
!            write(lu,*) 'node nr 3, x3,y3: ',ikles((j-1)*3+3),x3,y3
!            call plante2(-1)
!            stop
!          endif
!        end do
!
!----------------------------------------------------------------------
! read the boundary conditions file
!
!      write(lu,*) ' '
!      write(lu,*) '--------------------------'
!      write(lu,*) '  BC file: ',namecli
!      write(lu,*) '--------------------------'
!      write(lu,*) ' '
!
! but allocate first
!
      nptfrmax = npoin2   ! better idea ?
!
      allocate (lihbor(nptfrmax),stat=err)
      if (err.ne.0) call ALLOER (lu, 'lihbor')
      allocate (liubor(nptfrmax),stat=err)
      if (err.ne.0) call ALLOER (lu, 'liubor')
      allocate (livbor(nptfrmax),stat=err)
      if (err.ne.0) call ALLOER (lu, 'livbor')
      allocate (hbor(nptfrmax),stat=err)
      if (err.ne.0) call ALLOER (lu, 'hbor')
      allocate (ubor(nptfrmax),stat=err)
      if (err.ne.0) call ALLOER (lu, 'ubor')
      allocate (vbor(nptfrmax),stat=err)
      if (err.ne.0) call ALLOER (lu, 'vbor')
      allocate (aubor(nptfrmax),stat=err)
      if (err.ne.0) call ALLOER (lu, 'aubor')
      allocate (tbor(nptfrmax),stat=err)
      if (err.ne.0) call ALLOER (lu, 'tbor')
      allocate (atbor(nptfrmax),stat=err)
      if (err.ne.0) call ALLOER (lu, 'atbor')
      allocate (btbor(nptfrmax),stat=err)
      if (err.ne.0) call ALLOER (lu, 'btbor')
      allocate (litbor(nptfrmax),stat=err)
      if (err.ne.0) call ALLOER (lu, 'litbor')
      allocate (nbor(2*nptfrmax),stat=err)  ! for front2
      if (err.ne.0) call ALLOER (lu, 'nbor')
      allocate (numliq(nptfrmax),stat=err)
      if (err.ne.0) call ALLOER (lu, 'numliq')
      allocate (numsol(nptfrmax),stat=err)
      if (err.ne.0) call ALLOER (lu, 'numsol')
      allocate (check(nptfrmax),stat=err)
      if (err.ne.0) call ALLOER (lu, 'check')
!
! core name length
!
      i_s  = len(namecli)
      i_sp = i_s + 1
      do i=1,i_s
         if (namecli(i_sp-i:i_sp-i) .ne. ' ') exit
      enddo
      i_len=i_sp - i
      i_lencli = i_len
!
      if (i_leninp > maxlensoft) then
        write(lu,*) ' '
        write(lu,*) 'ATTENTION:'
        write(lu,*) 'The name of the boundary conditions file:'
        write(lu,*) namecli
        write(lu,*) 'is longer than ',maxlensoft,' characters' 
        write(lu,*) 'which is the longest applicable name for Telemac '
        write(lu,*) 'input and output files. Stopped. '
        call plante2(-1)
        stop
      endif
!
      open(ncli,file=namecli,status='OLD',form='FORMATTED')
      rewind ncli
!
!     reading boundary file and counting boundary points
!
      k=1
 900  continue
      read(ncli,*,end=901,err=901) lihbor(k),liubor(k),
     &                             livbor(k),
     &             hbor(k),ubor(k),vbor(k),aubor(k),litbor(k),
     &             tbor(k),atbor(k),btbor(k),nbor(k),check(k)
!
!     Now check is the boundary node colour
!     if(check(k).ne.k) then
!       write(lu,*) 'Error in boundary conditions file at line ',k
!       call plante2(-1)
!       stop
!     endif
      k=k+1
      goto 900
 901  continue
      nptfr = k-1
!      write (lu,*) ' '
!      write (lu,*) 'Number of boundary nodes in 2D mesh: ',nptfr
!      write (lu,*) ' '
      close(ncli)
!
!----------------------------------------------------------------------
! Numbering of open boundaries 
! numbering of liquid boundary, if 0 = solid
! opn: number of open boundary
! in order to do it in the same way as Telemac does, 
! it is best to call front2 here
!
! for calling BIEF mesh subroutines
! can be optimised / uses a lot of memory 
! the only reason is to obtain kp1bor and numliq
!
      allocate (dejavu(nptfr),stat=err)
      if (err.ne.0) call ALLOER (lu, 'dejavu')
      allocate (kp1bor(nptfr,2),stat=err)
      if (err.ne.0) call ALLOER (lu, 'kp1bor')
!jaj----------v ////
!     changed nelem to nelem2, ndp to 3 huh! 
!     causing errors when 3D restart/reference files are partitioned
!     and BC file is written again (what for, actually???) 
!     cause: calling voisin with nelem2 but ifabor(nelem=nelem3,ndp=6)
      allocate (ifabor(nelem2,3),stat=err)
      if (err.ne.0) call ALLOER (lu, 'ifabor')
      allocate (ifanum(nelem2,3),stat=err)
      if (err.ne.0) call ALLOER (lu, 'ifanum')
      allocate (iklbor(nptfr,2),stat=err)
      if (err.ne.0) call ALLOER (lu, 'iklbor')
      allocate (nelbor(nptfr),stat=err)
      if (err.ne.0) call ALLOER (lu, 'nelbor')
      allocate (nulone(nptfr,2),stat=err)
      if (err.ne.0) call ALLOER (lu, 'nulone')
      allocate (isegf(nptfr),stat=err)
      if (err.ne.0) call ALLOER (lu, 'isegf')
      allocate (ikle(nelem2,3),stat=err)
      if (err.ne.0) call ALLOER (lu, 'ikle')
      allocate (it1(npoin),stat=err)
      if (err.ne.0) call ALLOER (lu, 'it1')
      allocate (it2(npoin),stat=err)
      if (err.ne.0) call ALLOER (lu, 'it2')
      allocate (it3(npoin),stat=err)
      if (err.ne.0) call ALLOER (lu, 'it3')
!
! transform ikles--> ikle for 2D routines  (an old telemac disease) 
!
      do i = 1,3
        do j  = 1,nelem2
          ikle(j,i) = ikles((j-1)*3+i)
        enddo
      enddo
!
      call VOISIN_PARTEL(ifabor, nelem2, nelem2, 11, ikle, nelem2,
     &                   npoin2, it1, it2)
!
!      write(lu,'(/,'' Calling ELEBD'')')
!
      call ELEBD_PARTEL (nelbor, nulone, kp1bor, ifabor, nbor, ikle, 
     &                   nelem2, iklbor, nelem2, nelem2, nptfrmax,
     &                   npoin2, nptfr, 11, lihbor, 2, ifanum,
     &                   1, isegf, it1, it2, it3,npoin_tot )
!
!      write(lu,'(/,'' Boundary type numbering using FRONT2'')')
!      
      if (nameinp(1:3)== 'ART') THEN
         OPEN(UNIT=89,FILE='front_glob.dat')
         write(89,*) NPOIN_TOT
         write(89,*) NPTFR
         DO K=1,NPTFR
            write(89,*) NBOR(K)
         END DO 
         DO K=1,NPTFR
            write(89,*) KP1BOR(K,1)
         END DO
         DO K=1,NPTFR
            write(89,*) KP1BOR(K,2)
         END DO 
         call flush(89)
         close(89)
      end if
      call FRONT2_PARTEL (nfrliq,nfrsol,debliq,finliq,debsol,finsol,
     &             lihbor,liubor,f(1:npoin2,1),f(1:npoin2,2),
     &             nbor,kp1bor(1:nptfr,1),dejavu,npoin2,nptfr,
     &             2,.true.,numliq,numsol,nptfrmax)  
!
      deallocate (dejavu)
!jaj //// ifabor applied later for finding halo cell neighbourhoods 
!!!!      deallocate (ifabor)
      deallocate (ifanum)
      deallocate (iklbor)
!     deallocate (nelbor)
      deallocate (nulone)
      deallocate (isegf)
!      deallocate (ikle) !jaj #### we need it for sections 
      deallocate (it1)
      deallocate (it2)
      deallocate (it3)
!    COMMENTED BY CD 

!----------------------------------------------------------------------
! open and rewrite Metis software input files
! not necessary if visualisation or manual decompositions 
! are not required
!
c$$$      write(lu,*) ' '
c$$$      write(lu,*) '---------------------------'
c$$$      write(lu,*) ' Metis & PMVIS input files '
c$$$      write(lu,*) '---------------------------'
c$$$      write(lu,*) ' '
c$$$!
c$$$      open(nmet,file=namemet,status='UNKNOWN',form='FORMATTED')
c$$$      rewind nmet
c$$$      write(lu,*) 'Input file for partitioning: ', namemet
c$$$!
c$$$! the first line is not necessary in the latest version
c$$$! we write the files using C convention
c$$$!
c$$$! here the IKLE 2D is written, even in 3D (hence ndp considered to be 3)
c$$$!
c$$$      write(nmet,*) nelem2,'1'
c$$$      do k=1,nelem2
c$$$        write(nmet,'(3(i7,1x))') (ikles((k-1)*3+j)-1, j=1,3)
c$$$      end do
c$$$      close(nmet)
c$$$!
! write the node coordinates for visualisation 
! a check first...
!
c$$$      write (lu,'(/,'' As coordinates for visualisation taken: '')')
c$$$      write (*,'(1x,a20)') allvar(1:20)
c$$$      write (*,'(1x,a20)') allvar(22:41)
c$$$      write (*,'(1x,a20)') allvar(43:62)
!
!======================================================================
! partitioning
!
!

      !======================================================================
! Step 2 : partitioning the mesh 
!
! other partitioning methods should be used (SCOTCH for example)
!     all processors perform this task to avoid communication
!     The use of ParMetis or PTSCOTCH could be used for larger meshes
!     if there will be some memory allocation problem 
!======================================================================      
    
      allocate (epart(nelem2),stat=err)
      if (err.ne.0) call ALLOER (lu, 'epart')
      allocate (npart(npoin2),stat=err)
      if (err.ne.0) call ALLOER (lu, 'npart')
!
      if (ndp==3.or.ndp==6) then 
         etype = 1
      else
         write(lu,*) 'METIS: implemented for triangles or prisms only'
         call plante2(-1)
         stop
      endif 
      
! We only use METIS_PartMeshDual as only the finite elements partition
!     is relevant here.   
!     
!     IMPORTANT: we use fortran-like field elements numbering 1...n
!     in C version, 0...n-1 numbering is applied!!!
!     
      numflag = 1
!
      write(lu,*) 'Using only METIS_PartMeshDual subroutine'
      
      write(lu,*) ' The mesh partitioning step METIS starts'
      if (timecount) then 
         call system_clock (count=temps, count_rate=parsec)
         tdebp = temps
      endif
      call METIS_PartMeshDual 
     &     (nelem2, npoin2, ikles, etype, numflag, 
     &     nparts, edgecut, epart, npart)
     
      write(lu,*) ' The mesh partitioning step has finished'
      if (timecount) then
        call system_clock (count=temps, count_rate=parsec)
        tfinp = temps
        write(lu,*) ' Runtime of METIS ',
     &            (1.0*(tfinp-tdebp))/(1.0*parsec),' seconds'
      endif

      
!======================================================================
! Step 3 : allocate the global  arrays not depending of the partition
!     
!======================================================================   
 
!      write(lu,*) 'here '  
!     knogl(i) =>  global label of the local point i 
      allocate (knogl(npoin2,nparts),stat=err)
      if (err.ne.0) call ALLOER (lu, 'knogl')
      knogl(:,:)=0
      
!     nbre_ef(i) => number of finite element containing i
!     i is a global label 
      allocate (nbre_ef(npoin2),stat=err)
      if (err.ne.0) call ALLOER (lu, 'nbre_ef')
      
      if(nplan.eq.0) then
         allocate (f_p(npoin2,nvar+2,nparts),stat=err)
      else
         allocate (f_p(npoin2,nvar+2,nparts),stat=err)
      endif
      if (err.ne.0) call ALLOER (lu, 'f_p')
      
      allocate (part_p(npoin2,0:NBMAXNSHARE),stat=err)
      if (err.ne.0) call ALLOER (lu, 'part_p')
      part_p(:,:)=0
      
      allocate (cut_p(npoin2,nparts),stat=err)
      if (err.ne.0) call ALLOER (lu, 'cut_p')
      
      ALLOCATE (gelegl(nelem2,nparts),stat=err) 
      IF (err.ne.0) CALL ALLOER (lu, 'gelegl')
      
      allocate (sort(npoin2),stat=err)
      if (err.ne.0) call ALLOER (lu, 'cut_p')
      
      allocate (cut(npoin2),stat=err)
      if (err.ne.0) call ALLOER (lu, 'cut_p')
      
      allocate (nelem_p(nparts),stat=err)
      if (err.ne.0) call ALLOER (lu, 'nelem_p')
       
      allocate (npoin_p(nparts),stat=err)
      if (err.ne.0) call ALLOER (lu, 'npoin_p')
      
      allocate (nptfr_p(nparts),stat=err)
      if (err.ne.0) call ALLOER (lu, 'nptfr_p')
      
      allocate (nptir_p(nparts),stat=err)
      if (err.ne.0) call ALLOER (lu, 'nptir_p')
      
      allocate (nhalo(nparts),stat=err)
      if (err.ne.0) call ALLOER (lu, 'nhalo')

      allocate(tab_tmp( NBMAXNSHARE),stat=err)
      IF (err.ne.0) CALL ALLOER (lu, 'tab_tmp')
      
      allocate(ifapar(nparts,7,NBMAXHALO),stat=err)
         IF (err.ne.0) CALL ALLOER (lu, 'ifapar')
         ifapar(:,:,:)=0

!======================================================================
! Step 4 : Compute the number of finite elements and points
!     belonging to submesh i
!
!======================================================================   

      
!     Firstly, all MPI processes  work on the whole mesh
!     ----------------------------------------------      
!   
!     loop over the finite element of the mesh 
!     to compute the number of finite elements containing each point noeud
         if (nameinp(1:3) == 'ART') then     
            do ef=1,nelem2
               do k=1,ndp_2d
                  noeud=ikles((ef-1)*3+k)
                  if (irand(noeud) .ne. 0) then
                     epart(ef)=1
                  end if
               end do 
            end do
         end if
         
         nbre_ef(:)=0
      do ef=1,nelem2
         do k=1,ndp_2d
            noeud=ikles((ef-1)*3+k)
            nbre_ef(noeud)=nbre_ef(noeud)+1
         end do
      end do 
      do i=1,nparts
         
     

!     loop over the finite element of the mesh to compute 
!     the number of the finite element and points belonging 
!     to submesh i
   
         nelem_p(i)=0
         npoin_p(i)=0
         do ef=1,nelem2
            if (epart(ef) .eq. i) then
               nelem_p(i)=nelem_p(i)+1
               do k=1,ndp_2d
                  noeud=ikles((ef-1)*3+k)
                  if (knogl(noeud,i) .eq. 0) then
                     npoin_p(i)=npoin_p(i)+1
                     knogl(noeud,i)=npoin_p(i)
                  end if
               end do 
            end if
         end do
      end do  
    
!======================================================================
!     Step 4 : Allocation of local arrays needed by MPI Processus id
!              working on submesh id+1
!======================================================================   
 !     write(lu,*) 'after the first loop'
      max_nelem_p=maxval(nelem_p)
      max_npoin_p=maxval(npoin_p)


!     elegl(e) => global label of the finite element e
!     e is the local label on submesh i 
      allocate (elelg(max_nelem_p,nparts),stat=err)
      if (err.ne.0) call ALLOER (lu, 'elelg')
      elelg(:,:)=0
!     knolg(i) => global label of the point i
!     i is the local label on subdomain i
      if(nplan.eq.0) then
         allocate (knolg(max_nelem_p,nparts),stat=err)
      else
         allocate (knolg(max_npoin_p*nplan,nparts),stat=err)
      endif
      if (err.ne.0) call ALLOER (lu, 'knolg')
      knolg(:,:)=0
!     nbre_ef_loc(i) : number of finite elements containing the point i
!                      on submesh i  
!     i is the local label on submesh i
      allocate (nbre_ef_loc(max_nelem_p),stat=err)
      if (err.ne.0) call ALLOER (lu, 'nbre_ef_loc')

!     ef_i(e) is the global label of the interface finite element number e
      allocate (ef_i(max_nelem_p),stat=err)
      if (err.ne.0) call ALLOER (lu, 'ef_i')
!     ef_ii(e) is the local label of the interface finite element number e
      allocate (ef_ii(max_nelem_p),stat=err)
      if (err.ne.0) call ALLOER (lu, 'ef_ii')

      

!======================================================================
!     Step 5 : Initialisation  of local arrays 
!                  (gelelg and elelg, nbre_ef_loc)
!              
!======================================================================   
      do i=1,nparts
         nelem_p(i)=0
         do ef=1,nelem2
            if (epart(ef) .eq. i) then
               nelem_p(i)=nelem_p(i)+1
               elelg(nelem_p(i),i)=ef
               gelegl(ef,i)=nelem_p(i)
            end if
         end do
         do j=1,npoin_p(i)
            nbre_ef_loc(j)=0
         end do
        
!======================================================================
!     Step 5 : Compute the number of boundary and interface points
!              Initialisation of nbre_ef_loc and f_p 
!======================================================================   
         
         npoin_p(i)=0
         nptfr_p(i)=0
         nbre_noeud_interne=0
         nbre_noeud_interf=0
         
         do j=1,nelem_p(i)
            ef=elelg(j,i)
            do k=1,3
               noeud=ikles((ef-1)*3+k)
               nbre_ef_loc(knogl(noeud,i))=
     c              nbre_ef_loc(knogl(noeud,i))+1 
               if (nbre_ef_loc(knogl(noeud,i)) .eq. 1) then
!     the point noeud is encountered for the first time 
                  npoin_p(i)=npoin_p(i)+1    
!     is noeud a boundary point ?     
                  if (irand(noeud) .NE. 0) then
                     nptfr_p(i)= nptfr_p(i)+1
                  end if
!     modification of   knogl et f_p
                  knolg(npoin_p(i),i)=noeud
                  do l=1,nvar+2
                     f_p(npoin_p(i),l,i)=f(noeud,l)
                  end do
               end if
!    
!     noeud is a internal point if all finite elements
!     containing it belongs to the same submesh
               if (nbre_ef_loc(knogl(noeud,i)) .EQ. nbre_ef(noeud)) then 
                  nbre_noeud_interne=nbre_noeud_interne+1
               end if
            end do
         end do
        
         nbre_noeud_interf=npoin_p(i)-nbre_noeud_interne
         nptir_p(i)=0 
!     nombre de noeud interface du SDi
         nbre_ef_i=0            ! nombre d'elements finis interfaces du SDi
         do j=1,nelem_p(i)      ! On parcours a nouveau les elements finis du SDi
            interface=.false.
            ef=elelg(j,i)
            do k=1,ndp_2d
               noeud=ikles((ef-1)*3+k)
               if (abs(nbre_ef_loc(knogl(noeud,i))) .ne. nbre_ef(noeud))
     c          then
                  interface=.true.
               end if
               if (nbre_ef_loc(knogl(noeud,i)) .ne.  nbre_ef(noeud).and. 
     c              nbre_ef_loc(knogl(noeud,i)) .gt. 0) then
!     noeud est interface car il reste des elements finis hors de SDi qui le contient
                  interface=.true.
                  nptir_p(i)=nptir_p(i)+1
                  cut_p(nptir_p(i),i)=noeud
                   part_p(noeud,0)=part_p(noeud,0)+1
                   pos=part_p(noeud,0)
                   if (pos > NBMAXNSHARE-1) then
                     write(lu,*)  'Error : an interface node belongs to 
     &                     more than NBMAXNSHARE-1 subdomains'
                      call plante2(-1)
                      stop 
                   endif
                   part_p(noeud,pos)=i
                   nbre_ef_loc(knogl(noeud,i))=
     c                  -1*nbre_ef_loc(knogl(noeud,i))
               end if
            end do 
            if (interface .eqv. .true.) then 
               nbre_ef_i=nbre_ef_i+1 ! l'element fini est donc aussi interface
               ef_i(nbre_ef_i)=ef
               ef_ii(nbre_ef_i)=j
            end if
         end do
         
 ! first loop to compute the number of halo to allocate ifapar 
        

!     filling of  ifapar
         nhalo(i)=0
         do j=1,nbre_ef_i       ! on parcours juste les elements finis interfaces pour                             ! determiner des halo
            ef=ef_i(j)
            halo=.false.
            ifaloc(:)=ifabor(ef,:)
            
            

            where (ifaloc .gt. 0) 
               ifaloc=epart(ifaloc)
            end where
            halo=ANY(ifaloc .gt. 0 .and. ifaloc .ne. i)
            if (halo .EQV. .true.) then
               nhalo(i)=nhalo(i)+1
               if (nhalo(i) > NBMAXHALO) then
                  write(lu,*)  'Error : NBMAXHALO too small'
                  call plante2(-1)
                  stop 
               endif
               ifapar(i,1,nhalo(i))=ef_ii(j)
               ifapar(i,2:4,nhalo(i))=ifaloc(:)
               ifapar(i,5:7,nhalo(i))=ifabor(ef_i(j),:)
            end if
         end do 
         
         
                                !     
c       write(lu,*) 'Sous domaine ',i,'nbre points',npoin_p(i),
c     c        'nbre noeud inte',nbre_noeud_interne,'interface',
c     c        nptir_p(i),'nbre front',nptfr_p(i), 'halo',nhalo(i),
c     c        'nbre_efront',nbre_ef_i
       
c        if (.NOT. allocated(nbor_p)) then 
c            allocate(nbor_p(npoin2),stat=err)
c           if (err.ne.0) call ALLOER (lu, 'nbor_p')
c        end if 
      end do
c      deallocate(ifabor)
c      deallocate(nbre_ef)
c     deallocate(nbre_ef_loc)
c      deallocate(ef_i)
c      deallocate(ef_ii)
         
      max_n_neigh=maxval(part_p(:,0))
      if ( max_n_neigh > NBMAXNSHARE-1 ) then 
         write(lu,*) 'SERIOUS WARNING: ' 
         write(lu,*) 
     &        'An interface node belongs to ',
     &        'more than NBMAXNSHARE-1 subdomains'
         write(lu,*) 'Telemac may protest!'
      end if 
      if (max_n_neigh > maxnproc) then
         write (lu,*) 'There is a node which belongs to more than ',
     &        maxnproc,' processors, how come?'
         call plante2(-1)
          stop 
       endif
       if (max_n_neigh < NBMAXNSHARE-1) max_n_neigh = NBMAXNSHARE-1
       
      do i=1,nparts
!-----------------------------------------------------------------------
! the core names for the output bc files according to the number of parts
!
      nameclm = namecli    ! core name length is i_lencli
      nameout = nameinp    ! core name length is i_leninp
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!     work on the boundaries writing the BC files simultaneously...
!     
         nameclm(i_lencli+1:i_lencli+11) = extens(nparts-1,i-1)
!         write (lu,*) 'Working on BC file ',nameclm

         open(nclm,file=nameclm,
     c        status='UNKNOWN',form='FORMATTED')
         rewind(nclm)
!     
! file opened, now work on boundaries 
! -----------------------------------
!
! when the boundary node belongs to this subdomain it will be taken
! j is the running boundary node number
!
!        nptir = 0
         j = 0
!
         do k=1,nptfr
!
! boundary nodes belonging to this partition
!
            if ( knogl(nbor(k),i) /= 0) then
               j = j + 1
!               nbor_p(j) = nbor(k)
               iseg = 0
               xseg = 0.0
               yseg = 0.0
!     
!     if the original (global) boundary leads further into 
!     another partition then iseg is set not equal to zero
!     the next node along the global boundary has iptfr = m
!     (but check the case the circle closes)
!     
               m = kp1bor(k,1)
!
! nbor_p cannot be used, it is not fully filled with data
!
               iso = 0
!     modif JMH ON 10/03/2003 : checking if the adjacent element is not in the
!     sub-domain
               if (epart(nelbor(k)).NE.i) then
!     this was a test : if next boundary point not in the subdomain
!     but it can be in whereas the segment is not.
!     if ( knogl(nbor(m)) == 0 ) then
!                  write(lu,*) 
!     &                 'Global boundary leaves @node (#G,#L): ',
!     &                 nbor(k), knogl(nbor(k),i),
!     &                 ' --> (#G) ', nbor(m)
!     
                  iseg = nbor(m)
                  xseg = f(iseg,1)
                  yseg = f(iseg,2)
c                  nptir = nptir + 1
c                  cut(nptir) = irand_p(knogl(nbor(k)))
                  iso = iso + 1
            endif
!     
            m = kp1bor(k,2)
!     
!     modif JMH ON 10/03/2003 : same as above, but previous segment ,thus m, not k
            if (epart(nelbor(m)).NE.i) then
!     if ( knogl(nbor(m) ) == 0 ) then
!               write(lu,*) 
!     &              'Global boundary enters @node (#G,#L): ',
!     &              nbor(k), knogl(nbor(k),i),
!     &              ' <-- (#G) ', nbor(m)
!     
               iseg = -nbor(m)
               xseg = f(-iseg,1)
               yseg = f(-iseg,2)
               iso = iso + 1
c               nptir = nptir + 1
c               cut(nptir) = irand_p(knogl(nbor(k)))
            endif
!     
!     when both neighbours boundary nodes belong to another partition
!     
            if (iso == 2) then
               iseg = -9999
               iso = 0
               write(lu,*) 'Isolated boundary point', 
     &              nbor(k), knogl(nbor(k),i)
            endif
!     
!            nbor_p(j) = irand_p(knogl(nbor(k)))
!     
!     write a line of the first (classical) part of the boundary file
! concerning the node which has been researched
!
            write (nclm,4000) 
     &           lihbor(k), liubor(k), livbor(k),
     &           hbor(k), ubor(k), vbor(k), 
     &            aubor(k), litbor(k), tbor(k), atbor(k), btbor(k),
!     JMH 16/06/2008: initial line number or colour
     &           nbor(k),check(k), iseg, xseg, yseg, numliq(k)
!     &            nbor(k),    j   , iseg, xseg, yseg, numliq(k)
     
!     19/10/2007 ER+JMH SUR RECOMMANDATION CHARLES MOULINEC
!     MAIS XSEG ET YSEG NE SONT PLUS UTILISES
 4000       format (1x,i2,1x,2(i1,1x),3(f24.12,1x),1x,
     &           f24.12,3x,i1,1x,3(f24.12,1x),1i6,1x,1i6,
     &           1x,i7,1x,2(f27.15,1x),i6)
         endif
!     
      end do
      
      fmt4='(i6)'
      write (nclm,*) nptir_p(i)
       if (max_n_neigh < NBMAXNSHARE-1) max_n_neigh = NBMAXNSHARE-1
       fmt4='(   (i6,1x))'
       write (fmt4(2:4),'(i3)') max_n_neigh+1
  
       ! sorting node numbers to sort(j) so that cut_p(sort(j)) is ordered 
! cut is overwritten now
!
         do j=1,nptir_p(i)
           cut(j)=cut_p(j,i)
         end do
!
! if a node has been already found as min, cut(node) gets 0
!
         do j=1,nptir_p(i)
           idum = npoin2+1  ! largest possible node number + 1
           k=0
 401       continue
           k = k + 1
           if ( cut(k) /= 0 .and. cut_p(k,i) < idum ) then
             sort(j) = k
             idum = cut_p(k,i)
           endif
           if ( k < nptir_p(i) ) then
             goto 401
           else
             cut(sort(j)) = 0
           endif
         end do
!
         do j=1,nptir_p(i)
            tab_tmp=0
            l=0
            do k=1,max_n_neigh
              
               if (part_p(cut_p(sort(j),i),k) .NE. i .AND. 
     c        part_p(cut_p(sort(j),i),k) .NE. 0) then
                  l=l+1
               tab_tmp(l)=part_p(cut_p(sort(j),i),k)
            end if
         end do 
         write(nclm,fmt=fmt4) cut_p(sort(j),i),
     &                  (tab_tmp(k)-1, k=1,max_n_neigh)
         end do
                                 !     
         do j=1,nhalo(i)
            do m=0,2
               if (ifapar(i,2+m,j)>0) then
                  ifapar(i,5+m,j)=gelegl(ifapar(i,5+m,j),
     c                 ifapar(i,2+m,j))
               end if
            end do
         end do
          do j=1,nhalo(i)
           do m=0,2
              if (ifapar(i,2+m,j)>0) then
                 ifapar(i,2+m,j)=ifapar(i,2+m,j)-1
              end if
           end do
        end do
      
!
      WRITE(nclm,'(i9)') nhalo(i)
      DO k=1,nhalo(i)
         WRITE (nclm,'(7(i9,1x))') ifapar(i,:,k) 
      END DO 
     
      close(nclm)
      end do 

      deallocate(ifapar)
      deallocate(part_p)
      deallocate(lihbor)
      deallocate(liubor)
      deallocate(livbor)
      deallocate(hbor)
      deallocate(ubor)
      deallocate(vbor)
      deallocate(aubor)
      deallocate(litbor)
      deallocate(tbor)
      deallocate(atbor)
      deallocate(btbor)
      deallocate(nbor)
      deallocate(numliq)
      deallocate(tab_tmp)
      deallocate(numsol)
      deallocate(check)
      deallocate(gelegl)
      deallocate(cut)
      deallocate(cut_p)
      deallocate(sort)
     

      if (err.ne.0) call ALLOER (lu, 'f_p')
      allocate(ikles_p(max_nelem_p*3),stat=err)
      if(nplan.gt.1) then
         allocate(ikles3d_p(6,max_nelem_p,nplan-1),stat=err)
       endif
      if (err.ne.0) call ALLOER (lu, 'ikles3d_p')
!

      do i=1,nparts
!     ***************************************************************
!     writing geometry files for all parts/processors
!
      nameout(i_leninp+1:i_leninp+11) = extens(nparts-1,i-1)
!      write(lu,*) 'Writing geometry file: ',nameout      
      open(nout,file=nameout,form='unformatted'
     c     ,status='unknown')
      
      rewind(nout)
!     
!     title, the number of variables
!     
      write(nout) title
      write(nout) nvar,0
      do k=1,nvar
         write(nout) variable(k)
      end do
!     
!     10 integers...
! 1.  is the number of recordings in files
! 8.  is the number of boundary points (nptfr_p)
! 9.  is the number of interface points (nptir_p)
! 10. is 0 when no date passed; 1 if a date/time record follows
!
!       ib(7) = nplan   (already done)
        ib(8) = nptfr_p(i)
        ib(9) = nptir_p(i)
        write(nout) (ib(k), k=1,10)
        if (ib(10).eq.1) then 
           write(nout) date(1), date(2), date(3), 
     &                time(1), time(2), time(3)
           
        endif 

        if(nplan.le.1) then
          write(nout) nelem_p(i), npoin_p(i), ndp, ndum
        else
           write(nout) nelem_p(i)*(nplan-1),
     &          npoin_p(i)*nplan, ndp, ndum
        endif
!     
        do j=1,nelem_p(i)
           ef=elelg(j,i)
           do k=1,3
              ikles_p((j-1)*3+k) = knogl(ikles((ef-1)*3+k),i)
           end do
        end do
        if(nplan > 1) then
           do k = 1,nplan-1
                                  do j = 1,nelem_p(i)       
                ikles3d_p(1,j,k) = ikles_p(1+(j-1)*3) + (k-1)*npoin_p(i)
                ikles3d_p(2,j,k) = ikles_p(2+(j-1)*3) + (k-1)*npoin_p(i)
                ikles3d_p(3,j,k) = ikles_p(3+(j-1)*3) + (k-1)*npoin_p(i)
                ikles3d_p(4,j,k) = ikles_p(1+(j-1)*3) +  k   *npoin_p(i)
                ikles3d_p(5,j,k) = ikles_p(2+(j-1)*3) +  k   *npoin_p(i)
                ikles3d_p(6,j,k) = ikles_p(3+(j-1)*3) +  k   *npoin_p(i)
              enddo
           enddo
        endif
!
        if (nplan.eq.0) then
           write(nout) 
     &          ((ikles_p((j-1)*3+k),k=1,3),j=1,nelem_p(i))
        else
           
           
           write(nout)
     &         (((ikles3d_p(l,j,k),l=1,6),j=1,nelem_p(i)),k=1,nplan-1)
        endif
!     
! instead of irand, knolg is written !!!
! i.e. the table processor-local -> processor-global node numbers
!
        if (nplan.eq.0) then
          write(nout) (knolg(j,i), j=1,npoin_p(i))
        else
!     beyond npoin_p(i) : dummy values in knolg, never used
           write(nout) (knolg(j,i), j=1,npoin_p(i)*nplan)
        endif
!

! node coordinates x and y
!
        if (nplan.eq.0) then
          write(nout) (f_p(j,1,i),j=1,npoin_p(i))
          write(nout) (f_p(j,2,i),j=1,npoin_p(i))
        else
             
         
          write(nout) ((f(knolg(j,i)+(l-1)*npoin2,1),j=1,npoin_p(i)), 
     c          l=1,nplan)  
          write(nout) ((f(knolg(j,i)+(l-1)*npoin2,2),j=1,npoin_p(i)), 
     c          l=1,nplan)  
        endif
!
! time stamp (seconds) 
!
        write(nout) times
!
! now the time-dependent variables
!
        do k=3,nvar+2
          if(nplan.eq.0) then
             write(nout) (f_p(j,k,i),j=1,npoin_p(i))
          else
     
             write(nout) ((f(knolg(j,i)+(l-1)*npoin2,k),j=1,npoin_p(i)), 
     c            l=1,nplan) 
   

          
          endif
        end do
!     
        close (nout)
        write(lu,*) 'Finished subdomain ',i

      end do

! CD I HAVE COMMENTED THIS ... AVOIDING MULTIPLES FILES MAKING BUG ON SGI
!
!======================================================================
! writing epart and npart 
!
c$$$      chch='00000'
c$$$      if (nparts<10) then
c$$$        write (chch(5:5),'(i1)') nparts
c$$$        nameepart=nameinp(1:i_leninp) // '.epart.' // chch(5:5)
c$$$        namenpart=nameinp(1:i_leninp) // '.npart.' // chch(5:5)
c$$$      elseif (nparts<100) then
c$$$        write (chch(4:5),'(i2)') nparts
c$$$        nameepart=nameinp(1:i_leninp) // '.epart.' // chch(4:5)
c$$$        namenpart=nameinp(1:i_leninp) // '.npart.' // chch(4:5)
c$$$      elseif (nparts<1000) then
c$$$        write (chch(3:5),'(i3)') nparts
c$$$        nameepart=nameinp(1:i_leninp) // '.epart.' // chch(3:5)
c$$$        namenpart=nameinp(1:i_leninp) // '.npart.' // chch(3:5)
c$$$      elseif (nparts<10000) then
c$$$        write (chch(2:5),'(i4)') nparts
c$$$        nameepart=nameinp(1:i_leninp) // '.epart.' // chch(2:5)
c$$$        namenpart=nameinp(1:i_leninp) // '.npart.' // chch(2:5)
c$$$      else 
c$$$        write (chch(1:5),'(i5)') nparts
c$$$        nameepart=nameinp(1:i_leninp) // '.epart.' // chch(1:5)
c$$$        namenpart=nameinp(1:i_leninp) // '.npart.' // chch(1:5)
c$$$      endif
c$$$!
c$$$      write(lu,*) ' '
c$$$      write(lu,*) '------------------'
c$$$      write(lu,*) ' Partition files  '
c$$$      write(lu,*) '------------------'
c$$$      write(lu,*) ' '
c$$$!
c$$$      open(nepart,file=nameepart,status='UNKNOWN',form='FORMATTED')
c$$$      rewind nepart
c$$$      write(lu,*) 'Element partition file: ', nameepart
c$$$!
c$$$      open(nnpart,file=namenpart,status='UNKNOWN',form='FORMATTED')
c$$$      rewind nnpart
c$$$      write(lu,*) 'Node partition file: ', namenpart
c$$$!
c$$$! output absolutely the same as from partdnmesh (a C-program)
c$$$! that's why 1 substracted and the formats
c$$$!
c$$$      fmt1 = '(i1)'
c$$$      fmt2 = '(i2)'
c$$$      fmt3 = '(i3)'
c$$$!
c$$$      do j=1,nelem2
c$$$         k = epart(j) - 1
c$$$         if (k<10) then
c$$$           write (nepart,fmt=fmt1) k
c$$$         elseif (k<100) then
c$$$           write (nepart,fmt=fmt2) k
c$$$         else 
c$$$           write (nepart,fmt=fmt3) k
c$$$         endif
c$$$      end do
c$$$      close(nepart)
c$$$!
c$$$      do j=1,npoin2
c$$$         k = npart(j) - 1
c$$$         if (k<10) then
c$$$           write (nnpart,fmt=fmt1) k
c$$$         elseif (k<100) then
c$$$           write (nnpart,fmt=fmt2) k
c$$$         else
c$$$           write (nnpart,fmt=fmt3) k
c$$$         endif
c$$$      end do
c$$$      close(nnpart)

!
! //// jaj: la finita commedia for parallel characteristics, bye! 
!----------------------------------------------------------------------
! !jaj #### deal with sections 
!
      IF (nplan/=0) with_sections=.FALSE.
      IF (with_sections) THEN ! presently, for Telemac2D, ev. Sisyphe 

      WRITE(lu,*) 'Dealing with sections'
      OPEN (ninp,FILE=TRIM(namesec),FORM='formatted',STATUS='old') 
      READ (ninp,*) ! comment line
      READ (ninp,*) nsec, ihowsec
      IF (.NOT.ALLOCATED(chain)) ALLOCATE (chain(nsec))
      IF (ihowsec<0) THEN 
        DO isec=1,nsec
          READ (ninp,*) chain(isec)%descr
          READ (ninp,*) chain(isec)%npair(:)
          chain(isec)%xybeg(:)= (/f(chain(isec)%npair(1),1),
     &                            f(chain(isec)%npair(1),2)/)
          chain(isec)%xyend(:)= (/f(chain(isec)%npair(2),1),
     &                            f(chain(isec)%npair(2),2)/)
        END DO 
      ELSE
        DO isec=1,nsec
          READ (ninp,*) chain(isec)%descr
          READ (ninp,*) chain(isec)%xybeg(:), chain(isec)%xyend(:)
          chain(isec)%npair(:)=0
        END DO 
      ENDIF
      CLOSE(ninp) 
 
      ! if terminal points given by coordinates, find nearest nodes first

      WRITE(lu,*) 'npoin:',npoin
      IF (ihowsec>=0) THEN 
        DO isec=1,nsec
          xa=f(1,1) 
          ya=f(1,2)
          dminb = (chain(isec)%xybeg(1)-xa)**2 
     &          + (chain(isec)%xybeg(2)-ya)**2 
          dmine = (chain(isec)%xyend(1)-xa)**2 
     &          + (chain(isec)%xyend(2)-ya)**2 
          chain(isec)%npair(1)=1
          chain(isec)%npair(2)=1
          DO i=2,npoin ! computationally intensive 
            xa=f(i,1)
            ya=f(i,2)
            distb = (chain(isec)%xybeg(1)-xa)**2 
     &            + (chain(isec)%xybeg(2)-ya)**2 
            diste = (chain(isec)%xyend(1)-xa)**2 
     &            + (chain(isec)%xyend(2)-ya)**2 
            IF ( distb < dminb ) THEN 
              chain(isec)%npair(1)=i
              dminb=distb
            ENDIF
            IF ( diste < dmine ) THEN 
              chain(isec)%npair(2)=i
              dmine=diste 
            ENDIF 
          END DO
          WRITE(lu,'(a,3(1x,i9))') 
     &          ' -> section, terminal nodes: ', 
     &          isec, chain(isec)%npair(:)
        END DO  
      ELSE
        DO isec=1,nsec
          WRITE(lu,'(a,1x,i9,4(1x,1pg13.6))') 
     &          ' -> section, terminal coordinates: ', isec, 
     &          chain(isec)%xybeg, chain(isec)%xyend
        END DO 
      ENDIF 

      WRITE(lu,*) 'nsec,ihowsec: ',nsec,ihowsec
      WRITE(lu,*) 'anticipated sections summary:'
      DO isec=1,nsec
        WRITE(lu,*) chain(isec)%descr
        WRITE(lu,*) chain(isec)%xybeg(:), chain(isec)%xyend(:)
        WRITE(lu,*) chain(isec)%npair(:)
      END DO  

! now follow the flusec subroutine in bief to find sections 
! in the global mesh -> fill the field LISTE

      ncp = 2*nsec
      ALLOCATE(liste(nsemax,2),STAT=err) ! workhorse 
      IF (err.NE.0) CALL alloer (lu, 'liste')

      DO isec =1,nsec

        dep = chain(isec)%npair(1) 
        arr = chain(isec)%npair(2)

        pt = dep
        iseg = 0
        dist=(f(dep,1)-f(arr,1))**2+(f(dep,2)-f(arr,2))**2

 1010   CONTINUE ! a jump point 

        DO ielem =1,nelem
          i1 = ikle(ielem,1)
          i2 = ikle(ielem,2)
          i3 = ikle(ielem,3)
          IF (pt.EQ.i1.OR.pt.EQ.i2.OR.pt.EQ.i3) THEN
            dist1 = (f(i1,1)-f(arr,1))**2 + (f(i1,2)-f(arr,2))**2
            dist2 = (f(i2,1)-f(arr,1))**2 + (f(i2,2)-f(arr,2))**2
            dist3 = (f(i3,1)-f(arr,1))**2 + (f(i3,2)-f(arr,2))**2
            IF (dist1.LT.dist) THEN
              dist = dist1
              elbest = ielem
              igbest = i1
              ilbest = 1
              IF(i1.EQ.pt) ilprec = 1
              IF(i2.EQ.pt) ilprec = 2
              IF(i3.EQ.pt) ilprec = 3
            ENDIF
            IF (dist2.LT.dist) THEN 
              dist = dist2
              elbest = ielem
              igbest = i2
              ilbest = 2
              IF(i1.EQ.pt) ilprec = 1
              IF(i2.EQ.pt) ilprec = 2
              IF(i3.EQ.pt) ilprec = 3
            ENDIF
            IF(dist3.LT.dist) THEN
              dist = dist3
              elbest = ielem
              igbest = i3
              ilbest = 3
              IF(i1.EQ.pt) ilprec = 1
              IF(i2.EQ.pt) ilprec = 2
              IF(i3.EQ.pt) ilprec = 3
            ENDIF
          ENDIF

        END DO ! over elements 

        IF (igbest.EQ.pt) THEN
          WRITE(lu,*)'flusec : algorithm failed'
          CALL plante2(-1)
          STOP
        ELSE
          pt = igbest
          iseg = iseg + 1
          IF (iseg.GT.nsemax) THEN
            WRITE(lu,*) 'too many segments in a   '
            WRITE(lu,*) 'section. increase  nsemax'
            CALL plante2(-1)
            STOP
          ENDIF
          liste(iseg,1) = ikle(elbest,ilprec)
          liste(iseg,2) = ikle(elbest,ilbest)
          IF (igbest.NE.arr) GOTO 1010
        ENDIF
        chain(isec)%nseg = iseg
        ALLOCATE (chain(isec)%liste(chain(isec)%nseg,3), STAT=err)
        IF (err/=0) CALL alloer (lu, 'chain(isec)%liste') 
        DO iseg=1,chain(isec)%nseg
          chain(isec)%liste(iseg,1)=liste(iseg,1) 
          chain(isec)%liste(iseg,2)=liste(iseg,2) 
          chain(isec)%liste(iseg,3)=-1 ! initialise... for devel 
        END DO 
      END DO ! over sections 
      DEALLOCATE (liste) 

! now one can indicate the partitions the sections go through
! proceed segment-wise, usinf 2d knolg / knogl

!      DO i=1,npoin
!        WRITE(lu,*) i,knogl(i,:) 
!      END DO 

      ALLOCATE (anpbeg(NBMAXNSHARE), STAT=err)
      IF (err/=0) CALL alloer (lu, 'anpbeg') 
      ALLOCATE (anpend(NBMAXNSHARE), STAT=err) 
      IF (err/=0) CALL alloer (lu, 'anpend') 

      DO isec=1,nsec 
        DO iseg=1,chain(isec)%nseg 

          npbeg=COUNT( knogl(chain(isec)%liste(iseg,1),:)>0 )
          npend=COUNT( knogl(chain(isec)%liste(iseg,2),:)>0 )          

          IF (npbeg>NBMAXNSHARE .OR. npend>NBMAXNSHARE) THEN 
            WRITE(lu,*) 'npbeg or npend: ',npbeg,npend
            WRITE(lu,*) 'are larger than NBMAXNSHARE: ',NBMAXNSHARE
            CALL plante2(-1) 
            STOP
          ENDIF 

          ! the nice and usual case when both segment ends 
          ! belong to one subdomain - only one position in knogl 
          IF ( npbeg==1 .AND. npend==1) THEN  
             im(:) = MAXLOC ( knogl(chain(isec)%liste(iseg,1),:) ) 
             in(:) = MAXLOC ( knogl(chain(isec)%liste(iseg,2),:) )
             IF (im(1)==in(1)) THEN  
               chain(isec)%liste(iseg,3)=im(1) 
             ELSE ! they belong to different subdomains? how come?
               WRITE(lu,*) 'impossible case (1) by sections???'
               CALL plante2(-1)
               STOP
             ENDIF 
          ! at least one of the terminal nodes is on the interface
          ! take the largest common partition number they both belong to
          ELSE 
            IF (npbeg==1 .AND. npend>1) THEN ! the segment's end touches the interface 
              im(:) = MAXLOC ( knogl(chain(isec)%liste(iseg,1),:) )
              IF ( knogl(chain(isec)%liste(iseg,2),im(1))>0 ) THEN  
                chain(isec)%liste(iseg,3) = im(1) 
              ELSE 
                WRITE(lu,*) 'impossible case (2) by sections???'
                CALL plante2(-1)
                STOP
              ENDIF 
            ELSE IF (npbeg>1 .AND. npend==1) THEN ! the segment's beg. touches the interface
              in(:) = MAXLOC ( knogl(chain(isec)%liste(iseg,2),:) )
              IF ( knogl(chain(isec)%liste(iseg,1),in(1))>0 ) THEN  
                chain(isec)%liste(iseg,3) = in(1) 
              ELSE 
                WRITE(lu,*) 'impossible case (3) by sections???'
                CALL plante2(-1)
                STOP
              ENDIF 
            ELSE ! i.e. (npbeg>1 .AND. npend>1) - lies just on the interface or "a shortcut" 
              anpbeg=0
              anpend=0 
              i=0 
              DO n=1,nparts
                IF ( knogl(chain(isec)%liste(iseg,1),n)>0 ) THEN
                  i=i+1
                  anpbeg(i)=n 
                ENDIF  
              END DO 
              IF (i/=npbeg) WRITE(lu,*) 'oh! i/=npbeg'
              i=0 
              DO n=1,nparts
                IF ( knogl(chain(isec)%liste(iseg,2),n)>0 ) THEN
                  i=i+1
                  anpend(i)=n 
                ENDIF  
              END DO 
              IF (i/=npend) WRITE(lu,*) 'oh! i/=npend'

              WRITE(lu,*) 'anpbeg: ',anpbeg
              WRITE(lu,*) 'anpend: ',anpend

              found=.FALSE.
              DO i=npbeg,1,-1
                DO j=npend,1,-1
                  IF (anpbeg(i)==anpend(j)) THEN 
                     chain(isec)%liste(iseg,3) = anpbeg(i)
                    found=.TRUE.
                    EXIT
                  ENDIF 
                END DO 
                IF (found) EXIT 
              END DO 
              IF (.NOT.found) THEN 
                WRITE(lu,*) 'by section with nodes: ',
     &            chain(isec)%liste(iseg,1),chain(isec)%liste(iseg,2)
                WRITE(lu,*) 'impossible case (4) by sections???'
                CALL plante2(-1)
                STOP
              ENDIF 

            ENDIF
          ENDIF 

        END DO 
      END DO 

      DEALLOCATE (anpbeg,anpend) 

! devel printout 

!      WRITE(lu,*) 'summary of section chains partitioning'
!      DO isec=1,nsec
!        WRITE(lu,*) 'isec, nseg: ',isec,chain(isec)%nseg
!        WRITE(lu,*) 'descr: ',TRIM(chain(isec)%descr) 
!        DO iseg=1,chain(isec)%nseg
!          WRITE(lu,*) chain(isec)%liste(iseg,1), 
!     &                chain(isec)%liste(iseg,2),
!     &                chain(isec)%liste(iseg,3)
!        END DO 
!      END DO 

! write files 

      DO n=1,nparts
        nameout=TRIM(namesec)//extens(nparts-1,n-1)

        WRITE(lu,*) 'writing: ', TRIM(nameout) 

        OPEN (nout,FILE=TRIM(nameout),FORM='formatted',STATUS='unknown')
        REWIND(nout) 
        WRITE(nout,*) '# sections partitioned for ',extens(nparts-1,n-1)
        WRITE(nout,*) nsec, 1
        DO isec=1,nsec
          WRITE(nout,*) TRIM(chain(isec)%descr)
          i=COUNT(chain(isec)%liste(:,3)==n) 
          WRITE(nout,*) i
          DO iseg=1,chain(isec)%nseg
            IF (chain(isec)%liste(iseg,3)==n) THEN 
              WRITE(nout,*) 
     &          knogl(chain(isec)%liste(iseg,1),n),
     &          knogl(chain(isec)%liste(iseg,2),n)
            ENDIF
          END DO 
        END DO
        CLOSE(nout) 
      END DO 

      WRITE(lu,*) 'Finished dealing with sections'
      ENDIF ! nplan==0
!
!----------------------------------------------------------------------
!
!     note by j-m hervouet : deallocate causes errors on HP
!     (possible remaining bug ?)
!     note by jaj: deallocate(HP) ,^)
!
       deallocate (ikle) ! #### moved from far above 
       deallocate(npart)
       deallocate(epart)
       deallocate(npoin_p)
       deallocate(nelem_p)
       deallocate(nptfr_p)
       deallocate(nptir_p)
!
       deallocate(ikles)
      if(nplan.gt.1) then
         deallocate(ikles3d)
         deallocate(ikles3d_p)
      endif
      deallocate(ikles_p)
      deallocate(irand)
      deallocate(f)
      deallocate(f_p)

      deallocate(knolg)
      deallocate(knogl)
      deallocate(elelg)
      deallocate(kp1bor)

!
!----------------------------------------------------------------------      
!
 299  if (timecount) then 
        call system_clock (count=temps, count_rate=parsec)
        tfin = temps
        write(lu,*) 'Overall timing: ',
     &    (1.0*(tfin-tdeb))/(1.0*parsec),' seconds'
        write(lu,*) ' '
      endif
      write(lu,*) '+---- PARTEL: normal termination ----+'
      write(lu,*) ' '
!
      go to 999

 300  write(lu,*) 'Error by reading. '
      call plante2(-1)

 999  stop  
      end program PARTEL


      subroutine ALLOER (n, chfile)
      implicit none
      integer, intent(in) :: n
      character*(*), intent(in) :: chfile
      write(n,*) 'error by allocation of ',chfile
      call plante2(-1)
      stop
      end subroutine ALLOER
      
      subroutine ALLOER2(n,chfile)
      implicit none
      integer, intent(in) :: n
      character*(*), intent(in) :: chfile
      write(n,*)trim(chfile)
      call plante2(-1)
      stop
      end subroutine ALLOER2


      subroutine PLANTE2 (ival)
      implicit none
      integer, intent(in) :: ival
      integer icode      
      if (ival < 0) then      ! this indicates a controlled error
        icode = 1 
      else if (ival==0) then  ! this indicates a program failure
        icode = -1
      else                    ! this indicates a normal stop
        icode = 0
      endif 
      !!! write(*,*) 'Returning exit code: ', icode
      CALL EXIT(ICODE)
      stop    ! which is usually equivalent to call EXIT(0)
      end subroutine PLANTE2
C                       *********************************
                        CHARACTER(LEN=11) FUNCTION EXTENS
C                       *********************************
C
     *(N,IPID)
C
C***********************************************************************
C  PARA       VERSION 4.0         08/01/97        J-M HERVOUET (LNH)
C
C***********************************************************************
C
C      FONCTIONS: EXTENSION DES FICHIERS SUR CHAQUE PROCESSEUR.
C      ==========
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________|
C |     N          | -->| NOMBRE DE PROCESSEURS MOINS UN = NCSIZE-1
C |     IPID       | -->| NUMERO DU PROCESSEUR
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR :
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER IPID,N
C
C-----------------------------------------------------------------------
C
      IF(N.GT.0) THEN
C
        EXTENS='00000-00000'
C
        IF(N.LT.10) THEN
          WRITE(EXTENS(05:05),'(I1)') N
        ELSEIF(N.LT.100) THEN
          WRITE(EXTENS(04:05),'(I2)') N
        ELSEIF(N.LT.1000) THEN
          WRITE(EXTENS(03:05),'(I3)') N
        ELSEIF(N.LT.10000) THEN
          WRITE(EXTENS(02:05),'(I4)') N
        ELSE
          WRITE(EXTENS(01:05),'(I5)') N
        ENDIF
C
        IF(IPID.LT.10) THEN
          WRITE(EXTENS(11:11),'(I1)') IPID
        ELSEIF(IPID.LT.100) THEN
          WRITE(EXTENS(10:11),'(I2)') IPID
        ELSEIF(IPID.LT.1000) THEN
          WRITE(EXTENS(09:11),'(I3)') IPID
        ELSEIF(IPID.LT.10000) THEN
          WRITE(EXTENS(08:11),'(I4)') IPID
        ELSE
          WRITE(EXTENS(07:11),'(I5)') IPID
        ENDIF
C
      ELSE
C
        EXTENS='       '
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
C                       ************************
                        subroutine front2_partel
C                       ************************
C
     *(nfrliq,nfrsol,debliq,finliq,debsol,finsol,lihbor,liubor,
     * x,y,nbor,kp1bor,dejavu,npoin,nptfr,klog,listin,numliq,numsol,
     * nptfrmax)
C
C***********************************************************************
C bief version 5.5           04/05/04    j-m hervouet  01 30 87 80 18
C***********************************************************************
C
c  fonction  : reperage, numerotation des frontieres liquides et solides
c
c-----------------------------------------------------------------------
c                             arguments
c .________________.____.______________________________________________
c |      nom       |mode|                   role
c |________________|____|______________________________________________
c |   nfrliq       |<-- | nombre de frontieres liquides
c |   nfrsol       |<-- | nombre de frontieres solides
c |   debliq       |<-- | debuts des frontieres liquides
c |   finliq       |<-- | fins des frontieres liquides
c |   debsol       |<-- | debuts des frontieres solides
c |   finsol       |<-- | fins des frontieres solides
c |   lihbor       | -->| conditions aux limites sur h
c |   x , y        | -->| coordonnees du maillage.
c |   nbor         | -->| numeros globaux des points de bord
c |   kp1bor       | -->| numeros des extremites des segments de bord
c |                |    | dans la numerotation des points de bord
c |   dejavu       | -- | tableau de travail
c |   npoin        | -->| nombre de points du maillage
c |   nptfr        | -->| nombre de points frontiere
c |   klog         | -->| lihbor(k)=klog : frontiere solide
c |   listin       | -->| impressions sur listing (ou non)
c |________________|____|______________________________________________
c mode : -->(donnee non modifiee), <--(resultat), <-->(donnee modifiee)
c-----------------------------------------------------------------------
c
c  precautions d'emploi : les frontieres solides sont reperees par le
c                         fait que lihbor(k) = klog pour un point de
c                         bord de numero k.
c                         un segment compris entre un point liquide et
c                         un point solide est solide.
c
c***********************************************************************
c
      implicit none
      integer lng,lu
      common/info/lng,lu
c
c+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c
      integer, intent(in) :: npoin,nptfr,klog,nptfrmax
      integer, intent(out) :: nfrliq,nfrsol
c                                    *=maxfro (300 dans telemac-2d)
      integer, intent(out) :: debliq(*),finliq(*),debsol(*),finsol(*)
ccccccMOULINEC Begin
      integer , intent(in) :: lihbor(nptfrmax),liubor(nptfrmax)
c      integer , intent(in) :: lihbor(nptfr),liubor(nptfr)
ccccccMOULINEC End
      real, intent(in) :: x(npoin) , y(npoin)
ccccccMOULINEC Begin
      integer, intent(in) :: nbor(2*nptfrmax),kp1bor(nptfr)
c      integer, intent(in) :: nbor(nptfr),kp1bor(nptfr)
ccccccMOULINEC End
      integer, intent(out) :: dejavu(nptfr)
      logical, intent(in) :: listin
ccccccMOULINEC Begin
      integer, intent(out) :: numliq(nptfrmax)
      integer, intent(out) :: numsol(nptfrmax)
c      integer, intent(out) :: numliq(nptfr)
c      integer, intent(out) :: numsol(nptfr)
ccccccMOULINEC End
c
c+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c 
      integer k,kprev,idep,sol1,liq1,l1,l2,l3,nile
c
      logical solf,liqf,sold,liqd
c
      real minns,maxns,eps,ymin,ns
c
      intrinsic abs
c
c-----------------------------------------------------------------------
c
c  initialisations
c
c  dejavu : marque d'un 1 les points deja traites
c  nile   : nombre d'iles
c
      do 10 k=1,nptfr
        dejavu(k) = 0
        numliq(k) = 0
        numsol(k) = 0
10    continue
c
      nile = 0
      idep = 1
      nfrliq = 0
      nfrsol = 0
c
c-----------------------------------------------------------------------
c
c  on reviendra a l'etiquette 20 s'il y a au moins une ile
c
20    continue
c
c  recherche du point le plus sud-ouest (il peut y en avoir plusieurs)
c
      minns = x(nbor(idep)) + y(nbor(idep))
      maxns = minns
      ymin  = y(nbor(idep))
c
      do 30 k = 1 , nptfr
      if(dejavu(k).eq.0) then
        ns = x(nbor(k)) + y(nbor(k))
        if(ns.lt.minns) then
         idep = k
         minns = ns
         ymin = y(nbor(k))
        endif
        if(ns.gt.maxns) maxns = ns
      endif
30    continue
c
      eps = (maxns-minns) * 1.d-4
c
c  choix du point le plus sud parmi les candidats sud-ouest
c
      do 40 k = 1 , nptfr
      if(dejavu(k).eq.0) then
        ns = x(nbor(k)) + y(nbor(k))
        if(abs(minns-ns).lt.eps) then
          if(y(nbor(k)).lt.ymin) then
           idep = k
           ymin = y(nbor(k))
          endif
        endif
      endif
40    continue
c
c-----------------------------------------------------------------------
c
c  numerotation et reperage des frontieres du contour commencant
c  au point idep.
c
c  sold = .true. : la frontiere au depart de idep est solide
c  liqd = .true. : la frontiere au depart de idep est liquide
c  solf = .true. : la frontiere au retour a idep est solide
c  liqf = .true. : la frontiere au retour a idep est liquide
c  liq1 : numero de la premiere frontiere liquide du contour
c  sol1 : numero de la premiere frontiere solide du contour
c
      k = idep
c
      sol1 = 0
      liq1 = 0
      liqf = .false.
      solf = .false.
c
c nature du premier segment
c
c     loi de dominance du solide sur le liquide
      if(lihbor(k).eq.klog.or.lihbor(kp1bor(k)).eq.klog) then
c       le premier segment est solide
        nfrsol = nfrsol + 1
        sol1 = nfrsol
        sold = .true.
        liqd = .false.
      else
c       le premier segment est liquide
        nfrliq = nfrliq + 1
        liq1 = nfrliq
        liqd = .true.
        sold = .false.
      endif
c
      dejavu(k) = 1
      kprev = k
      k = kp1bor(k)
c
50    continue
c
c recherche des points de transition a partir du point suivant ideb
c
c on cherche aussi les cas de points isoles pour detecter les erreurs
c dans les donnees.
c
      l1 = lihbor(kprev)
      l2 = lihbor(k)
      l3 = lihbor(kp1bor(k))
c
      if(l1.eq.klog.and.l2.ne.klog.and.l3.ne.klog) then
c     transition solide-liquide au point k
        nfrliq = nfrliq + 1
        finsol(nfrsol) = k
        debliq(nfrliq) = k
        liqf = .true.
        solf = .false.
      elseif(l1.ne.klog.and.l2.ne.klog.and.l3.eq.klog) then
c     transition liquide-solide au point k
        nfrsol = nfrsol + 1
        finliq(nfrliq) = k
        debsol(nfrsol) = k
        liqf = .false.
        solf = .true.
      elseif(l1.ne.klog.and.l2.ne.klog.and.l3.ne.klog) then
c     recherche des transitions liquide-liquide au point k
        if(l2.ne.l3.or.liubor(k).ne.liubor(kp1bor(k))) then
          finliq(nfrliq) = k
          nfrliq = nfrliq + 1
          debliq(nfrliq) = kp1bor(k)
        endif
      elseif(l1.eq.klog.and.l2.ne.klog.and.l3.eq.klog) then
c     erreur dans les donnees
        if(lng.eq.1) write(lu,102) k
        if(lng.eq.2) write(lu,103) k
        call plante2(-1)
        stop
      elseif(l1.ne.klog.and.l2.eq.klog.and.l3.ne.klog) then
c     erreur dans les donnees
        if(lng.eq.1) write(lu,104) k
        if(lng.eq.2) write(lu,105) k
        call plante2(-1)
        stop
      endif
c
      dejavu(k) = 1
      kprev = k
      k = kp1bor(k)
      if(k.ne.idep) go to 50
c
c  cas d'un changement de frontiere au point de depart idep
c
      if(solf) then
c       la derniere frontiere etait solide
        if(sold) then
c         la premiere frontiere etait solide
          debsol(sol1) = debsol(nfrsol)
          nfrsol = nfrsol - 1
        elseif(liqd) then
c         la premiere frontiere etait liquide
          debliq(liq1) = idep
          finsol(nfrsol) = idep
        endif
c
      elseif(liqf) then
c       la derniere frontiere du contour etait liquide
        if(liqd) then
c         la premiere frontiere du contour etait liquide
          debliq(liq1) = debliq(nfrliq)
          nfrliq = nfrliq - 1
        elseif(sold) then
c         la premiere frontiere du contour etait solide
          debsol(sol1) = idep
          finliq(nfrliq) = idep
        endif
c
      else
c     cas ou tout le contour est du meme type
        if(sol1.ne.0) then
          debsol(sol1) = idep
          finsol(sol1) = idep
        elseif(liq1.ne.0) then
          debliq(liq1) = idep
          finliq(liq1) = idep
        else
          if(listin.and.lng.eq.1) then
           write(lu,'(1x,a)') 'cas impossible dans front2'
          endif
          if(listin.and.lng.eq.2) then
           write(lu,'(1x,a)') 'impossible case in front2'
          endif
          call plante2(-1)
          stop
        endif
      endif
c
c-----------------------------------------------------------------------
c
c  on regarde s'il reste des contours :
c
      do 60 k = 1 , nptfr
        if(dejavu(k).eq.0) then
          idep = k
          nile = nile + 1
          go to 20
        endif
60    continue
c
c-----------------------------------------------------------------------
c
      do 79 k=1,nptfr
        numliq(k)=0
79    continue
c
c  impression des resultats et calcul de numliq
c
      if(nile.ne.0.and.listin.and.lng.eq.1) write(lu,69) nile
      if(nile.ne.0.and.listin.and.lng.eq.2) write(lu,169) nile
c
      if(nfrliq.ne.0) then
        if(listin.and.lng.eq.1) write(lu,70) nfrliq
        if(listin.and.lng.eq.2) write(lu,170) nfrliq

        do 80 k = 1, nfrliq
c
c  marquage des numeros des frontieres liquides
c
          l1=debliq(k)
          numliq(l1)=k
707       l1=kp1bor(l1)
          numliq(l1)=k
          if(l1.ne.finliq(k)) go to 707
c
c  fin du marquage
c
          if(listin.and.lng.eq.1) write(lu,90)
     *                            k,debliq(k),nbor(debliq(k)),
     *                            x(nbor(debliq(k))),y(nbor(debliq(k))),
     *                            finliq(k),nbor(finliq(k)),
     *                            x(nbor(finliq(k))),y(nbor(finliq(k)))
          if(listin.and.lng.eq.2) write(lu,190)
     *                            k,debliq(k),nbor(debliq(k)),
     *                            x(nbor(debliq(k))),y(nbor(debliq(k))),
     *                            finliq(k),nbor(finliq(k)),
     *                            x(nbor(finliq(k))),y(nbor(finliq(k)))
80      continue
      endif
c
      if(nfrsol.ne.0) then
        if(listin.and.lng.eq.1) write(lu,100) nfrsol
        if(listin.and.lng.eq.2) write(lu,101) nfrsol

        do 110 k = 1, nfrsol
!
!  marking solid boundaries (why not?)
!  they get next boundary numbers 
!
          l1=debsol(k)
          numsol(l1)=k+nfrliq
708       l1=kp1bor(l1)
          numsol(l1)=k+nfrliq
          if(l1.ne.finsol(k)) go to 708
!
!  end od fmarking
!
          if(listin.and.lng.eq.1) write(lu,90)
     *                            k,debsol(k),nbor(debsol(k)),
     *                            x(nbor(debsol(k))),y(nbor(debsol(k))),
     *                            finsol(k),nbor(finsol(k)),
     *                            x(nbor(finsol(k))),y(nbor(finsol(k)))
          if(listin.and.lng.eq.2) write(lu,190)
     *                            k,debsol(k),nbor(debsol(k)),
     *                            x(nbor(debsol(k))),y(nbor(debsol(k))),
     *                            finsol(k),nbor(finsol(k)),
     *                            x(nbor(finsol(k))),y(nbor(finsol(k)))
110     continue
      endif
c
c-----------------------------------------------------------------------
c
c  formats
c
69    format(/,1x,'il y a ',1i3,' ile(s) dans le domaine')
169   format(/,1x,'there is ',1i3,' island(s) in the domain')
70    format(/,1x,'il y a ',1i3,' frontiere(s) liquide(s) :')
170   format(/,1x,'there is ',1i3,' liquid boundaries:')
100   format(/,1x,'il y a ',1i3,' frontiere(s) solide(s) :')
101   format(/,1x,'there is ',1i3,' solid boundaries:')
102   format(/,1x,'front2 : erreur au point de bord ',1i5,
     *       /,1x,'         point liquide entre deux points solides')
103   format(/,1x,'front2 : error at boundary point ',1i5,
     *       /,1x,'         liquid point between two solid points')
104   format(/,1x,'front2 : erreur au point de bord ',1i5,
     *       /,1x,'         point solide entre deux points liquides')
105   format(/,1x,'front2 : error at boundary point ',1i5,
     *       /,1x,'         solid point between two liquid points')
90    format(/,1x,'frontiere ',1i3,' : ',/,1x,
     *            ' debut au point de bord ',1i4,
     *            ' , de numero global ',1i6,/,1x,
     *            ' et de coordonnees : ',g16.7,3x,g16.7,
     *       /,1x,' fin au point de bord ',1i4,
     *            ' , de numero global ',1i6,/,1x,
     *            ' et de coordonnees : ',g16.7,3x,g16.7)
190   format(/,1x,'boundary ',1i3,' : ',/,1x,
     *            ' begins at boundary point: ',1i4,
     *            ' , with global number: ',1i6,/,1x,
     *            ' and coordinates: ',g16.7,3x,g16.7,
     *       /,1x,' ends at boundary point: ',1i4,
     *            ' , with global number: ',1i6,/,1x,
     *            ' and coordinates: ',g16.7,3x,g16.7)
c
c-----------------------------------------------------------------------
c
      if(nile.gt.300.or.nfrsol.gt.300.or.nfrliq.gt.300) then
        if(lng.eq.1) then
          write(lu,*) 'front2 : depassement de tableaux'
          write(lu,*) '         augmenter maxfro dans le code appelant' 
          write(lu,*) '         a la valeur ',max(nile,nfrsol,nfrliq)             
        endif
        if(lng.eq.2) then
          write(lu,*) 'front2: size of arrays exceeded'
          write(lu,*) '        increase maxfro in the calling program' 
          write(lu,*) '        up to the value ',max(nile,nfrsol,nfrliq)             
        endif
        call plante2(-1)
        stop
      endif
c
c-----------------------------------------------------------------------
c
      return
      end 
C                       ************************
                        SUBROUTINE VOISIN_PARTEL
C                       ************************
C
     *(IFABOR,NELEM,NELMAX,IELM,IKLE,SIZIKL,NPOIN,IADR,NVOIS)
C
C***********************************************************************
C BIEF VERSION 5.9         16/06/2008    J-M HERVOUET (LNHE) 30 87 80 18
C
C***********************************************************************
C
C    FONCTION : CONSTRUCTION DU TABLEAU IFABOR, OU IFABOR(IELEM,IFACE)
C               EST LE NUMERO GLOBAL DU VOISIN DE LA FACE IFACE DE
C               L'ELEMENT IELEM SI CE VOISIN EXISTE ET 0 SI LA FACE EST
C               SUR LA FRONTIERE DU DOMAINE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    IFABOR      |<-- | TABLEAU DES VOISINS DES FACES.
C |    NELEM       | -->| NOMBRE D'ELEMENTS DANS LE MAILLAGE.
C |    NELMAX      | -->| NOMBRE MAXIMUM D'ELEMENTS DANS LE MAILLAGE.
C |                |    | (CAS DES MAILLAGES ADAPTATIFS)
C |    IELM        | -->| 11: TRIANGLES
C |                |    | 21: QUADRILATERES
C |    IKLE        | -->| NUMEROS GLOBAUX DES POINTS DE CHAQUE ELEMENT
C |    NPOIN       | -->| NOMBRE TOTAL DE POINTS DU DOMAINE
C |________________|____|_______________________________________________
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT DANS TELEMAC 2D : PREDAT
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER IR1,IR2,IR3,IR4,SIZIKL
C
      INTEGER NELEM,NELMAX,IELM,IDIMAT,NPOIN,I,J,ERR,NDP
      INTEGER NFACE,KEL,I1,I2,IMAX,IFACE,IELEM,M1,M2,IV,IELEM2,IFACE2
      INTEGER IFABOR(NELMAX,*),IKLE(SIZIKL,*),NVOIS(NPOIN),IADR(NPOIN)
C
      INTEGER SOMFAC(2,4,2)
      DATA SOMFAC / 1,2 , 2,3 , 3,1 , 0,0   ,  1,2 , 2,3 , 3,4 , 4,1 /
C
C  TABLEAUX DE TRAVAIL ALLOUES DYNAMIQUEMENT
C
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAT1,MAT2,MAT3
C
C-----------------------------------------------------------------------
C
      IF(IELM.EQ.21) THEN
C       QUADRILATERES
        NFACE = 4
        NDP = 4
        KEL = 2
      ELSEIF(IELM.EQ.11.OR.IELM.EQ.41) THEN
C       TRIANGLES
        NFACE = 3
        NDP = 3
        KEL = 1
      ELSE
        IF(LNG.EQ.1) WRITE(LU,98) IELM
        IF(LNG.EQ.2) WRITE(LU,99) IELM
98      FORMAT(1X,'VOISIN: IELM=',1I6,' TYPE D''ELEMENT NON PREVU')
99      FORMAT(1X,'VOISIN: IELM=',1I6,' TYPE OF ELEMENT NOT AVAILABLE')
        CALL PLANTE2(-1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C     IDIMAT EST UNE MAJORATION DE LA SOMME DES NOMBRES DE VOISINS DE
C     TOUS LES POINTS.
C    
      IDIMAT = 2*NDP*NELEM
C
      ALLOCATE(MAT1(IDIMAT),STAT=ERR)
      ALLOCATE(MAT2(IDIMAT),STAT=ERR)
      ALLOCATE(MAT3(IDIMAT),STAT=ERR)
C
      IF(ERR.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,1000) ERR
        IF(LNG.EQ.2) WRITE(LU,2000) ERR
1000    FORMAT(1X,'VOISIN : ERREUR A L''ALLOCATION DE MEMOIRE : ',/,1X,
     *            'CODE D''ERREUR : ',1I6)
2000    FORMAT(1X,'VOISIN: ERROR DURING ALLOCATION OF MEMORY: ',/,1X,
     *            'ERROR CODE: ',1I6)
      ENDIF
C
C-----------------------------------------------------------------------
C
C  CALCUL DU TABLEAU NVOIS POUR CHAQUE POINT
C  ATTENTION : NVOIS N'EST PAS LE NOMBRE DE VOISINS MAIS PERMET DE
C              RESERVER ASSEZ DE PLACE DANS LES TABLEAUX MAT1,2,3.
C
      DO 10 I=1,NPOIN
        NVOIS(I) = 0
10    CONTINUE
C
      DO 20 IFACE = 1,NFACE
        DO 30 IELEM=1,NELEM
          I1 = IKLE( IELEM , SOMFAC(1,IFACE,KEL) )
          I2 = IKLE( IELEM , SOMFAC(2,IFACE,KEL) )
          NVOIS(I1) = NVOIS(I1) + 1
          NVOIS(I2) = NVOIS(I2) + 1
30      CONTINUE
20    CONTINUE
C
C-----------------------------------------------------------------------
C
C  CALCUL DES ADRESSES DE CHAQUE POINT DANS UNE STRUCTURE DE TYPE
C  MATRICE COMPACTE
C
      IADR(1) = 1
      DO 50 I= 2,NPOIN
        IADR(I) = IADR(I-1) + NVOIS(I-1)
50    CONTINUE
C
      IMAX = IADR(NPOIN) + NVOIS(NPOIN) - 1
      IF(IMAX.GT.IDIMAT) THEN
        IF(LNG.EQ.1) WRITE(LU,51) IDIMAT,IMAX
        IF(LNG.EQ.2) WRITE(LU,52) IDIMAT,IMAX
51      FORMAT(1X,'VOISIN: TAILLE DE MAT1,2,3 (',1I6,') INSUFFISANTE',/,
     *         1X,'IL FAUT AU MOINS : ',1I6)
52      FORMAT(1X,'VOISIN: SIZE OF MAT1,2,3 (',1I6,') TOO SHORT',/,
     *         1X,'MINIMUM SIZE: ',1I6)
        call plante2(-1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C  INITIALISATION A ZERO DE LA MATRICE COMPACTE
C
      DO 53 I=1,IMAX
        MAT1(I) = 0
53    CONTINUE
C
C-----------------------------------------------------------------------
C
C  BOUCLE SUR LES FACES DE CHAQUE ELEMENT :
C
      DO 60 IFACE = 1 , NFACE
      DO 70 IELEM = 1 , NELEM
C
      IFABOR(IELEM,IFACE) = -1
C
C        NUMEROS GLOBAUX DES POINTS DE LA FACE :
C
         I1 = IKLE( IELEM , SOMFAC(1,IFACE,KEL) )
         I2 = IKLE( IELEM , SOMFAC(2,IFACE,KEL) )
C
C        NUMEROS GLOBAUX ORDONNES :
C
         M1 = MIN0(I1,I2)
         M2 = MAX0(I1,I2)
C
         DO 80 IV = 1,NVOIS(M1)
C
           IF(MAT1(IADR(M1)+IV-1).EQ.0) THEN
              MAT1(IADR(M1)+IV-1)=M2
              MAT2(IADR(M1)+IV-1)=IELEM
              MAT3(IADR(M1)+IV-1)=IFACE
              GO TO 81
           ELSEIF(MAT1(IADR(M1)+IV-1).EQ.M2) THEN
              IELEM2 = MAT2(IADR(M1)+IV-1)
              IFACE2 = MAT3(IADR(M1)+IV-1)
              IFABOR(IELEM,IFACE) = IELEM2
              IFABOR(IELEM2,IFACE2) = IELEM
              GO TO 81
           ENDIF
C
80       CONTINUE
C
         IF(LNG.EQ.1) WRITE(LU,82)
         IF(LNG.EQ.2) WRITE(LU,83)
82       FORMAT(1X,'VOISIN : ERREUR DANS LE MAILLAGE       ',/,1X,
     *             '         PEUT-ETRE DES POINTS CONFONDUS')
83       FORMAT(1X,'VOISIN : ERROR IN THE MESH             ',/,1X,
     *             '         MAYBE SUPERIMPOSED POINTS     ')
         CALL PLANTE2(-1)
         STOP
C
81       CONTINUE
C
70    CONTINUE
60    CONTINUE
C
C-----------------------------------------------------------------------
C
      DEALLOCATE(MAT1)
      DEALLOCATE(MAT2)
      DEALLOCATE(MAT3)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
C                       ***********************
                        SUBROUTINE ELEBD_PARTEL
C                       ***********************
C
C
C      TAKEN FROM BIEF AND ADAPTED : USE BIEF REMOVED
C                                    CALL PLANTE REMOVED
C                                    ALL ACTIONS UNDER IF(NCSIZE.GT.1) REMOVED
C
C
     *(NELBOR,NULONE,KP1BOR,IFABOR,NBOR,IKLE,SIZIKL,
     * IKLBOR,NELEM,NELMAX,NPTFRMAX,
     * NPOIN,NPTFR,IELM,LIHBOR,KLOG,IFANUM,OPTASS,ISEG,T1,T2,T3,
     *     NPOIN_TOT)
C
C***********************************************************************
C BIEF VERSION 5.3           23/08/99    J-M HERVOUET (LNH) 30 87 80 18
C COPYRIGHT 1999
C***********************************************************************
C
C    PRISMES DECOUPES EN TETRAEDRES
C
C    FONCTION : 1) CONSTRUCTION DES TABLEAUX NELBOR ET NULONE
C               2) CONSTRUCTION DU TABLEAU KP1BOR
C               3) DISTINCTION DANS LE TABLEAU IFABOR ENTRE
C                  LES FACES DE BORD SOLIDES ET LES FACES LIQUIDES
C               4) COMPLEMENT DE NBOR. 
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |    NELBOR      |<-- | NUMERO DE L'ELEMENT ADJACENT AU KIEME SEGMENT|
C |    NULONE      |<-- | NUMERO LOCAL D'UN POINT DE BORD DANS         |
C |                |    | L'ELEMENT ADJACENT DONNE PAR NELBOR          |
C |    KP1BOR      |<-- | NUMERO DU POINT SUIVANT LE POINT DE BORD K.  |
C |    IFABOR      | -->| TABLEAU DES VOISINS DES FACES.
C |    NBOR        | -->| NUMERO GLOBAL DU POINT DE BORD K.
C |    IKLE        | -->| NUMEROS GLOBAUX DES POINTS DE CHAQUE ELEMENT.
C |    NELEM       | -->| NOMBRE TOTAL D'ELEMENTS DANS LE MAILLAGE.
C |    T1,2,3      | -->| TABLEAUX DE TRAVAIL ENTIERS.
C |    NPOIN       | -->| NOMBRE TOTAL DE POINTS DU DOMAINE.
C |    NPTFR       | -->| NOMBRE DE POINTS FRONTIERES.
C |    IELM        | -->| TYPE D'ELEMENT.
C |                |    | 11 : TRIANGLES.
C |                |    | 21 : QUADRILATERES.
C |    LIHBOR      | -->| TYPES DE CONDITIONS AUX LIMITES SUR H
C |    KLOG        | -->| CONVENTION POUR LA CONDITION LIMITE DE PAROI
C |    MXPTVS      | -->| NOMBRE MAXIMUM DE VOISINS D'UN POINT
C |    MXELVS      | -->| NOMBRE MAXIMUM D'ELEMENTS AUTOUR D'UN POINT
C |________________|____|______________________________________________|
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : 
C
C SOUS-PROGRAMME APPELE : NEANT
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: KLOG,NELMAX,NELEM,SIZIKL
      INTEGER, INTENT(IN)    :: NPOIN,NPTFR,IELM,OPTASS,NPTFRMAX
      INTEGER, INTENT(OUT)   :: NELBOR(NPTFR),NULONE(NPTFR,2)
      INTEGER, INTENT(OUT)   :: KP1BOR(NPTFR,2)
      INTEGER, INTENT(INOUT) :: NBOR(2*NPTFRMAX)
      INTEGER, INTENT(INOUT) :: IFABOR(NELMAX,*)
      INTEGER, INTENT(IN)    :: IKLE(SIZIKL,*)
      INTEGER, INTENT(IN)    :: LIHBOR(NPTFRMAX)
      INTEGER, INTENT(OUT)   :: IKLBOR(NPTFR,2)
      INTEGER, INTENT(INOUT) :: IFANUM(NELMAX,*)
      INTEGER, INTENT(IN)    :: ISEG(NPTFR)
      INTEGER, INTENT(OUT)   :: T1(NPOIN),T2(NPOIN),T3(NPOIN) 
      INTEGER, INTENT(IN)    :: NPOIN_TOT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,NFACE,NPT,KEL,IPOIN
      INTEGER K,IFACE,I1,I2,N1,N2,IPT,IEL,I,K1,K2
C
      INTEGER SOMFAC(2,4,2)
C
      DATA SOMFAC / 1,2 , 2,3 , 3,1 , 0,0   ,  1,2 , 2,3 , 3,4 , 4,1 /
C
C-----------------------------------------------------------------------
C
      IF(IELM.EQ.11.OR.IELM.EQ.41.OR.IELM.EQ.51) THEN
C       TRIANGLES
        NFACE = 3
        NPT = 3
        KEL = 1
      ELSE
        IF(LNG.EQ.1) WRITE(LU,900) IELM
        IF(LNG.EQ.2) WRITE(LU,901) IELM
900     FORMAT(1X,'ELEBD : IELM=',1I6,' TYPE D''ELEMENT INCONNU')
901     FORMAT(1X,'ELEBD: IELM=',1I6,' UNKNOWN TYPE OF ELEMENT')
        call plante2(-1)
        STOP
      ENDIF
C
C  INITIALISATION DE T1,2,3 A ZERO
C
      DO IPOIN=1,NPOIN
        T1(IPOIN) = 0
        T2(IPOIN) = 0
        T3(IPOIN) = 0
      ENDDO
C
C  ON STOCKE K DANS TRAV(*,3) A L'ADRESSE NBOR(K)
C  CE QUI PERMET DE PASSER DE NUMERO GLOBAL A NUMERO DE BORD
C
      DO K = 1, NPTFR
         T3(NBOR(K)) = K
      ENDDO
C
C  BOUCLE SUR TOUTES LES FACES DE TOUS LES ELEMENTS :
C
      DO 20 IFACE = 1 , NFACE
      DO 10 IELEM = 1 , NELEM
C
      IF(IFABOR(IELEM,IFACE).EQ.-1) THEN
C
C      C'EST UNE VRAIE FACE DE BORD (LES FACES INTERNES EN PARALLELISME
C                                    SONT MARQUEES AVEC DES -2).
C      NUMEROS GLOBAUX DES POINTS DE LA FACE :
C
       I1 = IKLE( IELEM , SOMFAC(1,IFACE,KEL) )
       I2 = IKLE( IELEM , SOMFAC(2,IFACE,KEL) )
C
C      ON STOCKE DANS T1 ET T2 A L'ADRESSE I1 : I2 ET IELEM
C
       T1(I1) = I2
       T2(I1) = IELEM
C
C      UNE FACE LIQUIDE EST RECONNUE AVEC LA CONDITION LIMITE SUR H
C
C      07/02/03 IF(NPTFR...  COURTESY OLIVER GOETHEL, HANNOVER UNIVERSITY
       IF(NPTFR.GT.0) THEN
       IF(LIHBOR(T3(I1)).NE.KLOG.AND.LIHBOR(T3(I2)).NE.KLOG) THEN
C        FACE LIQUIDE : IFABOR=0  FACE SOLIDE : IFABOR=-1
         IFABOR(IELEM,IFACE)=0
       ENDIF
       ENDIF
C
      ENDIF
C
10    CONTINUE
20    CONTINUE
C
C  BOUCLE SUR TOUS LES POINTS:
C
C     07/02/03 IF(NPTFR...  CORRECTION BY OLIVER GOETHELS, HANNOVER
      IF(NPTFR.GT.0) THEN
      DO I = 1 , NPOIN
         IF(T1(I).NE.0) THEN
C          POINT SUIVANT
           KP1BOR(T3(I),1)=T3(T1(I))
C          POINT PRECEDENT
           KP1BOR(T3(T1(I)),2)=T3(I)
           NELBOR(T3(I))=T2(I)
         ENDIF
      ENDDO
      ENDIF
C
C CALCUL DU TABLEAU NULONE
C
      DO 50 K1=1,NPTFR
C
      K2=KP1BOR(K1,1)
      IEL = NELBOR(K1)
      N1  = NBOR(K1)
      N2  = NBOR(K2)
C
      I1 = 0
      I2 = 0
C
      DO 60 IPT=1,NPT
C
        IF(IKLE(IEL,IPT).EQ.N1) THEN
          NULONE(K1,1) = IPT
          I1 = 1
        ENDIF
        IF(IKLE(IEL,IPT).EQ.N2) THEN
          NULONE(K1,2) = IPT
          I2 = 1
        ENDIF
C
60    CONTINUE
C
      IF(I1.EQ.0.OR.I2.EQ.0) THEN
        IF(LNG.EQ.1) WRITE(LU,810) IEL
        IF(LNG.EQ.2) WRITE(LU,811) IEL
810     FORMAT(1X,'ELEBD: ERREUR DE NUMEROTATION DANS L''ELEMENT:',I6,/,
     *         1X,'       CAUSE POSSIBLE :                       '   ,/,
     *         1X,'       LE FICHIER DES CONDITIONS AUX LIMITES NE'  ,/,
     *         1X,'       CORRESPOND PAS AU FICHIER DE GEOMETRIE  ')
811     FORMAT(1X,'ELEBD: ERROR OF NUMBERING IN THE ELEMENT:',I6,
     *         1X,'       POSSIBLE REASON:                       '   ,/,
     *         1X,'       THE BOUNDARY CONDITION FILE IS NOT      '  ,/,
     *         1X,'       RELEVANT TO THE GEOMETRY FILE           ')
        call plante2(-1)
        STOP
      ENDIF
C
50    CONTINUE
C
C  COMPLEMENT DU TABLEAU NBOR
C
C -----
C FD : BAW only ?
C
C      OPEN(UNIT=89,FILE='front_glob.dat')
C      write(89,*) NPOIN_TOT
C      write(89,*) NPTFR
C      DO K=1,NPTFR
C         write(89,*) NBOR(K)
C      END DO 
C -----
C
      DO 80 K=1,NPTFR
C
        NBOR(K+NPTFR) = NBOR(KP1BOR(K,1))
C
        IKLBOR(K,1) = K
        IKLBOR(K,2) = KP1BOR(K,1)
C
80    CONTINUE

      

C
C -----
C      DO K=1,NPTFR
C         write(89,*) KP1BOR(K,1)
C      END DO
C        DO K=1,NPTFR
C         write(89,*) KP1BOR(K,2)
C      END DO 
C ------
C
C ------
C      call flush(89)
C      CLOSE(89)
C ------
C-----------------------------------------------------------------------
C
      RETURN
      END
      
c                       ******************
!
! 29/01/2009:

c                       ******************
                        subroutine pares3d
c                       ******************

     *(nameinp,li)
c
c**********************************************************************
c  12/11/2009 Christophe Denis SINETICS/I23 
c  New version to decrease the pares3d computing time by improving 
c
c   - the tetra-tria connection
c   - the postprocessing 
c
c Comments on this new version ->  CD 
c *********************************************************************
c***********************************************************************
c partel version 5.6        08/06/06   o.boiteau/f.decung(sinetics/lnhe)
c partel version 5.8        02/07/07   f.decung(lnhe)
c version de developpement pour prise en compte pb decoupage
c f.decung/O.boiteau (janv 2008)
c copyright 2006
c***********************************************************************
c
c    constructions des fichiers pour alimenter le flot de donnees
c    parallele lors d'un calcul estel3d parallele en ecoulement
c 
c
c-----------------------------------------------------------------------
c                             arguments
c .________________.____.______________________________________________.
c |      nom       |mode|                   role                       |
c |________________|____|______________________________________________|
c |    nameinp     | -->| nom du fichier de geometrie estel3d
c |    li          | -->| unite logique d'ecriture pour monitoring 
c |________________|____|______________________________________________|
c  mode: -->(donnee non modifiee),<--(resultat),<-->(donnee modifiee)
c
c-----------------------------------------------------------------------
c
c appele par :  partel
c
c sous-programme appele :
c        alloer, alloer2 (gestion msgs)
c        METIS_PartMeshDual (from METIS library)
c***********************************************************************
C
      implicit none
!     integer, parameter :: maxnproc = 1000  ! max partition number [000..999]
      integer, parameter :: maxnproc = 100000 ! max partition number [00000..99999]
      integer, parameter :: maxlenhard = 250 ! hard max file name length
      integer, parameter :: maxlensoft = 144 ! soft max file name length
!
      integer lng,lu
      common/info/lng,lu
! parametres d'appel de la routine
      character(len=maxlenhard), intent(in) :: nameinp
      integer,                   intent(in) :: li
! variables locales
      character(len=maxlenhard) :: namelog,nameinp2,namelog2
      integer :: ninp=10,nlog=11,ninp2=12,nlog2=13
      integer :: nparts,i_s,i_sp,i,i_len,i_leninp,ierr,j,k,compt,
     &           n,numtet,numtri,numtrig,i_lenlog,l,ni,nf,nt,ibid,idd,
     &           compt1,compt2,compt3,nbtriidd,ibid1,m,ni1,nf1,color1,
     &           color2,pr1,pr2,iddbis,iddterce,nbtetj,iddnt,nit,nft,mt,
     &           numtrib,numtetb,ibidc,nbretouche
      logical :: is,isbis,isterce,linter
      character(len=300) :: texterror  ! texte msg d'erreur
      character(len=8)   :: str8       ! texte msg d'erreur
      character(len=300) :: str26      ! texte msg d'erreur
      character(len=80)  :: titre      ! mesh title in the file
      character(len=2)   :: moins1     ! "-1"
      character(len=4)   :: blanc      ! white space

      ! Addition JP renaud 15/02/2007
      character(len=200) :: line       ! One line, 200 characters maxaddch
      integer            :: pos        ! Position of a character in the line
      integer            :: ios        ! Status integer
      ! End addition JP Renaud
      character(len=72) :: theformat
      
      character(len=80), allocatable :: logfamily(:)  ! log informations
      integer            :: nsec       ! type of the section read
      integer, parameter :: nsec1=151  ! mesh title section id 
      integer, parameter :: nsec2=2411 ! nodes coordinates section id
      integer, parameter :: nsec3=2412 ! connectivity section id
      integer, parameter :: nsec4=2435 ! pour clore proprement la lecture
                                       ! du unv dans estel3d      
      logical            :: read_sec1  ! flag for reading section 1
      logical            :: read_sec2  ! flag for reading section 2
      logical            :: read_sec3  ! flag for reading section 3
      integer            :: nelemtotal ! total number of unv elements
      integer            :: npoint     ! total number of nodes 
      integer            :: nbfamily   ! total number of family
      integer            :: nelin      ! total number of inner triangles
      integer            :: size_flux  !  total number of inner surfaces
      integer, dimension(:), allocatable :: vectnb  ! vecteur aux pour nachb
!
      double precision, allocatable :: x1(:),y1(:),z1(:) ! coord nodes
      integer,          allocatable :: ncolor(:) ! nodes' colour
      integer,          allocatable :: ecolor(:) ! elements' colour
      integer            :: elem       ! type of the element 
      integer            :: ikle1,ikle2,ikle3,ikle4,ikleb   ! nodes
      integer, dimension(:), allocatable :: iklestet ! connectivite en
                   ! renumerotation global de la bief pour les tetraedres
      integer, dimension(:), allocatable :: iklestri ! connectivite en
                   ! renumerotation global de la bief pour les triangles
      integer, dimension(:,:), allocatable :: iklestrin ! connectivite en
                   ! renumerotation global de la bief pour les triangles
      integer, dimension(:,:), allocatable :: iklein ! copie ajustee de iklestrin
      integer, dimension(:,:), allocatable :: typelem ! type d'elt
      integer            :: nbtet,nbtri  ! nbre de tetra, triangle bord
      integer, dimension(:), allocatable :: tettri, tettri2 ! jointure
                                               !  tetra/triangle de bord
      integer            ::  edgecut ! var. auxiliaire pour METIS
      integer, dimension(:), allocatable :: epart ! numero de partition
                                                  ! par element
      integer, dimension(:), allocatable :: npart ! numero de partition
                                                  ! par noeuds
      integer, dimension(:), allocatable :: convtri,convtet ! convertisseur 
         ! numero local tria/tetra numero global; Inverse de typelem(:,2)
      integer            ::  tdeb,tfin,temps,parsec  ! Runtime
      character(len=11), external :: extens ! extension des noms fichier
      integer, dimension(:), allocatable :: npointsd, nelemsd ! nbre
                            ! de points et d'elements par sous-domaine
      integer, dimension(:), allocatable :: npointisd  ! nbre
                            ! de points d'interface par sous-domaine
                   ! Vecteurs lies aux connectivitees nodales inverses
      integer, dimension(:), allocatable :: nodes1,nodes2,nodes3,nodes4
      integer, dimension(:), allocatable :: nodes1t,nodes2t,nodes3t
      integer, dimension(:), allocatable :: triunv ! buffer pour ecrire
                 ! dans les .unv, d'abord les tetras puis les tria
! Pour traitement des Dirichlets confondus avec l'interface
      integer  :: nbcolor ! nbre de couleur de mailles externes
      integer, dimension(:), allocatable :: priority
      integer, dimension(:), allocatable :: ncolor2
! Pour traitement des Dirichlets sur les noeuds de tetra
      logical, dimension(:,:), allocatable :: tetcolor
      logical, dimension(:), allocatable :: deja_trouve
! Indispensable pour parallelisme TELEMAC
      integer, dimension(:), allocatable :: knolg 
      integer, dimension(:,:), allocatable :: nachb
      logical :: nachblog
!     MAXIMUM GEOMETRICAL MULTIPLICITY OF A NODE (variable aussi
!     presente dans la BIEF, ne pas changer l'une sans l'autre)
      INTEGER, PARAMETER :: NBMAXNSHARE =  10
! Cette variable est liee a la precedente et dimensionne differents
! vecteurs
! Note Size of NACHB will be here 2 more than in BIEF, but the extra 2 are
! local work arrays
      integer :: NBSDOMVOIS = NBMAXNSHARE + 2
!
      integer, parameter :: max_size_flux = 100
! number of inner surface (same as size_flux at the end)
      integer, dimension(max_size_flux) :: size_fluxin
! vecteur pour profiling
      integer  temps_sc(20)
!
!F.D
      INTEGER, DIMENSION(:  ), ALLOCATABLE  :: TEMPO,GLOB_2_LOC
      INTEGER, DIMENSION(:,:), ALLOCATABLE  :: IKLES,IKLE,IFABOR
      INTEGER, DIMENSION(:,:), ALLOCATABLE  :: NULONE,IKLBOR
      INTEGER                               :: N1,N2,N3,IKL,IFACE
      INTEGER                               :: NSOLS,NSOLS_OLD
      INTEGER                               :: IFACEBIS,IELEMBIS
      INTEGER                               :: IELEM,IPTFR,IELEB,IPOIN
      LOGICAL, DIMENSION(:), ALLOCATABLE    :: FACE_CHECK
      INTEGER, PARAMETER                    :: NCOL = 256
      INTEGER, DIMENSION(NCOL  )            :: COLOR_PRIO
      INTEGER                               :: PRIO_NEW,NPTFR
      INTEGER, DIMENSION(:), ALLOCATABLE    :: NBOR2,NBOR
      INTEGER, DIMENSION(:), ALLOCATABLE    :: NELBOR,IPOBO
CD******************************************************    ADDED BY Christophe Denis
      INTEGER, DIMENSION(:), ALLOCATABLE     :: nelem_p
C     size nparts, nelem_p(i) is the number of finite elements assigned to subdomain i
      INTEGER, DIMENSION(:), ALLOCATABLE     :: npoin_p
C     size nparts, npoin_p(i) is the number of nodes  assigned to subdomain i
      INTEGER :: node
c     one node ...
      INTEGER ::  pos_node 
c     position of one one node
      INTEGER :: max_nelem_p
c     maximum number of finite elements assigned among subdomains
      INTEGER :: max_npoin_p
c     maximum number of nodes assigned among subdomains
      INTEGER :: max_tria
c     maximum number of triangle sharing a node 
      INTEGER :: the_tri
c     one triangle 
      INTEGER :: jj
c     index counter
      INTEGER, DIMENSION(:), ALLOCATABLE :: number_tria
c     maximum number of triangle sharing a same node  
      INTEGER, DIMENSION(:,:), ALLOCATABLE  :: elegl
c     size max_nelem_p,nparts, elegl(j,i) is the global number of local finite element j in subdomain i
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: nodegl
c     size max_npoin_p,nparts, nodegl(j,i) is the global number of local node j in subdomain i
      INTEGER, DIMENSION(:), ALLOCATABLE :: nodelg
c     size npoint, nodelg(i)=j, j is the local number of global node i on one subdomain
      INTEGER,  DIMENSION(:,:), ALLOCATABLE :: tri_ref
c     size npoint*max_tria
CD********************************************************     
      INTEGER SOMFAC(3,4)
      DATA SOMFAC / 1,2,3 , 4,1,2 , 2,3,4 , 3,4,1  /


!-----------------------------------------------------------------------
! 1. Preambule
!---------------
      call system_clock(count=temps_sc(1),count_rate=parsec)
      allocate (vectnb(NBSDOMVOIS-3))
      write(lu,*)' '
      write(lu,*)'+-------------------------------------------------+'
      write(lu,*)'  Partel: telemac estel3d partitioner'
      write(lu,*)'+-------------------------------------------------+'
      write(lu,*)' Reading unv and log files'
! names of the input files:
      if (nameinp.eq.' ') then
        goto 149 
      else
        write(lu,89)nameinp
      endif
      inquire (file=nameinp,exist=is)
      if (.not.is) goto 140
      do
        read(li,'(a)')namelog
        if (namelog.eq.' ') then
          goto 150 
        else
          write(lu,90)namelog
          exit
        endif
      end do  
      inquire(file=namelog,exist=is)
      if (.not.is) goto 141   
      do 
        read(li,*)nparts
        if ( (nparts > maxnproc) .or. (nparts < 2) ) then
          write(lu,
     &    '('' Number of partitions must be in [2 -'',i6,'']'')') 
     &      maxnproc
        else
          write(lu,91)nparts
          exit
        endif 
      enddo
      
      
! find the input files core name length
      i_s  = len(nameinp)
      i_sp = i_s + 1
      do i=1,i_s
        if (nameinp(i_sp-i:i_sp-i) .ne. ' ') exit
      enddo
      i_len=i_sp - i
      i_leninp = i_len
      if (i_leninp > maxlensoft) goto 144
!
      i_s  = len(namelog)
      i_sp = i_s + 1
      do i=1,i_s
        if (namelog(i_sp-i:i_sp-i) .ne. ' ') exit
      enddo
      i_len=i_sp - i
      i_lenlog = i_len
      if (i_lenlog > maxlensoft) goto 145
!
      open(ninp,file=nameinp,status='old',form='formatted',err=131)
      rewind(ninp)
      open(nlog,file=namelog,status='old',form='formatted',err=130)
      rewind(nlog)
   
!----------------------------------------------------------------------
! 2a. Lecture du fichier .log
!---------------
      read(nlog,51,err=110,end=120)npoint
      read(nlog,52,err=110,end=120)nelemtotal
      read(nlog,53,err=110,end=120)nbfamily
      nbfamily=nbfamily+1            ! pour titre du bloc
      allocate(logfamily(nbfamily),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' logfamily')
      do i=1,nbfamily
        read(nlog,50,err=111,end=120)logfamily(i) 
      enddo
      nbcolor=0

!      read(nlog,531,err=110,end=120)nbcolor

      read(unit=nlog, fmt='(a200)', iostat=ios) line
      if (ios .ne. 0) then
         !         '!----------------------------------!'
         texterror='! Problem with the number of color !'
         call alloer2(lu,texterror)
         call plante2(-1)
      endif
      pos = index(line,':') + 1
      read(unit=line(pos:), fmt=*, iostat=ios) nbcolor
      if (ios .ne. 0) then
         !         '!-------------------------------!'
         texterror='! Problem with the number color !'
      endif

      allocate(priority(nbcolor),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' priority')
      write(lu,92) npoint
      write(lu,93) nelemtotal
      write(lu,94) nbcolor
      if (nbcolor.eq.0) then
        write(lu,*) 'Vous avez oublie de remplir le fichier log...'
        call plante2(-1)
        stop
      endif

      ! Modification JP Renaud 15/02/2007
      ! Some text has been added before the liost of priorities.
      ! Read a 200 character line, find the ':' and then
      ! read the values after the ':'
      read(unit=nlog, fmt='(a200)', iostat=ios) line
      if (ios .ne. 0) then
        !         '!------------------------------------------!'
        texterror='! Problem with the priority of color nodes !'
        call alloer2(lu,texterror)
        call plante2(-1)
      endif
      pos = index(line,':') + 1
      read(unit=line(pos:), fmt=*, iostat=ios) (priority(j),j=1,nbcolor)
      if (ios .ne. 0) then
        !         '!------------------------------------------!'
        texterror='! Problem with the priority of color nodes !'
        call alloer2(lu,texterror)
        call plante2(-1)
      endif
      ! End modification JP Renaud
      write(lu,*) (priority(j),j=1,nbcolor)
      close(nlog)
!
! 2b. Allocations memoires associees
!--------------- 

CD    ****************************** ALLOCATION MEMORY ADDED BY CD
      allocate(nelem_p(nparts),stat=ierr)
      if (ierr.ne.0) call alloer(lu,'nelem_p')
      allocate(npoin_p(nparts),stat=ierr)
      if (ierr.ne.0) call alloer(lu,'npoin_p')
      allocate(nodelg(npoint),stat=ierr)
      if (ierr.ne.0) call alloer(lu,'nodelg')
CD    *******************************      
      
      allocate(x1(npoint),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' x1')
      allocate(y1(npoint),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' y1')
      allocate(z1(npoint),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' z1')
      allocate(ncolor(npoint),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' ncolor')
      allocate(ncolor2(npoint),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' ncolor2')
      allocate(ecolor(nelemtotal),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' ecolor')
      allocate(iklestet(4*nelemtotal),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' iklestet')
      allocate(iklestri(3*nelemtotal),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' iklestri')
      allocate(iklestrin(nelemtotal,4),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' iklestrin')
      allocate(typelem(nelemtotal,2),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' typelem')
      allocate(convtri(nelemtotal),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' convtri')
      allocate(convtet(nelemtotal),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' convtet')
      allocate(npointsd(nparts),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' npointsd')
      allocate(nelemsd(nparts),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' nelemsd')
      allocate(npointisd(nparts),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' npointisd')

!F.D
      allocate(nbor2(npoint),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' nbor2')
      allocate(tempo(2*npoint),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' tempo')
      allocate(face_check(nbfamily),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' face_check')      
      allocate(glob_2_loc(npoint),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' glob_2_loc')
      allocate(ikles(nelemtotal,4),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' ikles')


      READ_SEC1 = .TRUE.
      READ_SEC2 = .TRUE.
      READ_SEC3 = .TRUE.

      DO WHILE ( READ_SEC1 .OR. READ_SEC2 .OR. READ_SEC3 )

          MOINS1 = '  '
          BLANC  = '1111'
          DO WHILE (MOINS1/='-1' .OR. BLANC/='    ')
              READ(NINP,2000, ERR=1100, END=1200) BLANC, MOINS1
          END DO

 2000     FORMAT(A4,A2)

          NSEC = -1

          DO WHILE (NSEC == -1)
              READ(NINP,*, ERR=1100, END=1200) NSEC
          END DO

 2100     FORMAT(I10)

          SELECT CASE (NSEC)

          CASE ( NSEC1 )

              READ_SEC1 = .FALSE.

              READ(NINP,25,ERR=1100, END=1200) TITRE

 25           FORMAT(A80)

          CASE ( NSEC2 )

              READ_SEC2 = .FALSE.
              NCOLOR(:) = -1
              TEMPO(:)  =  0

              DO IELEM = 1, NPOINT
                 READ(NINP,*,ERR=1100,END=1200) N, N1, N2, NCOLOR(IELEM)
          READ(NINP,*,ERR=1100,END=1200) X1(IELEM), Y1(IELEM), Z1(IELEM)
                 TEMPO(N) = IELEM
              END DO

          CASE (NSEC3 )

             READ_SEC3 = .FALSE.
                         
             NBTET         = 0  ! Number of tetra elements to 0
             NBTRI         = 0  ! Number of border elements to 0
             NPTFR         = 0  ! Number of border nodes to 0.
             NELIN         = 0  ! Number of inner surfaces to 0.
             size_flux     = 0  ! Number of user surfaces to 0.
             NBOR2(:)      = 0  ! local to global numbering
             GLOB_2_LOC(:) = 0  ! global to local numbering

!OB'S STUFF
             ECOLOR(:)    = -1
             IKLESTET(:)  = -1
             IKLESTRI(:)  = -1
             TYPELEM(:,:) = -1
             CONVTRI(:)   = -1
             CONVTET(:)   = -1
!
             iklestrin(:,:) = -1
                       
             FACE_CHECK(:) = .FALSE.
             !
             COLOR_PRIO(:)  = 0
             size_fluxin(:) = 0
             !
             DO k = 1, NBCOLOR
                COLOR_PRIO(PRIORITY(K)) = K
             END DO

             DO IELEM = 1, NELEMTOTAL

            READ(NINP,*,ERR=1100,END=1200) NSEC,ELEM,N1,N2,NSOLS,N3
            
                IF (NSEC == -1) EXIT

                SELECT CASE ( ELEM )

                CASE ( 111 )

                   NBTET        = NBTET + 1

                   ECOLOR(IELEM) = NSOLS

               READ(NINP,*, ERR=1100, END=1200) IKLE1,IKLE2,IKLE3,IKLE4

                   IKLES(IELEM, 1) = TEMPO(IKLE1)
                   IKLES(IELEM, 2) = TEMPO(IKLE2)
                   IKLES(IELEM, 3) = TEMPO(IKLE3)
                   IKLES(IELEM, 4) = TEMPO(IKLE4)

!OB'S STUFF
                   N=4*(NBTET-1)+1                       
                   IKLESTET(N)=IKLE1    ! vecteur de connectivite
                   IKLESTET(N+1)=IKLE2
                   IKLESTET(N+2)=IKLE3
                   IKLESTET(N+3)=IKLE4
                   TYPELEM(IELEM,1)=ELEM    ! pour typer les elts
                   TYPELEM(IELEM,2)=NBTET   ! pour conversion num elt> num tetra
                   CONVTET(NBTET)=IELEM     ! l'inverse
                   
                CASE ( 91 )
             
                   IF (NSOLS.GT.0) THEN

                      IF ( NSOLS > NCOL ) THEN
                         WRITE(LU,*) 'Color id pour surfaces externes ', 
     &                        ' trop grand. La limite est : ',NCOL
                      END IF

                      PRIO_NEW = color_prio(nsols)

                      IF ( PRIO_NEW .EQ. 0 ) THEN
                         WRITE(LU,*) ' Numero de face non declare',
     &                        'dans le tableau utilisateur LOGFAMILY ',      
     &                        'voir le fichier des parametres '
                         CALL PLANTE2(1)
                      END IF
                     
                      FACE_CHECK(PRIO_NEW) = .TRUE.

                      NBTRI = NBTRI + 1

                      ECOLOR(IELEM) = NSOLS

                     READ(NINP,*, ERR=1100, END=1200) IKLE1, IKLE2,IKLE3

                      IKLES(IELEM, 1) = TEMPO(IKLE1)
                      IKLES(IELEM, 2) = TEMPO(IKLE2)
                      IKLES(IELEM, 3) = TEMPO(IKLE3)

!OB'S STUFF
                      N=3*(NBTRI-1)+1
                      IKLESTRI(N)=IKLE1
                      IKLESTRI(N+1)=IKLE2
                      IKLESTRI(N+2)=IKLE3
                      TYPELEM(IELEM,1)=ELEM    ! Idem que pour tetra
                      TYPELEM(IELEM,2)=NBTRI
                      CONVTRI(NBTRI)=IELEM

                      DO J=1,3
                                               
                         IKL = IKLES(IELEM,J)
                      
                         IPTFR = GLOB_2_LOC(IKL)

                         IF ( IPTFR .EQ. 0 ) THEN
                            
                            NPTFR           = NPTFR+1
                            NBOR2(NPTFR)    = IKL
                            GLOB_2_LOC(IKL) = NPTFR        
                            IPTFR           = NPTFR

                         END IF                                                        
                         
                    ENDDO  ! Loop over the nodes of the element

                 else if (nsols.lt.0) then
                    !
                    ! User-defined surface for fluxes computation
                    !                      
                    ! NELIN is the counter for the internal elements.
                    ! Actually, we are reading the next internal element.

                    ! nsols_old is used for saving use of a new variable
                    nsols_old = -nsols
                    !
                    ! prio_new is used for saving use of a new variable
                    prio_new = size_fluxin(nsols_old)
                    !
                    if (prio_new.eq.0) then 
                       size_flux = size_flux + 1
                       size_fluxin(nsols_old) = 1
                    endif
                    !
                    nelin = nelin + 1
                    !
                    READ(NINP,*, ERR=1100, END=1200) IKLE1, IKLE2,IKLE3
                    !
                         iklestrin(nelin,1) = nsols
                         iklestrin(nelin,2) = tempo(ikle1)
                         iklestrin(nelin,3) = tempo(ikle2)
                         iklestrin(nelin,4) = tempo(ikle3)
                    !
                 ELSE           ! This is an inner surface, just read the line.

                    READ(NINP,*, ERR=1100, END=1200) IKLE1, IKLE2,IKLE3
                    
                 END IF
                                  
              CASE DEFAULT      ! This is an unknown element.
                 
                 WRITE(LU,*) 'Element inconnu dans le maillage'
                 
              END SELECT        ! The type of the mesh element
              
           END DO               ! Loop over elements to read.

           DO k=1,NBCOLOR
              IF ( .NOT. FACE_CHECK(K)) THEN
                 WRITE(LU,*) ' La couleur de face ',LOGFAMILY(K),
     &                ' n''apparait pas dans le maillage.'
              END IF
           END DO
           
!-----------------------------------------------------------------------

      END SELECT                ! Type of the Section
      
      END DO                    ! WHILE LOOP over sections to read

!------------------------------------------------------- Fin Version F.D

! Correction du nombre d'elements total car celui dans le .log est
! comporte des elements non pris en compte dans une etude ESTEL
      nelemtotal=nbtet+nbtri

      call system_clock(count=temps_sc(2),count_rate=parsec)
      write(lu,*)' TEMPS DE LECTURE FICHIERS LOG & UNV',
     &           (1.0*(temps_sc(2)-temps_sc(1)))/(1.0*parsec),' seconds'
!----------------------------------------------------------------------
! 3a. Construction de tettri/tettri2: correspondance tetra > tria
!---------------

      allocate(nelbor(nbtri),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' nelbor')
      allocate(nulone(nbtri,3),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' nulone')
      allocate(iklbor(nbtri,3),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' iklbor')
      allocate(ikle(nbtet,4),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' ikle')
      allocate(ifabor(nbtet,4),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' ifabor')
!
!------------------------------------------------------- Debut Version O. Boiteau
! Recherche de la correspondance tetraedre <--> triangles de bords
!      allocate(tettri(4*nbtet),stat=ierr)
!      if (ierr.ne.0) call alloer(lu,' tettri')
!      allocate(tettri2(nbtet),stat=ierr)
!      if (ierr.ne.0) call alloer(lu,' tettri2')
!      tettri(:)=-1
!      tettri2(:)=0
!      do i=1,nbtri
!        n=3*(i-1)+1
!        ikle1=iklestri(n)
!        ikle2=iklestri(n+1)
!        ikle3=iklestri(n+2)
!        do j=1,nbtet
!          compt=0
!          do k=1,4
!            ikleb=iklestet(4*(j-1)+k)
!            if (ikleb.eq.ikle1) compt=compt+1
!            if (ikleb.eq.ikle2) compt=compt+10
!            if (ikleb.eq.ikle3) compt=compt+100 
!          enddo ! boucle sur les noeuds du tetraedre j
!          if (compt.eq.111) then   ! tetraedre j associe au triangle i
!            ni=tettri2(j)
!            if (ni==4) then   ! tetra lie a plus de 4 tria, exit             
!              goto 153
!            else              ! prochain emplacement de libre
!              m=4*(j-1)+ni+1  ! dans tettri
!              tettri2(j)=ni+1
!              tettri(m)=i     ! en numerotation locale
!            endif
!            exit
!          endif
!          if (j.eq.nbtet) goto 143  ! erreur car triangle solitaire
!        enddo  ! sur la boucle en tetraedres
!      enddo ! sur la boucle en triangles de bords
!------------------------------------------------------- Fin Version O. Boiteau

!------------------------------------------------------- Debut Version F.D
      allocate(tettri(4*nbtet),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' tettri')
      allocate(tettri2(nbtet),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' tettri2')
!
      tettri (:) =-1
      tettri2(:) =0
!      
      DO IELEM = 1, NBTET
         DO I = 1,4
            IKLE(IELEM,I ) = IKLES (IELEM, I)
         END DO
      END DO
!
      deallocate(ikles)
!
      allocate(iklein(nelin,4),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' iklein')
!
      do ielem = 1, nelin
         do i = 1,4
            iklein(ielem,i ) = iklestrin (ielem, i)
         end do
      end do
!      
      deallocate(iklestrin)
!
      WRITE(LU,*) 'FIN DE LA COPIE DE LA CONNECTIVITE INITIALE'
!      
      allocate(nbor(nptfr),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' nbor')
!
      do ielem = 1, nptfr
            nbor(ielem) = nbor2(ielem)
      end do
!
      deallocate(nbor2)
!
      WRITE(LU,*) 'VOISIN31'

      CALL VOISIN31_PARTEL (IFABOR, NBTET, NBTET,
     &              IKLE,NBTET,NPOINT,NBOR,NPTFR)

      WRITE(LU,*) 'FIN DE VOISIN31'
      
      allocate(ipobo(npoint),stat=ierr)
      if (ierr.ne.0) call alloer(lu,'ipobo')

      CALL ELEBD31_PARTEL( NELBOR, NULONE, IKLBOR,    
     &              IFABOR, NBOR, IKLE,         
     &              NBTET, NBTRI, NBTET, NPOINT, 
     &              NPTFR,IPOBO)

      deallocate(ipobo)

      WRITE(LU,*) 'FIN DE ELEBD31'
      allocate(number_tria(npoint),stat=ierr)
      if (ierr.ne.0) call alloer(lu,'number_tria')
      number_tria = 0
!
      max_tria=0
      DO J = 1, NBTRI
         K = 3*(J-1)+1  
         ikle1 = IKLESTRI(K)
         ikle2 = IKLESTRI(K+1)
         ikle3 = IKLESTRI(K+2)
         the_tri=ikle1
         if (ikle2 < the_tri) the_tri=ikle2 
         if (ikle3< the_tri)  the_tri=ikle3
         number_tria(the_tri)=number_tria(the_tri)+1
      END DO
      max_tria=maxval(number_tria)
!    
      deallocate(number_tria)
!
      allocate(tri_ref(npoint,0:max_tria),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' tri_ref')
       TRI_REF=0 
       DO J = 1, NBTRI
         K = 3*(J-1)+1  
         ikle1 = IKLESTRI(K)
         ikle2 = IKLESTRI(K+1)
         ikle3 = IKLESTRI(K+2)
         the_tri=ikle1
         if (ikle2 < the_tri) the_tri=ikle2 
         if (ikle3< the_tri)  the_tri=ikle3
         TRI_REF(the_tri,0)=TRI_REF(the_tri,0)+1
         pos=TRI_REF(the_tri,0)
         TRI_REF(the_tri,pos)=J
      END DO
      

      DO IELEB = 1,NBTRI
         IELEM = NELBOR(IELEB)
         ikle1 = NBOR(IKLBOR(IELEB,1))
         ikle2 = NBOR(IKLBOR(IELEB,2))
         ikle3 = NBOR(IKLBOR(IELEB,3))
         the_tri=ikle1
         if (ikle2 < the_tri) the_tri=ikle2 
         if (ikle3<the_tri)  the_tri=ikle3
         pos=TRI_REF(the_tri,0)
         IS = .FALSE.
          M  = -1
         DO JJ = 1, pos
            J=TRI_REF(the_tri,JJ)
            K = 3*(J-1)+1            
            IF ((ikle1.EQ.IKLESTRI(K)).AND.
     &          (ikle2.EQ.IKLESTRI(K+1)).AND.
     &          (ikle3.EQ.IKLESTRI(K+2))) THEN
               IS = .TRUE.
            ELSE IF ((ikle1.EQ.IKLESTRI(K)).AND.
     &          (ikle3.EQ.IKLESTRI(K+1)).AND.
     &          (ikle2.EQ.IKLESTRI(K+2))) THEN
               IS = .TRUE.
            ELSE IF ((ikle2.EQ.IKLESTRI(K)).AND.
     &          (ikle1.EQ.IKLESTRI(K+1)).AND.
     &              (ikle3.EQ.IKLESTRI(K+2))) THEN
               IS = .TRUE.
            ELSE IF ((ikle2.EQ.IKLESTRI(K)).AND.
     &          (ikle3.EQ.IKLESTRI(K+1)).AND.
     &          (ikle1.EQ.IKLESTRI(K+2))) THEN
               IS = .TRUE.
            ELSE IF ((ikle3.EQ.IKLESTRI(K)).AND.
     &          (ikle1.EQ.IKLESTRI(K+1)).AND.
     &          (ikle2.EQ.IKLESTRI(K+2))) THEN
               IS = .TRUE.
            ELSE IF ((ikle3.EQ.IKLESTRI(K)).AND.
     &          (ikle2.EQ.IKLESTRI(K+1)).AND.
     &          (ikle1.EQ.IKLESTRI(K+2))) THEN
               IS = .TRUE.         
            ENDIF
            IF (IS) THEN
               M = J
               EXIT
            ENDIF
         ENDDO
         DO I = 1,4
            IF (IFABOR(IELEM,I).EQ.0) THEN
               IF ((IKLE1.EQ.(IKLE(NELBOR(IELEB),SOMFAC(1,I))))
     &         .AND.(IKLE2.EQ.(IKLE(NELBOR(IELEB),SOMFAC(2,I))))  
     &         .AND. (IKLE3.EQ.(IKLE(NELBOR(IELEB),SOMFAC(3,I)))))
     &         THEN
                  NI = TETTRI2(IELEM)
                  N  = 4*(IELEM-1)+NI+1
                  TETTRI(N) = M
c                  write(*,*) N, '---> ',M
                  TETTRI2(IELEM) = NI + 1
               ENDIF
            ENDIF
         END DO
      ENDDO
!
      deallocate(tri_ref)
!
!
!
!
      call system_clock(count=temps_sc(3),count_rate=parsec)
!------------------------------------------------------- Fin Version F.D




! 3b. Construction de nodes1/nodes2/nodes3: connectivite inverse noeud > tetra
!     pour l'ecriture a la volee des unv locaux
!---------------
! Parcours des mailles pour connaitre le nombre de mailles qui 
! les reference
      allocate(nodes1(npoint),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' nodes1')
      nodes1(:)=0
      do i=1,nbtet
        do k=1,4
          ikleb=iklestet(4*(i-1)+k)
          nodes1(ikleb)=nodes1(ikleb)+1
        enddo
      enddo
! nombre de referencement de points et pointeur nodes2 vers nodes3
! le ieme point a sa liste de tetra (en numerotation locale tetra)
! de nodes3(nodes2(i)) a nodes3(nodes2(i)+nodes1(i)-1)
      allocate(nodes2(npoint+1),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' nodes2')
      compt=0
      nodes2(1)=1
      do i=1,npoint
        compt=compt+nodes1(i)
        nodes2(i+1)=compt+1
      enddo
! Pour un noeuds donne, qu'elles sont les mailles qui le concernent
      allocate(nodes3(compt),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' nodes3')
      nodes3(:)=-1
      do i=1,nbtet
        do k=1,4
          ikleb=iklestet(4*(i-1)+k)
          ni=nodes2(ikleb)
          nf=ni+nodes1(ikleb)-1
          nt=-999
          do n=ni,nf ! on cherche le premier indice de libre de nodes3
            if (nodes3(n)==-1) then
              nt=n
              exit
            endif
          enddo ! en n
          if (nt==-999) then
            goto 146  ! pb de dimensionnement de vecteurs nodesi
          else
            nodes3(nt)=i  ! numero local du tetra i associe au noeud nt
          endif
        enddo
      enddo

! 3c. Construction de nodes1t/nodes2t/nodes3t: connectivite inverse noeud > tria
!     pour la couleur des noeuds (Dirichlet sur l'interface)
!---------------
      allocate(nodes1t(npoint),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' nodes1t')
      nodes1t(:)=0
      do i=1,nbtri
        do k=1,3
          ikleb=iklestri(3*(i-1)+k)
          nodes1t(ikleb)=nodes1t(ikleb)+1
        enddo
      enddo
      allocate(nodes2t(npoint+1),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' nodes2t')
      compt=0
      nodes2t(1)=1
      do i=1,npoint
        compt=compt+nodes1t(i)
        nodes2t(i+1)=compt+1
      enddo
      allocate(nodes3t(compt),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' nodes3t')
      nodes3t(:)=-1
      do i=1,nbtri
        do k=1,3
          ikleb=iklestri(3*(i-1)+k)
          ni=nodes2t(ikleb)
          nf=ni+nodes1t(ikleb)-1
          nt=-999
          do n=ni,nf ! on cherche le premier indice de libre de nodes3t
            if (nodes3t(n)==-1) then
              nt=n
              exit
            endif
          enddo ! en n
          if (nt==-999) then
            goto 146  ! pb de dimensionnement de vecteurs nodesi
          else
            nodes3t(nt)=i  ! numero local du tetra i associe au noeud nt
          endif
        enddo
      enddo
      call system_clock(count=temps_sc(4),count_rate=parsec)
      write(lu,*)' TEMPS CONNECTIVITE INVERSE PART1/ PART2',
     &          (1.0*(temps_sc(3)-temps_sc(2)))/(1.0*parsec),'/',
     &          (1.0*(temps_sc(4)-temps_sc(3)))/(1.0*parsec),' seconds'
         
!----------------------------------------------------------------------
! 4. Partitioning
!---------------

!        do i=1,4*nbtet
!                write(lu,*) 'tettrialpha',tettri(i)       
!        enddo
!        do i=1,nbtet
!                write(lu,*) 'tettribeta',tettri2(i)
!        enddo

      write(lu,*)' '
      write(lu,*)' Starting METIS mesh partitioning------------------+'
      allocate(epart(nbtet),stat=ierr)
      if (ierr.ne.0) call ALLOER (lu, 'epart')
      allocate (npart(npoint),stat=ierr)
      if (ierr.ne.0) call ALLOER (lu, 'npart')
      call system_clock(count=temps_sc(5),count_rate=parsec)
! Partitionnement des mailles
      call METIS_PartMeshDual(nbtet,npoint,iklestet,2,1,nparts,edgecut,
     &    epart,npart)
      call system_clock(count=temps_sc(6),count_rate=parsec)
      write(lu,*)' '
      write(lu,*)' End METIS mesh partitioning------------------+'
      write(lu,*)' TEMPS CONSOMME PAR  METIS ',
     &           (1.0*(temps_sc(6)-temps_sc(5)))/(1.0*parsec),' seconds'
      write(lu,80) nelemtotal,npoint
      write(lu,81) nbtet,nbtri
      write(lu,82) edgecut,nparts
      write(lu,*) 'Sortie de Metis correcte'
CD ******************************************************
CD     loop  over the tetra to computer the number and the label
CD     of finite elements assigned to  each subdomain
CD ******************************************************
CD     computation of the maximum number of finite elements assigned to one subdomain
      nelem_p(:)=0
      npoin_p(:)=0
       do i=1,nbtet
         nelem_p(epart(i))=nelem_p(epart(i))+1
      end do
      max_nelem_p=maxval(nelem_p)
      nelem_p(:)=0
CD     allocation of the elegl array 
      allocate(elegl(max_nelem_p,nparts),stat=ierr)
CD     elegl is the filled 
      if (ierr.ne.0) call alloer(lu,'elegl')
      do i=1,nbtet
         nelem_p(epart(i))=nelem_p(epart(i))+1
         elegl(nelem_p(epart(i)),epart(i))=i
       end do
CD     compute the maximum of nodes assigned to one subdomain
       nodelg(:)=0
CD     for each subdomain idd
       do idd=1,nparts  
          nodelg(:)=0
CD         loop on the finite elements ielem assigned to subdomain idd
          do pos=1,nelem_p(idd)
            ielem=elegl(pos,idd)
            n=4*(ielem-1)+1
CD          loop of the node contained in ielem            
            do k=0,3
               node=iklestet(n+k)
               if (nodelg(node) .eq. 0) then
                  npoin_p(idd)=npoin_p(idd)+1
                  nodelg(node)=npoin_p(idd)
               end if
            end do
         end do
      end do
CD    allocation and filling of  the nodegl array      
      max_npoin_p=maxval(npoin_p)
      npoin_p(:)=0
      nodelg(:)=0
!
      allocate(nodegl(max_npoin_p,nparts),stat=ierr)
       if (ierr.ne.0) call alloer(lu,'nodegl')
       do idd=1,nparts
          nodelg(:)=0
          do pos=1,nelem_p(idd)
             ielem=elegl(pos,idd)
             n=4*(ielem-1)+1
             do k=0,3
                node=iklestet(n+k)
                if (nodelg(node) .EQ. 0) then
                   npoin_p(idd)=npoin_p(idd)+1
                   nodelg(node)=npoin_p(idd)
                   nodegl(npoin_p(idd),idd)=node
                end if
             end do
          end do
       end do
!            
!----------------------------------------------------------------------
! 5a. Allocations pour ecriture des fichiers .unv/.log associant un sous-domaine
!     par proc
!------------

      nameinp2=nameinp
      namelog2=namelog
      blanc='    '
      moins1='-1'
      allocate(nodes4(npoint),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' nodes4')
c$$$      nodes4(:)=-1
      allocate(knolg(npoint),stat=ierr)      ! C'est sous-optimal en
      if (ierr.ne.0) call alloer(lu,' knolg')! terme de dimensionnement
      knolg(:)=-1      ! mais plus rapide pour le remplissage ulterieur
! 
! Parametre NBSDOMVOIS (Nombre de sous domaines voisins+2)
!      
      allocate(nachb(NBSDOMVOIS,npoint),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' nachb')
      nachb(1,:)=0
      do j=2,NBSDOMVOIS-1
        nachb(j,:)=-1
      enddo
      allocate(triunv(4*nbtri),stat=ierr)
      if (ierr.ne.0) call ALLOER (lu, 'triunv')
!
!
! 5b. Recherche de la vrai couleur aux noeuds pour eviter les pbs de Dirichlet
!     aux interfaces
!---------------
      ncolor2(:)=-1
      do j=1,npoint      ! boucle sur tous les points du maillages
        ni=nodes2t(j)
        nf=ni+nodes1t(j)-1
        do n=ni,nf       ! boucle sur les tetra contenant le point j
          numtet=nodes3t(n)   ! tria de numero local numtet
          numtrig=convtri(numtet)  ! numero global du triangle
          color1=ecolor(numtrig)   ! couleur du noeud avec ce tria
          color2=ncolor2(j)
          if (color2 > 0) then   ! On priorise les couleurs
            pr1=0
            pr2=0
            do l=1,nbcolor
              if (priority(l)==color1) then
                pr1=l
 1           endif
              if (priority(l)==color2) then
                pr2=l
              endif
            enddo
            if ((pr1==0).or.(pr2==0)) goto 154
            if (pr1<pr2) ncolor2(j)=color1  ! on change de couleur
          else        ! Premiere fois que ce noeud est traite
            ncolor2(j)=color1
          endif
        enddo
      enddo

      call system_clock(count=temps_sc(7),count_rate=parsec)

!      DO IELEM = 1, NPOINT
!         write(lu,*) 'ncolor2',ncolor2(ielem)
!      ENDDO

!      DO IELEM = 1, NBCOLOR
!         write(lu,*) 'prior',priority(ielem)
!      ENDDO
      
! ob d
!--------------
! Rajout pour tenir compte des couleurs des noeuds de tetras lies
! au tria de bord et situes dans d'autres sd
!--------------
      allocate(tetcolor(nbtet,4),stat=ierr)
      if (ierr.ne.0) call alloer(lu,' tetcolor')
      tetcolor(:,:)=.false.
      nbretouche=0
      do iptfr=1,nptfr      ! boucle sur tous les points de bord
        j=nbor(iptfr)
!       on ne fait qqe chose (eventuellement) que si il y a un tria
!       de bord (ecolor>0 et ncolor2 !=-1). Grace au traitement precedent
!       on s'en rend compte directement via ncolor2.
	linter=.false.
        nbtetj=nodes1(j) ! nbre de tetra rattaches a ce noeud
        ni=nodes2(j)     ! adresse dans nodes3 du premier
        nf=ni+nbtetj-1
	if (ncolor2(j) > 0) then
	  ! On cherche a savoir si le noeud est a l'interface linter=.true.
          do n=ni,nf       ! boucle sur les tetra contenant le point j
            nt=nodes3(n)   ! tetra de numero local nt
	    if (n == ni) then
	      iddnt=epart(nt)
	    else
	      if (epart(nt) /= iddnt) then
	        linter=.true.
		goto 20     ! on a le renseignement demande, on sort
	      endif
	    endif
	  enddo	          ! fin boucle sur les tetras
   20     continue
!         Le noeud j est un noeud d'interface. On va communiquer au noeud
!         correspondant des tetras (si un tria de bord n'est pas sur cette
!         face auxquel cas le pb est deja regle), la bonne couleur.
	  if (linter) then  
            do n=ni,nf       ! boucle sur les tetra contenant le point j
              nt=nodes3(n)   ! tetra de numero local nt
!         On va trier les cas non pathologiques et tres courant de tetra
!         dont une face coincide avec ce triangle
              if (tettri2(nt)>0) then   !tetra concerne par un tria
	        nit=4*(nt-1)+1
		nft=nit+tettri2(nt)-1
		do mt=nit,nft           ! boucle sur les tria du tetra
                  numtri=tettri(mt)     ! num local du tria
                  numtrib=3*(numtri-1)+1
                  ikle1=iklestri(numtrib) ! numero globaux des noeuds du tria
                  ikle2=iklestri(numtrib+1)
                  ikle3=iklestri(numtrib+2)
!                 Ce point j appartient deja a un tria acolle au tetra
!                 On saute le tetra nt
		  if ((ikle1==j).or.(ikle2==j).or.(ikle3==j)) then
! pour tests
!                    write(lu,*)'je saute le tetra ',nt,epart(nt),
!     &                          tettri2(nt),' nodes ',j
                    goto 21
                  endif		  
		enddo
	      endif            ! fin si tettri
!             le tetra nt est potentiellement oublie, on le traite au cas ou
!             le partage se fera dans estel3d/read_connectivity
	      numtetb=4*(nt-1)+1
	      do l=1,4
                ikle1=iklestet(numtetb+l-1) ! numero globaux des noeuds du tetra
		if (ikle1==j) then
                  tetcolor(nt,l)=(tetcolor(nt,l).or..true.)
                  nbretouche=nbretouche+1
                endif
	      enddo  ! en l
   21        continue
	    enddo	          ! fin boucle sur les tetras	    		    
	  endif            ! fin si linter   
	endif             ! fin si sur ncolor2
      enddo              ! fin boucle sur les points de bord     
! ob f
      call system_clock(count=temps_sc(8),count_rate=parsec)
      write(lu,*)' NOMBRE DE RETOUCHE DU PARTITIONNEMENT (PART2): ',
     &           nbretouche
      write(lu,*)' TEMPS DE RETOUCHE DU PARTITIONNEMENT PART1/PART2',
     &            (1.0*(temps_sc(7)-temps_sc(6)))/(1.0*parsec),'/',
     &           (1.0*(temps_sc(8)-temps_sc(7)))/(1.0*parsec),' seconds'
c$$$      write(lu,*)'idem version de reference'

! 5c. Remplissage effectif du unv par sd
!---------------
      ibid = 1
!
      allocate(deja_trouve(nelin),stat=ierr)
      if (ierr.ne.0) call alloer(lu,'deja_trouve')
      deja_trouve(:)=.false.
!      
      do idd=1,nparts  ! Boucle sur les sous-domaines

! nombre de triangles pour ce sous-domaine
        nbtriidd=0
! nom du fichier unv par sous-domaine
        nameinp2(i_leninp+1:i_leninp+11) = extens(nparts-1,idd-1)
        open(ninp2,file=nameinp2,status='unknown',form='formatted',
     &       err=132)
        rewind(ninp2)

! nom du fichier log par sous-domaine
        namelog2(i_lenlog+1:i_lenlog+11) = extens(nparts-1,idd-1)
        open(nlog2,file=namelog2,status='unknown',form='formatted',
     &       err=133)
        rewind(nlog2)

! titre (unv par sd)
        write(ninp2,60,err=112)blanc,moins1
        write(ninp2,61,err=112)nsec1
        write(ninp2,62,err=112)titre
        titre = ' ' 
        write(ninp2,62,err=112)titre
        write(ninp2,62,err=112)titre
        write(ninp2,62,err=112)titre
        write(ninp2,62,err=112)titre
        write(ninp2,62,err=112)titre
        write(ninp2,62,err=112)titre
        write(ninp2,60,err=112)blanc,moins1
!
! bloc sur les coordonnees/couleurs des noeuds (unv par sd)
        write(ninp2,60,err=112)blanc,moins1
        write(ninp2,61,err=112)nsec2
        compt=1
        nodes4(:)=-1
CD      new version of the loop to reduce the computing time
        do pos_node=1,npoin_p(idd) ! boucle sur tous les points du maillages
           j=nodegl(pos_node,idd)
CD       previous version of the loop
CD       ni=nodes2(j)
CD       nf=ni+nodes1(j)-1
CD       do n=ni,nf       ! boucle sur les tetra contenant le point j
CD           nt=nodes3(n)   ! tetra de numero local nt 
CD           if (epart(nt)==idd) then     ! C'est une maille du sous-domaine
              write(ninp2,63,err=112)compt,ibid,ibid,ncolor2(j)
              write(ninp2,64,err=112)x1(j),y1(j),z1(j)
              nodes4(j)=compt   ! le noeud j a le numero compt
                                ! pour le sous-domaine idd
! pour parallelisme TELEMAC
              knolg(compt)=j ! conversion sd (local)-->maillage entier (global)
              k=nachb(1,j)   ! nbre de sd contenant le noeud j
              nachblog=.true.       
              do l=1,k     ! noeud deja concerne par ce sd ?
                if (nachb(1+l,j)==idd) nachblog=.false.  ! oui      
              enddo
              if (nachblog) then                         ! non
                k=nachb(1,j)+1
                if (k.gt.NBSDOMVOIS-2) goto 151
                nachb(k+1,j)=idd  ! noeud global j concerne par idd
                nachb(1,j)=k      ! sa multiplicite
              endif 
              compt=compt+1
!              goto 10 ! on passe au noeud suivant           
!            endif  ! en epart
!          enddo ! en n
!   10     continue
        enddo   ! en j
! pour tests
!      do i=1,npoint
!        write(lu,*)'global numero point: ',i,' local: ',nodes4(i)
!      enddo
        npointsd(idd)=compt-1  ! nombre de noeuds du sous-domaine idd
        write(ninp2,60,err=112)blanc,moins1

! bloc sur les connectivites/couleurs des mailles (unv par sd)
        write(ninp2,60,err=112)blanc,moins1
        write(ninp2,61,err=112)nsec3
        compt=1
        ibid = 1
CD      previous version of the loop
CD      do j=1,nelemtotal
CD      if (typelem(j,1)==111) then ! C'est un tetraedre
CD        numtet=typelem(j,2) ! num local du tetra dans la liste des tetras
CD            if (epart(numtet)==idd) then 
        do pos=1,nelem_p(idd)
                                ! boucle sur tetra et tria pour ecolor
           j=elegl(pos,idd)
           numtet=typelem(j,2)  ! num local du tetra dans la liste des tetras 
           elem=111
C ob d
! pretraitement pour les eventuels pb de couleurs des noeuds de tetras
! a l'interface
              ibidc=0
	      if (tetcolor(numtet,1)) ibidc=ibidc+1000
	      if (tetcolor(numtet,2)) ibidc=ibidc+ 200
	      if (tetcolor(numtet,3)) ibidc=ibidc+  30
	      if (tetcolor(numtet,4)) ibidc=ibidc+   4
! pour monitoring
!              if (ibidc/=0) write(6,*)'idd',idd,'partel',j,compt,ibidc
! idem version de reference
!	      ibidc=0
! ob f
              write(ninp2,65,err=112)compt,elem,-ibidc,ibid,ecolor(j),4
              compt=compt+1
              n=4*(numtet-1)+1
              ikle1=nodes4(iklestet(n))
              ikle2=nodes4(iklestet(n+1))
              ikle3=nodes4(iklestet(n+2))
              ikle4=nodes4(iklestet(n+3))
              write(ninp2,66,err=112)ikle1,ikle2,ikle3,ikle4
       if ((ikle1.lt.0).or.(ikle2.lt.0).or.(ikle3.lt.0).or.(ikle4.lt.0))
     &          goto 147
              if (tettri2(numtet).ne.0) then
                ni=4*(numtet-1)+1
                nf=ni+tettri2(numtet)-1
                do m=ni,nf   ! on parcourt les triangles de bord associes
                  numtri=tettri(m)  ! au numtet tetraedre; num local du tria
                  numtrig=convtri(numtri)  ! numero global du triangle
                  elem=91
                  triunv(4*nbtriidd+1)=ecolor(numtrig)
                  n=3*(numtri-1)+1
                  ikle1=nodes4(iklestri(n))
                  ikle2=nodes4(iklestri(n+1))
                  ikle3=nodes4(iklestri(n+2))
                  triunv(4*nbtriidd+2)=ikle1
                  triunv(4*nbtriidd+3)=ikle2
                  triunv(4*nbtriidd+4)=ikle3
                  nbtriidd=nbtriidd+1
!
              if ((ikle1.lt.0).or.(ikle2.lt.0).or.(ikle3.lt.0)) goto 147
!
                enddo  ! en m
              endif  ! en tettri2
        !    endif  ! en epart
      !    endif  ! en typelem
        enddo ! en j

! Maintenant on peux recopier le bloc des triangles !
        elem=91
        do j=1,nbtriidd
          write(ninp2,65,err=112)compt,elem,ibid,ibid,
     &                           triunv(4*(j-1)+1),3
          ikle1=triunv(4*(j-1)+2)
          ikle2=triunv(4*(j-1)+3)
          ikle3=triunv(4*(j-1)+4)
          write(ninp2,67,err=112)ikle1,ikle2,ikle3        
          compt=compt+1
        enddo  ! en j
!
        elem=91
! Boucle surdimensionnee, on boucle sur le nombre de surface interne du maillage global...                                
        do j=1,nelin
           if (deja_trouve(j)) cycle
           ikle1=nodes4(iklein(j,2))
           ikle2=nodes4(iklein(j,3))
           ikle3=nodes4(iklein(j,4))
           if ((ikle1.eq.-1).or.(ikle2.eq.-1).or.(ikle3.eq.-1)) cycle
           write(ninp2,65,err=112) compt,elem,ibid,ibid,iklein(j,1),3
           write(ninp2,67,err=112) ikle1,ikle2,ikle3
           compt = compt+1
           deja_trouve(j) = .true.
        enddo ! en j
!
c$$$        write(lu,*) 'SubDomain',idd,'InnerTri',compt
!        
        write(ninp2,60,err=112)blanc,moins1
!        write(ninp2,60,err=112)blanc,moins1
!        write(ninp2,61,err=112)nsec4
!        write(ninp2,68,err=112) 1,0,0,0,0,0,0,0
        close(ninp2)
        nelemsd(idd)=compt-1  ! nombre de mailles du sous-domaine idd

! 5d. Remplissage effectif du log par sd
!---------------
! element standard du fichier log (log par sd)
        write(nlog2,51 ,err=113) npointsd(idd)      
        write(nlog2,52 ,err=113) nelemsd(idd)
        write(nlog2,523,err=113) size_flux
        write(nlog2,53 ,err=113) nbfamily-1
        do j=1,nbfamily
          write(nlog2,50,err=113)logfamily(j)
        enddo

        ! Addition by JP Renaud on 15/02/2007
        ! As the list of priorities has moved in ESTEL-3D from
        ! the steering file to the log file, we need to write "a"
        ! number of external faces + priority list here. As these
        ! are not used in parallel mode, we merely copy the list
        ! from the original log file.

        write(nlog2,531,err=113) nbcolor
        write(unit=theformat,fmt=1000) nbcolor
1000    format('(''Priority :'',',i3,'(x,i3,))')
        theformat=trim(theformat)
!        write(lu,*) 'FORMATT =',theformat
        write (nlog2,fmt=theformat(1:len(theformat)-1))
     &  (priority(i), i=1, nbcolor)

        ! End addition by JP Renaud

! knolg (log par sd)
        nt=npointsd(idd)
        ni=nt/6
        nf=nt-6*ni
        write(nlog2,54,err=113)ni,nf
        do j=1,ni
          write(nlog2,540,err=113)(knolg(6*(j-1)+k),k=1,6)
        enddo
        if (nf.eq.1) then
          write(nlog2,541,err=113)knolg(6*ni+1)
        else if (nf.eq.2) then
          write(nlog2,542,err=113)(knolg(6*ni+k),k=1,2)
        else if (nf.eq.3) then
          write(nlog2,543,err=113)(knolg(6*ni+k),k=1,3)
        else if (nf.eq.4) then
          write(nlog2,544,err=113)(knolg(6*ni+k),k=1,4)
        else if (nf.eq.5) then
          write(nlog2,545,err=113)(knolg(6*ni+k),k=1,5)
        endif
        write(nlog2,55,err=113)npoint  ! nombre de noeud du maillage
                    ! initial pour allocation KNOGL dans ESTEL
!
      enddo  ! boucle sur les sous-domaines

! 5e. Travaux supplementaires pour determiner le nachb avant de l'ecrire
!      dans le log
!---------------
      do idd=1,nparts  ! Boucle sur les sous-domaines
! Construction et dimensionnement du nachb propre a chaque sd
        compt=0
        nachb(NBSDOMVOIS,:)=-1
        do j=1,npoint      ! boucle sur tous les points du maillage
          n=nachb(1,j)
          if (n>1) then    ! point d'interface
            n=n+1
            do k=2,n
              if (nachb(k,j)==idd) then ! il concerne idd
                compt=compt+1   ! "compt"ieme point d'interface de idd
                nachb(NBSDOMVOIS,j)=compt  ! a retenir comme point d'interface
              endif
            enddo            ! fin boucle sur les sd du point j
          endif   
        enddo              ! fin boucle points
        npointisd(idd)=compt ! nombre de points d'interface de idd

! 5f. On continue l'ecriture du .log
!-------------
        namelog2(i_lenlog+1:i_lenlog+11) = extens(nparts-1,idd-1)
        open(nlog2,file=namelog2,status='old',form='formatted',
     &       position='append',err=133)
        write(nlog2,56,err=113) npointisd(idd)
        do j=1,npoint
          if (nachb(NBSDOMVOIS,j)>0) then  ! c'est un point d'interface de idd
            compt=0
            vectnb(:)=-1
            do k=1,NBSDOMVOIS-2    ! on prepare l'info pour le nachb telemac
              if (nachb(k+1,j)/= idd) then
                compt=compt+1
! Attention a celui-ci, surement lie au numero de points...
C OB D
                if (compt.gt.NBSDOMVOIS-3) goto 152
C OB F
                if (nachb(k+1,j)>0) then
! on stocke le numero de proc et non le numero de sous-domaine
! d'ou la contrainte, un proc par sous-domaine
                  vectnb(compt)=nachb(k+1,j)-1
                endif
              endif 
            enddo  ! en k
!            write(nlog2,561,err=113)j,(vectnb(k),k=1,NBSDOMVOIS-3)
             nt = NBSDOMVOIS-3	     
             ni=nt/6
             nf=nt-6*ni+1
	     write(nlog2,640,err=113)j,(vectnb(k),k=1,5)
	     do l=2,ni
	       write(nlog2,640,err=113)(vectnb(6*(l-1)+k),k=0,5)
	     enddo
	     if (nf.eq.1) then
	       write(nlog2,641,err=113)vectnb(6*ni)
	     else if (nf.eq.2) then
               write(nlog2,642,err=113)(vectnb(6*ni+k),k=0,1)
	     else if (nf.eq.3) then
	       write(nlog2,643,err=113)(vectnb(6*ni+k),k=0,2)
	     else if (nf.eq.4) then
	       write(nlog2,644,err=113)(vectnb(6*ni+k),k=0,3)
	     else if (nf.eq.5) then
               write(nlog2,645,err=113)(vectnb(6*ni+k),k=0,4)
             endif
          endif
        enddo  ! fin boucle en j
        write(nlog2,57,err=113)        
        close(nlog2)    
      enddo  ! boucle sur les sous-domaines
      call system_clock(count=temps_sc(9),count_rate=parsec)
      write(lu,*)' REMPLISSAGE DES FICHIERS UNV ET LOG',
     &           (1.0*(temps_sc(9)-temps_sc(8)))/(1.0*parsec),' seconds'
!----------------------------------------------------------------------
! 6. Affichages dans partel.log et test de completude du partitionnement
!------------
 
      write(lu,*)' '
      compt1=0
      compt2=0
      compt3=0
      do idd=1,nparts
        write(lu,86)idd,nelemsd(idd),npointsd(idd),npointisd(idd)
        compt3=compt3+npointisd(idd)
        compt2=compt2+npointsd(idd)
        compt1=compt1+nelemsd(idd)
      enddo
      write(lu,*)' ------------------------------------'
      write(lu,87)compt1,compt2,compt3
      write(lu,88)compt1/nparts,compt2/nparts,compt3/nparts
      write(lu,*)' '
      write(lu,83)(1.0*(temps_sc(9)-temps_sc(1)))/(1.0*parsec)
      write(lu,*)' Ending METIS mesh partitioning--------------------+'
      write(lu,*)' '
      write(lu,*)' Writing geometry file for each processor'
      write(lu,*)' Writing log file for each processor'

!----------------------------------------------------------------------
! 7. Divers
!---------------

! 7.a Format du log
!---------------
   50 format(a80)         ! les autres lignes
!             1234567890123456789012345678901234567890123456789
   51 format(' Total no. of nodes                   :    ',i10)
   52 format(' Total no. of elements                :    ',i10)
  523 format(' Total no. of user-flux               :    ',i10)
   53 format(' Total no. of families                :    ',i10)
  531 format(' Total number of external faces       :    ',i10)
   54 format(' Debut de KNOLG: ',i10,' ',i10)

 5401 format(6i5)              ! PRIORITY
 5411 format(i5)               ! 
 5421 format(2i5)              ! 
 5431 format(3i5)              ! 
 5441 format(4i5)              ! 
 5451 format(5i5)              ! 

  540 format(6i10)        ! ligne de bloc KNOLG et PRIORITY
  541 format(i10)         ! derniere ligne de bloc KNOLG
  542 format(2i10)        ! derniere ligne de bloc KNOLG
  543 format(3i10)        ! derniere ligne de bloc KNOLG
  544 format(4i10)        ! derniere ligne de bloc KNOLG
  545 format(5i10)        ! derniere ligne de bloc KNOLG

  641 format(i7)         ! derniere ligne de bloc NACHB
  642 format(2i7)        ! derniere ligne de bloc NACHB
  643 format(3i7)        ! derniere ligne de bloc NACHB
  644 format(4i7)        ! derniere ligne de bloc NACHB
  645 format(5i7)        ! derniere ligne de bloc NACHB
  640 format(6i7)        ! derniere ligne de bloc NACHB


  
   55 format(' Fin de KNOLG: ',i10)
   56 format(' Debut de NACHB: ',i10)
  561 format(10i10)        ! ligne de bloc NACHB
   57 format(' Fin de NACHB: ')

! 7b. Format du unv
!---------------
   60 format(a4,a2)       ! '    -1'   
   61 format(i6)          ! lecture nsec
   62 format(a80)         ! lecture titre      
   63 format(4i10)        ! ligne 1 bloc coord      
   64 format(3d25.16)     ! ligne 2 bloc coord      
   65 format(6i10)        ! ligne 1 bloc connectivite      
   66 format(4i10)        ! ligne 2 bloc connectivite si tetra      
   67 format(3i10)        ! ligne 2 bloc connectivite si triangle
   68 format(8i10)        ! bloc fantoche pour marquer la fin du bloc
                          ! connectivitee
      
! 7.c Affichages dans partel.log
!---------------
   80 format(' #Number total of elements: ',i8,
     &       ' #Nodes                 : ',i8)
   81 format(' #Tetrahedrons            : ',i8,
     &       ' #Triangle mesh border  : ',i8)
   82 format(' #Edgecuts                : ',i8,
     &       ' #Nparts                : ',i8)   
   83 format('  Runtime                 : ',f10.2,' s')
   86 format('  Domain: ',i3,' #Elements:   ',i8,' #Nodes:   ',i8,
     &       ' #InterfaceNodes:   ',i8)
   87 format('  Total values of Elements: ',i10,'  Nodes: ',i10,
     &       '  InterfaceNodes: ',i10)
   88 format('  Mean values of Elements :   ',i8,'  Nodes:   ',i8,
     &       '  InterfaceNodes:   ',i8)
   89 format('  Input unv file      :',a50)
   90 format('  Input log file      :',a50)
   91 format('  Number of partitions:',i5)
   92 format('  Number of nodes:',i10)
   93 format('  Number of elements:',i10)
   94 format('  Number of colors:',i5)

! 7.d Deallocate
!---------------
      deallocate(x1,y1,z1)
      deallocate(ncolor,ecolor)
      deallocate(iklestet,iklestri,typelem,convtri,tettri,tettri2)
      deallocate(epart,npart)
      deallocate(nelemsd,npointsd,npointisd)
      deallocate(nodes1,nodes2,nodes3,nodes4,triunv)
      deallocate(nodes1t,nodes2t,nodes3t)
      deallocate(knolg,nachb,priority,ncolor2)
      deallocate(elegl)
      deallocate(nodegl)
      deallocate(nodelg)
      deallocate(nelem_p)
      deallocate(npoin_p)
      return

! 7.e Messages d'erreurs
!---------------
  110 texterror='! Unexpected file format: '//namelog//' !'
      goto 999
  111 texterror='! Unexpected file format: '//nameinp//' !'
      goto 999
  112 texterror='! Unexpected file format: '//nameinp2//' !'
      goto 999
  113 texterror='! Unexpected file format: '//namelog2//' !'
      goto 999
  120 texterror='! Unexpected EOF while reading: '//namelog//' !'
      goto 999
  121 texterror='! Unexpected EOF while reading: '//nameinp//' !'
      goto 999
  130 texterror='! Problem while opening: '//namelog//' !'
      goto 999
  131 texterror='! Problem while opening: '//nameinp//' !'
      goto 999
  132 texterror='! Problem while opening: '//nameinp2//' !'
      goto 999
  133 texterror='! Problem while opening: '//namelog2//' !'
      goto 999
  140 texterror='! File does not exist: '//nameinp//' !'
      goto 999
  141 texterror='! File does not exist: '//namelog//' !'
      goto 999
  142 texterror='! Unknown type of element in the mesh !'
      goto 999
  143 do j = 1,nelemtotal
        if (typelem(j,2)==i) write(unit=str8,fmt='(i8)')j
      enddo
      write(unit=str26,fmt='(i8,x,i8,x,i8)')ikle1,ikle2,ikle3
      texterror='! Border surface of number '//str8//' and of nodes '//
     &          str26//' not link to a tetrahedron !'
      goto 999
  144 write(unit=str8,fmt='(i8)')maxlensoft
      texterror='! Name of input file '//nameinp//' is longer than '//
     &           str8(1:3)//' characters !'
      goto 999
  145 write(unit=str8,fmt='(i8)')maxlensoft
      texterror='! Name of input file '//namelog//' is longer than '//
     &           str8(1:3)//' characters !'
      goto 999
  146 texterror='! Problem with construction of inverse connectivity !'
      goto 999
  147 texterror='! Problem while writing: '//nameinp2//' !'
      goto 999
  148 texterror='! Several elements may be forgotten by partitionning !' 
      goto 999
  149 texterror='! No input unv file !' 
      goto 999
  150 texterror='! No input log file !' 
      goto 999
!  151 write(unit=str8,fmt='(i8)')j
!      write(unit=str26,fmt='(i3,x,i3,x,i3,x,i3,x,i3,x,i3)')
!     &                 (nachb(k,j),k=2,6),idd
  151 write(unit=str8,fmt='(i8)')j
      write(unit=str26,fmt='(i3,x,i3,x,i3,x,i3,x,i3,x,i3)')
     &                 (nachb(k,j),k=2,NBSDOMVOIS-1),idd
      texterror='! Node '//str8//' belongs to domains '//str26(1:23)
     &                 //' !' 
      goto 999
  152 texterror='! Problem with construction of vectnb for nachb !' 
      goto 999
  153 write(unit=str8,fmt='(i8)')convtet(j)
      texterror='! Tetrahedron '//str8//
     &          ' links to several border triangles !'
      goto 999
  154 texterror='! Problem with the priority of color nodes !'
      goto 999
! End of file and format errors :
 1100 texterror='ERREUR DE LECTURE DU FICHIER UNV '//
     &  'VIA MESH_CONNECTIVITY'
      goto 999
 1200 texterror='ERREUR DE FIN DE LECTURE DU FICHIER UNV '//
     &  'VIA MESH_CONNECTIVITY'
      goto 999

  999 call alloer2(lu,texterror)
!
      end subroutine pares3d
!
!                       *************************
                        SUBROUTINE ELEBD31_PARTEL
!                       *************************
!
     *(NELBOR,NULONE,IKLBOR,IFABOR,NBOR,IKLE,
     * NBTET,NBTRI,NELMAX,NPOINT,NPTFR,IPOBO)     
!    
!***********************************************************************
! BIEF VERSION 5.5           09/04/04    J-M HERVOUET (LNH) 30 87 80 18
! 
!***********************************************************************
!              
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________.
! |      NOM       |MODE|                   ROLE                       |
! |________________|____|______________________________________________|
! |    NELBOR      |<-- | NUMERO DE L'ELEMENT ADJACENT AU KIEME SEGMENT|
! |    NULONE      |<-- | NUMERO LOCAL D'UN POINT DE BORD DANS         |
! |                |    | L'ELEMENT ADJACENT DONNE PAR NELBOR          
! |    IKLBOR      |<-- | NUMERO LOCAL DES NOEUDS A PARTIR D'UN ELEMENT
! |                |    |  DE BORD
! |    IFABOR      | -->| TABLEAU DES VOISINS DES FACES.
! |    NBOR        | -->| NUMERO GLOBAL D'UN NOEUD A PARTIR DU NÃÂÃÂ° LOCAL
! |    IKLE        | -->| NUMEROS GLOBAUX DES POINTS DE CHAQUE ELEMENT.
! |    NBTET       | -->| NOMBRE TOTAL D'ELEMENTS DANS LE MAILLAGE.
! |    NPOINT      | -->| NOMBRE TOTAL DE POINTS DU DOMAINE.
! |    NPTFR       | -->| NOMBRE DE POINTS FRONTIERES.
C |    NBTRI       | -->| NOMBRE D'ELEMENTS DE BORD.
! |    31          | -->| TYPE D'ELEMENT.
! |________________|____|______________________________________________|
!  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
!
!-----------------------------------------------------------------------
!
! Subroutine written by Lam Minh-Phuong
!
!-----------------------------------------------------------------------
!
! Subroutine: elebd31_partel
!
! Function: Construction de NELBOR, NULONE, IKLBORD 
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)    :: NBTET,NBTRI,NELMAX
      INTEGER, INTENT(IN)    :: NPOINT,NPTFR
      INTEGER, INTENT(IN)    :: NBOR(NPTFR)
      INTEGER, INTENT(IN)    :: IFABOR(NELMAX,4)
      INTEGER, INTENT(IN)    :: IKLE(NBTET,4)
      INTEGER, INTENT(OUT)   :: NELBOR(NBTRI),NULONE(NBTRI,3)
      INTEGER, INTENT(OUT)   :: IKLBOR(NBTRI,3),IPOBO(NPOINT)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER   :: IELEM, IELEB, J,K,IPOIN
!      INTEGER   :: IPOBO(NPOINT)
!
      INTEGER SOMFAC(3,4)
      DATA SOMFAC / 1,2,3 , 4,1,2 , 2,3,4 , 3,4,1  /
!     face numero:    1       2       3       4
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      
! Creation de IPOBO qui permet de passer du numero global au numero local
      
      DO IPOIN=1,NPOINT
        IPOBO(IPOIN) = 0
      ENDDO
      
      DO K = 1, NPTFR
         IPOBO(NBOR(K)) = K
      ENDDO
           
       
! Construction de NELBOR, NULONE, IKLBORD 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IELEB = 0
      
      DO IELEM = 1,NBTET
         DO J = 1,4
            IF (IFABOR(IELEM,J)== 0) THEN
               IELEB           = IELEB +1
               IF ( IELEB .GT. NBTRI ) THEN
                 IF(LNG.EQ.1) WRITE(LU,101)
                 IF(LNG.EQ.2) WRITE(LU,102)
101              FORMAT(1X,'ELEBD31_PARTEL : Erreur dans le Maillage ')
102              FORMAT(1X,'ELEBD31_PARTEL : Error in Mesh. bye.')
                 CALL PLANTE2(1)
                 STOP
               END IF
               NELBOR(IELEB)   = IELEM
               NULONE(IELEB,1) = SOMFAC(1,J)
               NULONE(IELEB,2) = SOMFAC(2,J)
               NULONE(IELEB,3) = SOMFAC(3,J)
               IKLBOR(IELEB,1) = IPOBO(IKLE(NELBOR(IELEB),SOMFAC(1,J)))
               IKLBOR(IELEB,2) = IPOBO(IKLE(NELBOR(IELEB),SOMFAC(2,J)))
               IKLBOR(IELEB,3) = IPOBO(IKLE(NELBOR(IELEB),SOMFAC(3,J)))
               
            END IF
         END DO
      END DO
 
!
!-----------------------------------------------------------------------
!
      RETURN
      END SUBROUTINE ELEBD31_PARTEL   
! 
!                       **************************
                        SUBROUTINE VOISIN31_PARTEL
!                       **************************
!
     *(IFABOR,NBTET,NELMAX,IKLE,SIZIKL,
     * NPOIN,NBOR,NPTFR)
!
!***********************************************************************
! BIEF VERSION 5.6      02/03/06    REGINA NEBAUER (LNHE) 01 30 87 83 93
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
! |    NBTET       | -->| NOMBRE D'ELEMENTS DANS LE MAILLAGE.
! |    NELMAX      | -->| NOMBRE MAXIMUM D'ELEMENTS DANS LE MAILLAGE.
! |    NPOIN       | -->| NOMBRE TOTAL DE POINTS DU DOMAINE
! |    NPTFR       | -->| NOMBRE DE POINTS DE BORD
! |    IKLE        | -->| Table de connectivite domaine
! |    SIZIKLE     | -->| ??
! |    NBOR        | -->| Correspondance no noeud de bord/no global
! |    NACHB       | -->| TABLEAU DE VOISINAGE POUR PARALLELISME
! !  IKLETR,NBTRI | -->/ CONNECTIVITE DES TRIA DE BORD POUR ESTEL3D
! |________________|____|_______________________________________________
!  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
!***********************************************************************
!
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN   ) :: NPTFR
      INTEGER, INTENT(IN   ) :: NBTET
      INTEGER, INTENT(IN   ) :: NELMAX
      INTEGER, INTENT(IN   ) :: NPOIN
      INTEGER, INTENT(IN   ) :: SIZIKL
      INTEGER, INTENT(IN), DIMENSION(NPTFR) :: NBOR
      ! NOTE : on donne explicitement la deuxieme dimension de IFABOR et
      ! IKLE, car il s'agit ici toujours de tetraedres!
      INTEGER, INTENT(OUT), DIMENSION(NELMAX,4) :: IFABOR
      INTEGER, INTENT(IN), DIMENSION(SIZIKL,4)  :: IKLE
!
! VARIABLES LOCALES 
!-----------------------------------------------------------------------

      ! Le tableau qui est l'inverse de NBOR (ca donne pour chaque noeud
      ! du domaine le numero de noeud de bord, ou zero si le noeud est a
      ! l'interieur du domaine.
c$$$      INTEGER, DIMENSION(NPOIN)            :: NBOR_INV
 
      ! Le tableau definissant le nombre d'element (tetraedres) voisins
      ! d'un noeud.
      INTEGER, DIMENSION(:  ), ALLOCATABLE  :: NVOIS
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

      ! Le nombre maximal d'elements voisin d'un noeud.
      INTEGER :: NMXVOISIN
      INTEGER :: IMAX       ! Dimensionnement tableau IADR

      INTEGER :: NFACE      ! Nombre de faces de l'element (tetra : 4)
      INTEGER :: NB_TRI      ! Le nombre de triangles definis

      INTEGER :: IELEM      ! Compteur elements
      INTEGER :: IELEM2     ! Compteur elements
      INTEGER :: IPOIN      ! Compteur noeuds domaine
      INTEGER :: INOEUD     ! Compteur noeuds tetraedres/triangles
      INTEGER :: IFACE      ! Compteur face
      INTEGER :: IFACE2     ! Compteur face
      INTEGER :: ITRI       ! Compteur trianlges
      INTEGER :: IVOIS      ! Compteur voisins
      INTEGER :: NV         ! Nombre de voisins

      INTEGER :: ERR        ! Code d'erreur allocation memoire

      LOGICAL :: found      ! trouve ou pas ...

      INTEGER :: I1, I2, I3 ! Les trois noeuds d'un triangle   
      INTEGER :: M1, M2, M3 ! La meme chose en ordonne.

      INTEGER :: I,J,K      ! ca sert ...

      !????
      INTEGER :: IR1,IR2,IR3,IR4,IR5,IR6,COMPT
      LOGICAL :: BORD


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
!-----------------------------------------------------------------------
! ETAPE 1 : Comptage du nombre d'elements voisins d'un noeud.
!-----------------------------------------------------------------------
! Calcul du nombre d'elements voisins pour chaque noeud du maillage.
! RESULTAT : NVOIS(INOEUD) donne le nombre d'elements voisins pour le
! noeud INOEUD

      ! On commence avec l'initialisation a 0 du compteur des elements
      ! voisins.

      NFACE = 4
!
      ALLOCATE(NVOIS(NPOIN),STAT=ERR)
      IF(ERR.NE.0) GOTO 999
!
      DO I = 1, NPOIN
        NVOIS(I) = 0
      ENDDO   
      ! Puis on compte les elements voisins.
      ! En parcourant la table de connectivite, on incremente le
      ! compteur a chaque fois qu'un element reference le noeud IPOIN

      ! Boucle sur les 4 noeuds de l'element
      DO INOEUD = 1, 4   
        ! Boucle sur les elements
        DO IELEM = 1,NBTET
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

      ! Pour chaque noeud des elements, on recupere sont identifiant.
      DO inoeud = 1, 4  ! Boucle sur les noeuds de l'element
        DO IELEM=1,NBTET ! Boucle sur les elements
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
      NB_TRI = NMXVOISIN * 3
!
      allocate(IKLE_TRI(NB_TRI,3),STAT=ERR)
      IF(ERR.NE.0) GOTO 999
      allocate(VOIS_TRI(NB_TRI,2),STAT=ERR)
      IF(ERR.NE.0) GOTO 999

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
          NB_TRI         = 0
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
                  DO ITRI = 1, NB_TRI
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
                                  IF(LNG.EQ.1) WRITE(LU,908) 31
                                  IF(LNG.EQ.2) WRITE(LU,909) 31
908                               FORMAT(1X,'VOISIN: IELM=',1I6,', 
     &                            PROBLEME DE VOISIN')
909                               FORMAT(1X,'VOISIN: IELM=',1I6,',
     &                            NEIGHBOUR PROBLEM')
                                  CALL PLANTE2(1)
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
                                CALL PLANTE2(1)
                                STOP
                               END IF
                               ! l'element et son voisin : on note la 
                               ! correspondance dans IFABOR.
                               IFABOR(IELEM ,IFACE ) = IELEM2
                               IFABOR(IELEM2,IFACE2) = IELEM
                               FOUND = .TRUE.
                          END IF
                      END IF
                  END DO
                  ! Et non, on n'a pas encore trouve ce triangle, alors
                  ! on en cree un nouveau.
                  IF ( .NOT. FOUND) THEN
                      NB_TRI             = NB_TRI + 1
                      IKLE_TRI(NB_TRI,1) = M1
                      IKLE_TRI(NB_TRI,2) = M2
                      IKLE_TRI(NB_TRI,3) = M3
                      VOIS_TRI(NB_TRI,1) = IELEM
                      VOIS_TRI(NB_TRI,2) = IFACE
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
!
!-----------------------------------------------------------------------
!  
!  DISTINCTION DANS IFABOR ENTRE LES FACES DE BORD ET LES FACES LIQUIDES
!
!  INITIALISATION DE NBOR_INV A ZERO
!
!      DO IPOIN=1,NPOIN
!        NBOR_INV(IPOIN) = 0
!      ENDDO      
!
!  PERMET DE PASSER DU NUMERO GLOBAL AU NUMERO DE BORD
!
!      DO K = 1, NPTFR
!         NBOR_INV(NBOR(K)) = K
!      ENDDO           
!
!  BOUCLE SUR TOUTES LES FACES DE TOUS LES ELEMENTS :
!
!      DO 90 IFACE = 1 , NFACE
!      DO 100 IELEM = 1 , NBTET
!
!      IF(IFABOR(IELEM,IFACE).EQ.-1) THEN
!
!      C'EST UNE VRAIE FACE DE BORD (LES FACES INTERNES EN PARALLELISME
!                                    SONT MARQUEES AVEC DES -2).
!      NUMEROS GLOBAUX DES POINTS DE LA FACE :
!
!       I1 = IKLE( IELEM , SOMFAC(1,IFACE) )
!       I2 = IKLE( IELEM , SOMFAC(2,IFACE) )
!       I3 = IKLE( IELEM , SOMFAC(3,IFACE) )
!
!      ENDIF
!
!100    CONTINUE
!90    CONTINUE
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
      CALL PLANTE2(1)
      STOP
!
!-----------------------------------------------------------------------
!  
      END
