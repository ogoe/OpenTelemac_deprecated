
!                 **********************************
                  SUBROUTINE READ_SECTIONS_TELEMAC2D
!                 **********************************
!
!***********************************************************************
! TELEMAC 2D VERSION 6.0     15 Feb 2010                    jaj pinxit 
!                                               jacek.jankowski@baw.de
!***********************************************************************
!   
!  => reading SECTIONS INPUT FILE in serial and parallel cases
!  => defining the definition of control sections, or...
!  => ...re-defining the ones declared previously in the parameter file
!  => sections defined by global node numbers or, alternatively, 
!     by co-ordinetes of terminal points (then nearest nodes found)
!  => in parallel case, two options:
!     -> taking the "serial" file (as "previously")
!     -> taking a partitioned file - computing fluxes through sections 
!        crossing numerous mesh partitions possible  
!
!  => subroutine modifies CTRLSC and NCP 
!
!  subroutine called by point_telemac2d
!  see also flusec_telemac2d or flusec_sisyphe
!
!***********************************************************************
!
      USE BIEF, ONLY: ncsize
      USE DECLARATIONS_TELEMAC2D, ONLY: mesh, chain, ncp, ctrlsc, 
     &                                  t2d_files, t2dsec
      IMPLICIT NONE
      INTEGER lng,lu
      COMMON/INFO/lng,lu
!
      INTEGER :: nsec, ihowsec, i, n, err, inp
      DOUBLE PRECISION :: xa, ya, distb, diste, dminb, dmine
!
!-----------------------------------------------------------------------
!
!      WRITE(lu,*) '-> entering read_sections_telemac2d'
      inp=t2d_files(t2dsec)%lu
      READ (inp,*) ! the obligatory comment line
      READ (inp,*) nsec, ihowsec
      IF (.NOT.ALLOCATED(chain)) THEN 
        ALLOCATE (chain(nsec), STAT=err)
        IF (err/=0) THEN
          WRITE(lu,*)
     &      'READ_SECTIONS: Error by reallocating chain:',err
          CALL plante(1)
          STOP 
        ENDIF
      ENDIF 
      
      SELECT CASE (ihowsec) 
      CASE (:-1) ! section terminal points as provided global nodes 
        DO n=1,nsec
          READ (inp,*) chain(n)%descr
          READ (inp,*) chain(n)%npair(:)
          IF (ncsize>1) THEN 
            chain(n)%xybeg(:)=0.0d0 
            chain(n)%xyend(:)=0.0d0
          ELSE
            chain(n)%xybeg(:)= (/mesh%x%r(chain(n)%npair(1)),
     &                           mesh%y%r(chain(n)%npair(1))/)
            chain(n)%xyend(:)= (/mesh%x%r(chain(n)%npair(2)),
     &                           mesh%y%r(chain(n)%npair(2))/)
          ENDIF 
          chain(n)%nseg=-1
          NULLIFY(chain(n)%liste)
        END DO 
!        WRITE(lu,'(a)') ' -> section, terminal coordinates:'
!        DO n=1,nsec
!          WRITE(lu,'(i9,4(1x,1pg13.6))') n, 
!     &          chain(n)%xybeg, chain(n)%xyend
!        END DO 
      CASE (0) ! section terminal points provided in co-ordinantes
        DO n=1,nsec
          READ (inp,*) chain(n)%descr
          READ (inp,*) chain(n)%xybeg(:), chain(n)%xyend(:)
          chain(n)%npair(:)=0
          chain(n)%nseg=-1
          NULLIFY(chain(n)%liste)
        END DO 
        DO n=1,nsec         ! find nearest nodes 
          xa=mesh%x%r(1) 
          ya=mesh%y%r(1)
          dminb = SQRT( (chain(n)%xybeg(1)-xa)**2 
     &                + (chain(n)%xybeg(2)-ya)**2 )
          dmine = SQRT( (chain(n)%xyend(1)-xa)**2 
     &                + (chain(n)%xyend(2)-ya)**2 )
          chain(n)%npair(1)=1
          chain(n)%npair(2)=1
          DO i=2,mesh%npoin ! computationally intensive 
            xa=mesh%x%r(i)
            ya=mesh%y%r(i)
            distb = SQRT( (chain(n)%xybeg(1)-xa)**2 
     &                  + (chain(n)%xybeg(2)-ya)**2 )
            diste = SQRT( (chain(n)%xyend(1)-xa)**2 
     &                 + (chain(n)%xyend(2)-ya)**2 )
            IF ( distb < dminb ) THEN 
              chain(n)%npair(1)=i
              dminb=distb
            ENDIF
            IF ( diste < dmine ) THEN 
              chain(n)%npair(2)=i
              dmine=diste 
            ENDIF 
          END DO
!          WRITE(lu,'(a,3(1x,i9))') 
!     &          ' -> section, terminal nodes: ', n, chain(n)%npair(:)
        END DO  
      CASE (1:) ! partitioned, instead of terminal points, ready chains provided 
        DO n=1,nsec
          READ (inp,*) chain(n)%descr
          READ (inp,*) chain(n)%nseg
          IF (chain(n)%nseg>0) THEN
            ALLOCATE (chain(n)%liste(chain(n)%nseg,2), STAT=err)
            IF (err/=0) THEN
              WRITE(lu,*) 'READ_SECTIONS_TELEMAC2D: ',
     &         ' Error by reallocating chain(n)%liste, n, err:',n,err
              CALL plante(1)
              STOP
            ENDIF
            DO i=1,chain(n)%nseg
              READ(inp,*) chain(n)%liste(i,:) 
              chain(n)%npair=-1 ! hm...
              chain(n)%xybeg=0.0d0
              chain(n)%xyend=0.0d0
            END DO 
          ELSE
            NULLIFY(chain(n)%liste) 
          ENDIF 
        END DO 
      END SELECT 
!
!----------------------------------------------------------------------- 
!
!      WRITE(lu,*) 'sections summary:'
!      WRITE(lu,*) 'nsec,ihowsec: ',nsec,ihowsec
!      SELECT CASE (ihowsec) 
!      CASE(:0) ! serial case, or "classical case" in parallel (devel) 
!        DO n=1,nsec
!          WRITE(lu,*) chain(n)%descr
!          WRITE(lu,*) chain(n)%xybeg(:), chain(n)%xyend(:)
!          WRITE(lu,*) chain(n)%npair(:)
!        END DO  
!      CASE (1:) ! partitioned, ready segment chains given 
!        DO n=1,nsec
!          WRITE(lu,*) 'name: ', chain(n)%descr
!          WRITE(lu,*) 'nseg: ', chain(n)%nseg
!          DO i=1,chain(n)%nseg
!            WRITE(lu,*) chain(n)%liste(i,:)
!          END DO 
!        END DO 
!      END SELECT
!
!-----------------------------------------------------------------------
! transfer to the global telemac or sisyphe variables
! NCP is 2 * number of sections
! CTRLSC is the list of terminal nodes of the sections 
! we have to re-allocate CTRLSC carefully
!
!      WRITE (lu,*) 'arranging sections for telemac'
!      WRITE (lu,*) 'telemac ncp was: ',ncp
      ncp = 2*nsec 
      IF (ALLOCATED(ctrlsc)) THEN 
        DEALLOCATE(ctrlsc, STAT=err)
        IF (err/=0) THEN
          WRITE(lu,*) 
     &    'READ_SECTIONS_TELEMAC2D: Error by deallocating CTRLSC:',err
          CALL plante(1)
          STOP 
        ENDIF 
      ENDIF 
      ALLOCATE (ctrlsc(ncp), STAT=err)
      IF (err/=0) THEN
        WRITE(lu,*) 
     &  'READ_SECTIONS_TELEMAC2D: Error by reallocating CTRLSC:',err
        CALL plante(1)
        STOP 
      ENDIF
      i=1
      DO n=1,nsec
        ctrlsc(i)   = chain(n)%npair(1)
        ctrlsc(i+1) = chain(n)%npair(2)
        i=i+2
      END DO 
!      WRITE (lu,*) 'ncp@telemac: ',ncp
!      WRITE (lu,*) 'ctrlsc@telemac: ',ctrlsc
!
!-----------------------------------------------------------------------
!
!      WRITE(lu,*) '-> leaving read_sections_telemac2d'
      RETURN 
      END SUBROUTINE READ_SECTIONS_TELEMAC2D
!-----------------------------------------------------------------------
