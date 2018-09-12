! 
!  Program that efectuates Cluster analysis in fortran
!  please cite: 
!  
!  1) Facco, E., d’Errico, M., Rodriguez, A., & Laio, A. (2017). Estimating the
!  intrinsic dimension of datasets by a minimal neighborhood information.
!  Scientific reports, 7(1), 12140.
!
!  2) Rodriguez, A. d'Errico, M., Facco, E. & Laio, A. (2018) Computing the free
!  energy without collective variables. JCTC (Accepted for publication) 
!
!  3) d'Errico, M., Facco, E., Laio, A. & Rodriguez, A. (2018) Automatic
!  topography of high-dimensional data sets by non-parametric Density Peak
!  clustering (Submitted)
! 
 
program dp_advance
use m_rnkpar
implicit none
Integer, Parameter :: kdp = selected_real_kind(15)
Integer, Parameter :: idp = selected_int_kind(15)
!! Error control
integer (kind=idp) :: id_err
integer :: allocate_status
!!Global variables
integer (kind=idp) :: Nele                   ! Number of elements
integer (kind=idp) :: ND                     ! Number of distances
real (kind=kdp),allocatable :: Dist(:)       ! Distance matrix
real (kind=kdp) :: dimset                    ! dimension in which is embedded the dataset (real)
real (kind=kdp), Parameter :: minr = -huge(dimset) ! Minimum real value not -inf
real (kind=kdp), Parameter :: maxr = huge(dimset)  ! Maximum real value not inf
real (kind=kdp), Parameter :: nearzero = tiny(dimset)  ! Minimum positive value
!integer (kind=idp) :: dimint                 ! Not used anymore
integer (kind=idp),parameter :: maxknn=int(1001,idp)  ! maximum number of neighbours to explore+1
integer (kind=idp),parameter :: minknn=int(4,idp)    ! minimum number of neighbours to explore+1
integer (kind=idp),parameter :: zero=int(0,idp) ! some casting
integer (kind=idp),parameter :: uno=int(1,idp)  ! some casting
integer (kind=idp),parameter :: dos=int(2,idp)  ! some casting
real (kind=kdp),parameter :: rzero=real(0.,kdp) ! some casting
integer (kind=idp) :: limit
integer (kind=idp) :: kknn
integer (kind=idp) :: numlimit               !number of points that arrived to the limit when exploring
real (kind=kdp),allocatable :: Rho(:)        ! Density
real (kind=kdp),allocatable :: Rho_err(:)    ! Density error
real (kind=kdp),allocatable :: Rho_prob(:)   ! Probability of having maximum Rho
real (kind=kdp),allocatable :: dc(:)         ! Distance for density calculation
integer (kind=idp),allocatable :: Nlist(:,:) ! Neighbour list within dc
integer (kind=idp),allocatable :: Nstar(:)   ! Number of neighbours taken for computing density
integer (kind=idp),allocatable :: Centers(:) ! Centers of the peaks
logical,allocatable :: centquest(:) ! Answer to the question "Is that element the center of a cluster?"
integer (kind=idp),allocatable :: Cluster(:) ! Cluster ID for the element
integer (kind=idp) :: Nclus                  ! Number of Cluster
integer (kind=idp),allocatable :: candidates_B(:)     ! Cluster of the Border candidates
real (kind=kdp),allocatable :: Bord(:,:)     ! Border Densities
real (kind=kdp),allocatable :: Bord_err(:,:) ! Border Densities Error
real (kind=kdp),allocatable :: Bord_s(:,:)     ! Border Densities saved
real (kind=kdp),allocatable :: Bord_err_s(:,:) ! Border Densities Error saved
integer (kind=idp),allocatable :: eb(:,:)    ! Border elements
integer (kind=idp),allocatable :: Pop(:)     ! Cluster Population
real (kind=kdp),allocatable :: cent(:)       ! Center Density
real (kind=kdp),allocatable :: cent_err(:)   ! Center Error
! Underscore m implies data after automatic merging
integer (kind=idp) :: Nclus_m                    ! Number of Cluster merged
real (kind=kdp),allocatable :: Bord_m(:,:)     ! Border Densities
real (kind=kdp),allocatable :: Bord_err_m(:,:) ! Border Densities Error
integer (kind=idp),allocatable :: Pop_m(:)       ! Cluster Population
real (kind=kdp),allocatable :: cent_m(:)       ! Center Density
real (kind=kdp),allocatable :: cent_err_m(:)   ! Center Error
integer (kind=idp),allocatable :: Centers_m(:)   ! Centers of the peaks
integer (kind=idp),allocatable :: Cluster_m(:)   ! Cluster ID for the element
! Halo characteristics
integer (kind=idp),allocatable :: halo_m(:)   ! Cluster ID for the element with halo cut
real (kind=kdp),allocatable :: halo_cut(:)    ! density value cut for halo
integer (kind=idp),allocatable :: Pop_halo(:)       ! Cluster Population without halo points
real (kind=kdp) :: Zmerge                 ! Prefactor multiplying the errors for cluster merging
real (kind=kdp) :: Zinit                  ! First value of the scan
real (kind=kdp) :: Zend                   ! Last value of the scan
integer (kind=idp) :: numint              ! Number of intermediate values in the scan
integer (kind=idp) :: interval
real (kind=kdp) :: Rho_sum            ! Factor to sum to the logarithm of the densities for in such a way that all are positives
integer (kind=idp) :: Pop_cut        ! Ignore clusters with population under this threshold
real (kind=kdp) :: Average_k          ! Average value of exit k
real (kind=kdp) :: thrpval            ! Threshold value for D in the exit k process
! 
integer (kind=idp):: n
integer :: standard_int
! Modify if you are working with zero distances
logical :: ignze                      ! substitute zero distances by nearzero
! Output files

open (51,file="cluster.log")

! Initial set-up

thrpval=23.92812698 ! p-value 10E-6
ignze=.false.

! Start computing !!! 

call get_distances (id_err)                  ! read or compute distances
write (51,*) "Number of elements:",Nele
if (id_err.ne.zero) call error_man (id_err)
call get_densities_and_dc(id_err)            ! from the distances, compute densities and dc
if (id_err.ne.zero) call error_man (id_err)
call clustering(id_err)                      ! get Clusters
if (id_err.ne.zero) call error_man (id_err)
call scan_merging(id_err)                         
if (id_err.ne.zero) call error_man (id_err)

close (51)

stop
call error_man (id_err)
stop


contains

subroutine error_man (id_err)
implicit none
!! Error control
integer (kind=idp) :: id_err
select case (id_err)
  case (1)
    write (6,*) "Select case error for distance type"
    call exit (1)
  case (2)
    write (6,*)"Error opening distance file"
    call exit (2)
  case (3)
    write (6,*)"Error opening coordinates file"
    call exit (3)
  case (4)
    write (6,*)"Error on distance file format"
    call exit (4)
  case (5)
    write (6,*)"Error on coordinates file format"
    call exit (5)
  case (6)
    write (6,*)"Peridic Bounday Conditions option unrecognized"
    call exit (6)
  case (7)
    write (6,*)"Dimension calculation option unrecognized"
    call exit (7)
  case (8) 
    write (6,*)"Density calculation option unrecognized"
    call exit (8)
  case (9) 
    write (6,*)"Just one cluster"
    stop
  case (10) 
    write (6,*)"Just one cluster after merging, stop scan"
    stop
  case (11) 
    write (6,*)"Error in assignation"
    call exit (11)
  case (12) 
    write (6,*)"Error with different points with zero distance"
    call exit (12)
  case (15) 
    write (6,*)"STOP: please, check for problems in the distances"
    call exit (15)
  case (16) 
    write (6,*)"proposed limit greater than the maximum hard limit"
    call exit (16)
  case (17) 
    write (6,*) "No entries in the input coordinates/distances"
    call exit (17)
  case (18) 
    write (6,*) "Out of memory when allocating"
    call exit (18)
  case default
  stop "unrecognized error"
end select
endsubroutine error_man

subroutine get_distances(id_err)
implicit none
integer (kind=idp) :: id_err
!! Local variables
integer (kind=idp) :: i,j,k,l,m,n
integer (kind=idp) :: Nlines
integer (kind=idp) :: dis_type
integer (kind=idp) :: dim_type
integer (kind=idp) :: pbc_type
integer (kind=idp) :: nvar
real (kind=kdp)    :: d
real (kind=kdp),allocatable    :: period (:)
character (len=500) :: filename
character (len=5000000):: line
real (kind=kdp),allocatable :: CRD(:,:) ! Coordinates
logical*1 :: blanko
!!
id_err=0
write (6,*) "Distance calculation type:"
write (6,*) "(1) Read from file"
write (6,*) "Read a symmetric distance file written as e1,e2,d(e1,e2). Indeces must start from 1"
write (6,*) "(2) Compute from coordinates file"
write (6,*) "Read a coordinates file, each row corresponds to a data point, each column to a coordinate"
write (6,*) "(3) Read from file with only maxk neighbors"
write (6,*) "Same format as distance file (type 1), but only the nearest neighbors are given"
read (5,*) dis_type
write (51,*) "Distance type:",dis_type

select case (dis_type)
  case (1) ! Read a symetric distance file written as e1,e2,d(e1,e2). Element indexes must start by 1.
    write (6,*) "distance file"
    read (5,'(a500)') filename
    write (51,*) "distance file:", trim(filename)
    open (11,file=filename,status="old",err=111)
    Nele=0
    Nlines=zero
    do 
      read (11,*,err=113,end=1110) i,j,d
      Nlines=Nlines+uno
      if (i.gt.Nele) Nele=i
      if (j.gt.Nele) Nele=j
    enddo
1110 continue
    rewind (11)
    if (Nele.le.zero) then
      id_err=17
      return
    endif
    ND=(Nele*Nele-Nele)/dos
    allocate (Dist(ND),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    do n=1,Nlines
      read (11,*) i,j,d
      if (i.ne.j) then
        l=max(i,j)
        m=min(i,j)
        k=(m-uno)*Nele-(m*m+m)/dos+l
        Dist(k)=d
        if (Dist(k).le.nearzero) then
          if (ignze) then
            Dist(k)=nearzero*real(1.1,kdp)
          else
            write (6,*) i,j,d
            id_err=12
            return
          endif
        endif
      endif
    enddo
    close (11)
    write (6,*) "Dimension:"
    write (6,*) "(1) Compute it"
    write (6,*) "(2) Give it as input"
    read (5,*) dim_type
    select case (dim_type)
      case (1)
        call get_dimension (id_err)
        if (id_err.ne.zero) return
      case (2)
        write (6,*) "give the dimension"
        read (5,*) dimset
      case default
        id_err=7
        return
    end select
  case (2)
    write (6,*) "coordinates file"
    read (5,'(a500)') filename
    write (51,*) "coordinates file:",trim(filename)
    open (11,file=filename,status="old",err=112)
    read (11,'(a5000000)',end=1120) line
    nvar=0
    l = len (trim(line))
    blanko=.true.
    do i=1,l
      if ( line(i:i) == ' ' ) then
        blanko = .true.
      else if ( blanko ) then
        nvar = nvar + uno
        blanko = .false.
      end if
    enddo
! get number of files (elements,nele) by reading the file and rewind
    rewind (11)
    Nele=0
    do
      read (11,'(a5000)',end=1120,err=1120) line
      Nele=Nele+uno
    enddo
1120 continue
    if (Nele.le.zero) then
      id_err=17
      return
    endif
    rewind (11)
    write (51,*) "Number of columns:",nvar
    allocate (CRD(Nele,nvar),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    do i=uno,Nele
      read (11,*,err=114) (CRD(i,j),j=uno,nvar)
    enddo
    close (11)
    ND=(Nele*Nele-Nele)/dos
    allocate (Dist(ND),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    write (6,*) "Do you want to use PBC (Periodic Boundary Conditions, for example if your features are angles)?"
    write (6,*) "(1) No"
    write (6,*) "(2) Yes"
    read (5,*) pbc_type
    select case (pbc_type)
    case (1)
      k=0
      do i=uno,Nele-uno
        do j=i+uno,Nele
          k=k+uno
          Dist(k)=0.
          do l=uno,nvar
            d=CRD(i,l)-CRD(j,l)
            Dist(k)=Dist(k)+d*d
          enddo
          Dist(k)=dsqrt(Dist(k))
          if (Dist(k).le.nearzero) then
            if (ignze) then
              Dist(k)=nearzero*real(1.1,kdp)
            else
              write (6,*) i,j,d
              id_err=12
              return
            endif
          endif
        enddo
      enddo
    case (2)
      allocate (period(nvar),stat=allocate_status)
      if (allocate_status /= 0) then
        id_err=18
        return
      endif  
      write (6,*) "Specify the period for each variable"
      read (5,*) period(:)
      write (6,*) "------------------------------------"
      write (51,*) "PERIOD:", period(:)
      k=0
      do i=uno,Nele-uno
        do j=i+uno,Nele
          k=k+uno
          Dist(k)=0.
          do l=uno,nvar
            d=CRD(i,l)-CRD(j,l)
            if (d.gt.period(l)/real(2.,kdp)) d=d-period(l)
            if (d.lt.-period(l)/real(2.,kdp)) d=d+period(l)
            Dist(k)=Dist(k)+d*d
          enddo
          Dist(k)=dsqrt(Dist(k))
          if (Dist(k).le.nearzero) then
            if (ignze) then
              Dist(k)=nearzero*real(1.1,kdp)
            else
              write (6,*) i,j,d
              id_err=12
              return
            endif
          endif
        enddo
      enddo
    case default
      id_err=6
      return
    end select 
    write (6,*) "Dimension:"
    write (6,*) "(1) Compute it"
    write (6,*) "(2) Give it as input"
    read (5,*) dim_type
    select case (dim_type)
      case (1)
        call get_dimension (id_err)
        if (id_err.ne.zero) return
      case (2)
        write (6,*) "give the dimension"
        read (5,*) dimset
      case default
        id_err=int(7,idp)
        return
    end select
  case (3) ! Read a neighbor distance file written as e1,e2,d(e1,e2)
    write (6,*) "neighbor distance file"
    read (5,'(a500)') filename
    write (51,*) "neighbor distance file:",trim(filename)
    open (11,file=filename,status="old",err=111)
    Nele=zero
    Nlines=zero
    do
      read (11,*,err=113,end=3110) i,j,d
      Nlines=Nlines+uno
      if (i.gt.Nele) Nele=i
      if (j.gt.Nele) Nele=j
    enddo
3110 continue
    rewind (11)
    write (6,*) "Number of elements is:",Nele
    if (Nele.le.zero) then
      id_err=17
      return
    endif
    ND=(Nele*Nele-Nele)/dos
    allocate (Dist(ND),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    do j=1,ND
      Dist(j)=(0.9d0+0.1d0*dfloat(j)/dfloat(ND))*maxr
    enddo
    do n=1,Nlines
      read (11,*) i,j,d
      if (i.ne.j) then
        l=max(i,j)
        m=min(i,j)
        k=(m-uno)*Nele-(m*m+m)/dos+l
        Dist(k)=d
        if (Dist(k).le.nearzero) then
          if (ignze) then
            Dist(k)=nearzero*real(1.1,kdp)
          else
            write (6,*) i,j,d
            id_err=12
            return
          endif
        endif
      endif
    enddo
    close (11)
    write (6,*) "Dimension:"
    write (6,*) "(1) Compute it"
    write (6,*) "(2) Give it as input"
    read (5,*) dim_type
    select case (dim_type)
      case (1)
        call get_dimension (id_err)
        if (id_err.ne.zero) return
      case (2)
        write (6,*) "give the dimension"
        read (5,*) dimset
      case default
        id_err=7
        return
    end select
  case default
    id_err=1
    return
end select 
write (51,*) "DIMENSION EMPLOYED:", dimset
return
111 id_err=2
return
112 id_err=3
return
113 id_err=4
return
114 id_err=5
return

end subroutine get_distances

subroutine get_dimension(id_err)
implicit none
integer (kind=idp) :: id_err
!! Local variables
integer (kind=idp) :: i,j,k,jj
real (kind=kdp) :: a,b,c,r2,r1,x,y
real (kind=kdp) :: fdim   ! fraction of discarded points for dimension calculation
real (kind=kdp),allocatable :: nu(:) !Ordered distances
real (kind=kdp) :: d_mle 
integer (kind=idp),allocatable :: i_nu(:) !original position of the ordered distances
integer (kind=idp),parameter :: nfrac=int(17,idp)
integer (kind=idp) :: Ndec,block
integer (kind=idp),allocatable :: decset(:),set(:)
integer (kind=idp),allocatable :: boxes(:)
real (kind=kdp) :: nnu,cdmle,cdmle2
real (kind=kdp) :: id_MLE ! dimension estimated with MLE version of 2-NN


id_err=zero

write (6,*) "which fraction of distances do you want to discard for dimension calculation?"
read (5,*) fdim 
fdim=real(1.,kdp)-fdim
!
!! Compute nu=r2/r1
allocate (nu(Nele),i_nu(Nele),decset(Nele),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
d_mle=0.
do i=1,Nele
  a=gDist(i,uno)
  b=gDist(i,dos)
  if (i.eq.uno) a=maxr
  if (i.eq.dos) b=maxr
  r1=min(a,b)
  r2=max(a,b)
  do j=int(3,idp),Nele
    if (j.ne.i) then
      c=gDist(i,j)
      a=r1
      b=min(r2,c)
      r1=min(a,b)
      r2=max(a,b)
    endif
  enddo
  nu(i)=r2/r1
  d_mle=d_mle+log(nu(i))
enddo
d_mle=dfloat(Nele)/d_mle
id_MLE=d_mle
!
!! compute decimation plot with MLE
!
write (6,*) "Do you want to perform block analysis (time consuming)?"
write (6,*) "(0) No "
write (6,*) "(1) Yes"
read (5,*) block
if (block.eq.1) then
  allocate (boxes(nfrac),stat=allocate_status)
  if (allocate_status /= 0) then
    id_err=18
    return
  endif  
  boxes(1)=2
  boxes(2)=3
  boxes(3)=4
  boxes(4)=5
  boxes(5)=6
  boxes(6)=7
  boxes(7)=8
  boxes(8)=9
  boxes(9)=10
  boxes(10)=11
  boxes(11)=12
  boxes(12)=13
  boxes(13)=14
  boxes(14)=15
  boxes(15)=16
  boxes(16)=18
  boxes(17)=20

  open (999,file="decimation_plot.dat")
  write (999,*) Nele,d_mle,d_mle/sqrt(dfloat(Nele))
  allocate (set(Nele),stat=allocate_status)
  if (allocate_status /= 0) then
    id_err=18
    return
  endif  
  do jj=1,nfrac
    Ndec=Nint(Nele/float(boxes(jj)))
    cdmle=rzero
    cdmle2=rzero
! write (6,*) "Number of elements:",Ndec
    do k=uno,int(50,idp)
!   write (6,*) "dec:",jj,".",k
      do i=1,Nele
        set(i)=i
      enddo
      call Shuffle(set)
      d_mle=0.
      do i=1,Ndec
        a=gDist(set(i),set(1))
        b=gDist(set(i),set(2))
        if (i.eq.uno) a=maxr
        if (i.eq.dos) b=maxr
        r1=min(a,b)
        r2=max(a,b)
        do j=3,Ndec
          if (j.ne.i) then
            c=gDist(set(i),set(j))
             a=r1
             b=min(r2,c)
             r1=min(a,b)
             r2=max(a,b)
          endif
        enddo
        nnu=r2/r1
        d_mle=d_mle+log(nnu)
      enddo
      d_mle=dfloat(Ndec)/d_mle
      cdmle=cdmle+d_mle
      cdmle2=cdmle2+d_mle*d_mle
    enddo
    cdmle=cdmle/real(50.,kdp)
    cdmle2=sqrt(cdmle2/real(50.,kdp)-cdmle*cdmle)
    write (999,*) Ndec,cdmle,cdmle2
  enddo
  write (999,*) 1,0,0
  deallocate (set)
  close (999)
  open (999,file="decimation_plot.gpl")
  write (999,*) "pl '<sort -gk1 decimation_plot.dat' u 1:2:3 w e ps 2,'' u 1:2 w l"
  close (999)
!call execute_command_line("gnuplot -persists decimation_plot.gpl")
  deallocate (boxes)
endif
!! Sort nu and compute F(nu)=i/Nele
call HPSORT(Nele,nu,i_nu)
!! perform fit -log[1-F(nu)]=d*log(nu)

a=rzero
b=rzero
i=uno
open (21,file="xy_dim.dat")
open (23,file="xy_dim_all.dat")
open (22,file="xy_dim.gpl")
do while ((dfloat(i)/dfloat(Nele)).le.fdim)
  x=dlog(nu(i))
  y=-dlog(real(1.,kdp)-(dfloat(i)/dfloat(Nele)))
  write (21,*) x,y
  write (23,*) x,y
  a=a+x*y
  b=b+x*x
  i=i+uno
enddo
dimset=a/b
!write (6,*) real(Nele,kdp)
do j=i,Nele-1
  x=dlog(nu(j))
  y=-dlog(real(1.,kdp)-(real(j,kdp)/real(Nele,kdp)))
  write (23,*) x,y
enddo
write (22,'(a,f8.3,a)') "pl 'xy_dim_all.dat' pt 7 lc rgb 'black' t '','xy_dim.dat' pt 7 lc rgb 'blue' ps 2 t '',",dimset,"*x lc rgb 'blue' lw 2 t ''" 
close (22)
close (21)
!call execute_command_line("gnuplot -persists xy_dim.gpl")
write (6,*) "The relevant information for computing the dimension is in"
write (6,*) "files: xy_dim_all.dat &  decimation_plot.dat (if you performed block analysis)"
write (6,*) "We strongly recomend to plot them before choosing the dimension"
write (6,*) "If you use gnuplot, we suggest you to employ the scripts"
write (6,*) "xy_dim.gpl and decimation_plot.gpl generated by the code"
write (6,*) "Dimension estimated from the linear fit:",dimset
write (6,*) "MLE of d: ",id_MLE," with stdv ",id_MLE/sqrt(dfloat(Nele))
dimset=id_MLE
write (6,*) "Provide the dimension"
read (5,*) dimset
deallocate (nu,i_nu)
return
end subroutine get_dimension

subroutine get_densities_and_dc(id_err)
implicit none
integer (kind=idp) :: id_err
!! Local variables
integer (kind=idp) :: i,j,k,m
real (kind=kdp) :: ms,ns
integer (kind=idp) :: niter
real (kind=kdp),allocatable :: Vols(:)
real (kind=kdp),allocatable :: VV(:,:)
integer ,allocatable :: iVols(:)
real (kind=kdp), parameter :: pi=3.14159265359D0
real (kind=kdp) :: prefactor
real (kind=kdp) :: dL
real (kind=kdp),allocatable :: vi(:)
logical :: viol
real (kind=kdp) :: a,b
integer (kind=idp) :: nviol,ksel
real (kind=kdp) :: Cov2(2,2),Covinv2(2,2)
real (kind=kdp) :: L0,func,sigma,t,jf
real (kind=kdp) :: tt,gb,ga,sa,sb,s,a_err
real (kind=kdp) :: stepmax !! This variable controls the maximum variation on the log(rho) accepted during the optimization process
integer (kind=idp) :: den_type
! Variable that accounts for vi=0 --> then use kNN instead for density (error
! still PAk)
logical :: kNN

id_err=0

Rho_sum=0
  
write (6,*) "Which density estimator do you want to use?"
write (6,*) "(1) Pak"
write (6,*) "(2) k-NN (error 1./sqrt(k))"
read (5,*) den_type
write (51,*) "Density type:", den_type

allocate (Rho(Nele),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
if (allocate_status /= 0) then
  id_err=18
  return
endif  
allocate (Rho_err(Nele),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  

select case (den_type)
  case (1)
    write (51,*) "p-value for computing the density set to 1E-6, threshold set to 23.92812698"
    thrpval=23.92812698
    limit=min(maxknn,int(Nele/dos))
! get prefactor for Volume calculation
!   if (mod(dimint,dos).eq.zero) then
!     k=dimint/dos
!     m=uno
!     do i=uno,k
!       m=m*i
!     enddo
!     prefactor=pi**k/(dfloat(m))
!   else
!     k=(dimint-uno)/dos
!     ms=rzero
!     do i=uno,k
!       ms=ms+dlog(dfloat(i))
!     enddo
!     ns=ms
!     do i=k+uno,dimint
!       ns=ns+dlog(dfloat(i))
!     enddo
!     prefactor=real(2.,kdp)*dexp(ms-ns+k*dlog(real(4.,kdp)*pi))
!   endif
!   write (6,*) "old style W: ", prefactor
!   write (6,*) "new W :",
    prefactor=dexp(dimset/2.*dlog(pi)-log_gamma((dimset+2)/2.))
    allocate (Nstar(Nele),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    allocate (Nlist(Nele,limit),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    allocate (VV(Nele,limit),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    allocate (Vols(Nele),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    allocate (iVols(Nele),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    do i=uno,Nele
      Vols(:)=maxr
      do j=uno,Nele
        if (i.ne.j) then
          Vols(j)=gDist(i,j)
        endif
      enddo
      standard_int=limit
      call rnkpar(Vols,iVols,standard_int)
      do j=uno,limit
        Nlist(i,j)=iVols(j)
        VV(i,j)=prefactor*dexp(dimset*dlog(Vols(iVols(j))))
      enddo
    enddo
    deallocate (Vols,iVols)
    Rho_sum=maxr
    numlimit=zero
    Average_k=rzero
    do i=uno,Nele
      viol=.false.
      k=minknn
      nviol=0
      do while (.not.viol)
        ksel=k-uno
! Alternatively it can be defined as a function of r_i/r_j, the dimension and
! ksel: dL(x,ksel,d)=2*ksel*(log(x**d+x**(-d)+2)-log(4)) where x=r_i/r_j
        dL=real(-2.,kdp)*(float(ksel))*(log(VV(i,ksel))+log(VV(Nlist(i,k),ksel))-real(2.,kdp)*log(VV(i,ksel)+VV(Nlist(i,k),ksel))+log(real(4.,kdp)))
        k=k+uno
        if (dL.gt.thrpval) viol=.true. 
        if (k.gt.limit) then
          viol=.true.
          numlimit=numlimit+uno
        endif
      enddo
      Nstar(i)=k-dos
      Average_k=Average_k+dfloat(Nstar(i))
    enddo

    do i=uno,Nele
      allocate (vi(Nstar(i)),stat=allocate_status)
      if (allocate_status /= 0) then
        id_err=18
        return
      endif  
      Vi(uno)=VV(i,uno)
      kNN=.false.
      do j=dos,Nstar(i)
        Vi(j)=VV(i,j)-VV(i,j-uno)
        if ((Vi(j).lt.nearzero).and.(.not.kNN)) then
          write (51,*) "Warning: Pak cannot be employed at point",i
          write (51,*) "due to a shell volume equal to 0: Employing kNN"
          kNN=.true.
        endif
      enddo
      b=log(float(Nstar(i))/VV(i,Nstar(i)))
      L0=rzero
      a=rzero
      if (.not.kNN) then
        stepmax=real(0.1,kdp)*abs(b)
        gb=float(Nstar(i))
        ga=float(Nstar(i)+uno)*float(Nstar(i))/real(2.,kdp)
        Cov2(1,1)=rzero !gbb
        Cov2(1,2)=rzero !gab
        Cov2(2,2)=rzero !gaa
        do j=uno,Nstar(i)
          jf=float(j)
          t=b+a*jf
          s=exp(t)
          tt=vi(j)*s
          L0=L0+t-tt
          gb=gb-tt
          ga=ga-jf*tt
          Cov2(1,1)=Cov2(1,1)-tt
          Cov2(1,2)=Cov2(1,2)-jf*tt
          Cov2(2,2)=Cov2(2,2)-jf*jf*tt
        enddo
        Cov2(2,1)=Cov2(1,2)
        Covinv2=matinv2(Cov2)
        func=real(100.,kdp)
        niter=zero
        do while ((func>1D-3).and.(niter.lt.int(1000,idp)))
          sb=(Covinv2(1,1)*gb+Covinv2(1,2)*ga)
          sa=(Covinv2(2,1)*gb+Covinv2(2,2)*ga)
          niter=niter+uno
          sigma=real(0.1,kdp)
          if (abs(sigma*sb).gt.stepmax) then
            sigma=abs(stepmax/sb)
          endif
          b=b-sigma*sb
          a=a-sigma*sa
          L0=rzero
          gb=float(Nstar(i))
          ga=float(Nstar(i)+uno)*float(Nstar(i))/real(2.,kdp)
          Cov2(1,1)=rzero !gbb
          Cov2(1,2)=rzero !gab
          Cov2(2,2)=rzero !gaa
          do j=uno,Nstar(i)
            jf=float(j)
            t=b+a*jf
            s=exp(t)
            tt=vi(j)*s
            L0=L0+t-tt
            gb=gb-tt
            ga=ga-jf*tt
            Cov2(1,1)=Cov2(1,1)-tt
            Cov2(1,2)=Cov2(1,2)-jf*tt
            Cov2(2,2)=Cov2(2,2)-jf*jf*tt
          enddo
          Cov2(2,1)=Cov2(1,2)
          Covinv2=matinv2(Cov2)
          if ((abs(a).le.tiny(a)).or.(abs(b).le.tiny(b))) then
            func=max(gb,ga)
          else
            func=max(abs(gb/b),abs(ga/a))
          endif
        enddo
        Cov2(:,:)=-Cov2(:,:)
        Covinv2=matinv2(Cov2)
        a_err=dsqrt(Covinv2(2,2))
      endif
      Rho(i)=b
      if (ISNAN(Rho(i))) then
        write (6,*) "Density NaN at point",i 
        id_err=int(15,idp)
        return
      endif
      if ((Rho(i).gt.maxr).or.(Rho(i).lt.minr)) then
        write (6,*) "Density at point",i, Rho(i)
        id_err=int(15,idp)
        return
      endif
      Rho_err(i)=dsqrt(dfloat(int(4,idp)*Nstar(i)+dos)/dfloat(Nstar(i)*(Nstar(i)-uno)))
      if (Rho(i).lt.Rho_sum) Rho_sum=Rho(i)
      deallocate (vi)
    enddo
    write (6,*) "Number of points arrived to kmax:",numlimit
    write (51,*) "Number of points arrived to kmax:",numlimit
    if (float(numlimit)/float(Nele).gt.0.05) then
      write (6,*) "Warning, exceeds 5 %"
      write (6,*) "Consider changing maxknn"
      write (51,*) "Warning, exceeds 5 %"
      write (51,*) "Consider changing maxknn"
    endif
    deallocate (VV)
    write (51,*) "Added term to log(Rho)=",-Rho_sum+real(1,kdp)
    write (51,*) "Average exit k=",Average_k/dfloat(Nele)
    Rho(:)=Rho(:)-Rho_sum+real(1,kdp) !!! In this way the minimum Rho is equal to 1. (for hierarchies it does not matter so much)
! Once computed the densities we can emplot the Nstar for topography
    
!Get dc
    allocate (dc(Nele),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    do i=1,Nele
      j=Nlist(i,Nstar(i))
      dc(i)=gDist(i,j)
    enddo
  case (2)
    write (6,*) "which k?"
    read (5,*) kknn
    limit=kknn
    write (51,*) "k-NN, k set to:",limit
! get prefactor for Volume calculation
    prefactor=dexp(dimset/2.*dlog(pi)-log_gamma((dimset+2)/2.))

    allocate (Nstar(Nele),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    allocate (Nlist(Nele,limit),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    allocate (Vols(Nele),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    allocate (iVols(Nele),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    Rho_sum=maxr
    prefactor=log(float(kknn))-log(prefactor)
    do i=1,Nele
      Vols(:)=maxr
      do j=1,Nele
        if (i.ne.j) then
          Vols(j)=gDist(i,j)
        endif
      enddo
      standard_int=limit
      call rnkpar(Vols,iVols,standard_int)
      Nstar(i)=kknn
      do j=1,limit
        Nlist(i,j)=iVols(j)
      enddo
      Rho(i)=prefactor-dimset*log(Vols(iVols(kknn)))
      Rho_err(i)=real(1.,kdp)/sqrt(float(kknn))
      if (Rho(i).lt.Rho_sum) Rho_sum=Rho(i)
    enddo
    deallocate (Vols,iVols)
    Rho(:)=Rho(:)-Rho_sum+real(1.,kdp) !!! In this way the minimum Rho is equal to 1. (for hierarchies it does not matter so much)
    allocate (dc(Nele),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    do i=1,Nele
      j=Nlist(i,Nstar(i))
      dc(i)=gDist(i,j)
    enddo
  case default
    id_err=int(8,idp)
    return
end select
return
end subroutine get_densities_and_dc

pure function matinv2(A) result(B)
!! Performs a direct calculation of the inverse of a 2×2 matrix.
  real (kind=kdp), intent(in) :: A(2,2)   !! Matrix
  real (kind=kdp)             :: B(2,2)   !! Inverse matrix
  real (kind=kdp)             :: detinv
  ! Calculate the inverse determinant of the matrix
  detinv = real(1.,kdp)/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
  ! Calculate the inverse of the matrix
  B(1,1) = +detinv * A(2,2)
  B(2,1) = -detinv * A(2,1)
  B(1,2) = -detinv * A(1,2)
  B(2,2) = +detinv * A(1,1)
end function

subroutine clustering(id_err)
implicit none
integer (kind=idp) :: id_err
!! Local variables
integer (kind=idp) :: i,j,k
integer (kind=idp) :: ig,jg
integer (kind=idp) :: nassign
integer (kind=idp) :: iref
integer (kind=idp) :: l
logical :: idmax,extend
real (kind=kdp),allocatable :: Rho_copy(:)
integer (kind=idp),allocatable :: iRho(:)
integer (kind=idp),allocatable :: ordRho(:) 
real (kind=kdp) :: d,dmin
!
! Putative centers
integer (kind=idp),allocatable:: putative (:)
!!
id_err=0
allocate (putative(Nele),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
allocate (Cluster(Nele),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
allocate (candidates_B(Nele),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
! centers are points with maximum density in their neighborhood and do not
! belong to the neighborhood of any other point with higher density
putative(:)=0
Cluster(:)=0
Nclus=0
allocate (Rho_prob(Nele),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
Rho_prob(:)=0.
do i=1,Nele
  Rho_prob(i)=Rho(i)-Rho_err(i)
enddo
! copy of rho (not efficient, but it adds clarity to the code)

allocate (Rho_copy(Nele),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
allocate (iRho(Nele),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
Rho_copy(:)=-Rho_prob(:)
call HPSORT(Nele,Rho_copy,iRho) ! iRho contains the order in density (iRho(1) is the element with highest Rho...)
deallocate (Rho_copy)
allocate (ordRho(Nele),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
do i=1,Nele
  ordRho(iRho(i))=i                 ! ordRho is the complementary of iRho. Given an element, ordRho returns its order in density
enddo
do i=uno,Nele
  idmax=.true.
  j=uno
  do while (idmax .and. (j.le.Nstar(i)))
    if ((ordRho(i).gt.ordRho(Nlist(i,j))))  idmax=.false.
    j=j+uno
  enddo
  if (idmax) then
    putative(i)=uno
  endif
enddo
do i=uno,Nele
  if (putative(i).eq.uno)then
    idmax=.true. 
    j=uno
    do while  (idmax .and. (j.lt.ordRho(i))) 
      k=uno
      l=iRho(j)
      do while (idmax.and.(k.le.Nstar(l)))
        if (Nlist(l,k).eq.i) idmax=.false.
        k=k+uno
      enddo
      j=j+uno
    enddo
    if (idmax) then
      Nclus=Nclus+uno
      Cluster(i)=Nclus
    endif
  endif
enddo
write (51,*) "Number of clusters before automatic merging:",Nclus
allocate (Centers(Nclus),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
allocate (centquest(Nele),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
deallocate (putative)
centquest=.false.
do i=uno,Nele
  n=zero
  if (Cluster(i).ne.zero) then
    n=uno
    Centers(Cluster(i))=i
    centquest(i)=.true.
  endif
enddo 
nassign=Nclus
if (Nclus.gt.uno) then

write (51,*) "Performing assignation"
do i=uno,Nele
  j=iRho(i)
  if (Cluster(j).eq.zero) then
    dmin=maxr
    ig=int(-1,idp)
    do k=uno,i-uno
      l=iRho(k)
      if (gDist(j,l).le.dmin) then
        ig=l
        dmin=gDist(j,l)
      endif
    enddo
    if (ig.eq.int(-1,idp)) then
      do k=uno,i-uno
        l=iRho(k)
        if (gDist(j,l).le.dmin) then
          ig=l
          dmin=gDist(j,l)
        endif
      enddo
      stop
    endif
    Cluster(j)=Cluster(ig)
  endif
enddo
! count cluster population

allocate (Pop(Nclus),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
Pop(:)=zero
do i=uno,Nele
  Pop(Cluster(i))=Pop(Cluster(i))+uno
enddo

! find border densities

allocate (Bord(Nclus,Nclus),Bord_err(Nclus,Nclus),eb(Nclus,Nclus),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
Bord(:,:)=real(-1.,kdp)
Bord_err(:,:)=rzero
eb(:,:)=zero
candidates_B(:)=zero

do i=uno,Nele
  if (.not.centquest(i)) then
    dmin=maxr
    do j=uno,Nele
      if (cluster(j).ne.cluster(i)) then
        d=gDist(i,j)
        if (d.lt.dmin) then
          dmin=d
          jg=j
        endif
      endif
    enddo
    if (dmin.le.dc(i)) then
      iref=i
      k=zero
      extend=.true.
!    do while ( (k.lt.Nstar(i)).and.extend)
!      k=k+1
!      if (cluster(Nlist(i,k)).eq.cluster(i)) then
!        if (gDist(Nlist(i,k),jg).lt.dmin) extend=.false.
!      endif
!    enddo
      do while ( (k.lt.Nele).and.extend)
        k=k+uno
        if (cluster(k).eq.cluster(i)) then
          if (gDist(k,jg).lt.dmin) extend=.false.
        endif
      enddo
      if (extend) then
        candidates_B(i)=cluster(jg)
        if (Rho_prob(iref).gt. Bord(cluster(i),cluster(jg))) then
          Bord(cluster(i),cluster(jg))=Rho_prob(iref)
          Bord(cluster(jg),cluster(i))=Rho_prob(iref)
          Bord_err(cluster(i),cluster(jg))=Rho_err(iref)
          Bord_err(cluster(jg),cluster(i))=Rho_err(iref)
          eb(cluster(i),cluster(jg))=iref
          eb(cluster(jg),cluster(i))=iref
        endif
      endif
    endif
  endif
enddo
do i=uno,Nclus-uno
  do j=i+uno,Nclus
    if (eb(i,j).ne.zero) then
      Bord(i,j)=Rho(eb(i,j))
      Bord(j,i)=Rho(eb(i,j))
    else
      Bord(i,j)=rzero
      Bord(j,i)=rzero
    endif
  enddo
enddo

deallocate (Rho_prob)
deallocate (iRho)
deallocate (ordRho)
! Info per graph pre automatic merging
allocate (cent(Nclus),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
allocate (cent_err(Nclus),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
do i=uno,Nclus
  cent(i)=Rho(Centers(i))
  cent_err(i)=Rho_err(Centers(i))
enddo
else
  Cluster(:)=uno
  id_err=9
endif
return
end subroutine clustering

! This function allows to get the distances in a matrix like style

real (kind=kdp) function gDist(i,j)
integer (kind=idp) :: i,j,k,l,m
if (i.ne.j) then
  l=max(i,j)
  m=min(i,j)
  k=(m-uno)*Nele-(m*m+m)/dos+l
  gDist=Dist(k)
else
  gDist=rzero
endif
return
end function gDist
!
! 
subroutine scan_merging(id_err)
implicit none
integer (kind=idp) :: id_err
!! Local variables
integer (kind=idp) :: i,j,k,n,l,alive,dead,niter
logical :: Survive (Nclus)
logical :: change
integer (kind=idp) :: Nbarr
integer (kind=idp),allocatable :: Bcorr (:,:)
real (kind=kdp),allocatable :: Barrier (:)
real (kind=kdp),allocatable :: Barrier_err (:)
integer (kind=idp),allocatable :: iBarrier(:)
real (kind=kdp) :: c1,c2,b12
integer (kind=idp),allocatable :: M2O(:) !conversion from merged to original cluster number
integer (kind=idp) :: O2M(Nclus) ! Conversion from original cluster number to its equivalent in merged
real (kind=kdp) :: f1,f2,f12,b1,b2,sens_kin
character (len=3) :: cinter
!
!
!

id_err=zero
Nbarr=(Nclus*Nclus-Nclus)/dos

allocate (Bord_s(Nclus,Nclus),Bord_err_s(Nclus,Nclus),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
Bord_s=Bord
Bord_err_s=Bord_err

allocate (Barrier(Nbarr),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
allocate (Barrier_err(Nbarr),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
allocate (iBarrier(Nbarr),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
allocate (Bcorr(Nbarr,2),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
allocate (Cluster_m(Nele),stat=allocate_status)
if (allocate_status /= 0) then
  id_err=18
  return
endif  
write (6,*) "Provide the initial (min) Z parameter"
read (5,*) Zinit
write (6,*) "Provide the last (max) Z parameter"
read (5,*) Zend
write (6,*) "Provide the number of intermediate values"
read (5,*) numint
do interval=1,numint+2
  Bord=Bord_s
  Bord_err=Bord_err_s
  n=zero
  do i=uno,Nclus-uno
    do j=i+uno,Nclus
      n=n+uno
      Barrier(n)= Bord(i,j)
      Barrier_err(n)= Bord_err(i,j)
      Bcorr(n,1)=i
      Bcorr(n,2)=j
    enddo
  enddo
  if (Nbarr.gt.uno) then
    call HPSORT(Nbarr,Barrier,iBarrier)
    do n=uno,Nbarr
      k=iBarrier(n)
      i=Bcorr(k,1)
      j=Bcorr(k,2)
    enddo
  else
    iBarrier(1)=uno
  endif
  write (cinter,'(i3)') interval
  Zmerge=Zinit+(Zend-Zinit)/dfloat(numint+1)*dfloat(interval-1)
  write (51,*) "INTERVAL:",interval
  write (51,*) "Merging Z set to:",Zmerge
  sens_kin=rzero
  Survive(:)=.true.
  change=.true.
  Cluster_m(:)=Cluster(:)
  niter=0
  if (interval.gt.1) then
    deallocate (Pop_halo, halo_m, halo_cut, Centers_m, cent_m,cent_err_m, Bord_m,Bord_err_m, Pop_m, M2O)
  endif
  do while (change)
    niter=niter+uno
    change=.false.
  mdo:  do n=Nbarr,uno,int(-1,idp)
      k=iBarrier(n)
      i=Bcorr(k,uno)
      j=Bcorr(k,dos)
      if ((Bord(i,j).gt.rzero).and.(i.ne.j)) then
        if (Survive(i).and.Survive(j)) then
          c1=(cent(i)-Bord(i,j))/(Bord_err(i,j)+cent_err(i))
          c2=(cent(j)-Bord(i,j))/(Bord_err(i,j)+cent_err(j))
          b12=min(c1,c2)
          f1=cent(i)
          f2=cent(j)
          f12=Bord(i,j)
          b1=abs(f12-f1)
          b2=abs(f12-f2)
          if ((b12.lt.Zmerge).or.(b1.lt.sens_kin).or.(b2.lt.sens_kin)) then
            change=.true.
            Bord(i,j)=rzero
            Bord_err(i,j)=rzero
            if (c1.gt.c2) then
              alive=i
              dead=j
            else
              alive=j
              dead=i
            endif
            Bord(alive,alive)=rzero
            Bord_err(alive,alive)=rzero
            Survive(dead)=.false.
            do k=1,Nclus
              if (Survive(k)) then
                if (Bord(i,k).gt.Bord(j,k)) then
                  Bord(alive,k)=Bord(i,k)
                  Bord(k,alive)=Bord(k,i)
                  Bord_err(alive,k)=Bord_err(i,k)
                  Bord_err(k,alive)=Bord_err(k,i)
                else
                  Bord(alive,k)=Bord(j,k)
                  Bord(k,alive)=Bord(k,j)
                  Bord_err(alive,k)=Bord_err(j,k)
                  Bord_err(k,alive)=Bord_err(k,j)
                endif
              endif
            enddo
            do l=uno,Nele
              if (Cluster_m(l).eq.dead) Cluster_m(l)=alive
            enddo
            if (n.gt.uno) then
              do l=n-uno,uno,int(-1,idp)
                k=iBarrier(l)
                if (Bcorr(k,1).eq.dead) Bcorr(k,1)=alive
                if (Bcorr(k,2).eq.dead) Bcorr(k,2)=alive
              enddo
            endif
            exit mdo
          endif
        endif
      endif
    enddo mdo
! write (96,*) "MERGING ITER:",niter
  enddo


! get dictionary

  Nclus_m=zero
  do i=uno,Nclus
    if (Survive(i)) Nclus_m=Nclus_m+uno
  enddo
  write (51,*) "TOTAL CLUSTERS AFTER MERGING:",Nclus_m
  allocate (M2O(Nclus_m),stat=allocate_status)
  if (allocate_status /= 0) then
    id_err=18
    return
  endif  
  n=0
  O2M(:)=int(-1,idp)
  do i=uno,Nclus
    if (Survive(i)) then
      n=n+uno
      M2O(n)=i
      O2M(i)=n
    endif
  enddo

! get survival characteristics

  if (Nclus_m.gt.uno) then
    allocate (Pop_m(Nclus_m),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    allocate (Bord_m(Nclus_m,Nclus_m),Bord_err_m(Nclus_m,Nclus_m),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    allocate (cent_m(Nclus_m),cent_err_m(Nclus_m),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    allocate (Centers_m(Nclus_m),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
!   write (6,*) "Population threshold to ignore clusters in the dendrogram"
!   read (5,*) Pop_cut
!   write (51,*) "Clusters automatically ignored in the dendrogram if their population is under:",Pop_cut
!
    Pop_cut=zero
    Pop_m(:)=zero
    do i=uno,Nele
      Cluster_m(i)=O2M(Cluster_m(i))
      Pop_m(Cluster_m(i))=Pop_m(Cluster_m(i))+uno
    enddo
    do i=uno,Nclus_m
      do j=i+uno,Nclus_m
        Bord_m(i,j)=Bord(M2O(i),M2O(j)) 
        Bord_err_m(i,j)=Bord_err(M2O(i),M2O(j)) 
        Bord_m(j,i)=Bord(M2O(i),M2O(j)) 
        Bord_err_m(j,i)=Bord_err(M2O(i),M2O(j)) 
      enddo
      cent_m(i)=cent(M2O(i))
      cent_err_m(i)=cent_err(M2O(i))
      Centers_m(i)=Centers(M2O(i))
    enddo
! Compute halo
    allocate (halo_cut(Nclus_m),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    allocate (halo_m(Nele),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    allocate (Pop_halo(Nclus_m),stat=allocate_status)
    if (allocate_status /= 0) then
      id_err=18
      return
    endif  
    halo_cut(:)=maxval(Bord_m(:,:),dim=2)
    Pop_halo(:)=Pop_m(:)
    do i=uno,Nele
      halo_m(i)=Cluster_m(i)
      if (Rho(i).lt.halo_cut(Cluster_m(i))) then
        halo_m(i)=zero
        Pop_halo(Cluster_m(i))=Pop_halo(Cluster_m(i))-uno
      endif
    enddo
    open (21,file=trim("Point.info."//adjustl(cinter)))
    write (21,*) "# index kopt dc -F_i Err_{F_i} Cluster_i Halo_i"
    do i=uno,Nele
      write (21,'(i6,1x,i6,1x,3(es18.10,1x),i6,1x,i6)') i,Nstar(i),dc(i),Rho(i),Rho_err(i),Cluster_m(i),halo_m(i)
    enddo
    close (21)
    open (21,file=trim("Topography.info."//adjustl(cinter)))
    write (21,*) "# CENTERS"
    write (21,*) "# Cluster_index -F_center Error_{F_center} index_center Cluster_population Cluster_population_no_halo"
    do i=uno,Nclus_m
      write (21,'(i7,1x,2(es18.10,1x),i7,1x,i7,1x,i7)') i,cent_m(i),cent_err_m(i),Centers_m(i),Pop_m(i),Pop_halo(i)
    enddo
    write (21,*) "# SADDLES:"
    write (21,*) "# C1 C2 -F_bord_C1_C2 Err_F_bord_1_2"
    do i=uno,Nclus_m-uno
      do j=i+uno,Nclus_m
        write (21,'(i7,1x,i7,1x,2(es18.10,1x))') i,j,Bord_m(i,j),Bord_err_m(i,j)
      enddo
    enddo
    close (21)
  else
    id_err=int(10,idp)
    Cluster_m(:)=uno
    return
  endif
enddo

return
end subroutine 
!
! GRAPH:
!
! Build dendrogram and topography matrix
!
 subroutine graph2 (id_err)
 implicit none
! dummy
 integer (kind=idp) i,j
 integer (kind=idp) i1,i2
 integer (kind=idp) l1,l2,lp,lm1,lm2,l
 integer (kind=idp) id_err
 integer (kind=idp) Nlist
 integer (kind=idp) Nassign
 integer (kind=idp),dimension(Nclus_m,Nclus_m) :: peak_list
 integer (kind=idp),dimension(Nclus_m,Nclus_m,2) :: border_list
 integer (kind=idp),dimension(Nclus_m) :: N_peak_list
 integer (kind=idp),dimension(Nclus_m) :: list_peak
 integer (kind=idp),dimension(Nclus_m) :: Densort ! clusters sorted accordig with the dendrogram
 real (kind=kdp),dimension(Nclus_m,Nclus_m):: border_copy
 real (kind=kdp) :: xmax
 real (kind=kdp) :: b,eb
 real (kind=kdp),dimension(dos*Nclus_m-uno,2):: dendro_points
 logical,dimension(dos*Nclus_m-uno):: dendro_active_points
 logical,dimension(dos*Nclus_m-uno):: dendro_barrier
 real (kind=kdp),dimension(dos*Nclus_m-uno):: dendro_error
 real (kind=kdp),dimension(dos*Nclus_m-uno):: dendro_density
 real (kind=kdp),dimension(dos*Nclus_m-uno):: dendro_width
 character*3,dimension(dos*Nclus_m-uno):: dendro_label
 integer (kind=idp) :: dendro_active_barrier
 integer (kind=idp) :: dendro_n_active_points
 integer (kind=idp) :: dendro_left_nearest
 integer (kind=idp) :: dendro_right_nearest
 integer (kind=idp) :: Narrow
 integer (kind=idp) :: k
 real (kind=kdp) :: width
!
 real (kind=kdp) :: dminleft
 real (kind=kdp) :: dminright

 real (kind=kdp) :: last_barrier       ! Last node in the dendrogram
 real (kind=kdp) :: maxdens            ! Maximum density value in the dataset
 real (kind=kdp) :: hdendro


 id_err=zero
 border_copy(:,:)=Bord_m(:,:)
 peak_list(:,:)=zero
 border_list(:,:,:)=zero
 N_peak_list(:)=zero
 xmax=real(-1.,kdp)
 list_peak(:)=zero

! Find the highest border density and the center with higher density between the
! two defining this border
 do i=uno,Nclus_m-uno
   do j=i+uno,Nclus_m
     if (border_copy(i,j).ge.xmax) then
       xmax=border_copy(i,j)
       i1=i
       i2=j
     endif
   enddo
 enddo
 if (cent_m(i1).gt.cent_m(i2)) then
   l1=i1
   l2=i2
 else
   l1=i2
   l2=i1
 endif

! Generate a list of clusters
 Nlist=uno
 peak_list(Nlist,1)=l1
 peak_list(Nlist,2)=l2
 N_peak_list(Nlist)=dos
 border_list(Nlist,1,1)=l1
 border_list(Nlist,1,2)=l2
 border_copy(l1,l2)=real(-0.5,kdp)
 border_copy(l2,l1)=real(-0.5,kdp)
 list_peak(l1)=Nlist
 list_peak(l2)=Nlist

 Nassign=dos
 
 do while (N_peak_list(1).lt.Nclus_m)
! Find next barrier and highest max corresponding to this barrier
   xmax=real(-1.,kdp)
   do i=uno,Nclus_m-uno
     do j=i+uno,Nclus_m
       if (border_copy(i,j).ge.xmax) then
         xmax=border_copy(i,j)
         i1=i
         i2=j
       endif
     enddo
   enddo
   if (cent_m(i1).gt.cent_m(i2)) then
     l1=i1
     l2=i2
   else
     l1=i2
     l2=i1
   endif
! Control if none of the two maxima has been assigned to a list
   if (list_peak(l1).eq.zero) then
     if (list_peak(l2).eq.zero) then
!     in this case: New list
       border_copy(l1,l2)=real(-0.5,kdp)
       border_copy(l2,l1)=real(-0.5,kdp)
       Nlist=Nlist+uno
       peak_list(Nlist,1)=l1
       peak_list(Nlist,2)=l2
       N_peak_list(Nlist)=dos
       border_list(Nlist,1,1)=l1
       border_list(Nlist,1,2)=l2
       list_peak(l1)=Nlist
       list_peak(l2)=Nlist
     else
! If only l2 has been assigned to a list: Add l1 to list_peak(l2)==lp
       lp=list_peak(l2)
       do i=uno,N_peak_list(lp)
         border_copy(l1,peak_list(lp,i))=real(-0.5,kdp)
         border_copy(peak_list(lp,i),l1)=real(-0.5,kdp)
       enddo
       border_list(lp,N_peak_list(lp),1)=l1
       border_list(lp,N_peak_list(lp),2)=l2
       N_peak_list(lp)=N_peak_list(lp)+uno
       peak_list(lp,N_peak_list(lp))=l1
       list_peak(l1)=lp
     endif
   else
     if (list_peak(l2).eq.zero) then
! If only l1 has been assigned to a list: Add l2 to list_peak(l1)
       lp=list_peak(l1)
       do i=uno,N_peak_list(lp)
         border_copy(l2,peak_list(lp,i))=real(-0.5,kdp)
         border_copy(peak_list(lp,i),l2)=real(-0.5,kdp)
       enddo
       border_list(lp,N_peak_list(lp),1)=l1
       border_list(lp,N_peak_list(lp),2)=l2
       N_peak_list(lp)=N_peak_list(lp)+uno
       peak_list(lp,N_peak_list(lp))=l2
       list_peak(l2)=lp
     else
! if l1 belongs to a list and l2 too, merge the two lists in one (lm1)
       lm1=min(list_peak(l1),list_peak(l2))
       lm2=max(list_peak(l1),list_peak(l2))
       border_list(lm1,N_peak_list(lm1),1)=l1
       border_list(lm1,N_peak_list(lm1),2)=l2
       do i=uno,N_peak_list(lm2)-uno
         border_list(lm1,N_peak_list(lm1)+i,1)=border_list(lm2,i,1) 
         border_list(lm1,N_peak_list(lm1)+i,2)=border_list(lm2,i,2) 
       enddo
       do i=uno,N_peak_list(lm2)
         N_peak_list(lm1)=N_peak_list(lm1)+uno
         peak_list(lm1,N_peak_list(lm1))=peak_list(lm2,i)
         list_peak(peak_list(lm2,i))=lm1
       enddo
       N_peak_list(lm2)=zero
       border_list(lm2,:,:)=zero
       peak_list(lm2,:)=zero
       do i=uno,N_peak_list(lm1)-uno
         do j=uno,N_peak_list(lm1)
           border_copy(peak_list(lm1,i),peak_list(lm1,j))=real(-0.5,kdp)
           border_copy(peak_list(lm1,j),peak_list(lm1,i))=real(-0.5,kdp)
         enddo
       enddo
     endif
   endif
 enddo
! Write graph and recover information for dendrogram
! dendro_points: coordinates of peaks and barriers 
! dendro_active_points: are this point still alive ?
! dendro_barrier :: is this point a barrier ?
! dendro_active_barrier :: index of the active barrier
 dendro_active_points(:)=.false.
 dendro_barrier(:)=.false.
 dendro_label(:)="   "
 dendro_width(:)=rzero
 l=Pop_m(peak_list(1,1))/dos
 k=uno
 write (dendro_label(k),'(i3)') peak_list(1,1)
 dendro_points(k,1)=dfloat(l)
 dendro_points(k,2)=cent_m(peak_list(1,1))
 dendro_density(k)=cent_m(peak_list(1,1))
 dendro_error(k)=cent_err_m(peak_list(1,1))
 dendro_active_points(k)=.true.
 dendro_barrier(k)=.false.
 do i=dos,N_peak_list(1)
   j=i-uno
   b=Bord_m(border_list(1,j,1),border_list(1,j,2))
   eb=Bord_err_m(border_list(1,j,1),border_list(1,j,2))
   l=l+Pop_m(peak_list(1,j))/dos
   k=k+uno
   dendro_points(k,1)=dfloat(l)
   dendro_points(k,2)=b
   dendro_density(k)=b
   dendro_error(k)=eb
   dendro_active_points(k)=.false.
   dendro_barrier(k)=.true.
   l=l+Pop_m(peak_list(1,i))/dos
   k=k+uno
   write (dendro_label(k),'(i3)') peak_list(1,i)
   dendro_points(k,1)=dfloat(l)
   dendro_points(k,2)=cent_m(peak_list(1,i))
   dendro_density(k)=cent_m(peak_list(1,i))
   dendro_error(k)=cent_err_m(peak_list(1,i))
   dendro_active_points(k)=.true.
   dendro_barrier(k)=.false.
 enddo
 l=l+Pop_m(peak_list(1,N_peak_list(1)))/dos

 open (21,file="Dendrogram_graph_2.gpl")
 open (22,file="Dendrogram_labels_2.dat")
 dendro_n_active_points=zero
 Narrow=zero
 3333 FORMAT ("set arrow ",i4," from ",e14.7,",",e14.7," to ",e14.7,",",e14.7," nohead lw ",f14.7)
  xmax=real(-1.,kdp)
 do i=uno,dos*Nclus_m-uno
   if (dendro_barrier(i).and. (dendro_points(i,2).gt.xmax)) then
      dendro_active_barrier=i
      xmax=dendro_points(i,2)
   endif
 enddo
 do i=uno,dos*Nclus_m-uno
   if (dendro_active_points(i)) then
      dendro_n_active_points=dendro_n_active_points+uno
   endif
 enddo
 do while (dendro_n_active_points.gt.uno)
! Look for the higher barrier
   xmax=real(-1.,kdp)
   do i=uno,dos*Nclus_m-uno
     if (dendro_barrier(i).and. (dendro_points(i,2).gt.xmax)) then
        dendro_active_barrier=i
        xmax=dendro_points(i,2)
     endif
   enddo 
! Look for the nearest active points at left and right
   dminleft=maxr
   dminright=maxr
   do i=uno,dos*Nclus_m-uno
     if (dendro_active_points(i)) then
        if ((dendro_points(i,1).lt.dendro_points(dendro_active_barrier,1))) then
          if ((dendro_points(dendro_active_barrier,1)-dendro_points(i,1)).lt.dminleft) then
            dminleft=dendro_points(dendro_active_barrier,1)-dendro_points(i,1) 
            dendro_left_nearest=i
          endif
       else
         if ((dendro_points(i,1)-dendro_points(dendro_active_barrier,1)).lt.dminright) then
            dminright=dendro_points(i,1)-dendro_points(dendro_active_barrier,1)
            dendro_right_nearest=i
          endif
       endif
     endif
   enddo
   width=real(2.,kdp)
   dendro_width(:)=width
   
! Write coordinates

   Narrow=Narrow+uno
   write (21,3333) Narrow,dendro_points(dendro_left_nearest,1),dendro_points(dendro_left_nearest,2),&
   & dendro_points(dendro_left_nearest,1),dendro_points(dendro_active_barrier,2),dendro_width(dendro_left_nearest)
   Narrow=Narrow+uno
   write (21,3333) Narrow,dendro_points(dendro_left_nearest,1),dendro_points(dendro_active_barrier,2),&
   & dendro_points(dendro_right_nearest,1),dendro_points(dendro_active_barrier,2),width
   Narrow=Narrow+uno
   write (21,3333) Narrow,dendro_points(dendro_right_nearest,1),dendro_points(dendro_active_barrier,2),&
   & dendro_points(dendro_right_nearest,1),dendro_points(dendro_right_nearest,2),dendro_width(dendro_right_nearest)

   if (dendro_label(dendro_left_nearest).ne."   ") write (22,*) dendro_points(dendro_left_nearest,1),dendro_points(dendro_left_nearest,2),dendro_label(dendro_left_nearest)
   if(dendro_label(dendro_right_nearest).ne."   ") write (22,*) dendro_points(dendro_right_nearest,1),dendro_points(dendro_right_nearest,2),dendro_label(dendro_right_nearest)

! update barrier list and dendro_active_points and densities and errors
   dendro_active_points(dendro_left_nearest)=.false.
   dendro_active_points(dendro_right_nearest)=.false.
   dendro_active_points(dendro_active_barrier)=.true.
   dendro_width(dendro_active_barrier)=width
   if ((dendro_density(dendro_left_nearest)-dendro_error(dendro_left_nearest)).gt.(dendro_density(dendro_right_nearest)-dendro_error(dendro_right_nearest))) then
     dendro_density(dendro_active_barrier)=dendro_density(dendro_left_nearest)
     dendro_error(dendro_active_barrier)=dendro_error(dendro_left_nearest)
   else
     dendro_density(dendro_active_barrier)=dendro_density(dendro_right_nearest)
     dendro_error(dendro_active_barrier)=dendro_error(dendro_right_nearest)
   endif
   dendro_barrier(dendro_active_barrier)=.false.
! count of active points
   dendro_n_active_points=0
   do i=uno,dos*Nclus_m-uno
     if (dendro_active_points(i)) dendro_n_active_points=dendro_n_active_points+uno
   enddo
 enddo
 do i=uno,dos*Nclus_m-uno
   if (dendro_active_points(i)) k=i
 enddo
 Narrow=Narrow+uno
 last_barrier=dendro_points(k,2)
 maxdens=maxval(Rho(:))
 hdendro=(maxdens-last_barrier)/real(0.9,kdp)
 write (21,3333) Narrow,dendro_points(k,1),dendro_points(k,2),dendro_points(k,1),dendro_points(k,2)-real(0.05,kdp)*hdendro,2.
! write (21,'(a,i6,a)') "pl [0:",Nele,"] [0:] 'Dendrogram_labels_2.dat' u 1:2:3 w p pt 7 ps 2 lc palette z t '','' u 1:2:3 w labels t ''"
 write (21,'(a,i6,a,es17.10,a,es17.10,a)') "pl [0:",Nele,"] [",dendro_points(k,2)-real(0.05,kdp)*hdendro,":",maxdens+real(0.05,kdp)*hdendro,"] 'Dendrogram_labels_2.dat' u 1:2:3 w p pt 7 ps 2 lc palette z t '','' u 1:2:3 w labels t ''"
!write (23,3333) Narrow,dendro_points(k,1),dendro_points(k,2),dendro_points(k,1),-1.,width_max+1.
 close (21)
 close (22)
! call execute_command_line("gnuplot -persists Dendrogram_graph_2.gpl")
!
 open (21,file="mat_labels.dat")
 open (22,file="Topography.mat")
 do i=uno,N_peak_list(1)
   Bord_m(i,i)=cent_m(i)
 enddo
 do i=uno,N_peak_list(1)
   write (22,'(500(f15.5,1x))') (Bord_m(peak_list(1,i),peak_list(1,j)),j=1,N_peak_list(1))
 enddo
 
 do i=uno,N_peak_list(1)
    write (21,*) float(i)-1.0,-1.0,peak_list(1,i)
 enddo
 return
 end subroutine graph2

!
SUBROUTINE HPSORT(N,RA,s_order)
  implicit none
  integer (kind=idp) N,s_order(N)
  real (kind=kdp) RA(N)
  integer (kind=idp) L,IR,I,J,sord
  real (kind=kdp) RRA
  do i=1,N
     s_order(i)=i
  enddo
  L=N/2+1
  IR=N
    !The index L will be decremented from its initial value during the
    !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
    !will be decremented from its initial value down to 1 during the
    !"retirement-and-promotion" (heap selection) phase.
10  continue
  if(L > 1)then
     L=L-1
     RRA=RA(L)
     sord=s_order(L)
  else
     RRA=RA(IR)
     sord=s_order(IR)
     RA(IR)=RA(1)
     s_order(ir)=s_order(1)
     IR=IR-1
     if(IR.eq.1)then
        RA(1)=RRA
        s_order(1)=sord
        return
     end if
   end if
   I=L
   J=L+L
20 if(J.le.IR)then
     if(J < IR)then
        if(RA(J) < RA(J+1))  J=J+1
     end if
     if(RRA < RA(J))then
        RA(I)=RA(J)
        s_order(i)=s_order(j)
        I=J; J=J+J
     else
        J=IR+1
     end if
     goto 20
  end if
  RA(I)=RRA
  s_order(I)=sord
  goto 10
END subroutine HPSORT
SUBROUTINE init_random_seed()
 INTEGER i, clock
 INTEGER seed (100)
 integer N
 n=100
 CALL RANDOM_SEED(size = n)
 CALL SYSTEM_CLOCK(COUNT=clock)
 do i=1,100
   seed (i) = clock + 37 * (i-1)
 enddo
 CALL RANDOM_SEED(PUT = seed)
END SUBROUTINE
subroutine Shuffle(a)
  integer(kind=idp), intent(inout) :: a(:)
  integer(kind=idp) :: i, randpos, temp
  real :: r
  do i = int(size(a),idp), dos, int(-1,idp)
    call random_number(r)
    randpos = int(r * real(i),idp) + uno
    temp = a(randpos)
    a(randpos) = a(i)
    a(i) = temp
  end do
 
end subroutine Shuffle

end program dp_advance
