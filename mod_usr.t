module mod_usr
  use mod_ffhd
  implicit none
  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav,trelax,tstop,temax
  double precision :: heatunit,gzone,SRadius,dya
  double precision, allocatable :: pa(:),ra(:)
  integer, parameter :: jmax=50000

  double precision :: lQ,lQ0,lambdah

  double precision :: bQ0,Qalfa,Qgama
  double precision :: dh,Bh,dv,Bv,xv,cv,ch
  integer :: Btot_
  logical :: ffhd_Btot

  !> parameters for special amr over specific points
  character(len=std_len) :: csvfile
  double precision, allocatable :: pos(:,:)

  !> parameters for set b1,b2,b3 for visualization
  integer :: i_b1,i_b2,i_b3
contains

 subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ Qalfa,Qgama,bQ0,tstop,lQ,lambdah,&
                        dh,Bh,dv,Bv,xv,ffhd_Btot,csvfile
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do
  end subroutine usr_params_read

  subroutine usr_init()
    call set_coordinate_system("Cartesian_3D")

    unit_length        = 1.d9 !< cm
    unit_temperature   = 1.d6 !< K
    unit_numberdensity = 1.d9 !< cm^-3

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_source          => special_source
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_set_B0          => specialset_B0
    usr_modify_output   => set_output_vars

!    usr_process_global  => special_global
!    usr_set_field_w     => special_field_w
    call usr_params_read(par_files)
    call ffhd_activate()

    if (ffhd_Btot) then
      Btot_ = var_set_wextra()
    end if

    i_b1 = var_set_extravar("b1", "b1")
    i_b2 = var_set_extravar("b2", "b2")
    i_b3 = var_set_extravar("b3", "b3")
  end subroutine usr_init

  subroutine initglobaldata_usr()
    heatunit=unit_pressure/unit_time !< 3.697693390805347E-003 erg*cm^-3/s
    usr_grav=-2.74d4*unit_length/unit_velocity**2 !< solar gravity
    gzone=0.5d0 !< thickness of a ghostzone below the bottom boundary
    dya=(2.d0*gzone+xprobmax3-xprobmin3)/dble(jmax) !< cells size of high-resolution 1D solar atmosphere
    SRadius=69.61d0 !< Solar radius
    trelax=100.d0 !< time for relaxation
    temax=2.d0 !< artificial maximum temperature for the top boundary

    lQ0=lQ/heatunit
    call inithdstatic

    ! for fan-spine field use in wyper2016a_field
    cv = (Bv*dv**3.d0)*half
    ch = (Bh*dh**3.d0)*half

    ! initialize possible points for special amr
    ! call read_csv(pos, csvfile)
  end subroutine initglobaldata_usr

  subroutine inithdstatic
    ! use mod_ionization
    integer :: j,na,nb,ibc,flag
    double precision:: rpho,Tcor,Tpho,wtra,res,rhob,pb,htra,Ttr,Fc,invT,kappa,hflag
    double precision, allocatable :: ya(:),Ta(:),gg(:)

    Tcor=1.d0
    Tpho=1.3d-2
    htra=0.6d0
    wtra=0.05d0
    rpho=0.5d0
    hflag=2.d0
    Fc=2.d5/heatunit/unit_length
    kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
    allocate(ya(jmax),Ta(jmax),gg(jmax),pa(jmax),ra(jmax))
    do j=1,jmax
      ya(j)=(dble(j)-0.5d0)*dya-gzone
      if(ya(j)<hflag) then
        flag=j
      end if
      if(ya(j)>=htra) then
        Ta(j)=(3.5d0*Fc/kappa*(ya(j)-htra)+Tcor**3.5d0)**(2.d0/7.d0)
      else
        Ta(j)=Tpho+0.5d0*(Tcor-Tpho)*(1.d0+tanh((ya(j)-htra)/wtra))
      end if
      gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
    end do
    nb=int((gzone+hflag)/dya)
    ra(nb)=rpho
    pa(nb)=rpho*Ta(nb)
    ! if(ffhd_ionization) call ion_get_density(pa(nb),Ta(nb),ra(nb))
    invT=0.d0
    do j=nb+1,jmax
      invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
      pa(j)=pa(nb)*dexp(invT*dya)
      ra(j)=pa(j)/Ta(j)
      ! if(ffhd_ionization) call ion_get_density(pa(j),Ta(j),ra(j))
    enddo
    invT=0.d0
    do j=nb-1,1,-1
      invT=invT-(gg(j)/Ta(j)+gg(j+1)/Ta(j+1))*0.5d0
      pa(j)=pa(nb)*dexp(invT*dya)
      ra(j)=pa(j)/Ta(j)
      ! if(ffhd_ionization) call ion_get_density(pa(j),Ta(j),ra(j))
    enddo
    na=floor(gzone/dya+0.5d0)
    res=gzone-(dble(na)-0.5d0)*dya
    rhob=ra(na)+res/dya*(ra(na+1)-ra(na))
    pb=pa(na)+res/dya*(pa(na+1)-pa(na))
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))
    do ibc=nghostcells,1,-1
      na=floor((gzone-dx(3,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0))/dya+0.5d0)
      res=gzone-dx(3,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dya
      rbc(ibc)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
    end do
    deallocate(ya,gg,Ta)
  end subroutine inithdstatic

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    !> initialize one grid
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: res
    integer :: ix^D,na

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      na=floor((x(ix^D,3)-xprobmin3+gzone)/dya+0.5d0)
      res=x(ix^D,3)-xprobmin3+gzone-(dble(na)-0.5d0)*dya
      w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
      w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
    {end do\}
    w(ixO^S,mom(:))=zero

    call ffhd_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    !> special boundary types, user defined
    integer, intent(in) :: ixO^L,iB,ixI^L
    double precision, intent(in) :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth(ixI^S),tmp(ixI^S),ggrid(ixI^S),invT(ixI^S)
    double precision :: Q(ixI^S),Qp(ixI^S)
    integer :: ix^D,ixOs^L,ixC^L,hxC^L,jxO^L

    select case(iB)
    case(5)
      w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(1))&
                      /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
      do ix3=ixOmin3,ixOmax3
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,rho_)=rbc(ix3)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,p_)=pbc(ix3)
      enddo
      call ffhd_to_conserved(ixI^L,ixO^L,w,x)
    case(6)
      ixOs^L=ixO^L;
      ixOsmin3=ixOmin3-1;ixOsmax3=ixOmin3-1;
      call ffhd_get_pthermal(w,x,ixI^L,ixOs^L,pth)
      ixOsmin3=ixOmin3-1;ixOsmax3=ixOmax3;
      call getggrav(ggrid,ixI^L,ixOs^L,x)
      invT(ixOmin3-1^%3ixO^S)=w(ixOmin3-1^%3ixO^S,rho_)/pth(ixOmin3-1^%3ixO^S)
      tmp=0.d0
      do ix3=ixOmin3,ixOmax3
        tmp(ixOmin3-1^%3ixO^S)=tmp(ixOmin3-1^%3ixO^S)+0.5d0*&
            (ggrid(ix3^%3ixO^S)+ggrid(ix3-1^%3ixO^S))*invT(ixOmin3-1^%3ixO^S)
        w(ix3^%3ixO^S,p_)=pth(ixOmin3-1^%3ixO^S)*dexp(tmp(ixOmin3-1^%3ixO^S)*dxlevel(3))
        w(ix3^%3ixO^S,rho_)=w(ix3^%3ixO^S,p_)*invT(ixOmin3-1^%3ixO^S)
      enddo
      w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,mom(1))&
                      /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
      call ffhd_to_conserved(ixI^L,ixO^L,w,x)
    case default
      call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine specialbound_usr

  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,3)=ggrid(ixO^S)
  end subroutine gravity

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(SRadius/(SRadius+x(ixO^S,3)))**2
  end subroutine

  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    integer, intent(in) :: ixI^L,ixO^L,iw^LIM
    double precision, intent(in) :: qdt,qtC,qt
    double precision, intent(in) :: x(ixI^S,1:ndim),wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: lQgrid(ixI^S),bQgrid(ixI^S)

    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)
    call getlQ(lQgrid,ixI^L,ixO^L,qt,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*lQgrid(ixO^S)
  end subroutine special_source

  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
    !> calculate background heating bQ
    use mod_radiative_cooling
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qt,x(ixI^S,1:ndim),w(ixI^S,1:nw)
    double precision, intent(out) :: bQgrid(ixI^S)
    double precision :: bQ1,bQ2,lambda

    select case(iprob)
    case(1,2,3)
      bQ1=max((-0.04d0*qt+3.d0)*1.d-4,1.d-4)/heatunit
      lambda=5.d0
      bQgrid(ixO^S)=bQ1*dexp(-x(ixO^S,3)/lambda)
    case default
      call mpistop("iprob wrong in getbQ")
    end select
    !if(ffhd_Btot) then
    !  bQ1=1.d-3/heatunit
    !  bQ2=3.d-6/heatunit
    !  bQgrid(ixO^S)=bQ0*block%wextra(ixO^S,Btot_)**Qalfa*w(ixO^S,rho_)**Qgama
    !  where(bQgrid(ixO^S) .lt. bQ2)
    !    bQgrid(ixO^S)=bQ2
    !  !elsewhere(bQgrid(ixO^S) .gt. bQ1)
    !  !  bQgrid(ixO^S)=bQ1
    !  endwhere
    !else
    !endif
  end subroutine getbQ

  subroutine getlQ(lQgrid,ixI^L,ixO^L,qt,w,x)
  !> calculate localized heating lQ
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)

    double precision :: lQgrid(ixI^S)
    double precision :: tramp,tr,theat,tres

    tramp=500.d0/unit_time !< time for localized heating to increase from 0 to 1
    theat=qt-trelax
    tr=zero
    if (theat .lt. zero) then
      tr=zero
    elseif(theat .lt. tramp) then
      tr=theat/tramp
    elseif(theat .lt. tstop-trelax-tramp) then
      tr=one
    elseif(theat .lt. tstop-trelax) then
      tr=(tstop-trelax-theat)/tramp
    endif
    lQgrid(ixO^S)=lQ0*tr*dexp(-(x(ixO^S,3)-0.15d0)**2/lambdah**2)
    select case(iprob)
    case(1)
      lQgrid(ixO^S)=lQgrid(ixO^S)*(block%wextra(ixO^S,Btot_)/4.d0)**2
    case default
      call mpistop("iprob wrong xxxxxxxxxxxxxx")
    end select
  end subroutine getlQ

  subroutine wyper2016a_field(ixI^L,ixO^L,x,Bvec)

    ! from wyper to here: x --> z; y --> x; z --> y
    ! from here to wyper: x --> y; y --> z; z --> x

    integer, intent(in)  :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out) :: Bvec(ixI^S,1:ndir)

    double precision :: xh, yh, zh, yv, zv

    xh = 0.d0
    yh = 0.d0
    zh = -1 * dh

    ! xv already defined in amrvac.par
    yv = 0.d0
    zv = -1 * dv

    Bvec(ixO^S,1) = -3.d0*ch*((x(ixO^S,3)-zh)**2+(x(ixO^S,2)-yh)**2) / ((x(ixO^S,1)-xh)**2+(x(ixO^S,2)-yh)**2+(x(ixO^S,3)-zh)**2)**(5./2.) + &
                    3.d0*cv*(x(ixO^S,1)-xv)*(x(ixO^S,3)-zv) / ((x(ixO^S,1)-xv)**2+(x(ixO^S,2)-yv)**2+(x(ixO^S,3)-zv)**2)**(5./2.) + &
                    2*ch / ((x(ixO^S,1)-xh)**2+(x(ixO^S,2)-yh)**2+(x(ixO^S,3)-zh)**2)**(3./2.)

    Bvec(ixO^S,2) = 3.d0*ch*(x(ixO^S,1)-xh)*(x(ixO^S,2)-yh) / ((x(ixO^S,1)-xh)**2+(x(ixO^S,2)-yh)**2+(x(ixO^S,3)-zh)**2)**(5./2.) + &
                    3.d0*cv*(x(ixO^S,2)-yv)*(x(ixO^S,3)-zv) / ((x(ixO^S,1)-xv)**2+(x(ixO^S,2)-yv)**2+(x(ixO^S,3)-zv)**2)**(5./2.)

    Bvec(ixO^S,3) = 3.d0*ch*(x(ixO^S,1)-xh)*(x(ixO^S,3)-zh) / ((x(ixO^S,1)-xh)**2+(x(ixO^S,2)-yh)**2+(x(ixO^S,3)-zh)**2)**(5./2.) + & 
                    (-3.d0)*cv*((x(ixO^S,2)-yv)**2+(x(ixO^S,1)-xv)**2) / ((x(ixO^S,1)-xv)**2+(x(ixO^S,2)-yv)**2+(x(ixO^S,3)-zv)**2)**(5./2.) + & 
                    2.d0*cv/((x(ixO^S,1)-xv)**2+(x(ixO^S,2)-yv)**2+(x(ixO^S,3)-zv)**2)**(3./2.)

  end subroutine wyper2016a_field

  subroutine curl2(qvec,ixI^L,ixO^L,xL,curlvec)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, dimension(ixI^S,1:3), intent(in) :: xL,qvec
    double precision, dimension(ixI^S,1:3), intent(inout) :: curlvec
    !> local
    integer :: idir,jdir,kdir,hxO^L,jxO^L
    double precision, dimension(ixI^S) :: tmp

    curlvec=zero
    do idir=1,3;do jdir=1,3;do kdir=1,3
      if(lvc(idir,jdir,kdir)/=0) then
        tmp(ixI^S)=qvec(ixI^S,kdir)
        hxO^L=ixO^L-kr(jdir,^D);
        jxO^L=ixO^L+kr(jdir,^D);
        tmp(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(xL(jxO^S,jdir)-xL(hxO^S,jdir))
        if(lvc(idir,jdir,kdir)==1) then
          curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp(ixO^S)
        else
         curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp(ixO^S)
        end if
      end if
    end do;end do;end do
  end subroutine curl2

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
    !> output: wB0 must be a unit field vector
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)
    !> local
    integer :: ixK^L
    ! double precision :: B1(ixI^S,1:ndir) !> background field
    ! double precision :: B2(ixI^S,1:ndir) !> flux rope
    double precision :: B3(ixI^S,1:ndir) !> wyper field
    double precision :: Btot(ixI^S)

    ! B1=zero
    ! B2=zero
    ixK^L=ixO^L^LADD1;
    call wyper2016a_field(ixI^L,ixO^L,x,B3)
    wB0=B3
    Btot(ixO^S)=dsqrt(wB0(ixO^S,1)**2+wB0(ixO^S,2)**2+wB0(ixO^S,3)**2)
    if(ffhd_Btot) then
      block%wextra(ixO^S,Btot_)=Btot(ixO^S)
    endif
    wB0(ixO^S,1)=wB0(ixO^S,1)/Btot(ixO^S)
    wB0(ixO^S,2)=wB0(ixO^S,2)/Btot(ixO^S)
    wB0(ixO^S,3)=wB0(ixO^S,3)/Btot(ixO^S)
  end subroutine specialset_B0

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    integer, intent(in) :: igrid,level,ixI^L,ixO^L
    double precision, intent(in) :: qt
    double precision, dimension(ixI^S,1:nw), intent(in) :: w
    double precision, dimension(ixI^S,1:ndim), intent(in) :: x
    integer, intent(inout) :: refine,coarsen

    integer :: i

      ! do i=1,npoints
      !   if (pos(1,i) >= minval(x(ixO^S,1)) .and. pos(1,i) <= maxval(x(ixO^S,1)) .and. &
      !       pos(2,i) >= minval(x(ixO^S,2)) .and. pos(2,i) <= maxval(x(ixO^S,2)) .and. &
      !       pos(3,i) >= minval(x(ixO^S,3)) .and. pos(3,i) <= maxval(x(ixO^S,3))) then
      !     refine=1
      !     coarsen=-1
      !     exit
      !   else 
      !     refine=-1
      !     coarsen=-1
      !   end if
      ! end do

    ! if(any(abs(x(ixO^S,3)) .lt. 0.3)) then
    !   refine=1
    !   coarsen=-1
    ! end if

  end subroutine special_refine_grid

  subroutine set_output_vars(ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,nw)
    double precision                :: divb(ixI^S)
    double precision                :: v0(ndim)
    double precision                :: bfield(ixI^S, ndir)
    double precision                :: ws(ixGs^T, nws)
    double precision :: wlocal(ixI^S,1:nw),pth(ixI^S)
    double precision :: wion
    double precision :: cor(ixI^S)
    integer :: ix^D

    ! wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! call ffhd_get_pthermal(wlocal,x,ixI^L,ixI^L,pth)

    w(ixO^S,i_b1)=block%B0(ixO^S,1,0)*block%wextra(ixO^S,Btot_)
    w(ixO^S,i_b2)=block%B0(ixO^S,2,0)*block%wextra(ixO^S,Btot_)
    w(ixO^S,i_b3)=block%B0(ixO^S,3,0)*block%wextra(ixO^S,Btot_)

  end subroutine set_output_vars

end module mod_usr
