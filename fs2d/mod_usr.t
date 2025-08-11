module mod_usr
  use mod_ffhd
  implicit none
  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav,trelax,tstop
  double precision :: heatunit,gzone,SRadius,dya
  double precision, allocatable :: pa(:),ra(:)
  integer, parameter :: jmax=10000

  double precision :: lQ,lQ0,sigma,asym

  double precision :: htra, xl, xr

  !> parameters for background heating
  double precision :: bQ0,Qalfa,Qgama
  !> parameters for modified output
  integer :: i_Te,i_b1,i_b2
  !> parameters for fan-spine field in wyper2016
  double precision :: dh,Bh,dv,Bv,xv,cv,ch
  !> parameters for Btot as extra variable
  integer :: Btot_
contains

 subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ Qalfa,Qgama,bQ0,& !> for background heating
                        trelax,tstop,lQ,& !> for localized heating (temporal profile)
                        htra,asym,xr,xl,sigma,& !> for localized heating (spaitial profile)
                        dh,Bh,dv,Bv,xv !> for fan-spine field in wyper2016
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do
  end subroutine usr_params_read

  subroutine usr_init()
    call set_coordinate_system("Cartesian_2D")

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

    Btot_ = var_set_wextra()

    i_b1 = var_set_extravar("b1", "b1")
    i_b2 = var_set_extravar("b2", "b2")
    i_Te = var_set_extravar("Te", "Te")
  end subroutine usr_init

  subroutine initglobaldata_usr()
    heatunit=unit_pressure/unit_time !< 3.697693390805347E-003 erg*cm^-3/s
    usr_grav=-2.74d4*unit_length/unit_velocity**2 !< solar gravity
    gzone=0.5d0 !< thickness of a ghostzone below the bottom boundary
    dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax) !< cells size of high-resolution 1D solar atmosphere
    SRadius=69.61d0 !< Solar radius

    lQ0=lQ/heatunit

    ! for fan-spine field use in wyper2016a_field
    cv = (Bv*dv**3.d0)*half
    ch = (Bh*dh**3.d0)*half
    call inithdstatic
  end subroutine initglobaldata_usr

  subroutine inithdstatic
    ! use mod_ionization
    integer :: j,na,nb,ibc,flag
    double precision:: rpho,Tcor,Tpho,wtra,res,rhob,pb,Fc,invT,kappa,hflag
    double precision, allocatable :: ya(:),Ta(:),gg(:)

    Tcor=1.d0
    Tpho=1.3d-2
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
      na=floor((gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0))/dya+0.5d0)
      res=gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dya
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
      na=floor((x(ix^D,ndim)-xprobmin2+gzone)/dya+0.5d0)
      res=x(ix^D,ndim)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
      w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
      w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
    {end do\}
    w(ixO^S,mom(:))=zero
    if (ffhd_hyperbolic_thermal_conduction) w(ixO^S,q_)=zero
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
    case(3)
      w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(1))&
                   /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
        w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
        if(ffhd_hyperbolic_thermal_conduction) w(ix2^%2ixO^S,q_)=zero
      enddo
      call ffhd_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      ixOs^L=ixO^L;
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmin2-1;
      call ffhd_get_pthermal(w,x,ixI^L,ixOs^L,pth)
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmax2;
      call getggrav(ggrid,ixI^L,ixOs^L,x)
      invT(ixOmin2-1^%2ixO^S)=w(ixOmin2-1^%2ixO^S,rho_)/pth(ixOmin2-1^%2ixO^S)
      tmp=0.d0
      do ix2=ixOmin2,ixOmax2
        tmp(ixOmin2-1^%2ixO^S)=tmp(ixOmin2-1^%2ixO^S)+0.5d0*&
            (ggrid(ix2^%2ixO^S)+ggrid(ix2-1^%2ixO^S))*invT(ixOmin2-1^%2ixO^S)
        w(ix2^%2ixO^S,p_)=pth(ixOmin2-1^%2ixO^S)*dexp(tmp(ixOmin2-1^%2ixO^S)*dxlevel(2))
        w(ix2^%2ixO^S,rho_)=w(ix2^%2ixO^S,p_)*invT(ixOmin2-1^%2ixO^S)
        if(ffhd_hyperbolic_thermal_conduction) w(ix2^%2ixO^S,q_)=zero
      enddo
      w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,mom(1))&
                      /w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,rho_)
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
    gravity_field(ixO^S,ndim)=ggrid(ixO^S)
  end subroutine gravity

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(SRadius/(SRadius+x(ixO^S,ndim)))**2
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

    bQ1=max((-0.04d0*qt+3.d0)*1.d-4,1.d-4)/heatunit
    lambda=5.d0
    bQgrid(ixO^S)=bQ1*dexp(-x(ixO^S,ndim)/lambda)
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

    double precision :: lQgrid(ixI^S),factors(ixI^S),tr
    double precision :: al,ar
    double precision :: tramp, theat
    double precision :: heatx_left(100),heaty_left(100),heatx_right(100),heaty_right(100)
    double precision :: heatx_center(100),heaty_center(100)
    integer :: tloc
    double precision :: xd1,yd1,xd2,yd2,xd3,yd3,tres

    !> asymmetric footpoint local heating
    al=two*asym/(1.0+asym)
    ar=two*1.0/(1.0+asym)

    !> time profile for localized heating (relax->ramp->peak->ramp->relax; ramp is 500s)
    tramp=500.d0/unit_time !< time for localized heating to increase from 0 to 1, about 6 ta
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

    !> heating profile depending on random point at the separatrix layer
    heatx_left = (/ &
    -7.466381, -7.459293, -7.459293, -7.459293, -7.468071, -7.459293, -7.466381, -7.468071, -7.459569, -7.468071, &
    -7.459293, -7.464707, -7.466381, -7.468071, -7.468071, -7.466381, -7.466381, -7.459293, -7.466381, -7.464707, &
    -7.459569, -7.459569, -7.464707, -7.459293, -7.468071, -7.459293, -7.464707, -7.459569, -7.468071, -7.459569, &
    -7.468071, -7.459569, -7.466381, -7.464707, -7.466381, -7.464707, -7.459293, -7.464707, -7.468071, -7.459293, &
    -7.459569, -7.468071, -7.459293, -7.459569, -7.459569, -7.459293, -7.459293, -7.459293, -7.459569, -7.459569, &
    -7.459569, -7.459293, -7.468071, -7.466381, -7.459293, -7.459293, -7.468071, -7.464707, -7.466381, -7.459569, &
    -7.459569, -7.468071, -7.459569, -7.459569, -7.468071, -7.468071, -7.464707, -7.466381, -7.466381, -7.464707, &
    -7.459569, -7.468071, -7.464707, -7.464707, -7.459569, -7.459569, -7.468071, -7.464707, -7.459293, -7.459293, &
    -7.468071, -7.459293, -7.459293, -7.466381, -7.468071, -7.468071, -7.459293, -7.464707, -7.459569, -7.466381, &
    -7.466381, -7.466381, -7.468071, -7.459293, -7.459293, -7.459293, -7.459293, -7.466381, -7.466381, -7.459569 /)

    heaty_left = (/ &
    0.044294, 0.111672, 0.111672, 0.111672, 0.016427, 0.111672, 0.044294, 0.016427, 0.122934, 0.016427, &
    0.111672, 0.048332, 0.044294, 0.016427, 0.016427, 0.044294, 0.044294, 0.111672, 0.044294, 0.048332, &
    0.122934, 0.122934, 0.048332, 0.111672, 0.016427, 0.111672, 0.048332, 0.122934, 0.016427, 0.122934, &
    0.016427, 0.122934, 0.044294, 0.048332, 0.044294, 0.048332, 0.111672, 0.048332, 0.016427, 0.111672, &
    0.122934, 0.016427, 0.111672, 0.122934, 0.122934, 0.111672, 0.111672, 0.111672, 0.122934, 0.122934, &
    0.122934, 0.111672, 0.016427, 0.044294, 0.111672, 0.111672, 0.016427, 0.048332, 0.044294, 0.122934, &
    0.122934, 0.016427, 0.122934, 0.122934, 0.016427, 0.016427, 0.048332, 0.044294, 0.044294, 0.048332, &
    0.122934, 0.016427, 0.048332, 0.048332, 0.122934, 0.122934, 0.016427, 0.048332, 0.111672, 0.111672, &
    0.016427, 0.111672, 0.111672, 0.044294, 0.016427, 0.016427, 0.111672, 0.048332, 0.122934, 0.044294, &
    0.044294, 0.044294, 0.016427, 0.111672, 0.111672, 0.111672, 0.111672, 0.044294, 0.044294, 0.122934 /)

    heatx_right = (/ &
    -2.909872, -2.915581, -2.938322, -2.915089, -2.938322, -2.937916, -2.909872, -2.913944, -2.911274, -2.916208, &
    -2.913395, -2.915089, -2.938322, -2.938322, -2.914738, -2.915089, -2.915089, -2.913944, -2.938322, -2.997004, &
    -2.921004, -2.940426, -2.914033, -2.953556, -2.938322, -3.000093, -2.929754, -2.914738, -2.938322, -2.916908, &
    -2.920487, -2.915089, -2.954845, -2.938322, -2.929754, -2.915581, -2.938322, -2.938322, -2.938322, -2.911274, &
    -2.915089, -2.914033, -2.938322, -2.969098, -2.945798, -2.938322, -2.909574, -2.915089, -2.938322, -2.938322, &
    -2.938322, -2.915089, -2.914738, -2.944577, -2.969098, -2.997004, -2.938322, -2.938322, -2.938322, -2.909872, &
    -2.938322, -2.938322, -2.915089, -2.916208, -2.940426, -2.914089, -2.937916, -2.914089, -2.917553, -2.938322, &
    -2.997004, -2.915089, -2.921262, -2.944577, -2.915089, -2.913240, -2.920487, -2.915089, -2.914033, -2.938322, &
    -2.938322, -2.909872, -2.938322, -2.915581, -2.915089, -2.938322, -2.919700, -2.915089, -2.915089, -2.938322, &
    -2.921262, -2.942626, -2.915089, -2.918642, -2.945798, -2.938322, -2.917553, -3.000093, -2.919700, -2.921004 /)

    heaty_right = (/ &
    0.012592, 0.022363, 0.102452, 0.021142, 0.102452, 0.103026, 0.012592, 0.024032, 0.017367, 0.035996, &
    0.021334, 0.021142, 0.102452, 0.102452, 0.012573, 0.021142, 0.021142, 0.024032, 0.102452, 0.286774, &
    0.047379, 0.106613, 0.018814, 0.158901, 0.102452, 0.283096, 0.079969, 0.012573, 0.102452, 0.020182, &
    0.042344, 0.021142, 0.149076, 0.102452, 0.079969, 0.022363, 0.102452, 0.102452, 0.102452, 0.017367, &
    0.021142, 0.018814, 0.102452, 0.202137, 0.134486, 0.102452, 0.016159, 0.021142, 0.102452, 0.102452, &
    0.102452, 0.021142, 0.012573, 0.123228, 0.202137, 0.286774, 0.102452, 0.102452, 0.102452, 0.012592, &
    0.102452, 0.102452, 0.021142, 0.035996, 0.106613, 0.025727, 0.103026, 0.025727, 0.040434, 0.102452, &
    0.286774, 0.021142, 0.050604, 0.123228, 0.021142, 0.024370, 0.042344, 0.021142, 0.018814, 0.102452, &
    0.102452, 0.012592, 0.102452, 0.022363, 0.021142, 0.102452, 0.046022, 0.021142, 0.021142, 0.102452, &
    0.050604, 0.125117, 0.021142, 0.026890, 0.134486, 0.102452, 0.040434, 0.283096, 0.046022, 0.047379 /)

    heatx_center = (/ &
    -5.133157, -5.147233, -5.145180, -5.127937, -5.123587, -5.127937, -5.119565, -5.139627, -5.124581, -5.124714, &
    -5.124714, -5.115056, -5.124714, -5.145180, -5.113275, -5.118808, -5.124714, -5.139787, -5.116294, -5.114500, &
    -5.124714, -5.145180, -5.120397, -5.124714, -5.111620, -5.145180, -5.118808, -5.145180, -5.110535, -5.115056, &
    -5.147646, -5.114500, -5.137662, -5.113388, -5.114466, -5.112173, -5.123265, -5.147646, -5.118808, -5.114500, &
    -5.114500, -5.139388, -5.125105, -5.115056, -5.124714, -5.112500, -5.145180, -5.121586, -5.145180, -5.132793, &
    -5.135142, -5.112345, -5.145180, -5.130021, -5.124714, -5.139627, -5.149610, -5.124714, -5.110430, -5.118543, &
    -5.114466, -5.133157, -5.124714, -5.113275, -5.145180, -5.125105, -5.145180, -5.116294, -5.110153, -5.116294, &
    -5.113027, -5.124714, -5.124714, -5.124714, -5.113388, -5.124714, -5.146618, -5.146618, -5.145180, -5.114466, &
    -5.114500, -5.138992, -5.114500, -5.123754, -5.110623, -5.114466, -5.114500, -5.114500, -5.123587, -5.114500, &
    -5.113275, -5.145180, -5.124208, -5.114500, -5.110430, -5.114500, -5.139627, -5.147646, -5.147646, -5.114466 /)

    heaty_center = (/ &
    0.175589, 0.255934, 0.258732, 0.119222, 0.125996, 0.119222, 0.065171, 0.227625, 0.086821, 0.111871, &
    0.111871, 0.062969, 0.111871, 0.258732, 0.017041, 0.083755, 0.111871, 0.235881, 0.054100, 0.040152, &
    0.111871, 0.258732, 0.083858, 0.111871, 0.020337, 0.258732, 0.083755, 0.258732, 0.015056, 0.062969, &
    0.282310, 0.040152, 0.203655, 0.026352, 0.043160, 0.031607, 0.111483, 0.282310, 0.083755, 0.040152, &
    0.040152, 0.197361, 0.104900, 0.062969, 0.111871, 0.023812, 0.258732, 0.084391, 0.258732, 0.163893, &
    0.180031, 0.019344, 0.258732, 0.138734, 0.111871, 0.227625, 0.284771, 0.111871, 0.028939, 0.063028, &
    0.043160, 0.175589, 0.111871, 0.017041, 0.258732, 0.104900, 0.258732, 0.054100, 0.029090, 0.054100, &
    0.019814, 0.111871, 0.111871, 0.111871, 0.026352, 0.111871, 0.256369, 0.256369, 0.258732, 0.043160, &
    0.040152, 0.203111, 0.040152, 0.129971, 0.013183, 0.043160, 0.040152, 0.040152, 0.125996, 0.040152, &
    0.017041, 0.258732, 0.110782, 0.040152, 0.028939, 0.040152, 0.227625, 0.282310, 0.282310, 0.043160 /)

    factors(ixO^S) = one

    if ((theat > 0) .and. (theat < tstop-trelax)) then
      tloc = floor(theat/(300.d0/unit_time)) ! about five minutes?
      tres = (theat-300.d0/unit_time*dble(tloc))/(300.d0/unit_time)
      xd1 = heatx_left(tloc+1)*(1.0-tres)+heatx_left(tloc+2)*tres
      yd1 = heaty_left(tloc+1)*(1.0-tres)+heaty_left(tloc+2)*tres
      xd2 = heatx_right(tloc+1)*(1.0-tres)+heatx_right(tloc+2)*tres
      yd2 = heaty_right(tloc+1)*(1.0-tres)+heaty_right(tloc+2)*tres
      xd3 = heatx_center(tloc+1)*(1.0-tres)+heatx_center(tloc+2)*tres
      yd3 = heaty_center(tloc+1)*(1.0-tres)+heaty_center(tloc+2)*tres
      factors(ixO^S) = dexp(-(x(ixO^S,1)-xd1)**2/sigma**2-(x(ixO^S,2)-yd1)**2/sigma**2)+&
                      dexp(-(x(ixO^S,1)-xd2)**2/sigma**2-(x(ixO^S,2)-yd2)**2/sigma**2)
                      ! dexp(-(x(ixO^S,1)-xd3)**2/sigma**2-(x(ixO^S,2)-yd3)**2/sigma**2)
    end if


    lQgrid(ixO^S) = lQ0*tr*factors(ixO^S)

    ! lQgrid(ixO^S) = lQgrid(ixO^S)*(dexp(-(x(ixO^S,1)-xr)**2/sigma**2)*ar+&
    !                                dexp(-(x(ixO^S,1)-xl)**2/sigma**2)*al)

    ! if (ffhd_Btot) then
    !   lQgrid(ixO^S)=lQgrid(ixO^S)*(block%wextra(ixO^S,Btot_)/9.25d0)**2
    !   lQgrid(ixO^S)=lQgrid(ixO^S)*((ar-al)/(xr-xl)*(x(ixO^S,1)-xl)+al)
    ! end if
  end subroutine getlQ

  subroutine wyper2016a_field(ixI^L,ixO^L,x,Bvec)

    ! from wyper to here: x --> z; y --> x; z --> y
    ! from here to wyper: x --> y; y --> z; z --> x

    integer, intent(in)  :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out) :: Bvec(ixI^S,1:ndir)

    double precision :: xh, zh, zv

    xh = 0.d0
    zh = -1 * dh

    ! xv already defined in amrvac.par
    zv = -1 * dv

    Bvec(ixO^S,1) = -3.d0*ch*((x(ixO^S,2)-zh)**2) / ((x(ixO^S,1)-xh)**2+(x(ixO^S,2)-zh)**2)**(5./2.) + &
                    3.d0*cv*(x(ixO^S,1)-xv)*(x(ixO^S,2)-zv) / ((x(ixO^S,1)-xv)**2+(x(ixO^S,2)-zv)**2)**(5./2.) + &
                    2*ch / ((x(ixO^S,1)-xh)**2+(x(ixO^S,2)-zh)**2)**(3./2.)

    Bvec(ixO^S,2) = 3.d0*ch*(x(ixO^S,1)-xh)*(x(ixO^S,2)-zh) / ((x(ixO^S,1)-xh)**2+(x(ixO^S,2)-zh)**2)**(5./2.) + & 
                    (-3.d0)*cv*((x(ixO^S,1)-xv)**2) / ((x(ixO^S,1)-xv)**2+(x(ixO^S,2)-zv)**2)**(5./2.) + & 
                    2.d0*cv/((x(ixO^S,1)-xv)**2+(x(ixO^S,2)-zv)**2)**(3./2.)

  end subroutine wyper2016a_field

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
    !> output: wB0 must be a unit field vector
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    double precision :: Btot(ixI^S)

    call wyper2016a_field(ixI^L,ixO^L,x,wB0)
    Btot(ixO^S)=dsqrt(wB0(ixO^S,1)**2+wB0(ixO^S,2)**2)
    block%wextra(ixO^S,Btot_)=Btot(ixO^S)
    wB0(ixO^S,1)=wB0(ixO^S,1)/Btot(ixO^S)
    wB0(ixO^S,2)=wB0(ixO^S,2)/Btot(ixO^S)
  end subroutine specialset_B0

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    integer, intent(in) :: igrid,level,ixI^L,ixO^L
    double precision, intent(in) :: qt
    double precision, dimension(ixI^S,1:nw), intent(in) :: w
    double precision, dimension(ixI^S,1:ndim), intent(in) :: x
    integer, intent(inout) :: refine,coarsen

    if(any(abs(x(ixO^S,ndim)) .lt. htra)) then
      refine=1
      coarsen=-1
    end if

  end subroutine special_refine_grid

  subroutine set_output_vars(ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    ! use mod_ionization

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,nw)
    double precision :: wlocal(ixI^S,1:nw),pth(ixI^S)

    w(ixO^S,i_b1)=block%B0(ixO^S,1,0)*block%wextra(ixO^S,Btot_)
    w(ixO^S,i_b2)=block%B0(ixO^S,2,0)*block%wextra(ixO^S,Btot_)

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    call ffhd_get_pthermal(wlocal,x,ixI^L,ixI^L,pth)
    w(ixO^S,i_Te)=pth(ixO^S)/w(ixO^S,rho_)
  end subroutine set_output_vars

end module mod_usr
