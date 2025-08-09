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
    integer :: tloc
    double precision :: xd1,yd1,xd2,yd2,tres

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
    -7.464633, -7.450113, -7.472507, -7.462049, -7.437082, -7.461177, -7.468804, -7.472507, -7.464633, -7.465919, &
    -7.461177, -7.453475, -7.464633, -7.470363, -7.470363, -7.466514, -7.465569, -7.470269, -7.465569, -7.470075, &
    -7.457613, -7.467940, -7.467389, -7.467389, -7.470363, -7.465569, -7.441896, -7.467940, -7.467940, -7.464633, &
    -7.466514, -7.453475, -7.470075, -7.465919, -7.441896, -7.468849, -7.461177, -7.466514, -7.468849, -7.470075, &
    -7.457613, -7.470269, -7.466514, -7.468804, -7.468804, -7.468804, -7.472507, -7.470269, -7.465919, -7.453475, &
    -7.467389, -7.455415, -7.470269, -7.472507, -7.441896, -7.470363, -7.453723, -7.441896, -7.468804, -7.453475, &
    -7.461177, -7.453475, -7.468804, -7.470075, -7.453723, -7.437082, -7.466514, -7.462049, -7.472507, -7.468804, &
    -7.457613, -7.457613, -7.470269, -7.472507, -7.466514, -7.470269, -7.470363, -7.467940, -7.472507, -7.465569, &
    -7.472507, -7.462049, -7.450113, -7.465919, -7.462049, -7.467940, -7.466514, -7.465919, -7.465569, -7.453723, &
    -7.465569, -7.470269, -7.466237, -7.465569, -7.468849, -7.465919, -7.468849, -7.465569, -7.472507, -7.441896 /)

    heaty_left = (/ &
    0.035742, 0.195737, 0.023934, 0.091613, 0.299468, 0.138087, 0.014741, 0.023934, 0.035742, 0.042336, &
    0.138087, 0.213476, 0.035742, 0.025319, 0.025319, 0.015777, 0.032359, 0.059909, 0.032359, 0.034679, &
    0.146637, 0.029906, 0.016739, 0.016739, 0.025319, 0.032359, 0.270388, 0.029906, 0.029906, 0.035742, &
    0.015777, 0.213476, 0.034679, 0.042336, 0.270388, 0.049563, 0.138087, 0.015777, 0.049563, 0.034679, &
    0.146637, 0.059909, 0.015777, 0.014741, 0.014741, 0.014741, 0.023934, 0.059909, 0.042336, 0.213476, &
    0.016739, 0.139936, 0.059909, 0.023934, 0.270388, 0.025319, 0.162861, 0.270388, 0.014741, 0.213476, &
    0.138087, 0.213476, 0.014741, 0.034679, 0.162861, 0.299468, 0.015777, 0.091613, 0.023934, 0.014741, &
    0.146637, 0.146637, 0.059909, 0.023934, 0.015777, 0.059909, 0.025319, 0.029906, 0.023934, 0.032359, &
    0.023934, 0.091613, 0.195737, 0.042336, 0.091613, 0.029906, 0.015777, 0.042336, 0.032359, 0.162861, &
    0.032359, 0.059909, 0.088776, 0.032359, 0.049563, 0.042336, 0.049563, 0.032359, 0.023934, 0.270388 /)

    heatx_right = (/ &
    -2.925009, -2.987860, -2.951121, -2.951121, -2.925009, -2.921696, -2.925009, -2.970336, -2.912734, -2.964204, &
    -2.912399, -2.925009, -2.925009, -2.913982, -2.925009, -2.970949, -2.914512, -2.976825, -2.928897, -2.964204, &
    -2.925009, -2.925009, -2.940842, -2.917113, -2.911652, -2.925009, -2.925009, -2.921696, -2.918969, -2.925009, &
    -2.912399, -2.933827, -2.925009, -2.925009, -2.964204, -2.940842, -2.948856, -2.976825, -2.925009, -2.918969, &
    -2.925009, -2.925009, -2.925009, -2.912734, -2.925009, -2.949335, -2.925009, -2.915770, -2.941692, -2.949335, &
    -2.949335, -2.925009, -2.925009, -2.925009, -2.928897, -2.913982, -2.925009, -2.917056, -2.925009, -2.925009, &
    -2.970949, -2.925009, -2.923806, -2.909286, -2.912399, -2.925009, -2.917113, -2.918969, -2.914512, -2.951121, &
    -2.923806, -2.917056, -2.912734, -2.925009, -2.925009, -2.933827, -2.949335, -2.925009, -2.910203, -2.925009, &
    -2.964204, -2.912399, -2.950107, -2.917056, -2.915770, -2.911652, -2.951121, -2.970336, -2.923806, -2.951121, &
    -2.964204, -2.933827, -2.928897, -2.921696, -2.915770, -2.976825, -2.933827, -2.925009, -2.917056, -2.949335 /)

    heaty_right = (/ &
    0.059409, 0.252597, 0.148177, 0.148177, 0.059409, 0.050034, 0.059409, 0.220226, 0.012594, 0.192700, &
    0.012779, 0.059409, 0.059409, 0.022130, 0.059409, 0.207788, 0.035292, 0.227743, 0.084341, 0.192700, &
    0.059409, 0.059409, 0.106762, 0.023678, 0.017364, 0.059409, 0.059409, 0.050034, 0.036767, 0.059409, &
    0.012779, 0.095822, 0.059409, 0.059409, 0.192700, 0.106762, 0.139459, 0.227743, 0.059409, 0.036767, &
    0.059409, 0.059409, 0.059409, 0.012594, 0.059409, 0.140586, 0.059409, 0.032204, 0.109909, 0.140586, &
    0.140586, 0.059409, 0.059409, 0.059409, 0.084341, 0.022130, 0.059409, 0.028388, 0.059409, 0.059409, &
    0.207788, 0.059409, 0.054089, 0.015998, 0.012779, 0.059409, 0.023678, 0.036767, 0.035292, 0.148177, &
    0.054089, 0.028388, 0.012594, 0.059409, 0.059409, 0.095822, 0.140586, 0.059409, 0.017167, 0.059409, &
    0.192700, 0.012779, 0.141574, 0.028388, 0.032204, 0.017364, 0.148177, 0.220226, 0.054089, 0.148177, &
    0.192700, 0.095822, 0.084341, 0.050034, 0.032204, 0.227743, 0.095822, 0.059409, 0.028388, 0.140586 /)

    factors(ixO^S) = one

    if (theat > 0) then
      tloc = floor(theat/(300.d0/unit_time)) ! about five minutes?
      tres = (theat-300.d0/unit_time*dble(tloc))/(300.d0/unit_time)
      xd1 = heatx_left(tloc+1)*(1.0-tres)+heatx_left(tloc+2)*tres
      yd1 = heaty_left(tloc+1)*(1.0-tres)+heaty_left(tloc+2)*tres
      xd2 = heatx_right(tloc+1)*(1.0-tres)+heatx_right(tloc+2)*tres
      yd2 = heaty_right(tloc+1)*(1.0-tres)+heaty_right(tloc+2)*tres
      factors(ixO^S) = dexp(-(x(ixO^S,1)-xd1)**2/sigma**2-(x(ixO^S,2)-yd1)**2/sigma**2)+&
                      dexp(-(x(ixO^S,1)-xd2)**2/sigma**2-(x(ixO^S,2)-yd2)**2/sigma**2)
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
