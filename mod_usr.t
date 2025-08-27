module mod_usr
  use mod_ffhd
  implicit none
  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav,trelax,tstop,temax
  double precision :: heatunit,gzone,SRadius,dya
  double precision, allocatable :: pa(:),ra(:)
  integer, parameter :: jmax=50000

  integer, parameter :: np=10 !> num of points along the ring
  double precision :: a0 !> radius of the cross section
  double precision :: F0,F1,rad,qc,Bzp,x_cir,z_cir
  double precision :: L_cha0,d_cha0,q_cha0,zh,lQ,lQ0,lambdah,sigma,asym
  double precision :: xs0(np,3)

  integer :: numL1,numL2,numFL,numLP,nwP,nwL
  integer :: ngrid_ground,ngrid_leng
  integer, allocatable :: ground_grid(:),leng_grid(:)
  double precision :: dL,max_len,finegrid
  double precision :: bQ0,Qbeta,Qalfa,Qgama
  integer :: i_Te,i_b1,i_b2,i_b3
  double precision :: dh,Bh,dv,Bv,xv,cv,ch
  integer :: Btot_
  logical :: ffhd_Btot

  !> parameters for special amr over specific points
  character(len=std_len) :: csvfile
  integer :: npoints=0
  double precision, allocatable :: pos(:,:)
contains

 subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ a0,Qbeta,Qalfa,Qgama,bQ0,tstop,zh,lQ,lambdah,sigma,asym,F1,&
                        dh,Bh,dv,Bv,xv,ffhd_Btot,csvfile,npoints
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
    i_Te = var_set_extravar("Te", "Te")
  end subroutine usr_init

  subroutine initglobaldata_usr()
    heatunit=unit_pressure/unit_time !< 3.697693390805347E-003 erg*cm^-3/s
    usr_grav=-2.74d4*unit_length/unit_velocity**2 !< solar gravity
    gzone=0.5d0 !< thickness of a ghostzone below the bottom boundary
    dya=(2.d0*gzone+xprobmax3-xprobmin3)/dble(jmax) !< cells size of high-resolution 1D solar atmosphere
    SRadius=69.61d0 !< Solar radius
    trelax=100.d0 !< time for relaxation
    temax=2.d0 !< artificial maximum temperature for the top boundary

    !> parameters for magnetic flux rope
    rad=8.d0          !> radius of the ring
    x_cir=zero
    z_cir=-4.5d0
    Bzp=1.63d0        !> B_perp along the isocontour
    F0=F1*6.0*sqrt(2.0)/5.0*dpi*a0*rad*Bzp/(log(8.0*rad/a0)-25.0/24.0)
    xs0=0.d0
    ! call calc_geom(xs0)
    !> parameters for background potential field
    L_cha0=3.d0  !> distance between p/m
    d_cha0=4.5d0 !> depth of the dipoles
    qc=600.d0    !> strength of the dipoles
    q_cha0=qc/sqrt(4*dpi)
    lQ0=lQ/heatunit
    !>
    finegrid=1
    dL=min((xprobmax3-xprobmin3)/domain_nx3,&
           (xprobmax2-xprobmin2)/domain_nx2,&
           (xprobmax1-xprobmin1)/domain_nx1)/2**(refine_max_level-1)
    max_len=sqrt((xprobmax1-xprobmin1)**2+(xprobmax2-xprobmin2)**2+(xprobmax3-xprobmin3)**2)*3.d0
    numLP=floor(max_len/dL)
    numL1=floor((xprobmax1-xprobmin1)/dL/finegrid)
    numL2=floor((xprobmax2-xprobmin2)/dL/finegrid)
    numFL=numL1*numL2
    call inithdstatic

    ! for fan-spine field use in wyper2016a_field
    cv = (Bv*dv**3.d0)*half
    ch = (Bh*dh**3.d0)*half

    ! initialize possible points for special amr
    if (npoints > 0) then
      allocate(pos(3, npoints))
      call read_csv(npoints, pos, csvfile)
    end if
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

    double precision :: lQgrid(ixI^S),factors(ixI^S),tr
    double precision :: sigmax,xr,xl,yr,yl
    double precision :: tramp
    double precision :: heatx(100),heaty(100),theat
    double precision :: r(200),theta(200),heatz(200),rd,thetad
    integer :: tloc, i_tloc
    double precision :: tres,xd,yd,zd,al,ar

    factors=zero
    al=two*asym/(1.0+asym)
    ar=two*1.0/(1.0+asym)
    xr=-2.9
    xl=-7.5
    yr=zero
    yl=zero
    sigmax=a0*sigma
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
    lQgrid(ixO^S)=lQ0*tr*dexp(-(x(ixO^S,3)-zh)**2/lambdah**2)
    select case(iprob)
    case(1)
      lQgrid(ixO^S)=lQgrid(ixO^S)*(dexp(-(x(ixO^S,1)-xr)**2/sigmax**2-(x(ixO^S,2)-yr)**2/sigmax**2)*ar+&
                                   dexp(-(x(ixO^S,1)-xl)**2/sigmax**2-(x(ixO^S,2)-yl)**2/sigmax**2)*al)
    case(2)
      ! 9.25 = |B|max
      lQgrid(ixO^S)=lQgrid(ixO^S)*(block%wextra(ixO^S,Btot_)/9.25d0)**2
      lQgrid(ixO^S)=lQgrid(ixO^S)*((ar-al)/(xr-xl)*(x(ixO^S,1)-xl)+al)
    case(3)

      r = (/ & 
        1.874540,  2.450714,  2.231994,  2.098658,  1.656019, &
        1.655995,  1.558084,  2.366176,  2.101115,  2.208073, &
        1.520584,  2.469910,  2.332443,  1.712339,  1.681825, &
        1.683405,  1.804242,  2.024756,  1.931945,  1.791229, &
        2.111853,  1.639494,  1.792145,  1.866362,  1.956070, &
        2.285176,  1.699674,  2.014234,  2.092415,  1.546450, &
        2.107545,  1.670524,  1.565052,  2.448886,  2.465632, &
        2.308397,  1.804614,  1.597672,  2.184233,  1.940152, &
        1.622038,  1.995177,  1.534389,  2.409320,  1.758780, &
        2.162522,  1.811711,  2.020068,  2.046710,  1.684854, &
        2.469585,  2.275133,  2.439499,  2.394827,  2.097900, &
        2.421874,  1.588493,  1.695983,  1.545227,  1.825330, &
        1.888677,  1.771349,  2.328738,  1.856753,  1.780935, &
        2.042696,  1.640924,  2.302197,  1.574551,  2.486887, &
        2.272245,  1.698716,  1.505522,  2.315461,  2.206857, &
        2.229007,  2.271270,  1.574045,  1.858466,  1.615869, &
        2.363103,  2.123298,  1.830898,  1.563558,  1.810982, &
        1.825183,  2.229606,  2.137557,  2.387213,  1.972215, &
        1.619594,  2.213245,  2.260785,  2.061277,  2.270967, &
        1.993796,  2.022733,  1.927541,  1.525419,  1.607891, &
        1.531429,  2.136410,  1.814356,  2.008571,  2.407566, &
        1.749292,  1.910383,  2.255551,  1.728798,  1.576980, &
        1.789751,  1.661221,  2.429698,  2.308120,  2.133404, &
        2.371461,  2.303672,  1.686570,  2.392559,  2.039342, &
        2.307440,  2.396091,  1.818003,  1.610052,  1.727935, &
        1.927108,  2.318015,  2.360731,  1.506952,  2.010747, &
        1.917411,  1.722108,  1.619865,  1.837615,  2.442910, &
        1.823203,  2.018791,  2.203019,  1.863630,  2.471782, &
        2.462447,  1.751782,  1.997249,  1.800878,  1.784840, &
        1.536887,  2.109564,  2.002679,  1.551479,  1.778646, &
        2.408266,  1.739562,  1.644895,  1.989453,  2.485650, &
        1.742055,  2.172136,  2.261620,  1.737638,  2.228216, &
        1.867783,  2.132306,  2.133530,  2.035775,  1.590290, &
        2.335302,  1.820780,  1.686519,  1.540775,  2.090893, &
        2.177564,  1.516588,  2.012093,  1.726496,  2.145173, &
        1.674366,  2.190938,  1.886735,  2.436730,  1.637521, &
        1.841066,  1.613474,  2.424694,  2.377339,  1.757942, &
        2.159984,  2.317222,  2.055201,  2.029651,  1.741852, &
        1.593103,  2.397216,  2.400418,  2.133101,  1.839030, &
        1.849210,  2.225956,  2.397110,  2.387086,  2.279876 /)

        theta = (/ & 
        4.034004,  0.528667,  1.015543,  5.645782,  3.810306, &
        0.057787,  0.637565,  4.168905,  0.031803,  1.010387, &
        3.447796,  4.347306,  4.096393,  1.409126,  4.474754, &
        1.490680,  2.044547,  4.690344,  4.081764,  5.335828, &
        4.131904,  3.570788,  0.588576,  2.310427,  1.666316, &
        1.533032,  6.113606,  2.469906,  5.604894,  3.965561, &
        4.993947,  3.158162,  3.624794,  3.094580,  1.226748, &
        4.539301,  1.764145,  0.152782,  4.055622,  1.112819, &
        5.909076,  5.993710,  5.748262,  2.325776,  0.097117, &
        5.832798,  2.690360,  6.073671,  6.054603,  5.359616, &
        1.850077,  2.419640,  5.347849,  1.991280,  1.064954, &
        3.498486,  5.882034,  4.373284,  3.581800,  0.610578, &
        3.864204,  6.220692,  0.880174,  3.256761,  5.512698, &
        4.654386,  4.379479,  4.413838,  2.258750,  1.844692, &
        5.085366,  5.090093,  5.447976,  5.738060,  3.212859, &
        3.151120,  5.015837,  4.083844,  4.410588,  5.000113, &
        5.592068,  2.123686,  2.359857,  0.590506,  3.633441, &
        0.225832,  2.925439,  3.409537,  1.800392,  3.712315, &
        0.191639,  0.234666,  5.168552,  2.263145,  0.798345, &
        3.281351,  4.838012,  1.356044,  3.913736,  0.536254, &
        0.324726,  3.338600,  3.396911,  4.005090,  4.562166, &
        6.131459,  3.244011,  2.029195,  4.996302,  1.701689, &
        2.758139,  0.492956,  0.159283,  6.048498,  5.252618, &
        4.372935,  2.569527,  1.088840,  0.982923,  1.572323, &
        3.450893,  4.489939,  4.148142,  1.758877,  5.999596, &
        4.636343,  3.483109,  3.843555,  2.636425,  1.556540, &
        2.236642,  4.761688,  0.090437,  0.729306,  0.289043, &
        0.255907,  5.375017,  4.421213,  2.979322,  0.614710, &
        3.088914,  2.974911,  1.088259,  2.725970,  2.503879, &
        3.869500,  3.990411,  0.284653,  2.353760,  3.932394, &
        3.161298,  5.381484,  4.138694,  1.023747,  0.443397, &
        4.036439,  0.166575,  3.680537,  5.907641,  3.615811, &
        2.438944,  4.041899,  2.879288,  3.428211,  5.915398, &
        2.425954,  6.039338,  5.688486,  1.230192,  0.435810, &
        0.633207,  0.114491,  0.593403,  4.291458,  0.447291, &
        2.004183,  5.308508,  0.146222,  5.117456,  1.770946, &
        0.742452,  4.377729,  3.951764,  5.513319,  4.618588, &
        5.048420,  1.772075,  1.114886,  4.716252,  5.069492, &
        6.223527,  2.592553,  2.337459,  4.878347,  2.141332, &
        5.848121,  5.393566,  2.695449,  4.717862,  4.740933 /)

        heatz = (/ & 
        0.103124,  0.902553,  0.505252,  0.826457,  0.320050, &
        0.895523,  0.389202,  0.010838,  0.905382,  0.091287, &
        0.319314,  0.950062,  0.950607,  0.573438,  0.631837, &
        0.448446,  0.293211,  0.328665,  0.672518,  0.752375, &
        0.791579,  0.789618,  0.091206,  0.494420,  0.057559, &
        0.549529,  0.441531,  0.887704,  0.350915,  0.117067, &
        0.142992,  0.761511,  0.618218,  0.101123,  0.084107, &
        0.700969,  0.072763,  0.821860,  0.706242,  0.081349, &
        0.084838,  0.986640,  0.374271,  0.370642,  0.812800, &
        0.947249,  0.986001,  0.753378,  0.376260,  0.083501, &
        0.777147,  0.558404,  0.424222,  0.906354,  0.111197, &
        0.492625,  0.011354,  0.468661,  0.056303,  0.118818, &
        0.117526,  0.649210,  0.746045,  0.583369,  0.962173, &
        0.374871,  0.285712,  0.868599,  0.223596,  0.963223, &
        0.012154,  0.969879,  0.043160,  0.891143,  0.527701, &
        0.992965,  0.073797,  0.553854,  0.969303,  0.523098, &
        0.629399,  0.695749,  0.454541,  0.627558,  0.584314, &
        0.901158,  0.045446,  0.280963,  0.950411,  0.890264, &
        0.455657,  0.620133,  0.277381,  0.188121,  0.463698, &
        0.353352,  0.583656,  0.077735,  0.974395,  0.986211, &
        0.698162,  0.536096,  0.309528,  0.813795,  0.684731, &
        0.162617,  0.910927,  0.822537,  0.949800,  0.725720, &
        0.613415,  0.418243,  0.932728,  0.866064,  0.045219, &
        0.026367,  0.376463,  0.810553,  0.987276,  0.150417, &
        0.594131,  0.380891,  0.969914,  0.842119,  0.838329, &
        0.468693,  0.414820,  0.273407,  0.056375,  0.864722, &
        0.812901,  0.999718,  0.996637,  0.555432,  0.768987, &
        0.944766,  0.849647,  0.247348,  0.450544,  0.129159, &
        0.954051,  0.606175,  0.228643,  0.671701,  0.618128, &
        0.358163,  0.113558,  0.671573,  0.520308,  0.772318, &
        0.520164,  0.852182,  0.551907,  0.560938,  0.876654, &
        0.403483,  0.134015,  0.028783,  0.755137,  0.620310, &
        0.704080,  0.212964,  0.136371,  0.014545,  0.350588, &
        0.589918,  0.392244,  0.437475,  0.904159,  0.348255, &
        0.513989,  0.783653,  0.396543,  0.622087,  0.862364, &
        0.949521,  0.147073,  0.926588,  0.492116,  0.258244, &
        0.459136,  0.980033,  0.492618,  0.328752,  0.633401, &
        0.240146,  0.075863,  0.128880,  0.128046,  0.151903, &
        0.138827,  0.640875,  0.181880,  0.345667,  0.896788, &
        0.473962,  0.667558,  0.172320,  0.192289,  0.040869 /)

      if(theat .gt. zero) then
        tloc=floor(theat/(300.d0/unit_time))
        tres=(theat-300.d0/unit_time*dble(tloc))/(300.d0/unit_time)
        if(tres .lt. zero .or. tres .gt. one) call mpistop("tres wrong")
        ! xd=heatx(tloc+1)*(1.0-tres)+heatx(tloc+2)*tres
        ! yd=heaty(tloc+1)*(1.0-tres)+heaty(tloc+2)*tres
        do i_tloc=1,50
          rd=r(tloc+i_tloc)*(1.0-tres)+r(tloc+i_tloc+1)*tres
          thetad=theta(tloc+i_tloc)*(1.0-tres)+theta(tloc+i_tloc+1)*tres
          zd=heatz(tloc+i_tloc)*(1.0-tres)+heatz(tloc+i_tloc+1)*tres+half
          factors(ixO^S)=factors(ixO^S)+zd*dexp(-(x(ixO^S,1)-xv-rd*cos(thetad))**2/sigmax**2-(x(ixO^S,2)-rd*sin(thetad))**2/sigmax**2)
        end do
        lQgrid(ixO^S)=lQgrid(ixO^S)*factors(ixO^S)*(block%wextra(ixO^S,Btot_)/20)**2
        ! lQgrid(ixO^S)=lQgrid(ixO^S)*zd*(dexp(-(x(ixO^S,1)-xr-xd*a0*0.25d0)**2/sigmax**2-(x(ixO^S,2)-yd*a0*0.25d0)**2/sigmax**2))+&
        !               lQgrid(ixO^S)*(dexp(-(x(ixO^S,1)-xl)**2/sigmax**2-(x(ixO^S,2)-zero)**2/sigmax**2))
      endif
    case default
      call mpistop("iprob wrong xxxxxxxxxxxxxx")
    end select
  end subroutine getlQ

  subroutine calc_geom(xs)
    double precision, dimension(np,3), intent(inout) :: xs

    integer :: i
    double precision :: theta

    xs=zero 
    do i = 1,np
      theta = (i-0.5)*2*dpi/np
      xs(i,1) = x_cir + rad * cos(theta)
      xs(i,3) = z_cir + rad * sin(theta)
    end do
  end subroutine calc_geom

  subroutine rbslB(ixI^L,ixO^L,ixK^L,x,a,F,xs,Bout)
    integer, intent(in) :: ixI^L,ixO^L,ixK^L
    double precision, intent(in) :: a,F
    double precision, dimension(np,3), intent(in) :: xs
    double precision, dimension(ixI^S,3), intent(in) :: x
    double precision, dimension(ixI^S,3), intent(inout) :: Bout
    !> local
    integer :: ix^D,ixp,idirmin
    double precision :: re_pi,I_cur,r_mag,dl,KIr,KFr
    double precision, dimension(3) :: Rpl,r_vec,Rcr
    double precision, dimension(ixK^S,3) :: xL,BfrI,BfrF,AIx,AFx

    re_pi=1.d0/dpi
    !> (eq.12)
    I_cur=-1.0*F*5.0d0*sqrt(2.0d0)/3.0d0/a
    AIx=0.d0
    AFx=0.d0

    xL=zero
    xL(ixO^S,1:3)=x(ixO^S,1:3)
    ix1=ixKmin1
    do ix2=ixOmin2,ixOmax2
    do ix3=ixOmin3,ixOmax3
      xL(ix^D,1)=2.d0*xL(ixOmin1,ix2,ix3,1)-xL(ixOmin1+1,ix2,ix3,1)
      xL(ix^D,2)=xL(ixOmin1,ix2,ix3,2)
      xL(ix^D,3)=xL(ixOmin1,ix2,ix3,3)
    end do
    end do
    ix1=ixKmax1
    do ix2=ixOmin2,ixOmax2
    do ix3=ixOmin3,ixOmax3
      xL(ix^D,1)=2.d0*xL(ixOmax1,ix2,ix3,1)-xL(ixOmax1-1,ix2,ix3,1)
      xL(ix^D,2)=xL(ixOmax1,ix2,ix3,2)
      xL(ix^D,3)=xL(ixOmax1,ix2,ix3,3)
    end do
    end do
    ix2=ixKmin2
    do ix1=ixKmin1,ixKmax1
    do ix3=ixOmin3,ixOmax3
      xL(ix^D,1)=xL(ix1,ixOmin2,ix3,1)
      xL(ix^D,2)=2.d0*xL(ix1,ixOmin2,ix3,2)-xL(ix1,ixOmin2+1,ix3,2)
      xL(ix^D,3)=xL(ix1,ixOmin2,ix3,3)
    end do
    end do
    ix2=ixKmax2
    do ix1=ixKmin1,ixKmax1
    do ix3=ixOmin3,ixOmax3
      xL(ix^D,1)=xL(ix1,ixOmax2,ix3,1)
      xL(ix^D,2)=2.d0*xL(ix1,ixOmax2,ix3,2)-xL(ix1,ixOmax2-1,ix3,2)
      xL(ix^D,3)=xL(ix1,ixOmax2,ix3,3)
    end do
    end do
    ix3=ixKmin3
    do ix1=ixKmin1,ixKmax1
    do ix2=ixKmin2,ixKmax2
      xL(ix^D,1)=xL(ix1,ix2,ixOmin3,1)
      xL(ix^D,2)=xL(ix1,ix2,ixOmin3,2)
      xL(ix^D,3)=2.d0*xL(ix1,ix2,ixOmin3,3)-xL(ix1,ix2,ixOmin3+1,3)
    end do
    end do
    ix3=ixKmax3
    do ix1=ixKmin1,ixKmax1
    do ix2=ixKmin2,ixKmax2
      xL(ix^D,1)=xL(ix1,ix2,ixOmax3,1)
      xL(ix^D,2)=xL(ix1,ix2,ixOmax3,2)
      xL(ix^D,3)=2.d0*xL(ix1,ix2,ixOmax3,3)-xL(ix1,ix2,ixOmax3-1,3)
    end do
    end do

    {do ix^DB=ixKmin^DB,ixKmax^DB\}
    do ixp=1,np
      !> position vector from source point to field point
      r_vec(1)=(xL(ix^D,1)-xs(ixp,1))/a
      r_vec(2)=(xL(ix^D,2)-xs(ixp,2))/a
      r_vec(3)=(xL(ix^D,3)-xs(ixp,3))/a
      r_mag=sqrt(r_vec(1)**2+r_vec(2)**2+r_vec(3)**2)
      !> calculate tangential vector of the axis
      if(ixp .eq. 1) then
        Rpl(1)=0.5d0*(xs(ixp+1,1)-xs(np,1))
        Rpl(2)=0.5d0*(xs(ixp+1,2)-xs(np,2))
        Rpl(3)=0.5d0*(xs(ixp+1,3)-xs(np,3))
      else if(ixp .eq. np) then
        Rpl(1)=0.5d0*(xs(1,1)-xs(ixp-1,1))
        Rpl(2)=0.5d0*(xs(1,2)-xs(ixp-1,2))
        Rpl(3)=0.5d0*(xs(1,3)-xs(ixp-1,3))
      else
        Rpl(1)=0.5d0*(xs(ixp+1,1)-xs(ixp-1,1))
        Rpl(2)=0.5d0*(xs(ixp+1,2)-xs(ixp-1,2))
        Rpl(3)=0.5d0*(xs(ixp+1,3)-xs(ixp-1,3))
      end if
      !> dl=sqrt(Rpl(1)**2+Rpl(2)**2+Rpl(3)**2)
      !> Rpl corss r_vec
      Rcr(1)=Rpl(2)*r_vec(3)-Rpl(3)*r_vec(2)
      Rcr(2)=Rpl(3)*r_vec(1)-Rpl(1)*r_vec(3)
      Rcr(3)=Rpl(1)*r_vec(2)-Rpl(2)*r_vec(1)
      if(r_mag .le. 1.d0) then
        !> (eq.13)
        KIr=2.d0*re_pi*(asin(r_mag)/r_mag+(5.d0-2.d0*r_mag**2)/3.d0*sqrt(1-r_mag**2))
        !> (eq.14)
        KFr=2.d0*re_pi/r_mag**2*(asin(r_mag)/r_mag-sqrt(1-r_mag**2))&
           +2.d0*re_pi*sqrt(1.d0-r_mag**2)+(5.d0-2.d0*(r_mag**2))*0.5d0/sqrt(6.d0)*(1.d0&
           -2.d0*re_pi*asin((1.d0+2.d0*r_mag**2)/(5.d0-2.d0*r_mag**2)))
        !> (eq.2)
        AIx(ix^D,:)=AIx(ix^D,:)+I_cur*0.25d0*re_pi*KIr*Rpl(:)/a
        !> (eq.3)
        AFx(ix^D,:)=AFx(ix^D,:)+F*0.25d0*re_pi*KFr*Rcr(:)/a**2
      else
        !> Jackson 1962
        KIr=1.d0/r_mag
        KFr=KIr**3
        !> (eq.2)
        AIx(ix^D,:)=AIx(ix^D,:)+I_cur*0.25d0*re_pi*KIr*Rpl(:)/a
        !> (eq.3)
        AFx(ix^D,:)=AFx(ix^D,:)+F*0.25d0*re_pi*KFr*Rcr(:)/a**2
      end if
    end do
    {end do\}
    call curl2(AIx,ixK^L,ixO^L,xL,BfrI)
    call curl2(AFx,ixK^L,ixO^L,xL,BfrF)
    Bout(ixO^S,1:3)=BfrI(ixO^S,1:3)+BfrF(ixO^S,1:3) 
  end subroutine rbslB

  subroutine bipoB(ixI^L,ixO^L,x,L_cha,d_cha,q_cha,Bout)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: L_cha,d_cha,q_cha
    double precision, dimension(ixI^S,3), intent(in) :: x
    double precision, dimension(ixI^S,3), intent(inout) :: Bout
    !> local
    double precision, dimension(ixI^S) :: rpv,rmv
    double precision, dimension(ixI^S,1:3) :: rplus,rminu

    Bout=zero
    rplus(ixO^S,1)=x(ixO^S,1)
    rminu(ixO^S,1)=x(ixO^S,1)
    rplus(ixO^S,2)=x(ixO^S,2)-L_cha
    rminu(ixO^S,2)=x(ixO^S,2)+L_cha
    rplus(ixO^S,3)=x(ixO^S,3)+d_cha
    rminu(ixO^S,3)=x(ixO^S,3)+d_cha
    rpv(ixO^S)=dsqrt(^D&rplus(ixO^S,^D)**2+)
    rmv(ixO^S)=dsqrt(^D&rminu(ixO^S,^D)**2+)
    Bout(ixO^S,1)=Bout(ixO^S,1)+q_cha&
      *(rplus(ixO^S,1)/rpv(ixO^S)**3-rminu(ixO^S,1)/rmv(ixO^S)**3)
    Bout(ixO^S,2)=Bout(ixO^S,2)+q_cha&
      *(rplus(ixO^S,2)/rpv(ixO^S)**3-rminu(ixO^S,2)/rmv(ixO^S)**3)
    Bout(ixO^S,3)=Bout(ixO^S,3)+q_cha&
      *(rplus(ixO^S,3)/rpv(ixO^S)**3-rminu(ixO^S,3)/rmv(ixO^S)**3)
  endsubroutine bipoB

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

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)!,ffBtot)
    !> output: wB0 must be a unit field vector
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)
    ! logical, intent(in) :: ffBtot
    !> local
    integer :: ixK^L
    ! double precision :: B1(ixI^S,1:ndir) !> background field
    ! double precision :: B2(ixI^S,1:ndir) !> flux rope
    double precision :: B3(ixI^S,1:ndir) !> wyper field
    double precision :: Btot(ixI^S)

    ! B1=zero
    ! B2=zero
    ixK^L=ixO^L^LADD1;
    ! call bipoB(ixI^L,ixO^L,x,L_cha0,d_cha0,q_cha0,B1)
    ! call rbslB(ixI^L,ixO^L,ixK^L,x,a0,F0,xs0,B2)
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

    if (npoints > 0) then
      do i=1,npoints
        if (pos(1,i) >= minval(x(ixO^S,1)) .and. pos(1,i) <= maxval(x(ixO^S,1)) .and. &
            pos(2,i) >= minval(x(ixO^S,2)) .and. pos(2,i) <= maxval(x(ixO^S,2)) .and. &
            pos(3,i) >= minval(x(ixO^S,3)) .and. pos(3,i) <= maxval(x(ixO^S,3))) then
          refine=1
          coarsen=-1
          exit
        else 
          refine=-1
          coarsen=-1
        end if
      end do
    else 
      if(level .ge. 2) then
        refine=-1
      endif
      if(level .ge. 3) then
        coarsen=1
      endif
    end if

    if(any(abs(x(ixO^S,3)) .lt. 0.3)) then
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
    double precision                :: divb(ixI^S)
    double precision                :: v0(ndim)
    double precision                :: bfield(ixI^S, ndir)
    double precision                :: ws(ixGs^T, nws)
    double precision :: wlocal(ixI^S,1:nw),pth(ixI^S)
    double precision :: wion
    double precision :: cor(ixI^S)
    integer :: ix^D

    w(ixO^S,i_b1)=block%B0(ixO^S,1,0)*block%wextra(ixO^S,Btot_)
    w(ixO^S,i_b2)=block%B0(ixO^S,2,0)*block%wextra(ixO^S,Btot_)
    w(ixO^S,i_b3)=block%B0(ixO^S,3,0)*block%wextra(ixO^S,Btot_)
    cor=zero
    where(x(ixO^S,3)>0.75d0)
      cor(ixO^S)=1.d0
    endwhere
    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    call ffhd_get_pthermal(wlocal,x,ixI^L,ixI^L,pth)
    ! if(ffhd_ionization) then
    !   {do ix^DB=ixOmin^DB,ixOmax^DB\}
    !     call ion_get_temp(wlocal(ix^D,rho_),pth(ix^D),w(ix^D,i_Te),w(ix^D,i_ion))
    !   {end do\}
    ! else
    !   w(ixO^S,i_Te)=pth(ixO^S)/w(ixO^S,rho_)
    !   w(ixO^S,i_ion)=1.0
    ! endif
    w(ixO^S,i_Te)=pth(ixO^S)/w(ixO^S,rho_)
  end subroutine set_output_vars

  subroutine read_csv(npoints, pos, ifile)
    character(len=std_len), intent(in) :: ifile
    integer, intent(in) :: npoints
    double precision, intent(out) :: pos(3, npoints)
    integer :: iunit, i, ierr
    character(len=std_len) :: line
    double precision :: x, y, z

    iunit = 10
    open(unit=iunit, file=ifile, status='old', action='read')
    read(iunit, *, end=100) pos(:,:)
    100 close(iunit)
  end subroutine read_csv

end module mod_usr
