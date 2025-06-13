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
  integer :: i_Te,i_v2,i_v3,i_ion,i_T,i_b
  double precision :: dh,Bh,dv,Bv,xv,cv,ch
  integer :: Btot_
  logical :: ffhd_Btot
contains

 subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ a0,Qbeta,Qalfa,Qgama,bQ0,tstop,zh,lQ,lambdah,sigma,asym,F1,dh,Bh,dv,Bv,xv,ffhd_Btot
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

    i_v2 = var_set_extravar("m2", "m2")
    i_v3 = var_set_extravar("m3", "m3")
    i_Te = var_set_extravar("Te", "Te")
    i_ion = var_set_extravar("ion", "ion")
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
    case(1,2)
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

    double precision :: lQgrid(ixI^S),tr
    double precision :: sigmax,xr,xl,yr,yl
    double precision :: tramp
    double precision :: heatx(100),heaty(100),heatz(100),theat
    integer :: tloc
    double precision :: tres,xd,yd,zd,al,ar

    al=two*asym/(1.0+asym)
    ar=two*1.0/(1.0+asym)
    if(rad-zh+z_cir<0) call mpistop("error in getlQ")
    xr=sqrt(rad**2-(zh-z_cir)**2)
    xl=-xr
    yr=zero
    yl=zero
    sigmax=a0*sigma
    tramp=500.d0/unit_time !< time for localized heating to increase from 0 to 1
    theat=qt-trelax
    tr=zero
    if(theat .lt. tramp) then
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
      heatx = (/ &
       -0.65454899,  0.28458194, -0.22744813, -0.11025966,  0.12486398,& 
        1.67986545,  0.62547323, -0.78672243, -0.64661351,  0.57655719,&
       -0.05092373, -1.29963115, -0.24798119,  0.04078717, -0.39009043,&
        0.45614569, -0.2929759 ,  0.37962491,  0.0503413 ,  1.67795815,&
       -0.01763759, -0.47948665, -0.2771821 ,  0.5771912 , -1.91192736,&
       -0.50666706,  1.38760602,  0.39948909,  1.48355513,  1.00806203,&
        0.55092868,  1.4254353 ,  0.43811228, -2.14842539, -0.3999542 ,&
        1.38132305, -0.17225423, -0.48298547, -0.9487324 , -2.76969865,&
        0.86106731, -0.0675009 , -0.51425981,  0.64578342, -0.43522053,&
        0.48742831,  0.01157438,  0.29103381,  0.12189464, -0.47781916,&
       -0.79409682,  1.60402155,  0.17687554,  0.60651735,  1.35122346,&
       -0.70989718,  0.66557593,  1.20654423, -0.525087  , -2.11950442,&
        0.0492665 ,  0.54802776, -0.27578264,  0.54574691,  1.10856832,&
       -0.08182947, -0.68410033, -1.21231785, -1.34515078,  0.77129213,&
       -2.38824713,  0.47587433,  0.33316896, -0.97460803, -2.40390393,&
       -0.00940567, -0.77949011, -0.26671203,  0.85338401,  0.55620533,&
       -1.33480893, -0.0424146 ,  1.43181972,  0.39396895, -0.05333105,&
        0.79357678, -2.40365088, -0.41255822,  0.03111715,  0.87507241,&
        0.68403418, -0.03348593,  0.67557926,  1.02679693, -2.06213923,&
       -1.41086825, -2.7759568 , -0.38200447, -0.98147183, -0.26921884 /)
      heaty = (/ &
        0.20620114,  1.57571623,  0.99273184, -0.81739949, -0.09028757,& 
        0.26852579, -0.16850966,  1.50546915, -0.720076  , -0.31272231,&
        1.25941132,  0.26672533, -0.7914178 , -0.33946165, -1.22503327,&
        1.43180163,  0.6337601 , -2.71806572, -0.07684867, -0.07674896,&
       -1.08964026, -1.04476894, -0.27188359,  1.63424976,  0.77022016,&
        0.96704038, -0.26433368, -0.10366884, -0.27514909,  1.07046762,&
        0.17439194, -0.33976762, -1.77717898,  0.10208177,  1.1474956 ,&
        0.00932566,  0.15594628, -0.92199679,  1.82042009, -0.61959455,&
       -0.33204945, -0.62671587,  1.48692591, -0.29138805,  0.07887798,&
       -0.03801963, -0.35773417, -0.24258521, -0.54199331,  1.14443916,&
       -0.20816858,  0.93649641, -0.25490412,  0.33826083, -0.46534694,&
        1.35811135, -1.79826145,  0.16643203, -0.27804007, -0.25578936,&
        1.52704993,  0.04188611, -1.06840646, -0.27185478, -0.82089466,&
       -0.923261  ,  0.87314006, -0.80073293,  1.49824823, -0.33047469,&
       -0.74968787,  0.50871257,  2.00046146, -0.85082335,  0.13259392,&
       -0.41552142, -0.17040082, -0.39202753, -0.13708549,  0.4035436 ,&
        0.83317583, -1.22049835, -0.81478837,  1.55131072,  0.62158102,&
       -1.38449911, -0.37726387,  0.03014308, -1.49904103, -0.60172849,&
       -0.37537371, -0.67739449,  1.22787828, -1.61610727,  0.18564605,&
        0.20644934,  0.0350815 , -1.60732156,  0.7852108 ,  0.87321395 /)
      heatz = (/ &
       0.73067589, 0.35482033, 0.26675653, 0.93757457, 0.08839678,& 
       0.63908825, 0.2893326 , 0.65893077, 0.97322695, 0.27295517,&
       0.36166088, 0.68992796, 0.42371476, 0.5685879 , 0.53048824,&
       0.94806443, 0.73030684, 0.77324391, 0.16685977, 0.18081555,&
       0.50025676, 0.43265277, 0.3716569 , 0.09890832, 0.37196836,&
       0.81442104, 0.14057131, 0.39876124, 0.58151049, 0.2567612 ,&
       0.43387366, 0.62050069, 0.66587827, 0.05283403, 0.87520463,&
       0.69839553, 0.96411771, 0.91098996, 0.38926923, 0.45287646,&
       0.2983996 , 0.06681691, 0.64508937, 0.48900514, 0.19094137,&
       0.60272712, 0.88487611, 0.28261415, 0.57149667, 0.85670345,&
       0.70335735, 0.5563727 , 0.12679228, 0.29209735, 0.14028817,&
       0.79833874, 0.00730277, 0.93513725, 0.20295754, 0.62338279,&
       0.18882757, 0.73370351, 0.78643077, 0.74508176, 0.00613495,&
       0.16170847, 0.11445703, 0.13424633, 0.77818848, 0.36289674,&
       0.8606341 , 0.00463685, 0.63241199, 0.6198887 , 0.20494183,&
       0.95472318, 0.30358861, 0.91527117, 0.7390046 , 0.98781702,&
       0.43516623, 0.44490055, 0.5485382 , 0.93965329, 0.38029733,&
       0.83875486, 0.39907584, 0.81057131, 0.03176455, 0.93005177,&
       0.53639446, 0.45443253, 0.21837997, 0.69120038, 0.94786863,&
       0.86016832, 0.31628997, 0.68816025, 0.24262838, 0.11878956 /)

      if(theat .gt. zero) then
        tloc=floor(theat/(300.d0/unit_time))
        tres=(theat-300.d0/unit_time*dble(tloc))/(300.d0/unit_time)
        if(tres .lt. zero .or. tres .gt. one) call mpistop("tres wrong")
        xd=heatx(tloc+1)*(1.0-tres)+heatx(tloc+2)*tres
        yd=heaty(tloc+1)*(1.0-tres)+heaty(tloc+2)*tres
        zd=heatz(tloc+1)*(1.0-tres)+heatz(tloc+2)*tres+half
        lQgrid(ixO^S)=lQgrid(ixO^S)*zd*(dexp(-(x(ixO^S,1)-xr-xd*a0*0.25d0)**2/sigmax**2-(x(ixO^S,2)-yd*a0*0.25d0)**2/sigmax**2))+&
                      lQgrid(ixO^S)*(dexp(-(x(ixO^S,1)-xl)**2/sigmax**2-(x(ixO^S,2)-zero)**2/sigmax**2))
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

    double precision :: r1(ixI^S)
    
    ! r1(ixO^S)=sqrt(x(ixO^S,1)**2+(x(ixO^S,3)+d_cha0)**2)

    select case(refine_max_level)
    case(1,2)
      if(any(abs(x(ixO^S,3)) .lt. 0.2d0)) then
        refine=1
        coarsen=-1
      end if
    case(3)
      if(level .ge. 2) then
        refine=-1
      endif
      if(any(abs(x(ixO^S,3)) .lt. 0.2d0)) then
        refine=1
        coarsen=-1
      end if
    case(4)
      if(level .ge. 2) then
        refine=-1
      endif
      if(level .eq. 4) then
        coarsen=1
      endif
      if(any(abs(x(ixO^S,3)) .lt. 0.2d0)) then
        refine=1
        coarsen=-1
      end if
    case(5)
      if(level .ge. 2) then
        refine=-1
      endif
      if(level .ge. 3) then
        coarsen=1
      endif
      if(any(abs(x(ixO^S,3)) .lt. 0.2d0)) then
        refine=1
        coarsen=-1
      end if
    case default
      call mpistop("refine_max_level too high, not defined")
    end select

    ! if(any(abs(x(ixO^S,2)) .lt. a0 .and. abs(r1(ixO^S)-rad) .lt. a0)) then
    !   refine=1
    !   coarsen=-1
    ! endif
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

    w(ixO^S,i_v2)=zero
    w(ixO^S,i_v3)=zero
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
    w(ixO^S,i_ion)=1.0
  end subroutine set_output_vars

end module mod_usr
