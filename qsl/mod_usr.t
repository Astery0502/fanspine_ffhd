module mod_usr
  use mod_mf
  use mod_qsl
  implicit none
  double precision :: dh,Bh,dv,Bv,xv,cv,ch
contains

 subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ dh,Bh,dv,Bv,xv
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
    usr_special_convert => custom_processing

    call usr_params_read(par_files)
    call mf_activate()
    call qsl_init()
  end subroutine usr_init

  subroutine initglobaldata_usr()

    ! for fan-spine field use in wyper2016a_field
    cv = (Bv*dv**3.d0)*half
    ch = (Bh*dh**3.d0)*half
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    !> initialize one grid
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: res
    integer :: ix^D,na

    call wyper2016a_field(ixI^L,ixO^L,x,w(ixI^S,mag(1):mag(3)))
    w(ixO^S,mom(:))=zero

  end subroutine initonegrid_usr

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

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    integer, intent(in) :: ixO^L,iB,ixI^L
    double precision, intent(in) :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: ix^D

    select case(iB)
    case(1)
      w(ixO^S,mom(:))=zero
      do ix3=ixOmin3,ixOmax3
      do ix2=ixOmin2,ixOmax2
      do ix1=ixOmax1,ixOmin1,-1
        w(ix^D,mag(:))=(w(ix1+3,ix2,ix3,mag(:)) &
                  -5.d0*w(ix1+2,ix2,ix3,mag(:)) &
                  +7.d0*w(ix1+1,ix2,ix3,mag(:)))/3.d0
      enddo
      enddo
      enddo
      if(mf_glm) w(ixO^S,psi_)=0.d0
    case(2)
      w(ixO^S,mom(:))=zero
      do ix3=ixOmin3,ixOmax3
      do ix2=ixOmin2,ixOmax2
      do ix1=ixOmin1,ixOmax1
        w(ix^D,mag(:))=(w(ix1-3,ix2,ix3,mag(:)) &
                  -5.d0*w(ix1-2,ix2,ix3,mag(:)) &
                  +7.d0*w(ix1-1,ix2,ix3,mag(:)))/3.d0
      enddo
      enddo
      enddo
      if(mf_glm) w(ixO^S,psi_)=0.d0
    case(3)
      w(ixO^S,mom(:))=zero
      do ix3=ixOmin3,ixOmax3
      do ix2=ixOmax2,ixOmin2,-1
      do ix1=ixOmin1,ixOmax1
        w(ix^D,mag(:))=(w(ix1,ix2+3,ix3,mag(:)) &
                  -5.d0*w(ix1,ix2+2,ix3,mag(:)) &
                  +7.d0*w(ix1,ix2+1,ix3,mag(:)))/3.d0
      enddo
      enddo
      enddo
      if(mf_glm) w(ixO^S,psi_)=0.d0
    case(4)
      w(ixO^S,mom(:))=zero
      do ix3=ixOmin3,ixOmax3
      do ix2=ixOmin2,ixOmax2
      do ix1=ixOmin1,ixOmax1
        w(ix^D,mag(:))=(w(ix1,ix2-3,ix3,mag(:)) &
                  -5.d0*w(ix1,ix2-2,ix3,mag(:)) &
                  +7.d0*w(ix1,ix2-1,ix3,mag(:)))/3.d0
      enddo
      enddo
      enddo
      if(mf_glm) w(ixO^S,psi_)=0.d0
    case(5)
      do ix2=ixOmin2,ixOmax2
      do ix1=ixOmin1,ixOmax1
        w(ix1,ix2,ixOmin3:ixOmax3,mom(1))=-w(ix1,ix2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(1))
        w(ix1,ix2,ixOmin3:ixOmax3,mom(2))=-w(ix1,ix2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(2))
        w(ix1,ix2,ixOmin3:ixOmax3,mom(3))=-w(ix1,ix2,ixOmax3+nghostcells:ixOmax3+1:-1,mom(3))
      enddo
      enddo
      do ix3=ixOmax3,ixOmin3,-1
      do ix2=ixOmin2,ixOmax2
      do ix1=ixOmin1,ixOmax1
        w(ix^D,mag(:))=0.12d0*w(ix1,ix2,ix3+5,mag(:)) &
                      -0.76d0*w(ix1,ix2,ix3+4,mag(:)) &
                      +2.08d0*w(ix1,ix2,ix3+3,mag(:)) &
                      -3.36d0*w(ix1,ix2,ix3+2,mag(:)) &
                      +2.92d0*w(ix1,ix2,ix3+1,mag(:))
      enddo
      enddo
      enddo
      if(mf_glm) w(ixO^S,psi_)=0.d0
    case(6)
      w(ixO^S,mom(:))=zero
      do ix3=ixOmin3,ixOmax3
      do ix2=ixOmin2,ixOmax2
      do ix1=ixOmin1,ixOmax1
        w(ix^D,mag(:))=(w(ix1,ix2,ix3-3,mag(:)) &
                  -5.d0*w(ix1,ix2,ix3-2,mag(:)) &
                  +7.d0*w(ix1,ix2,ix3-1,mag(:)))/3.d0
      enddo
      enddo
      enddo
      if(mf_glm) w(ixO^S,psi_)=0.d0
    case default
      call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine specialbound_usr

  subroutine custom_processing(qunitconvert)
    use mod_global_parameters
    integer, intent(in) :: qunitconvert
    character(len=20)   :: userconvert_type

    call qsl_preprocessing()
    call trace_Bfield()
  end subroutine custom_processing
end module mod_usr
