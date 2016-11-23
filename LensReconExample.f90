!Simple example code to simulate a lensed CMB map, then use it to do temperature quadratic lensing reconstruction
!Not at all realistic (no masks, not even any simulated noise, only isotropic noise in the filter weights)
!AL Apr 2014

    module Recon
    use HealpixObj
    use HealpixVis
    use Random
    use spinalm_tools
    use IniFile
    use AMLUtils
    implicit none

    contains

    subroutine NoiseInit(Noise, NoiseVar, noise_fwhm_deg, lmax)
    Type(HealpixPower):: Noise
    real(dp) NoiseVar, noise_fwhm_deg
    real(sp) amp, xlc, sigma2
    integer lmax, L

    call healpixPower_Init(Noise,lmax,.false.)

    xlc= 180*sqrt(8.*log(2.))/3.14159
    sigma2 = (noise_fwhm_deg/xlc)**2
    do l=0, lmax
        Noise%Cl(l,1) = NoiseVar*exp(l*(l+1)*sigma2)
    end do

    end subroutine NoiseInit

    subroutine Getg(P, Noise, g, lmax_phi, lmax_est, est_type)
    Type(HealpixPower) :: P, Noise
    integer lmax_est, lmax_phi
    integer l1,l2,L, est_type
    complex(dp), allocatable :: threejs3(:), threejs5(:)
    complex(dp), dimension(:, :, :):: g
    allocate(threejs3(0:lmax_est + lmax_phi+1))

!    open (1, file = 'g.txt')
    if (est_type == 1)
    do L=1, lmax_phi
        do l1=2, lmax_est
            call GetThreeJs(threejs3(abs(l1-L)),l1,L,0,0)
            do l2 = max(2,abs(l1-L)), min(l1+L,lmax_est)
                if (threejs3(l2)==0 .or. mod(L+l1+l2, 2)) cycle
                g(L, l1, l2) = threejs3(l2)*sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(4*HO_fourpi))* &
                 (( real((L*(L+1) + l2*(l2+1)-l1*(l1+1)), dp)*P%Cl(l2,1)  ) + real((L*(L+1) - l2*(l2+1)+l1*(l1+1)), dp)*P%Cl(l1,1)) / &
                 ( 2*(P%Cl(l1,1) +Noise%Cl(l1,1))*(P%Cl(l2,1)+Noise%Cl(l2,1)))

 
!                write(1, *) L, l1, l2, g(L, l1, l2)
            end do
        end do    
    end do
    else if (est_type == 2) !EE
    do L=1, lmax_phi
        do l1=2, lmax_est
            call GetThreeJs(threejs3(abs(l1-L)),l1,L,-2,0)
            do l2 = max(2,abs(l1-L)), min(l1+L,lmax_est)
                if (threejs3(l2)==0 .or. mod(L+l1+l2, 2)) cycle
                g(L, l1, l2) =conjg(threejs3(l2))*sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(4*HO_fourpi))* &
                 (( real((L*(L+1) + l2*(l2+1)-l1*(l1+1)), dp)*P%Cl(l2,2)  ) + real((L*(L+1) - l2*(l2+1)+l1*(l1+1)), dp)*P%Cl(l1,2)) / &
                 ( 2*(P%Cl(l1,2) +Noise%Cl(l1,2))*(P%Cl(l2,2)+Noise%Cl(l2,2)))


!                write(1, *) L, l1, l2, g(L, l1, l2)
            end do
        end do
    end do
    else if (est_type == 3) !TE
    do L=1, lmax_phi
        do l1=2, lmax_est
            call GetThreeJs(threejs3(abs(l1-L)),l1,L,0,0)
            call GetThreeJs(threejs5(abs(l1-L)),l1,L,-2,0)
            do l2 = max(2,abs(l1-L)), min(l1+L,lmax_est)
                if (mod(L+l1+l2, 2)) cycle
                g(L, l1, l2) = sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(4*HO_fourpi))* &
                 (real(1-2*mod(L+l1+l2, 2) , dp) *conjg(threejs5(l2))* (real((L*(L+1) - l2*(l2+1)+l1*(l1+1)), dp)*P_tilde%Cl(l2,4)*P%(l2,1)*P%Cl(l1,2)- &
                                                                        real((L*(L+1) + l2*(l2+1)-l1*(l1+1)), dp)*P%Cl(l1,4)*P%Cl(l2,4)*P_tilde%Cl(l2, 4)) + &
                  threejs3(l2)*(real((L*(L+1) + l2*(l2+1)-l1*(l1+1)), dp)*P%Cl(l1,2)*P%Cl(l2,1)*P_tilde%Cl(l2, 4)+ &
                                real((L*(L+1) - l2*(l2+1)+l1*(l1+1)), dp)*P%Cl(l1,4)*P%Cl(l2,4)*P_tilde%Cl(l1, 4))) / &
                 ( (P%Cl(l1,1) +Noise%Cl(l1,1))*(P%Cl(l2,1)+Noise%Cl(l2,1))* (P%Cl(l1,2) +Noise%Cl(l1,2))*(P%Cl(l2,2)+Noise%Cl(l2,2))- &
                    (P%Cl(l1,4)*P%Cl(l2,4))**2)


!                write(1, *) L, l1, l2, g(L, l1, l2)
            end do
        end do
    end do
    else if (est_type == 4) ! TB 
    do L=1, lmax_phi
        do l1=2, lmax_est
            call GetThreeJs(threejs3(abs(l1-L)),l1,L,-2,0)
            do l2 = max(2,abs(l1-L)), min(l1+L,lmax_est)
                if (threejs3(l2)==0 .or. (1-mod(L+l1+l2, 2))) cycle
                g(L, l1, l2) = threejs3(l2)*sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(4*HO_fourpi)) &
                 (( cmplx(0, -1)*real((L*(L+1)+l2*(l2+1)-l1*(l1+1)), dp)*P_tilde%Cl(l1,4))) / &
                 ( (P%Cl(l1,1) +Noise%Cl(l1,1))*(P%Cl(l2,3)+Noise%Cl(l2,3)))


!                write(1, *) L, l1, l2, g(L, l1, l2)
            end do
        end do
    end do
    else if (est_type == 5) !EB
    do L=1, lmax_phi
        do l1=2, lmax_est
            call GetThreeJs(threejs3(abs(l1-L)),l1,L,-2,0)
            do l2 = max(2,abs(l1-L)), min(l1+L,lmax_est)
                if (threejs3(l2)==0 .or. (1-mod(L+l1+l2, 2))) cycle
                g(L, l1, l2) = cmplx(0, -1)*conjg(threejs3(l2))*sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(4*HO_fourpi))* &
                 (-( real((L*(L+1) + l2*(l2+1)-l1*(l1+1)), dp)*P_tilde%Cl(l2,3)) - real((L*(L+1) - l2*(l2+1)+l1*(l1+1)), dp)*P_tilde%Cl(l1,2) ) / &
                 ( (P%Cl(l1,2) +Noise%Cl(l1,2))-*(P%Cl(l2,3)+Noise%Cl(l2,3)))


!                write(1, *) L, l1, l2, g(L, l1, l2)
            end do
        end do
    end do

    end subroutine Getg

    subroutine GetA_L(P, Noise, AL, g, lmax_phi, lmax_est, est_type)
    Type(HealpixPower) :: P, Noise
    real(dp) AL(:), g1, g2, g3
    integer lmax_est, lmax_phi
    integer l1,l2,L, est_type
    real(dp), allocatable :: threejs4(:), threejs6(:)
    real(dp), dimension(:, :, :):: g
    open (2, file = 'ALL.txt')
    allocate(threejs4(0:lmax_est + lmax_phi+1))
    allocate(threejs6(0:lmax_est + lmax_phi+1))
    AL = 0
    if (est_type == 1)
    do L=1, lmax_phi
        do l1=2, lmax_est
            call GetThreeJs(threejs4(abs(l1-L)),l1,L,0,0)
            do l2 = max(2,abs(l1-L)), min(l1+L,lmax_est)
                if (mod(L+l1+l2, 2) .or. threejs4(l2)==0) cycle

                AL(L) = AL(L)+ g(L, l1, l2)* (threejs4(l2)*sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(4*HO_fourpi))* &
                (real((L*(L+1) + l2*(l2+1)-l1*(l1+1)), dp)*P%Cl(l2,1) + real((L*(L+1) - l2*(l2+1)+l1*(l1+1)), dp)*P%Cl(l1,1)))
            end do
        end do
        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L)
!        write (2,*) L, AL(L)
    end do
    else if (est_type == 2) !EE
    do L=1, lmax_phi
        do l1=2, lmax_est
            call GetThreeJs(threejs4(abs(l1-L)),l1,L,-2,0)
            do l2 = max(2,abs(l1-L)), min(l1+L,lmax_est)
                if (mod(L+l1+l2, 2) .or threejs4(l2)==0) cycle
                AL(L) = AL(L)+ g(L, l1, l2)* (sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(4*HO_fourpi)) * &
                threejs4(l2)*(real((L*(L+1) + l2*(l2+1)-l1*(l1+1)), dp)*P_tilde%Cl(l2,2) + real((L*(L+1) - l2*(l2+1)+l1*(l1+1)), dp)*P%Cl(l1,2)))
            end do
        end do
        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L)
!        write (2,*) L, AL(L)
    end do
    else if (est_type == 3) !TE
    do L=1, lmax_phi
        do l1=2, lmax_est
            call GetThreeJs(threejs4(abs(l1-L)),l1,L,0,0)
            call GetThreeJs(threejs4(abs(l1-L)),l1,L,-2,0)
            do l2 = max(2,abs(l1-L)), min(l1+L,lmax_est)
                if (mod(L+l1+l2, 2) .or. threejs4(l2)==0) cycle

                AL(L) = AL(L)+ g(L, l1, l2)* (*sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(4*HO_fourpi))* &
                (threejs4(l2)*real((L*(L+1)+l2*(l2+1)-l1*(l1+1)), dp)*P_tilde%Cl(l2,4) + threejs6(l2)*real((L*(L+1) - l2*(l2+1)+l1*(l1+1)), dp)*P_tilde%Cl(l1,4)))
            end do
        end do
        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L)
!        write (2,*) L, AL(L)
    end do
    else if (est_type == 4) !TB
    do L=1, lmax_phi
        do l1=2, lmax_est
            call GetThreeJs(threejs4(abs(l1-L)),l1,L,-2,0)
            do l2 = max(2,abs(l1-L)), min(l1+L,lmax_est)
                if ((mod(L+l1+l2, 2)==0) .or. (threejs4(l2)==0)) cycle

                AL(L) = AL(L)- g(L, l1, l2)* (cmplx(0,1) * threejs4(l2)*sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(4*HO_fourpi))* &
                (real((L*(L+1) - l2*(l2+1)+l1*(l1+1)), dp)*P_tilde%Cl(l1,4)))
            end do
        end do
        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L)
!        write (2,*) L, AL(L)
    end do
    else if (est_type == 5)
    do L=1, lmax_phi
        do l1=2, lmax_est
            call GetThreeJs(threejs4(abs(l1-L)),l1,L,-2,0)
            do l2 = max(2,abs(l1-L)), min(l1+L,lmax_est)
                if (mod(L+l1+l2, 2)==0 .or. threejs4(l2)==0) cycle

                AL(L) = AL(L)- g(L, l1, l2)* (cmplx(0,1)*threejs4(l2)*sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(4*HO_fourpi))* &
                (real((L*(L+1) + l2*(l2+1)-l1*(l1+1)), dp)*P_tilde%Cl(l2,3) + real((L*(L+1) - l2*(l2+1)+l1*(l1+1)), dp)*P_tilde%Cl(l1,2)))
            end do
        end do
        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L)
!        write (2,*) L, AL(L)
    end do

    end subroutine GetA_L

    


    subroutine PhiReconQuadra(P, A, Noise, P_tilde, phi_power, lmax_phi, lmax_est, AL, g, L_up, L_down, est_type)
    Type(HealpixPower) :: P, Noise, P_tilde, phi_power
    !real(dp) AL(:)
    integer lmax_est, lmax_phi, L_up, L_down, est_type
    integer l1,l2,L, m, m1, m2
    real(dp), allocatable :: AL(:), threejs(:), threejs2(:)
    real(dp), dimension(:,:,:), allocatable :: g
    Type(HealpixAlm) :: A, d_alm
    complex(dp)  C
    allocate(threejs(0:lmax_est + lmax_phi+1))
    !allocate(F(0:lmax_est + lmax_phi+1), 0:lmax_est + lmax_phi+1), 0:lmax_est + lmax_phi+1))
    allocate(threejs2(0:lmax_est + lmax_phi+1))
    !allocate (AL(1:lmax_phi))
    !allocate(g(1:lmax_phi, 1:lmax_phi, 1:lmax_phi))
   ! call HealpixAlm_nullify(d_alm)
    call HealpixPower_nullify(phi_power)


    call HealpixAlm_Init(d_alm, lmax_est , npol = 0, HasPhi = .true.)
    !call HealpixPower_Init(phi_power, lmax_est )


    !g = 0
    !d_alm%spin = 1
  !  d_alm%HasPhi = .true.
  !  d_alm%TEB(1, 1, 1) = (0.0, 0.0)    
 !        d_alm%Phi(0, 1, 1) = (0.0, 0.0)
     !    d_alm%Phi(0, 1, -1) = (0.0, 0.0)
     !    d_alm%Phi(0, 0, 0) = (0.0, 0.0)


    open (unit = 11, file = 'strange.txt')
    open (unit = 10, file = 'record.txt') 
!   do L = 1, lmax_phi
   do L = L_up, L_down
      do m = 1, L
!         m = 10
         d_alm%Phi(1, L, m) = (0.0, 0.0)
         do l2=2, lmax_est
           ! call GetThreeJs(threejs2(abs(l2-L)),l2,L,0,0)
            if ( P%Cl(l2,1) == 0 ) cycle
            do m2 = 1, l2
               call GetThreeJs(threejs(abs(l2-L)),l2,L,m2,-m)

               do l1 = max(2,abs(l2-L)), min(l2+L,lmax_est)
                    m1 = m - m2
                    if (g(L, l1, l2) == 0 .or.threejs(l1)==0 .or. P%Cl(l1,1) == 0 .or. ISNAN(real( A%TEB(1, l1, abs(m1)))) .or. isnan(real( A%TEB(1, l2, m2)))) cycle   
                     
                    !write (10,*) C, L, l1, m, m1
                    !g = ((L*(L+1)+l1*(l1+1)-l2*(l2+1))*HO_twopi* P_tilde%Cl(l1,1)/( l1*(l1+1)) + (-1)**(L+l1+l2)*(L*(L+1)+l2*(l2+1)-l1*(l1+1))*HO_twopi* P_tilde%Cl(l2,1)  / ( l2*(l2+1))  ) / (2*HO_twopi**2 * P%Cl(l1,1) * P%Cl(l2,1) / ( l1*(l1+1) * l2*(l2+1) ))
!                    write (10,*) l1, l2, m1, m2, C, g(L, l1, l2), AL(L), A%TEB(1, l1, abs(m1)), A%TEB(1, l2, m2)   ! d_alm%Phi(1, L, m), AL(L)
                    C= AL(L) *  ((-1)**m) *threejs(l1) * ( A%TEB(1, l1, abs(m1)) * A%TEB(1, l2, m2)) * g(L, l1, l2)/(real(L*(L+1), dp))
                   !C = (1E-08, 1E-08)
                   !  if (abs(C) .GE. 0.0001) write (11,*) l1, l2, m1, m2, C, g(L, l1, l2), AL(L), A%TEB(1, l1, abs(m1)), A%TEB(1, l2, m2)
                    d_alm%Phi(1, L, m) = d_alm%Phi(1, L, m) + C! AL(L) *  ((-1)**m) *threejs(l1) * ( A%TEB(1, l1, abs(m1)) * A%TEB(1, l2, m2)) * g(L, l1, l2)

!                    if ((abs(C) .GT. 0.1) .or. (abs(g(L, l1, l2)) .GT. 100) .or. (abs(AL(L)) .GT. 0.001) .or. (abs(A%TEB(1, l1, m1)) .GT. 100))  write (11,*) l1, l2, m1, m2, C, g(L, l1, l2), AL(L), A%TEB(1, l1, m1), A%TEB(1, l2, m2), HO_twopi * AL(L) *  ((-1)**m) *C * ( A%TEB(1, l1, m1) * A%TEB(1, l2, m2)) * g(L, l1, l2)  ! d_alm%Phi(1, L, m), AL(L)

               end do
               
            end do
         end do
      end do
   end do



    
    call HealpixAlm2Power(d_alm, phi_power)
    
  !  do L=1, lmax_phi
!	phi_power%PhiCl(L, 1) = phi_power%PhiCl(L, 1)/(L*(L+1)  !transfer d to phi now
!    end do


    call HealpixAlm_Free(d_alm)
    



    end subroutine PhiReconQuadra



    subroutine QuadraticPhi(H,A, MGradT, LensedCl, Noise, lminest,lmaxest)
    Type(HealpixInfo)  :: H
    Type(HealpixPower) :: LensedCl, Noise
    Type(HealpixAlm) :: A, AOut , PhiSol
    Type(HealpixMap) :: MFilt, MGradT
    integer l,lmaxest, lminest

    print *,'get quadratic', lminest, lmaxest
    call HealpixAlm_Init(AOut, lmaxest,1)
    AOut%TEB=0
    do l=max(2,lminest),lmaxest
        AOut%TEB(1,l,:) = A%TEB(1,l,:) / (LensedCl%Cl(l,1) + Noise%Cl(l,1))
    end do
    call HealpixAlm2Map(H, AOut,MFilt,H%npix)

    do l=max(2,lminest),lmaxest
        AOut%TEB(1,l,:) = AOut%TEB(1,l,:)  * LensedCl%Cl(l,1)! initially Aout for second one 
    end do
    call HealpixAlm2GradientMap(H, AOut,MGradT,H%npix,'T')

    MGradT%SpinField =  MGradT%SpinField * MFilt%TQU(:,1)
    call HealpixMap_Free(Mfilt)
    call HealpixAlm_Free(AOut)

    end subroutine QuadraticPhi





    end module Recon

    program SimLensReconTT
    use HealpixObj
    use HealpixVis
    use Random
    use spinalm_tools
    use IniFile
    use AMLUtils
    use Recon
    implicit none
    Type(HealpixInfo)  :: H
    Type(HealpixMap)   :: M, GradPhi, MGradT
    Type(HealpixPower) :: UnlensedCl, LensedCl, Noise, P
    Type(HealpixAlm)   :: A,SimAlm,PhiRecon

    integer            :: nside, lmax
    integer(I_NPIX)    :: npix, est_type
    character(LEN=1024)  :: w8name = '../Healpix_2.00/data/'
    character(LEN=1024)  :: file_stem, cls_file, out_file_root, cls_lensed_file
    character(LEN=1024) :: healpixloc, aname, in_map
    integer, parameter :: lens_interp =1, lens_exact = 2
    integer :: lens_method = lens_interp
    integer :: mpi_division_method = division_equalrows
    integer ::   rand_seed
    logical :: err, want_pol
    real :: interp_factor
    integer status, i, L
    integer lmax_phi, lmax_est, L_up, L_down
    real(dp), allocatable :: AL(:)
    real(dp) Noisevar, noise_fwhm_deg
    real(dp), dimension(:,:,:), allocatable :: g



#ifdef MPIPIX

    call mpi_init(i)
#endif

    Ini_Fail_On_Not_Found = .true.
    call Ini_Open(GetParam(1), 3,err)
    if (err) then
#ifdef MPIPIX
        call mpi_finalize(i)
#endif
        stop 'No ini'
    end if
    nside  = Ini_Read_Int('nside')
    npix = nside2npix(nside)

    lmax   = Ini_Read_Int('lmax')
    cls_file = Ini_Read_String('cls_file')
    cls_lensed_file = Ini_Read_String('cls_lensed_file')
    out_file_root = Ini_Read_String('out_file_root')

    want_pol = Ini_Read_Logical('want_pol')
    rand_seed = Ini_Read_Int('rand_seed')

    noiseVar = Ini_read_real('noise')/mK**2
    noise_fwhm_deg = Ini_read_real('noise_fwhm')/60
    call NoiseInit(Noise, NoiseVar, noise_fwhm_deg, lmax)

    lmax_phi = Ini_read_int('lmax_phi')
    lmax_est = Ini_read_int('lmax_est',lmax_phi)

    est_type = Ini_read_int('est_type')

    L_up = Ini_read_int('L_up')
    L_down = Ini_read_int('L_down')

    Ini_Fail_On_Not_Found = .false.

    in_map = Ini_read_String('input_map')
    if (in_map=='') in_map = 'lensed_sim_cache.fits'  


    w8name = Ini_Read_String('w8dir')
    interp_factor=0
    if (lens_method == lens_interp) interp_factor = Ini_Read_Real('interp_factor',3.)
#ifdef MPIPIX
    mpi_division_method = Ini_Read_Int('mpi_division_method',division_balanced);
#endif 

    call Ini_Close

    file_stem =  trim(out_file_root)//'_lmax'//trim(IntToStr(lmax))//'_nside'//trim(IntTOStr(nside))// &
    '_interp'//trim(RealToStr(interp_factor,3))

    if (want_pol) file_stem=trim(file_stem)//'pol_'
    file_stem = trim(file_stem)//trim(IntToStr(lens_method)) 

    call SetIdlePriority()

    if (w8name=='') then
        call get_environment_variable('HEALPIX', healpixloc, status=status)
        if (status==0) then
            w8name = trim(healpixloc)//'/data/'
        end if
    end if

    if (w8name=='') then
        write (*,*) 'Warning: using unit weights as no w8dir found'
        call HealpixInit(H,nside, max(lmax, lmax_phi*2),.true., w8dir='', method= mpi_division_method) 
    else
        call HealpixInit(H,nside, max(lmax, lmax_phi*2),.true., w8dir=w8name,method=mpi_division_method) 
    end if 

    if (H%MpiID ==0) then !if we are main thread
        !All but main thread stay in HealpixInit

        call HealpixAlm_nullify(A)
        call HealpixMap_nullify(GradPhi)
        call HealpixMap_nullify(M)

        call HealpixPower_ReadFromTextFile(UnlensedCl,cls_file,lmax,pol=.true.,dolens = .true.)
        !Reads in unlensed C_l text files as produced by CAMB (or CMBFAST if you aren't doing lensing)

        call HealpixPower_ReadFromTextFile(LensedCl,cls_lensed_file,lmax,pol=.true.)

        print *,'at L=1500, Lens, Unlens, noise C_l =', LensedCl%Cl(1500,1), UnlensedCl%Cl(1500,1), Noise%Cl(1500,1)

        !allocate(AL(lmax_phi))
        allocate (AL(1:lmax_phi))
        allocate(g(1:lmax_phi, 1:lmax_phi, 1:lmax_phi))

        aname=trim(out_file_root)//'_AL.txt'
        call Getg(LensedCl, Noise, g, lmax_phi, lmax, est_type)
        if (.not. FileExists(aname)) then
            call GetA_L(LensedCl, Noise, AL, g, lmax_phi, lmax, est_type)
            call CreateTxtFile(aname,1)
            do i=1, lmax_phi
                write(1,*) i, AL(i)
            end do
            close(1)
        else
            call OpenTxtFile(aname,1)
            do i=1, lmax_phi
                read(1,*) L, AL(i)
            end do
            close(1)
        end if




        if (FileExists(in_map)) then   !not sure here
            print *, 'reading input map', trim(in_map)
            call HealpixMap_Read(M,in_map)
            call HealpixAlm_Nullify(SimAlm)
        else
            call HealpixAlm_Sim(SimAlm, UnlensedCl, rand_seed,HasPhi=.true., dopol = want_pol)!
            call HealpixAlm2Power(SimAlm,P)
            call HealpixPower_Write(P,trim(file_stem)//'_unlensed_simulated.dat')

            call HealpixAlm2GradientMap(H,SimAlm, GradPhi,H%npix,'PHI')
            call HealpixInterpLensedMap_GradPhi(H,SimAlm,GradPhi, M, interp_factor, interp_cyl)

            call HealpixMap_Write(M, in_map)
        end if



        !Have simulated lensed map M
        !Test very simple TT reconstruction
        !Note no noise added to map, N0=A_L will not be correct noise bias

        call HealpixMap2Alm(H,M, A,lmax)
        call PhiReconQuadra(LensedCl, A, Noise, UnlensedCl, P, lmax_phi, lmax, AL, g, L_up, L_down, est_type)
       ! call QuadraticPhi(H,A, MGradT, LensedCl, Noise, 2,lmax)
       ! call HealpixMap2Alm(H, MGradT, A,lmax_est)
       ! call HealpixAlm_Init(PhiRecon,lmax_phi,npol=0,HasPhi=.true.)
        !do i=1,lmax_phi
       !     PhiRecon%Phi(1,i,:) =  A%SpinEB(1,i,:) * AL(i) / sqrt(i*(i+1.))
       ! end do
       ! call HealpixAlm2Power(PhiRecon, P)
        call HealpixPower_write(P, trim(file_stem)//'recon_power.dat')
       ! if (associated(SimAlm%Phi)) then
       !     call HealpixAlm2CrossPhi(PhiRecon, SimAlm, P)
       !     call HealpixPower_write(P, trim(file_stem)//'recon_cross_power.dat')
       ! end if
    end if

#ifdef MPIPIX
    call HealpixFree(H)
    call mpi_finalize(i)
#endif

#ifdef DEBUG
    write (*,*) 'End of program'
    pause
#endif
    end program SimLensReconTT

