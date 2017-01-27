!routines to add systematics into simulation

      module Contamination
      use MPIstuff
      use healpix_types, ONLY: SP,DP,I4B,I8B,SPC
      use healpix_fft, only: real_fft
      USE head_fits, ONLY : add_card, get_card
      USE pix_tools, ONLY :  npix2nside, nside2npix, query_disc
      !USE fitstools, ONLY : getsize_fits, input_map, read_par, read_dbint
      USE spinalm_tools
      USE AMLutils
      USE HealpixObj
      use HealpixVis
      use Random
      

      implicit none


      type noisespinMap
        integer(I_NPIX) npix
        integer(I4B) nmaps, ordering, nside, type
        complex(SPC), dimension(:), pointer :: Map => NULL()
        integer :: spin
      end type noisespinMap



 
      type noisescalMap
        integer(I_NPIX) npix
        integer(I4B) nmaps, ordering, nside, type
        real(SPC), dimension(:), pointer :: Map => NULL()
        integer :: spin
      end type noisescalMap




      type noiseAlm
             integer(I4B) lmax
             COMPLEX(KIND=SPC), DIMENSION(:,:,:), pointer :: Alm  => NULL()
             integer spin
      end type noiseAlm





      type noisePower
             integer(I4B) lmax
             COMPLEX(KIND=SPC), DIMENSION(:), pointer :: Power  => NULL()
             integer spin
      end type noisePower


      contains

      subroutine noiseAlmPack(almin, almout, nlmax, spin)

      integer, intent(in) :: spin, nlmax
      COMPLEX(SPC), INTENT(in), dimension(spin:nlmax,0:nlmax) :: almin
      COMPLEX(SPC), INTENT(out), dimension((int(nlmax+1,I_NPIX)*(nlmax+2))/2) :: almout
      integer m, l
      integer(I_NPIX) a_ix

      a_ix = 0
      do m = 0, nlmax
        do l=m, max(m, spin)
            a_ix = a_ix+1
            almout(a_ix) = 0
        end do

        do l=max(m, spin), nlmax
            a_ix = a_ix + 1
            almout(a_ix) = almin(l,m)
        end do
      end do

      end subroutine noiseAlmPack


      subroutine noisescalMap_init(M, npix, spin, nside)

      type(noisescalMap) :: M
      integer(I_NPIX), intent(in):: npix
      integer, intent(in) :: spin
      integer, intent(in), optional :: nside
      integer status


!      if (present(npix)) then
          M%npix = npix
!          if (present(nside)) then
!              if (npix /= nside2npix(nside)) then 
!               call MpiStop('HealpixMap_Init: nside and npix specified')
!              end if
!              M%nside =nside
!          else
              M%nside = npix2nside(npix)
!          end if
!      else
!          if (present(nside)) then
!              M%npix = nside2npix(nside)
!              M%nside = nside
!          else
!              call MpiStop('HealpixMap_Init: must specifc nside or npix')
!          end if
!      end if
      M%spin = spin
      call noisescalMap_free(M)
 
      deallocate(M%Map, stat = status)
      ALLOCATE(M%Map(0:M%npix-1), stat = status)
      if (status /= 0) call MpiStop('NoiseMap_ALLOCATE : allocate')
      M%Map = 0

      end subroutine noisescalMap_init
   


      subroutine noisespinMap_init(M, npix, spin, nside)

      type(noisespinMap) :: M
      integer(I_NPIX), intent(in):: npix
      integer, intent(in) :: spin
      integer, intent(in), optional :: nside
      integer status


!      if (present(npix)) then
          M%npix = npix
!          if (present(nside)) then
!              if (npix /= nside2npix(nside)) then 
!               call MpiStop('HealpixMap_Init: nside and npix specified')
!              end if
!              M%nside =nside
!          else
              M%nside = npix2nside(npix)
!          end if
!      else
!          if (present(nside)) then
!              M%npix = nside2npix(nside)
!              M%nside = nside
!          else
!              call MpiStop('HealpixMap_Init: must specifc nside or npix')
!          end if
!      end if
      M%spin = spin
      call noisespinMap_free(M)
 
      deallocate(M%Map, stat = status)
      ALLOCATE(M%Map(0:M%npix-1), stat = status)
      if (status /= 0) call MpiStop('NoiseMap_ALLOCATE : allocate')
      M%Map = 0

      end subroutine noisespinMap_init
  
  

      subroutine noisespinMap_free(M)
      Type(noisespinMap) :: M
      integer status

      deallocate(M%Map, stat = status)
      nullify(M%Map)


      end subroutine noisespinMap_free


      subroutine noisescalMap_free(M)
      Type(noisescalMap) :: M
      integer status

      deallocate(M%Map, stat = status)
      nullify(M%Map)


      end subroutine noisescalMap_free
 
     
 
     
      subroutine noiseAlm2spinMap(H, map, A, npix, inlmax, spin)
      use MPIstuff
      type(HealpixInfo) :: H
      Type(noisespinMap) :: map
      Type(noiseAlm), intent(in) :: A
      integer, intent(in) :: spin
      integer(I_NPIX) :: npix
     
      INTEGER(I4B) :: nsmax
      INTEGER(I4B), INTENT(IN) :: inlmax
  
  
      call noisespinMap_init(map, npix, spin) 
      
!======================================================

!calculation part. modified from healpix routine

!======================================================


      if (spin>3 ) then 
           write(*,*)   ! need rewrite spin>= 4 part
      else 
           call spinalm2map(H, inlmax, A%Alm, map%Map, spin)
      end if  
  
  
      end subroutine noiseAlm2spinMap

      subroutine noiseAlm2scalMap(H, map, A, npix, inlmax, spin)
      use MPIstuff
      type(HealpixInfo) :: H
      Type(noisescalMap) :: map
      Type(noiseAlm), intent(in) :: A
      integer, intent(in) :: spin
      integer(I_NPIX) :: npix
     
      INTEGER(I4B) :: nsmax
      INTEGER(I4B), INTENT(IN) :: inlmax
  
  
      call noisescalMap_init(map, npix, spin) 
      
      if (spin /= 0) then
           write(*,*) 'error: please check if the field is scalar'
           stop
      else
           call scalalm2map(H, inlmax, A%Alm, map%Map)
      end if
  
      end subroutine noiseAlm2scalMap

 

      subroutine Normalization(noise_fwhm, rms,lmax, C_02, spin)
      
      real(dp) :: sigma2, rms, C_02, K, xlc, noise_fwhm     !noise_fwhm is noise coherance sclae or alpha_S in our paper. rms is A_S in our paper   FWHM = (8 ln2)^(1/2) noise_fwhm
      integer lmax, l, spin
       
      xlc= 180*60*sqrt(8.*log(2.))/3.14159
      sigma2 = (noise_fwhm/xlc)**2

                 !need a part to correct unit here
      K = 0

      do l=spin, 2*lmax
               K = K+ real(l,dp)*exp(-real(l*(l+1.),dp)*sigma2)/HO_twopi  !approx for integral, not precise
      end do
      C_02 = rms**2/K

      end subroutine Normalization




      subroutine noisealm_init(A, lmax, spin)
      Type(noiseAlm) A
      integer, intent(in) :: lmax, spin
      integer status

 
      deallocate(A%alm, stat=status)
      nullify(A%alm)
      A%lmax = lmax
      
      if(spin == 0) then
            ALLOCATE(A%alm(1:1, spin:lmax, 0:lmax), stat = status)
      else 
            ALLOCATE(A%alm(1:2, spin:lmax, 0:lmax), stat = status)       
      end if


      if (status /= 0) call MpiStop('No Mem: noisealm_init lmax = '//IntToStr(lmax))
      A%spin = spin 
      A%alm = 0    

      end subroutine noisealm_init


      subroutine noisepower_init(P, lmax, spin)
      Type(noisePower) P
      integer, intent(in) :: lmax, spin
      integer status


      deallocate(P%Power, stat=status)
      nullify(P%Power)
      P%lmax = lmax
      P%spin = spin
      ALLOCATE(P%Power(spin:lmax), stat = status)
      if (status /= 0) call MpiStop('No Mem: noisepower_init lmax = '//IntToStr(lmax))

      P%Power = 0

      end subroutine noisepower_init



      subroutine noise_scalalm_sim(A, P, spin, seed)

      use random
      use alm_tools
      use ran_tools
      Type(noiseAlm) :: A          
      Type(noisePower) :: P      
      integer, intent(in) :: spin
      integer, intent(in), optional ::seed
  
      integer l,m

      real(sp) xamp, corr, tamp, Bamp, Examp, sqrt2


!random seed module here
      if (present(seed)) then
             call InitRandom(seed)
      else
             if (.not. RandInited) call InitRandom
             RandInited = .true.
      end if
      

      call noisealm_init(A, P%lmax, spin)            
      

      sqrt2 = sqrt(2.)

      do l=max(spin, 1), P%lmax
             A%alm(1,l,0) = Gaussian1()* sqrt(P%Power(l))
             tamp = sqrt(P%Power(l)/2.)
             do m = 1,l
                      A%alm(1,l,m) = cmplx(Gaussian1(),Gaussian1())*tamp
             end do
      end do


      end subroutine noise_scalalm_sim




      subroutine noise_spinalm_sim(A, P1, P2, spin, seed)

      use random
      use alm_tools
      use ran_tools
      Type(noiseAlm) :: A          
      Type(noisePower) :: P1, P2      
      integer, intent(in) :: spin
      integer, intent(in), optional ::seed
  
      integer l,m

      real(sp) xamp, corr, tamp1, tamp2, sqrt2


!random seed module here
      if (present(seed)) then
             call InitRandom(seed)
      else
             if (.not. RandInited) call InitRandom
             RandInited = .true.
      end if
      

      call noisealm_init(A, P1%lmax, spin)            
      

      sqrt2 = sqrt(2.)

      do l=max(spin, 1), P1%lmax
             A%alm(1,l,0) = Gaussian1()* sqrt(P1%Power(l))
             A%alm(2,l,0) = Gaussian1()*sqrt(P2%Power(l))
             tamp1 = sqrt(P1%Power(l)/2.)
             tamp2 = sqrt(P2%Power(l)/2.)
             do m = 1,l
                      A%alm(1,l,m) = Gaussian1()*tamp1
                      A%alm(2,l,m) = Gaussian1()*tamp2
             end do
      end do


      end subroutine noise_spinalm_sim






      subroutine GetscalNoise(Noise_Cl, Noise_Alm, lmax, noise_fwhm, rms, spin, con_type)
     
      real(dp) :: C_02
      Type(noisePower) :: Noise_Cl
      Type(noiseAlm) :: Noise_Alm
      real(dp) :: noise_fwhm, rms
      integer, optional :: spin

      integer  l, lmax 
      character(LEN=*), intent(IN), optional :: con_type

!      if (.not. present(spin)) spin =0
!write something for spin \=0

      
      if (rms == 0 .or. noise_fwhm ==0) then
          write(*,*) 'ignore noise', con_type
      else
          call Normalization(noise_fwhm, rms,lmax, C_02, spin) 
          call noisepower_init(Noise_Cl, lmax, spin)
          do l = 1, 2*lmax
              Noise_Cl%Power(l) = C_02 * exp(- l*(l+1)* noise_fwhm**2)
          end do
          call noise_scalalm_Sim(Noise_Alm, Noise_Cl, spin)  
      end if



      end subroutine GetscalNoise

      subroutine GetNoise(Noise_Cl1, Noise_Cl2, Noise_Alm, lmax, noise_fwhm1, noise_fwhm2, rms1, rms2, spin, con_type)
 
      real(dp) :: C_01, C_02
      Type(noisePower) :: Noise_Cl1, Noise_Cl2
      Type(noiseAlm) :: Noise_Alm
      real(dp) :: noise_fwhm1, rms1, noise_fwhm2, rms2

      integer, optional :: spin

      integer  l, lmax 
      character(LEN=*), intent(IN), optional :: con_type

      if (.not. present(spin)) spin =0


      
      if ((rms1 == 0 .or. noise_fwhm1 ==0).and. (rms2 ==0 .or. noise_fwhm2 == 0)) then
          write(*,*) 'ignore noise', con_type
      else
          call Normalization(noise_fwhm1, rms1, lmax, C_01, spin) 
          call noisepower_init(Noise_Cl1, lmax, spin)
          call Normalization(noise_fwhm2, rms2, lmax, C_02, spin) 
          call noisepower_init(Noise_Cl2, lmax, spin)
          do l = max(spin, 1), 2*lmax
              Noise_Cl1%Power(l) = C_01 * exp(- l*(l+1)* noise_fwhm1**2)
              Noise_Cl2%Power(l) = C_02 * exp(- l*(l+1)* noise_fwhm2**2)
          end do
          call noise_spinalm_Sim(Noise_Alm, Noise_Cl1, Noise_Cl2, spin) 
      end if



      end subroutine GetNoise




!write a file read noise from file in the future




      end module Contamination


    
 
      Program Noise

      use HealpixObj
      use HealpixVis
      use Random
      use spinalm_tools
      use IniFile
      use AMLUtils
      use Recon
      use Contamination 
      implicit none
     
      Type(noisePower) :: Noise_Cl1, Noise_Cl2,Noise_Cl3,Noise_Cl4,Noise_Cl5,Noise_Cl6,Noise_Cl7,Noise_Cl8,Noise_Cl9,Noise_Cl10,Noise_Cl11
!      Type(noiseAlm) :: Noise_Alm1, Noise_Alm2,Noise_Alm3,Noise_Alm4,Noise_Alm5,Noise_Alm6,Noise_Alm7,Noise_Alm8,Noise_Alm9,Noise_Alm10,Noise_Alm11
       Type(noiseAlm) :: Noise_Alm1,Noise_Alm2, Noise_Alm3,Noise_Alm5,Noise_Alm7,Noise_Alm9,Noise_Alm11
      
      Type(noisespinMap) :: Noise_Map3, Noise_Map5, Noise_Map7, Noise_Map9
      Type(noisescalMap) :: Noise_Map1,Noise_Map2,Noise_Map11
     
      Type(HealpixInfo) :: H 

      Type(HealpixMap) :: Map, Map_con
      Type(HealpixAlm) :: Alm, Alm_con
      Type(HealpixPower) :: Cl, Cl_con
      integer :: nside, m, m1, m2, noise_input_type, L_up, L_down, l, l1, l2, lmax_est, lmax_phi, lmax, rand_seed
      integer(I_NPIX) :: npix
      character(LEN=32) :: noise_file, k, cls_lensed_file, Interp_Method
      character(LEN=1024) :: w8name = '../Healpix_3.31/data/'    ! adjust if lenspix and Healpix are not in same directory
      integer i
      complex(dp) :: C, D
      real(dp) beam_width
      real(dp), dimension(11) :: noise_fwhm, rms
      real(dp), allocatable :: threejs(:), threejs1(:), threejs2(:), threejs3(:)
      logical :: err
     !1: calibration a, 2: rotation omega, 3,4: spin flip f_a & f_b, 5,6: pointing p_a &p_b 7,8: monopole leakage gamma_a&b 9,10: dipole leakage d_a&b 11: quadrupole leakage q
      allocate(threejs(0:lmax_est +lmax_phi+ 1))
      allocate(threejs1(0:lmax_phi +lmax_est+ 1))
      allocate(threejs2(0:lmax_phi +lmax_est+ 1))
      allocate(threejs3(0:lmax_phi +lmax_est+ 1))

      call mpi_init(i)
!endif

      Ini_Fail_On_Not_Found = .true.
      call Ini_Open(GetParam(1), 3,err)
      if (err) then
!ifdef MPIPIX
        call mpi_finalize(i)
!endif
        stop 'No ini'
      end if


      Interp_Method = Ini_Read_string('Interp_Method')   !SH or Map only

      lmax_phi =Ini_Read_Int('lmax_phi')
      lmax_est =Ini_Read_Int('lmax_est', lmax_phi)

      noise_file = Ini_Read_string('noise_file')

      cls_lensed_file = Ini_Read_String('cls_lensed_file')
      
      rand_seed = Ini_Read_Int('rand_seed')

      nside = Ini_Read_Int('nside')


      L_up = Ini_read_int('L_up')
      L_down = Ini_read_int('L_down')

      lmax   = Ini_Read_Int('lmax')

    
      call Ini_Close



      call SetIdlePriority()



    

      npix = nside2npix(nside)

      !noise initiallization block
      call HealpixInit(H, nside, lmax, .true., w8dir=w8name, method=3) ! remember to check w8, and change mpi division method if you want


!
      call HealpixAlm_Init(alm, lmax_est , npol = 3, HasPhi = .false.)
      call HealpixAlm_Init(alm_con, lmax_est , npol = 3, HasPhi = .false.)
      call HealpixPower_Init(Cl, lmax_est , pol = .true. )
      call HealpixPower_Init(Cl_con, lmax_est , pol = .true. )

      call HealpixMap_Init(Map, npix, pol = .true.)
      call HealpixMap_Init(Map_con, npix, pol = .true.)


      call HealpixPower_ReadFromTextFile(Cl,cls_lensed_file,lmax,pol=.true.,dolens = .false.)  ! should not be real lensed file, since we need to do Map simulation and will be not accurate to compare with the original map. but phi are not intended to be included

      call HealpixAlm_Sim(alm, Cl, rand_seed, HasPhi = .false.)

      call polalm2map(H, lmax, alm%TEB, Map%TQU)


!=======================part for generating gradient map, need write new subroutines===============
!      call HealpixAlm2GradientMap()


!================================================================



!      do i =1, 11
!           noise_fwhm(i) = Ini_read_Double('noise_fwhm'//trim(IntToStr(i)))
!           rms(i) = Ini_read_Double('rms'//trim(IntToStr(i)))
!      end do
!this part does not work due to unknown error,shit

           noise_fwhm(1) = Ini_read_Double('noise_fwhm1')
           rms(1) = Ini_read_Double('rms1')
           noise_fwhm(2) = Ini_read_Double('noise_fwhm2')
           rms(2) = Ini_read_Double('rms2')
           noise_fwhm(3) = Ini_read_Double('noise_fwhm3')
           rms(3) = Ini_read_Double('rms3')
           noise_fwhm(4) = Ini_read_Double('noise_fwhm4')
           rms(4) = Ini_read_Double('rms4')
           noise_fwhm(5) = Ini_read_Double('noise_fwhm5')
           rms(5) = Ini_read_Double('rms5')
           noise_fwhm(6) = Ini_read_Double('noise_fwhm6')
           rms(6) = Ini_read_Double('rms6')
           noise_fwhm(7) = Ini_read_Double('noise_fwhm7')
           rms(7) = Ini_read_Double('rms7')
           noise_fwhm(8) = Ini_read_Double('noise_fwhm8')
           rms(8) = Ini_read_Double('rms8')
           noise_fwhm(9) = Ini_read_Double('noise_fwhm9')
           rms(9) = Ini_read_Double('rms9')
           noise_fwhm(10) = Ini_read_Double('noise_fwhm10')
           rms(10) = Ini_read_Double('rms10')
           noise_fwhm(11) = Ini_read_Double('noise_fwhm11')
           rms(11) = Ini_read_Double('rms11')






      if (.not. FileExists(noise_file)) then
          if (noise_fwhm(1) /= 0) then
            call Getscalnoise(Noise_Cl1, Noise_Alm1, lmax, noise_fwhm(1), rms(1), spin = 0)
            call noiseAlm2scalMap(H, Noise_Map1, Noise_Alm1, npix, lmax, 0)
          end if


          if (noise_fwhm(2) /= 0) then
            call Getscalnoise(Noise_Cl2, Noise_Alm2, lmax, noise_fwhm(2), rms(2), spin = 0)
            call noiseAlm2scalMap(H, Noise_Map2, Noise_Alm2, npix, lmax, 0)
          end if


          if (noise_fwhm(3) /= 0 .or. noise_fwhm(4) /= 0) then
            call GetNoise(Noise_Cl3, Noise_Cl4, Noise_Alm3, lmax, noise_fwhm(3), noise_fwhm(4), rms(3), rms(4), spin = 4)
            call noiseAlm2spinMap(H, Noise_Map3, Noise_Alm3, npix, lmax, 4)
          end if


          if (noise_fwhm(5) /= 0 .or. noise_fwhm(6) /= 0) then
            call GetNoise(Noise_Cl5, Noise_Cl6, Noise_Alm5, lmax, noise_fwhm(5),noise_fwhm(6), rms(5), rms(6), spin = 1)
            call noiseAlm2spinMap(H, Noise_Map5, Noise_Alm5, npix, lmax, 1)
          end if


          if (noise_fwhm(7) /= 0 .or. noise_fwhm(8) /= 0) then
            call GetNoise(Noise_Cl7, Noise_Cl8, Noise_Alm7, lmax, noise_fwhm(7),noise_fwhm(8), rms(7), rms(8), spin = 2)
            call noiseAlm2spinMap(H, Noise_Map7, Noise_Alm7, npix, lmax, 2)
          end if



          if (noise_fwhm(9) /= 0 .or. noise_fwhm(10) /= 0) then
            call GetNoise(Noise_Cl9, Noise_Cl10, Noise_Alm9, lmax, noise_fwhm(9),noise_fwhm(10), rms(9), rms(10), spin = 1)
            call noiseAlm2spinMap(H, Noise_Map9, Noise_Alm9, npix, lmax, 1)
          end if


          if (.not. noise_fwhm(11) == 0) then
            call GetscalNoise(Noise_Cl11, Noise_Alm11, lmax, noise_fwhm(11), rms(11), spin = 0)
            call noiseAlm2scalMap(H, Noise_Map11, Noise_Alm11, npix, lmax, 0)
          end if
          call CreateTxtFile(noise_file,1)
!            do i=1, lmax_phi
!                write(1,*) i, Noise_Cl1%Cl(1,i), Noise_Cl2%Cl(1,i), Noise_Cl3%Cl(1,i), Noise_Cl4%Cl(1,i), Noise_Cl5%Cl(1,i), Noise_Cl6%Cl(1,i), Noise_Cl7%Cl(1,i), Noise_Cl8%Cl(1,i), Noise_Cl9%Cl(1,i), Noise_Cl10%Cl(1,i), Noise_Cl11%Cl(1,i)
!            end do
            close(1)
      else
            call OpenTxtFile(noise_file,1)
!            do i=1, lmax_phi
!               read(1,*) L, Noise_Cl(1)%Cl(i), Noise_Cl(2)%Cl(i), Noise_Cl(3)%Cl(i), Noise_Cl(4)%Cl(i), Noise_Cl(5)%Cl(i), Noise_Cl(6)%Cl(i), Noise_Cl(7)%Cl(i), Noise_Cl(8)%Cl(i), Noise_Cl(9)%Cl(i), Noise_Cl(10)%Cl(i), Noise_Cl(11)%Cl(i)
!                call HealpixAlm_Sim(Noise_Alm(i), Noise_Cl(i), HasPhi=.false.)  !may occupy too much space, can consider rewrite the simulation for noise , also need to check if sim apply for different spin field
!           end do
            close(1)
      end if     
        !add input map in the future
!      do i = 1, 11
!             call HealpixAlm_Sim(Noise_Alm(i), Noise_Cl(i), HasPhi=.false.)  !may occupy too much space, can consider rewrite the simulation for noise 
!      end do

     ! Alm_con%TEB(2,l,m) = 0
      !calculation part.
      !calibration and rotation
      

      if(Interp_Method == 'SH' .and. noise_fwhm(1)+noise_fwhm(2) /= 0) then
      do L=L_down, L_up
        do l2=2, lmax_est
           call GetThreeJs(threejs(abs(l2-L)),l2,L,2,-2)  !(L, l1, l2)(-2, 0, 2)
           do m = 1, L
              do m2 = 1, l2
                call GetThreeJs(threejs1(abs(l2-L)),l2,L,m2,-m)
                do l1 = max(2,abs(l2-L)), min(l2+L,lmax_est)
                    if (threejs(l1)==0 .or. threejs1(l1)==0) cycle
                    if (mod(L+l1+l2, 2)) then
                       
                       C = Alm%TEB(2,l2,m2)*Noise_Alm1%Alm(1,l1,m1)
                       D = Alm%TEB(2,l2,m2)*Noise_Alm2%Alm(1,l1,m1)*cmplx(0,-2)
                    else
                       C = Alm%TEB(2,l2,m2)*Noise_Alm2%Alm(1,l1,m1)*cmplx(0,-2)
                       D = Alm%TEB(2,l2,m2)*Noise_Alm1%Alm(1,l1,m1)
                    
                    Alm_con%TEB(2,L,m) = Alm_con%TEB(2,l,m)+(-1)**m * threejs(l1) &
                                     * threejs1(l1) * C*&
                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
                    Alm_con%TEB(3,L,m) = Alm_con%TEB(3,l,m)+cmplx(0,-1)*(-1)**m &
                                         * threejs(l1)* threejs1(l1)*D*&
                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
                    end if
                enddo
              enddo
           enddo
        enddo
!        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L) 
!       write (2,*) L, AL(L)
      enddo
      end if
     !spin flip 3,4 f, spin 4

!=================== map part========================


      if(Interp_Method == 'Map' .and. noise_fwhm(1) /= 0) then
         do L = 0, npix-1
            
            Map_con%TQU(L, 2) = Map%TQU(L,2) * Noise_Map1%Map(L)
            Map_con%TQU(L, 3) = Map%TQU(L,3) * Noise_Map1%Map(L)

         end do
      end if 


 
      if(Interp_Method == 'Map' .and. noise_fwhm(2) /= 0) then
         do L = 0, npix-1
            
            Map_con%TQU(L, 2) = -Map%TQU(L,3) *2* Noise_Map2%Map(L)
            Map_con%TQU(L, 3) = Map%TQU(L,2) *2* Noise_Map2%Map(L)

         end do
      end if 
 


!      if(Interp_Method == 'SH' .and. noise_fwhm(3) /= 0) then
!      do L=L_down, L_up
!        do l2=2, lmax_est
!           call GetThreeJs(threejs(abs(l2-L)),l2,L,-2,-2)  !(L, l1, l2)(-2, 0, 2)
!           do m = 1, L
!              do m2 = 1, l2
!                call GetThreeJs(threejs1(abs(l2-L)),l2,L,m2,-m)
!                do l1 = max(2,abs(l2-L)), min(l2+L,lmax_est)
!                    if (threejs(l1)==0 .or. threejs1(l1)==0) cycle
!                    if (mod(L+l1+l2, 2)) then
!                       C = Alm%TEB(2,l2,m2)*Noise_Alm3%Alm(l1,m1)
!                       D = Alm%TEB(2,l2,m2)*Noise_Alm4%Alm(l1,m1)*cmplx(0,-1)
!                    else
!                       C = Alm%TEB(2,l2,m2)*Noise_Alm4%Alm(l1,m1)*cmplx(0,-1)
!                       D = Alm%TEB(2,l2,m2)*Noise_Alm3%Alm(l1,m1)
!
!                    Alm_con%TEB(2,L,m) = Alm_con%TEB(2,l,m)+(-1)**m * threejs(l1) &
!                                     * threejs1(l1) * C*&
!                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
!                    Alm_con%TEB(3,L,m) = Alm_con%TEB(3,l,m)+cmplx(0,-1)*(-1)**m &
!                                         * threejs(l1)* threejs1(l1)*D*&
!                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
!
!                    end if
!                enddo
!                enddo
!              enddo
!           enddo
!        enddo
!!        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L) 
!!       write (2,*) L, AL(L)
!      end if
!
!     !poingting, 5.6. pa, pb   (waiting to be finished)
!      if(.not. noise_fwhm(5)==0) then
!      do L=L_down, L_up
!        do l2=2, lmax_est
!           call GetThreeJs(threejs(abs(l2-L)),l2,L,-1,2)  !(L, l1, l2)(-2, 0, 2)
!           call GetThreeJs(threejs3(abs(l2-L)),l2,L,-3,2)  
!           do m = 1, L
!              do m2 = 1, l2
!                call GetThreeJs(threejs1(abs(l2-L)),l2,L,m2,-m)
!                do l1 = max(2,abs(l2-L)), min(l2+L,lmax_est)
!                    if (threejs1(l1)==0) cycle
!                    if (mod(L+l1+l2, 2)) then
!                       C = Alm%TEB(2,l2,m2)*Noise_Alm5%Alm(l1,m1)*(sqrt(real((l2+2)*(l2-1),dp))*threejs(l1)+sqrt(real((l2-2)*(l2+3),dp))*threejs3(l1))
!                       D = Alm%TEB(2,l2,m2)*Noise_Alm6%Alm(l1,m1)*(sqrt(real((l2+2)*(l2-1),dp))*threejs(l1)-sqrt(real((l2-2)*(l2+3),dp))*threejs3(l1))*cmplx(0,1)
!                    else
!                       C = Alm%TEB(2,l2,m2)*Noise_Alm6%Alm(l1,m1)*cmplx(0,-1)*(sqrt(real((l2+2)*(l2-1),dp))*threejs(l1)-sqrt(real((l2-2)*(l2+3),dp))*threejs3(l1))
!                       D = Alm%TEB(2,l2,m2)*Noise_Alm5%Alm(l1,m1)*cmplx(0,-1)*(sqrt(real((l2+2)*(l2-1),dp))*threejs(l1)+sqrt(real((l2-2)*(l2+3),dp))*threejs3(l1))*cmplx(0,1)
!
!
!                    Alm_con%TEB(2,L,m) = Alm_con%TEB(2,l,m)-(-1)**m * threejs1(l1) &
!                                     * C*beam_width/2 *&
!                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
!                    Alm_con%TEB(3,L,m) = Alm_con%TEB(3,l,m)-cmplx(0,-1)*(-1)**m &
!                                         *  threejs1(l1)*D*beam_width/2 *&
!                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
!
!                    end if
!                enddo
!                enddo
!              enddo
!           enddo
!        enddo
!!        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L) 
!!       write (2,*) L, AL(L)
!      end if
!
!
!
!      !monopole leakage 7,8 gamma_a,b spin 2
!
!      if(.not. noise_fwhm(7)==0) then
!      do L=L_down, L_up
!        do l2=2, lmax_est
!           call GetThreeJs(threejs(abs(l2-L)),l2,L,0,-2)  !(L, l1, l2)(-2, 0, 2)
!           do m = 1, L
!              do m2 = 1, l2
!                call GetThreeJs(threejs1(abs(l2-L)),l2,L,m2,-m)
!                do l1 = max(2,abs(l2-L)), min(l2+L,lmax_est)
!                    if (threejs(l1)==0 .or. threejs1(l1)==0) cycle
!                    if (mod(L+l1+l2, 2)) then
!                       C = -Alm%TEB(1,l2,m2)*Noise_Alm7%Alm(l1,m1)
!                       D = -Alm%TEB(2,l2,m2)*Noise_Alm8%Alm(l1,m1)*cmplx(0,-1)
!                    else
!                       C = -Alm%TEB(1,l2,m2)*Noise_Alm8%Alm(l1,m1)*cmplx(0,-1)
!                       D = -Alm%TEB(2,l2,m2)*Noise_Alm7%Alm(l1,m1)
!
!                    Alm_con%TEB(2,L,m) = Alm_con%TEB(2,l,m)+(-1)**m * threejs(l1) &
!                                     * threejs1(l1) * C*&
!                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
!                    Alm_con%TEB(3,L,m) = Alm_con%TEB(3,l,m)+cmplx(0,-1)*(-1)**m &
!                                         * threejs(l1)* threejs1(l1)*D*&
!                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
!
!                    end if
!                enddo
!                enddo
!              enddo
!           enddo
!        enddo
!!        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L) 
!!       write (2,*) L, AL(L)
!      end if
!
!
!
!
!     !dipole leakage 9, 10, d_a, d_b spin 1
!
!      if(.not. noise_fwhm(9)==0) then
!
!      do L=L_down, L_up
!        do l2=2, lmax_est
!           call GetThreeJs(threejs(abs(l2-L)),l2,L,1,-2)  !(L, l1, l2)(-2, 0, 2)
!           do m = 1, L
!              do m2 = 1, l2
!                call GetThreeJs(threejs1(abs(l2-L)),l2,L,m2,-m)
!                do l1 = max(2,abs(l2-L)), min(l2+L,lmax_est)
!                    if (threejs(l1)==0 .or. threejs1(l1)==0) cycle
!                    if (mod(L+l1+l2, 2)) then
!                       D = -beam_width*sqrt(real(l2*(l2+1), dp))*Alm%TEB(2,l2,m2)*Noise_Alm9%Alm(l1,m1)*cmplx(0,1)
!                       C = -beam_width*sqrt(real(l2*(l2+1), dp))*Alm%TEB(2,l2,m2)*Noise_Alm10%Alm(l1,m1)
!                    else
!                       C = -beam_width*sqrt(real(l2*(l2+1), dp))*Alm%TEB(2,l2,m2)*Noise_Alm9%Alm(l1,m1)*cmplx(0,-1)
!                       D = -beam_width*sqrt(real(l2*(l2+1), dp))*Alm%TEB(2,l2,m2)*Noise_Alm10%Alm(l1,m1)
!
!
!                    Alm_con%TEB(2,L,m) = Alm_con%TEB(2,l,m)+(-1)**m * threejs(l1) &
!                                     * threejs1(l1) * C*&
!                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
!                    Alm_con%TEB(3,L,m) = Alm_con%TEB(3,l,m)+cmplx(0,-1)*(-1)**m &
!                                         * threejs(l1)* threejs1(l1)*D*&
!                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
!
!                    end if
!                enddo
!                enddo
!              enddo
!           enddo
!        enddo
!!        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L) 
!!       write (2,*) L, AL(L)
!      end if
!
!
!     !quadrupole leakage 11 p
!
!
!      if(.not. noise_fwhm(11)==0) then
!
!      do L=L_down, L_up
!        do l2=2, lmax_est
!           call GetThreeJs(threejs(abs(l2-L)),l2,L,2,-2)  !(L, l1, l2)(-2, 0, 2)
!           do m = 1, L
!              do m2 = 1, l2
!                call GetThreeJs(threejs1(abs(l2-L)),l2,L,m2,-m)
!                do l1 = max(2,abs(l2-L)), min(l2+L,lmax_est)
!                    if (threejs(l1)==0 .or. threejs1(l1)==0) cycle
!                    if (mod(L+l1+l2, 2)) then
!                         C = beam_width**2*sqrt(real((l2+2)*(l2+1)*l2*(l2-1))) *Alm%TEB(2,l2,m2)*Noise_Alm11%Alm(l1,m1)*cmplx(0,1)
!                    else
!                         D = beam_width**2*sqrt(real((l2+2)*(l2+1)*l2*(l2-1))) *Alm%TEB(2,l2,m2)*Noise_Alm11%Alm(l1,m1)
!
!
!                    Alm_con%TEB(2,L,m) = Alm_con%TEB(2,l,m)+(-1)**m * threejs(l1) &
!                                     * threejs1(l1) * C*&
!                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
!                    Alm_con%TEB(3,L,m) = Alm_con%TEB(3,l,m)+cmplx(0,-1)*(-1)**m &
!                                         * threejs(l1)* threejs1(l1)*D*&
!                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
!
!                    end if
!                enddo
!                enddo
!              enddo
!           enddo
!        enddo
!!        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L) 
!!       write (2,*) L, AL(L)
!      end if
!

      call HealpixAlm2Power(Alm_con,Cl_con)
      call HealpixPower_write(Cl_con, 'Cl_con.dat')
       !remember to modify here to change output file
      call HealpixPower_nullify(Cl_con)
 
      call HealpixMap_Write(Map_con, 'Map_con.fits')



      end program Noise
