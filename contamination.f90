!routines to add systematics into simulation

      module Contamination
      use healpix_types, ONLY: SP,DP,I4B,I8B,SPC
      USE head_fits, ONLY : add_card, get_card
      USE pix_tools, ONLY :  npix2nside, nside2npix, query_disc
      !USE fitstools, ONLY : getsize_fits, input_map, read_par, read_dbint
      USE spinalm_tools
      USE AMLutils
      USE HealpixObj
      use HealpixVis
      use Random
      

      implicit none

      contains
   
     !type Noise
        ! integer(I4B) lmax
        ! COMPLEX(KIND=SPC), DIMENSION(:,:,:), pointer :: TEB  => NULL()    !T,E,B index, l,m
 
     !end type
 
      subroutine GetNoise(Noise_Cl, Noise_Alm, lmax, nscale, rms)
     
      real(dp) :: C_0
      Type(HealpixPower), dimension(11) :: Noise_Cl
      Type(HealpixAlm), dimension(11) :: Noise_Alm
      real(dp), dimension(11) :: nscale, rms
      integer status l, i

      do i=1, 13
            C_0 = 
            do l = 1, lmax
               Noise_Cl(i)%Cl(l,1) = C_0 * exp(- l*(l+1)* nscale(i)**2)
            end do
            call HealpixAlm_Sim(Noise_Alm(i), Noise_Cl(i), HasPhi=.false.)  !may occupy too much space, can consider rewrite the simulation for noise 

      end do

      end subroutine GetNoise

      end module Contamination


    
 
      Program Noise()

      use HealpixObj
      use HealpixVis
      use Random
      use spinalm_tools
      use IniFile
      use AMLUtils
      use Recon
      
      implicit none
     
      Type(HealpixPower), dimension(11) :: Noise_Cl
      Type(HealpixAlm), dimension(11) :: Noise_Alm
      Type(HealpixAlm) :: Alm, Alm_con
      Type(HealpixPower) :: Cl, Cl_con
      !HealpixMap, dimension(11) :: Noise_map
      integer :: noise_input_type, L_up, L_down, l, l1, l2, lmax_est, lmax_phi, lmax, rand_seed
      character(LEN=32) :: noise_file, k
      integer status i, L
      complex(dp) :: C, D
      real(dp): beam_width
      real(dp), dimension(11) :: noise_scale, rms
      real(dp), allocatable :: threejs(:), threejs1(:), threejs2(:), threejs3(:)
     
     !1: calibration a, 2: rotation omega, 3,4: spin flip f_a & f_b, 5,6: pointing p_a &p_b 7,8: monopole leakage gamma_a&b 9,10: dipole leakage d_a&b 11: quadrupole leakage q
      allocate(threejs(0:lmax_est +lmax_phi+ 1))
      allocate(threejs1(0:lmax_phi +lmax_est+ 1))
      allocate(threejs2(0:lmax_phi +lmax_est+ 1))
      allocate(threejs3(0:lmax_phi +lmax_est+ 1))



      lmax_phi =Ini_Read_Int('lmax_phi')
      lmax_est =Ini_Read_Int('lmax_est', lmax_phi)

      noise_file = Ini_Read_string('noise_file')

      cls_lensed_file = Ini_Read_String('cls_lensed_file')
      
      rand_seed = Ini_Read_Int('rand_seed')



      L_up = Ini_read_int('L_up')
      L_down = Ini_read_int('L_down')

      lmax   = Ini_Read_Int('lmax')

    
      !noise initiallization block
      

      do i =1, 11
           write(k,*),i
           noise_scale(i) = Ini_read_Double('noise_scale'//trim(k))
           rms(i) = Ini_read_Double('rms'//trim(k))
      end do





      call Ini_Close


      call SetIdlePriority()




!mpi initiation










    
      if (.not. FileExists(noise_file)) then
            call GetNoise(Noise_Cl, Noise_Alm, lmax, noise_scale, rms)
            call CreateTxtFile(noise_file,1)
            do i=1, lmax_phi
                write(1,*) i, Noise_Cl(1)%Cl(i), Noise_Cl(2)%Cl(i), Noise_Cl(3)%Cl(i), Noise_Cl(4)%Cl(i), Noise_Cl(5)%Cl(i), Noise_Cl(6)%Cl(i), Noise_Cl(7)%Cl(i), Noise_Cl(8)%Cl(i), Noise_Cl(9)%Cl(i), Noise_Cl(10)%Cl(i), Noise_Cl(11)%Cl(i)
            end do
            close(1)
      else
            call OpenTxtFile(noise_file,1)
            do i=1, lmax_phi
                read(1,*) L, Noise_Cl(1)%Cl(i), Noise_Cl(2)%Cl(i), Noise_Cl(3)%Cl(i), Noise_Cl(4)%Cl(i), Noise_Cl(5)%Cl(i), Noise_Cl(6)%Cl(i), Noise_Cl(7)%Cl(i), Noise_Cl(8)%Cl(i), Noise_Cl(9)%Cl(i), Noise_Cl(10)%Cl(i), Noise_Cl(11)%Cl(i)
                call HealpixAlm_Sim(Noise_Alm(i), Noise_Cl(i), HasPhi=.false.)  !may occupy too much space, can consider rewrite the simulation for noise , also need to check if sim apply for different spin field
            end do
            close(1)
      end if     
        

      

        !add input map in the future
!      do i = 1, 11
!             call HealpixAlm_Sim(Noise_Alm(i), Noise_Cl(i), HasPhi=.false.)  !may occupy too much space, can consider rewrite the simulation for noise 
!      end do

      Alm_con%TEB(2,l,m) = 0
      !calculation part.
      !calibration and rotation


      do L=1, lmax_phi
        do l2=2, lmax_est
           call GetThreeJs(threejs(abs(l2-L)),l2,L,2,-2)  !(L, l1, l2)(-2, 0, 2)
           do m = 1, L
              do m2 = 1, l2
                call GetThreeJs(threejs1(abs(l2-L)),l2,L,m2,-m)
                do l1 = max(2,abs(l2-L)), min(l2+L,lmax_est)
                    if (threejs(l1)==0 .or. threejs1(l1)==0) cycle
                    if (mod(L+l1+l2, 2)) then
                       C = Alm%TEB(2,l2,m2)*Noise_Alm(1)%TEB(1,l1,m1)
                       D = Alm%TEB(2,l2,m2)*Noise_Alm(2)%TEB(1,l1,m1)*cmplx(0,-2)
                    else
                       C = Alm%TEB(2,l2,m2)*Noise_Alm(2)%TEB(1,l1,m1)*cmplx(0,-2)
                       D = Alm%TEB(2,l2,m2)*Noise_Alm(1)%TEB(1,l1,m1)
                    
                    Alm_con%TEB(2,L,m) = Alm_con%TEB(2,l,m)+(-1)**m * threejs(l1) &
                                     * threejs1(l1) * C*&
                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
                    Alm_con%TEB(3,L,m) = Alm_con%TEB(3,l,m)+cmplx(0,-1)*(-1)**m &
                                         * threejs(l1)* threejs1(l1)*D*&
                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))

                end do
              end do
           end do
        end do
!        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L) 
!       write (2,*) L, AL(L)
      end do

     !spin flip 3,4 f, spin 4
      do L=1, lmax_phi
        do l2=2, lmax_est
           call GetThreeJs(threejs(abs(l2-L)),l2,L,-2,-2)  !(L, l1, l2)(-2, 0, 2)
           do m = 1, L
              do m2 = 1, l2
                call GetThreeJs(threejs1(abs(l2-L)),l2,L,m2,-m)
                do l1 = max(2,abs(l2-L)), min(l2+L,lmax_est)
                    if (threejs(l1)==0 .or. threejs1(l1)==0) cycle
                    if (mod(L+l1+l2, 2)) then
                       C = Alm%TEB(2,l2,m2)*Noise_Alm(3)%TEB(1,l1,m1)
                       D = Alm%TEB(2,l2,m2)*Noise_Alm(4)%TEB(1,l1,m1)*cmplx(0,-1)
                    else
                       C = Alm%TEB(2,l2,m2)*Noise_Alm(4)%TEB(1,l1,m1)*cmplx(0,-1)
                       D = Alm%TEB(2,l2,m2)*Noise_Alm(3)%TEB(1,l1,m1)

                    Alm_con%TEB(2,L,m) = Alm_con%TEB(2,l,m)+(-1)**m * threejs(l1) &
                                     * threejs1(l1) * C*&
                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
                    Alm_con%TEB(3,L,m) = Alm_con%TEB(3,l,m)+cmplx(0,-1)*(-1)**m &
                                         * threejs(l1)* threejs1(l1)*D*&
                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))

                end do
              end do
           end do
        end do
!        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L) 
!       write (2,*) L, AL(L)
      end do





     !poingting, 5.6. pa, pb   (waiting to be finished)

      do L=1, lmax_phi
        do l2=2, lmax_est
           call GetThreeJs(threejs(abs(l2-L)),l2,L,-1,2)  !(L, l1, l2)(-2, 0, 2)
           call GetThreeJs(threejs3(abs(l2-L)),l2,L,-3,2)  
           do m = 1, L
              do m2 = 1, l2
                call GetThreeJs(threejs1(abs(l2-L)),l2,L,m2,-m)
                do l1 = max(2,abs(l2-L)), min(l2+L,lmax_est)
                    if (threejs1(l1)==0) cycle
                    if (mod(L+l1+l2, 2)) then
                       C = Alm%TEB(2,l2,m2)*Noise_Alm(5)%TEB(1,l1,m1)*(sqrt(real((l2+2)*(l2-1),dp))*threejs(l1)+sqrt(real((l2-2)*(l2+3),dp))*threejs3(l1))
                       D = Alm%TEB(2,l2,m2)*Noise_Alm(6)%TEB(1,l1,m1)*(sqrt(real((l2+2)*(l2-1),dp))*threejs(l1)-sqrt(real((l2-2)*(l2+3),dp))*threejs3(l1))*cmplx(0,1)
                    else
                       C = Alm%TEB(2,l2,m2)*Noise_Alm(6)%TEB(1,l1,m1)*cmplx(0,-1)*(sqrt(real((l2+2)*(l2-1),dp))*threejs(l1)-sqrt(real((l2-2)*(l2+3),dp))*threejs3(l1))
                       D = Alm%TEB(2,l2,m2)*Noise_Alm(5)%TEB(1,l1,m1)*cmplx(0,-1)*(sqrt(real((l2+2)*(l2-1),dp))*threejs(l1)+sqrt(real((l2-2)*(l2+3),dp))*threejs3(l1))*cmplx(0,1)


                    Alm_con%TEB(2,L,m) = Alm_con%TEB(2,l,m)-(-1)**m * threejs1(l1) &
                                     * C*beam_width/2 *&
                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
                    Alm_con%TEB(3,L,m) = Alm_con%TEB(3,l,m)-cmplx(0,-1)*(-1)**m &
                                         *  threejs1(l1)*D*beam_width/2 *&
                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))

                end do
              end do
           end do
        end do
!        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L) 
!       write (2,*) L, AL(L)
      end do

      !monopole leakage 7,8 gamma_a,b spin 2

      do L=1, lmax_phi
        do l2=2, lmax_est
           call GetThreeJs(threejs(abs(l2-L)),l2,L,0,-2)  !(L, l1, l2)(-2, 0, 2)
           do m = 1, L
              do m2 = 1, l2
                call GetThreeJs(threejs1(abs(l2-L)),l2,L,m2,-m)
                do l1 = max(2,abs(l2-L)), min(l2+L,lmax_est)
                    if (threejs(l1)==0 .or. threejs1(l1)==0) cycle
                    if (mod(L+l1+l2, 2)) then
                       C = -Alm%TEB(1,l2,m2)*Noise_Alm(7)%TEB(1,l1,m1)
                       D = -Alm%TEB(2,l2,m2)*Noise_Alm(8)%TEB(1,l1,m1)*cmplx(0,-1)
                    else
                       C = -Alm%TEB(1,l2,m2)*Noise_Alm(8)%TEB(1,l1,m1)*cmplx(0,-1)
                       D = -Alm%TEB(2,l2,m2)*Noise_Alm(7)%TEB(1,l1,m1)

                    Alm_con%TEB(2,L,m) = Alm_con%TEB(2,l,m)+(-1)**m * threejs(l1) &
                                     * threejs1(l1) * C*&
                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
                    Alm_con%TEB(3,L,m) = Alm_con%TEB(3,l,m)+cmplx(0,-1)*(-1)**m &
                                         * threejs(l1)* threejs1(l1)*D*&
                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))

                end do
              end do
           end do
        end do
!        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L) 
!       write (2,*) L, AL(L)
      end do





     !dipole leakage 9, 10, d_a, d_b spin 1


      do L=1, lmax_phi
        do l2=2, lmax_est
           call GetThreeJs(threejs(abs(l2-L)),l2,L,1,-2)  !(L, l1, l2)(-2, 0, 2)
           do m = 1, L
              do m2 = 1, l2
                call GetThreeJs(threejs1(abs(l2-L)),l2,L,m2,-m)
                do l1 = max(2,abs(l2-L)), min(l2+L,lmax_est)
                    if (threejs(l1)==0 .or. threejs1(l1)==0) cycle
                    if (mod(L+l1+l2, 2)) then
                       D = -beam_width*sqrt(real(l2*(l2+1), dp))*Alm%TEB(2,l2,m2)*Noise_Alm(9)%TEB(1,l1,m1)*cmplx(0,1)
                       C = -beam_width*sqrt(real(l2*(l2+1), dp))*Alm%TEB(2,l2,m2)*Noise_Alm(10)%TEB(1,l1,m1)
                    else
                       C = -beam_width*sqrt(real(l2*(l2+1), dp))*Alm%TEB(2,l2,m2)*Noise_Alm(9)%TEB(1,l1,m1)*cmplx(0,-1)
                       D = -beam_width*sqrt(real(l2*(l2+1), dp))*Alm%TEB(2,l2,m2)*Noise_Alm(10)%TEB(1,l1,m1)


                    Alm_con%TEB(2,L,m) = Alm_con%TEB(2,l,m)+(-1)**m * threejs(l1) &
                                     * threejs1(l1) * C*&
                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
                    Alm_con%TEB(3,L,m) = Alm_con%TEB(3,l,m)+cmplx(0,-1)*(-1)**m &
                                         * threejs(l1)* threejs1(l1)*D*&
                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))

                end do
              end do
           end do
        end do
!        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L) 
!       write (2,*) L, AL(L)
      end do



     !quadrupole leakage 11 p

      do L=1, lmax_phi
        do l2=2, lmax_est
           call GetThreeJs(threejs(abs(l2-L)),l2,L,2,-2)  !(L, l1, l2)(-2, 0, 2)
           do m = 1, L
              do m2 = 1, l2
                call GetThreeJs(threejs1(abs(l2-L)),l2,L,m2,-m)
                do l1 = max(2,abs(l2-L)), min(l2+L,lmax_est)
                    if (threejs(l1)==0 .or. threejs1(l1)==0) cycle
                    if (mod(L+l1+l2, 2)) then
                         C = beam_width**2*sqrt(real((l2+2)*(l2+1)*l2*(l2-1))) *Alm%TEB(2,l2,m2)*Noise_Alm(11)%TEB(1,l1,m1)*cmplx(0,1)
                    else
                         D = beam_width**2*sqrt(real((l2+2)*(l2+1)*l2*(l2-1))) *Alm%TEB(2,l2,m2)*Noise_Alm(11)%TEB(1,l1,m1)


                    Alm_con%TEB(2,L,m) = Alm_con%TEB(2,l,m)+(-1)**m * threejs(l1) &
                                     * threejs1(l1) * C*&
                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))
                    Alm_con%TEB(3,L,m) = Alm_con%TEB(3,l,m)+cmplx(0,-1)*(-1)**m &
                                         * threejs(l1)* threejs1(l1)*D*&
                                     sqrt(real((2*l1+1)*(2*l2+1)*(2*L+1),dp)/(HO_fourpi))

                end do
              end do
           end do
        end do
!        AL(L) =real( L*(L+1)*(2*L+1), dp)/AL(L) 
!       write (2,*) L, AL(L)
      end do






 




      end subroutine Noise
