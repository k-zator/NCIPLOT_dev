! Julia Conteras-Garcia <julia.contreras.garcia@gmail.com>,
! Erin R. Johnson <ejohnson29@ucmerced.edu>,
! A. Otero-de-la-Roza <aoterodelaroza@ucmerced.edu>
! Weitao Yang <weitao.yang@duke.edu>,
! Roberto A. Boto <robalboto@gmail.com>,
! Chaoyou Quan  <quanchaoyu@gmail.com>,
! Ruben Laplaza <rlaplaza@lct.jussieu.fr>
! Erna Wieduwilt <erna@sdu.dk>
!
! nciplot is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! This is NCIPLOT Ver. 4.2.1 alpha
! Integration additivity problems. Solved. 
! Reproducibity proplems: Connected with the parallelization of the 
! calcprops_wfn routine. Lines 91-93/252 of routine calcprops_wfn commented. 

program nciplot
   use param
   use tools_io
   use tools_math
   use reader
   use props
!   USE IEEE_ARITHMETIC
   use, intrinsic :: iso_c_binding
   implicit none

   ! python ones
   integer(c_int) :: py_status
   character(kind=c_char, len=1000) :: command_ncicluster
   integer :: iounit_p1
   logical :: doclustering

   integer, parameter :: mfiles = 100 ! max number of geometry files

   integer :: aux
   real*8 :: value

   character*(mline) :: argv(2), oname
   integer :: argc, nfiles, ifile, idx, istat, ntotal
   integer :: i, j, k, nn0, nnf, lp, n1, molid, it1,it2,it3
   character*(mline) :: filein, line, oline, word
   logical :: ok, ispromol,isnotcube
   real*8 :: rdum, deltag
   ! the molecular info
   type(molecule), allocatable :: m(:)
   ! logical units
   integer :: lugc, ludc, luvmd, ludat
   ! cubes
   real*8, allocatable, dimension(:, :, :) :: crho, cgrad
   ! ligand, intermolecular keyword
   logical :: ligand, inter, intra
   real*8 :: rthres
   integer :: udat0
   ! radius and cube keywords
   logical :: autor, flag_dens_neg
   real*8 :: x(3), xinit(3), xmax(3), xinc(3), xcom(3)
   integer :: nstep(3)
   ! noutput keyword
   integer :: noutput
   ! cutoffs
   real*8 :: rhocut, dimcut, maxrho, minrho
   ! cutplot
   real*8 :: rhoplot, isordg
   ! discarding rho parameter
   real*8 :: rhoparam, rhoparam2
   ! properties of rho
   real*8 :: rho, grad(3), dimgrad, grad2, hess(3, 3)
   integer, parameter :: mfrag = 100 ! max number of fragments
   real*8 :: rhom(mfrag)
   ! eispack
   real*8 :: wk1(3), wk2(3), heigs(3), hvecs(3, 3)
   ! Fragments
   integer :: nfrag
   logical :: autofrag
   ! chk file
   real*8 :: xinc_init(3)
   ! Modification
   logical ::  dointeg
   ! Atcube
   integer :: natommax, igroup, natom
   integer, allocatable, dimension(:, :) :: group
   real*8, allocatable, dimension(:, :) :: xinitat, xmaxat, nstepat

   ! Initializing multilevel grids
   integer :: ind_g, ng
   integer :: indx(3), nstep_coarse(3)
   real*8, allocatable :: fginc(:)
   real*8 :: xinc_coarse(3)
   real*8, allocatable, dimension(:, :, :) :: tmp_crho, tmp_cgrad
   logical :: flag, firstgrid
   logical, allocatable :: rmbox_coarse(:, :, :), tmp_rmbox(:, :, :), rmpoint_coarse(:,:,:)
   integer :: cr, c0, c1, c2, c3, c4, c5, c6
   integer :: i0, j0, k0
   integer :: lumesh, lusol
   logical, allocatable :: vert_use(:, :, :)
   real*8 :: sum_rhon_vol(9), sum_signrhon_vol(7)
   real*8 :: sum_rhon_area(9)
   real*8, allocatable ::  rho_n(:)
   real*8, allocatable :: crho_n(:, :, :, :), cheig(:, :, :)
   ! crho_n is a variable created for multilevel grids. It is equivalen to rhom (density per fragment)
   real*8, allocatable :: tmp_crho_n(:, :, :, :), tmp_cheigs(:, :, :)
   real*8 :: percent

   ! Range integration variables
   logical :: dorange, flag_range(2, 2, 2), IsInter(2,2,2)
   integer :: nranges, l, l1(2), j1(2), k1(2),i1(2), ll, jj, kk
   real*8, allocatable :: srhorange(:, :),rho_range(:)
   real*8 :: upperbound, lowerbound,total_rho,sumrangedensity, densitydifference
   logical, allocatable:: tmp_rmbox_range(:,:,:,:), tmp_rmbox_range_tmp(:,:,:), box_in_range(:,:,:,:)
   !===============================================================================!
   ! System clock to measure run time.
   !===============================================================================!
   call system_clock(count_rate=cr)
   call system_clock(count=c0)

   call param_init()
   !===============================================================================!
   ! Reading positional parameters/arguments given to code. Input file to read.
   !===============================================================================!
   call getargs(argc, argv)
   if (argc >= 1) then
      open (uin, file=argv(1), status='old')
      if (argc == 2) then
         open (uout, file=argv(2), status='unknown')
      endif
   endif

   !open(unit=99, file='debug.log', status='unknown')
   !===============================================================================!
   ! Internal clock starts! Header drops.
   !===============================================================================!
   call header()
   call tictac(' # Start')

   !===============================================================================!
   ! Check number of structure/wavefunction files, read them, check them, group them
   !===============================================================================!
   read (uin, *) nfiles   ! number of files to read
   if (nfiles > mfiles) call error('nciplot', 'too many files, increase mfiles', faterr)
   allocate (m(nfiles), stat=istat)
   if (istat /= 0) call error('nciplot', 'could not allocate memory for molecules', faterr)
   do ifile = 1, nfiles ! read files
      read (uin, '(a)') filein ! filein is read
      filein = trim(adjustl(filein))
      inquire (file=filein, exist=ok) ! check for existence
      if (.not. ok) &
         call error('nciplot', 'requested file does not exist: '//trim(filein), faterr)
      m(ifile) = readfile(filein)  ! reading wave function file (wfn or wfx) or structure (xyz file) or cube (cube file)
      if (ifile == 1) oname = filein(1:index(filein, '.', .true.) - 1)
!        write(6,*) oname
      do i = 1, m(ifile)%n       ! assigning atoms to fragments, every atom in ifile to fragment ifile
         m(ifile)%ifrag(i) = ifile
      end do
   enddo
   nfrag = nfiles ! by default each file defines a fragment
   autofrag = .true.
   if (nfrag > mfrag) then ! too many fragments?
      call error('nciplot', 'too many fragments. Increase mfrag', faterr)
   end if

   !===============================================================================!
   ! Stop if user tries to mix computed and promolecular densities.
   !===============================================================================!
   if (any(m(:)%ifile == ifile_xyz) .and. any(m(:)%ifile == ifile_wfn)) then
      call error('nciplot', 'mixing xyz and wfn  not allowed', faterr)
   end if
   if (any(m(:)%ifile == ifile_xyz) .and. any(m(:)%ifile == ifile_wfx)) then
      call error('nciplot', 'mixing xyz and wfx not allowed', faterr)
   end if
   if (any(m(:)%ifile == ifile_wfn) .and. any(m(:)%ifile == ifile_wfx)) then
      call error('nciplot', 'mixing wfn and wfx is not advised', warning)
   end if
   if (any(m(:)%ifile == ifile_wfn) .and. any(m(:)%ifile == ifile_cube)) then
      call error('nciplot', 'mixing wfn and cube is not allowed', faterr)
   end if
   if (any(m(:)%ifile == ifile_xyz) .and. any(m(:)%ifile == ifile_cube)) then
      call error('nciplot', 'mixing xyz and cube is not allowed', faterr)
   end if

   !===============================================================================!
   ! Implicitely checking if the run mode is promolecular or not.
   !===============================================================================!
   ispromol = .not. (all(m(:)%ifile == ifile_wfn) .or. (all(m(:)%ifile == ifile_wfx)) .or. (all(m(:)%ifile == ifile_cube)))

   ! read density grids (props.f90)
   ! by default, use density grids for heavier or charged atoms if promolecularity is on
   if (ispromol) then
   call init_rhogrid(m, nfiles)
! this call was initially out of the if loop, but seems best inside
      do i = 1, nfiles
         if (m(i)%ifile == ifile_xyz .and. (any(m(i)%z > atomic_zmax) .or. any(m(i)%q > 0))) then
            m(i)%ifile = ifile_grd
            do j = 1, m(i)%n
               if (m(i)%z(j) > atomic_zmax .and. .not. grd(iztype(m(i)%z(j), m(i)%q(j)))%init) then
                  call error('nciplot', 'Some atomic density grids for heavy atoms are needed but not initialized', faterr)
               end if
            end do
         end if
      end do
   end if
    isnotcube = .not. (all(m(:)%ifile == ifile_cube))
   ! This checks whether the calculation uses cube files
   !===============================================================================!
   ! Input files read and processed. Set defaults for running NCIPLOT now.
   !===============================================================================!
   rhocut = 0.5d0 ! density cutoff
   dimcut = 1.0d0  ! RDG cutoff
   if (.not.isnotcube) then
      rhoplot=0.05d0 ! density cutoff for printing, in cube  case
      isordg = 0.5d0
   endif
   if (isnotcube) then
      xinc = 0.1d0/bohrtoa   ! grid step
   if (any(m(:)%ifile == ifile_wfn) .or. any(m(:)%ifile == ifile_wfx)) then
      isordg = 0.5d0  ! RDG isosurface
      rhoplot = 0.05d0 ! Density isosurface
   else
      isordg = 0.3d0  ! RDG isosurface
      rhoplot = 0.07d0 ! Density isosurface
   end if
   rhoparam = 0.95d0  ! cutoffs for inter or intramolecularity
   rhoparam2 = 0.75d0 !
   noutput = 3        ! number of outputs
   udat0 = 1
   autor = .true.     ! build the cube automatically
   ligand = .false.   ! ligand keyword
   inter = .false.    ! intermolecular keyword
   rthres = 0.75d0/bohrtoa     ! box limits around the molecule
   dointeg = .false.  ! integrating properties or not
   doclustering = .false. ! clustering python script NCICLUSTER
   dorange = .false.  ! do not integrate range
   firstgrid = .true. ! flag for the first adaptive grid run
   if (.not. allocated(fginc)) then ! default setting CG2FG 3 4 2 1
      ng = 4
      allocate (fginc(ng))
      fginc = (/8, 4, 2, 1/)
   end if
     !===============================================================================!
   ! Estimating box around the molecule using xinit and xmax for the main system.
   !===============================================================================!
   xinit = m(1)%x(:, 1)
   xmax = m(1)%x(:, 1)
   do i = 1, nfiles
      do j = 1, m(i)%n
         xinit = min(xinit, m(i)%x(:, j))
         xmax = max(xmax, m(i)%x(:, j))
      enddo
   enddo
   ntotal = 0
   do i = 1, nfiles
      ntotal = ntotal + m(i)%n    ! compute the total number of atoms
   enddo
 else
    xinc = m(1)%xinc0 ! if cube file, then only one molecule (file), and increment is contained in molecule type
    xinit= m(1)%xinit0 ! same for initial point in grid
    xmax= m(1)%xmax0 ! same for last
    xcom = m(1)%xcom0 ! Number of points given by the input cube
    noutput = 3
   do i = 1, nfiles
      ntotal = ntotal + m(i)%n    ! compute the total number of atoms
   enddo
end if !isnotcube
   !===============================================================================!
   ! Read optional keywords.
   ! Accepted keywords in this version:
   ! - RTHRES
   ! - LIGAND
   ! - RADIUS
   ! - INTERMOLECULAR
   ! - ONAME
   ! - INCREMENTS
   ! - OUTPUT
   ! - CUBEparam
   ! - ATCUBE
   ! - FRAGMENT
   ! - CUTOFFS
   ! - CUTPLOT
   ! - ISORDG
   ! - INTERCUT
   ! - DGRID
   ! - RANGE
   ! - CG2FG
   ! - INTEGRATE
   ! - FINE
   ! - ULTRAFINE
   ! - COARSE
   ! - CLUSTERING
   !===============================================================================!

do while (.true.) 
      read (uin, '(a)', end=11) line
      line = trim(adjustl(line))
      oline = line
      call upper(line)
      if (line(1:1) == "#") cycle ! skip comments
      if (len(trim(line)) < 1) cycle ! skip blank lines
      idx = index(line, ' ')
      word = line(1:idx - 1)
      line = line(idx:)
      oline = oline(idx:)
      select case (trim(word))
      case ("ONAME")         ! output name
         read (oline, '(a)') oname
         oname = trim(adjustl(oname))
         oname = oname(1:index(oname, ' '))

      case ("OUTPUT")
         read (line, *) noutput          ! number of output files

      case ("CUTOFFS")          ! density and RDG cutoffs
         read (line, *) rhocut, dimcut

      case ("CUTPLOT")          ! density cutoff used in the VMD script
         read (line, *) rhoplot, isordg

      case ("ISORDG")           !RDG isosurface used in the RDG script
         read (line, *) isordg

      case ("CLUSTERING")  ! integration
         doclustering = .true.              ! python script cluster
         write (uout, 138) 
      !case default ! something else is read
      !   call error('nciplot', 'Don''t know what to do with '//trim(word)//' keyword', faterr)
      end select

   enddo

   
11   continue 
   rewind(uin)
   
   do while (isnotcube)
      read (uin, '(a)', end=13) line
      line = trim(adjustl(line))
      oline = line
      call upper(line)
      if (line(1:1) == "#") cycle ! skip comments
      if (len(trim(line)) < 1) cycle ! skip blank lines
      idx = index(line, ' ')
      word = line(1:idx - 1)
      line = line(idx:)
      oline = oline(idx:)
      select case (trim(word))

      case ("RTHRES")       ! extra box limits
         read (line, *) rthres
         rthres = max(rthres, 1d-3)
         rthres = rthres/bohrtoa ! change to bohr

      case ("LIGAND")
         ligand = .true.              ! ligand option
         inter = .true.               ! intermolecular option automatically on
         read (line, *) udat0, rthres ! system in (udat0) as center
         rthres = max(rthres, 1d-3)
         rthres = rthres/bohrtoa      ! change to bohr

      case ("INTERMOLECULAR")
         inter = .true.              ! intermolecular option

      case ("RADIUS")
         autor = .false.
         read (line, *) x, rdum      ! center of the box
         rdum = max(rdum, 1d-3)
         xinit = (x - rdum)/bohrtoa  ! box limits
         xmax = (x + rdum)/bohrtoa

      case ("CUBEparam")         ! defining cube limits from coordinates in angstroms. Example:
         autor = .false.    !CUBE x0,y0,z0,x1,y1,z1 format
         read (line, *) xinit, xmax

      case ("ATCUBE")          ! defining cube limits from atoms. Example:
         autor = .false.       !ATCUBE
         xinit = 1d40          !ifile atom1,atom2,...,atomn
         xmax = -1d40          !END
         natommax = 1
         do i = 1, nfiles
            natommax = max(natommax, m(i)%n)
         enddo

         allocate (group(ntotal, natommax))
         allocate (xmaxat(ntotal, 3))
         allocate (xinitat(ntotal, 3))
         allocate (nstepat(ntotal, 3))

         igroup = 0
         do while (.true.)
            igroup = igroup + 1
            read (uin, '(a)') line
            line = trim(adjustl(line))
            call upper(line)
            if (line(1:1) == "#") cycle ! skip comments
            if (len(trim(line)) < 1) cycle ! skip blank lines
            lp = 1
            idx = index(line, ' ')
            word = line(1:idx - 1)
            if (trim(word) /= "END" .and. trim(word) /= "ENDATCUBE") then
               ok = isinteger(ifile, line, lp)
               if (.not. ok) call error('nciplot', 'bad atcube syntax', faterr)
               if (ifile < 1 .or. ifile > nfiles) call error('nciplot', 'atcube: wrong file number', faterr)
               ok = isinteger(n1, line, lp)
               natom = 0
               do while (ok)
                  if (n1 < 1 .or. n1 > m(ifile)%n) call error('nciplot', 'atcube: wrong atom number', faterr)
                  natom = natom + 1
                  group(igroup, natom) = n1
                  xinit = min(xinit, m(ifile)%x(:, n1))
                  xmax = max(xmax, m(ifile)%x(:, n1))
                  ok = isinteger(n1, line, lp)
               end do
            else
               exit
            end if

            xinitat(igroup, :) = xinit - rthres
            xmaxat(igroup, :) = xmax + rthres
            nstepat(igroup, :) = abs(ceiling((xmax - xinit)/xinc))
         end do

      case ("FRAGMENT")          !defining fragments. Example:
         if (autofrag) then      !FRAGMENT
            nfrag = 0            !ifile atom1, atom2,...,atomn
            do ifile = 1, nfiles !END
               do i = 1, m(ifile)%n
                  m(ifile)%ifrag(i) = 0
               end do
            end do
         end if
         autofrag = .false.
         inter = .true.

         nfrag = nfrag + 1
         do while (.true.)
            read (uin, '(a)') line
            line = trim(adjustl(line))
            call upper(line)
            if (line(1:1) == "#") cycle ! skip comments
            if (len(trim(line)) < 1) cycle ! skip blank lines
            lp = 1
            idx = index(line, ' ')
            word = line(1:idx - 1)
            if (trim(word) /= "END" .and. trim(word) /= "ENDFRAGMENT") then
               ok = isinteger(ifile, line, lp)
               if (.not. ok) call error('nciplot', 'bad fragment syntax', faterr)
               if (ifile < 1 .or. ifile > nfiles) call error('nciplot', 'fragment: wrong file number', faterr)
               ok = isinteger(n1, line, lp)
               do while (ok)
                  if (n1 < 1 .or. n1 > m(ifile)%n) call error('nciplot', 'fragment: wrong atom number', faterr)
                  m(ifile)%ifrag(n1) = nfrag
                  ok = isinteger(n1, line, lp)
               end do
            else
               exit
            end if
         end do

      case ("INCREMENTS")   ! grid increments
         read (line, *) xinc
         xinc = max(xinc, 1d-4)
         xinc = xinc/bohrtoa ! transforming to angstrom to bohr

      case ("FINE")          ! FINE defaults
         xinc = 0.05d0/bohrtoa   ! grid step
         if (allocated(fginc)) then
            deallocate (fginc)
         end if
         ng = 4
         allocate (fginc(ng))
         fginc = (/12, 8, 4, 1/)

      case ("ULTRAFINE")          ! ULTRAFINE defaults
         xinc = 0.025d0/bohrtoa   ! grid step
         if (allocated(fginc)) then
            deallocate (fginc)
         end if
         ng = 4
         allocate (fginc(ng))
         fginc = (/24, 16, 8, 1/)

      case ("COARSE")          ! COARSE defaults
         xinc = 0.15d0/bohrtoa   ! grid step
         if (allocated(fginc)) then
            deallocate (fginc)
         end if
         allocate (fginc(4))
         fginc = (/8, 4, 2, 1/)
         if (.not. inter) then
            rthres = 0.5d0/bohrtoa     ! box limits around the molecule
         end if

      case ("INTERCUT")          ! cutoffs for intermolecularity definition
         read (line, *) rhoparam, rhoparam2

      case ("DGRID")            ! using grids for promolecular densities
         do i = 1, nfiles
            if (m(i)%ifile == ifile_xyz) m(i)%ifile = ifile_grd
         end do

      case ("CG2FG") ! coarse grid to fine grid multi-level
         read (line, *) ng ! number of multi-level grids
         if (allocated(fginc)) then
            deallocate (fginc)
         end if
         allocate (fginc(ng)) ! factors of grid increments. Example:
         read (line(:), *) ng, fginc ! CG2FG 4 8 4 2 1

      case ("RANGE")  ! range integration
         dointeg = .true.              ! integration is required as well to apply the s=dimcut cutoff
         read (line, *) nranges ! number of ranges
         if (nranges .le. 0) then
            call error('nciplot', 'No ranges were given', faterr)
         else
            dorange = .true.
            allocate (srhorange(nranges, 2), stat=istat)
            if (istat /= 0) call error('nciplot', 'could not allocate memory for range intervals', faterr)
            do i = 1, nranges
               read (uin, *) srhorange(i, :)
               do j = 1, 2
                  if (abs(srhorange(i, j)) .lt. 1d-30) then
                     srhorange(i, j) = srhorange(i, j) + 1d-30
                  endif
               enddo
            enddo
         endif
      case ("INTEGRATE")  ! integration
         dointeg = .true.              ! integrate
!        end if
        end select 
end do
    ! all the previous options are not compatible with cube file inputs
    ! following are

13 continue


   !===============================================================================!
   ! Defining box in detail now.
   !===============================================================================!
    if (.not. isnotcube) then
         autor=.false.
         nstep= abs(xcom) !dimensions given on principal cube
    end if
    ! if we have a cube file, parameters are already known

   if (autor) then ! automatically build grid, it has not been done
      if (ligand) then ! ligand mode is enabled
         xinit = m(udat0)%x(:, 1)
         xmax = m(udat0)%x(:, 1)
         do j = 1, m(udat0)%n
            xinit = min(xinit, m(udat0)%x(:, j))
            xmax = max(xmax, m(udat0)%x(:, j))
         end do
      end if
      xinit = xinit - rthres
      xmax = xmax + rthres
      nstep = abs(ceiling((xmax - xinit)/xinc)) !number of grid steps case not cube
   end if

   !nstep = abs(ceiling(xmax - xinit)/(xinc)) !number of grid steps
   !nstep = x
   !write(*,*) nstep
   ! Arregle los valores de las dimensiones del cubo, ceiling lo tiraba siempre hacia arriba.
   !===============================================================================!
   ! Information for user and output files. Default logical units first.
   !===============================================================================!
   lugc = -1 ! RDG logical unit
   ludc = -1 ! Density logical unit
   luvmd = -1 ! VMD logical unit
   write (uout, 131)
   if (inter) then
      write (uout, 132) ! tell user intermolecular mode is on
      if (nfrag .eq. 1) then ! not enough fragments
         call error('nciplot', 'not enough fragments for intermolecular', faterr)
      end if
   end if
   if (ligand) write (uout, 130) trim(m(udat0)%name)
   if (ispromol) write (uout, 133) ! tell user promolecular mode is on
   write (uout, 120)
   write (uout, 110) 'RHO  THRESHOLD   (au):', rhocut
   write (uout, 110) 'RDG  THRESHOLD   (au):', dimcut
   if (inter) write (uout, 110) 'DISCARDING RHO PARAM :', rhoparam
   if (ligand) write (uout, 110) 'RADIAL THRESHOLD  (A):', rthres*bohrtoa
   write (uout, *)
!  write(uout,121) xinit, xmax, xinc, nstep ! this is currently not used because it will be adapted
   if (noutput >= 2) then    ! number of outputs --> 2: Only .CUBE files
      lugc = 9
      ludc = 10
      luvmd = 11
      open (lugc, file=trim(oname)//"-grad.cube")    ! RDG cube file
      open (ludc, file=trim(oname)//"-dens.cube")    ! Density cube file
      open (luvmd, file=trim(oname)//".vmd")         ! VMD script
   endif

   if (noutput == 1 .or. noutput == 3) then
      ludat = 16
      open (ludat, file=trim(oname)//".dat")         ! RDG vs sign(lambda2) file
   else
      ludat = -1
   endif

   !===============================================================================!
   ! Write output files.
   !===============================================================================!
   write (uout, 122)
   if (noutput == 1 .or. noutput == 3) then
      write (uout, 123) trim(oname)//".dat"
   end if

   if (noutput >= 2) then
      write (uout, 124) trim(oname)//"-grad.cube", &
         trim(oname)//"-dens.cube", &
         trim(oname)//".vmd"
   end if
   if (lugc > 0) call write_cube_header(lugc, 'grad_cube', '3d plot, reduced density gradient')
   if (ludc > 0) call write_cube_header(ludc, 'dens_cube', '3d plot, density')

   !===============================================================================!
   ! Start run, using multi-level grids.
   !===============================================================================!
    if (isnotcube) then
   ind_g = 1  ! index of the multi-level grids, starts at 1 always
   xinc_init = xinc ! initial coarse grid
   allocate (rho_n(1:nfiles))

12 continue    ! set grids from coarse to fine
   xinc = fginc(ind_g)*xinc_init
   nstep = ceiling((xmax - xinit)/xinc)

   write (uout, *)    ! punch info
   write (uout, 121) ind_g, xinit, xmax, xinc, nstep
   !===============================================================================!
   ! Allocate memory for density and gradient.
   !===============================================================================!
   if (ind_g .eq. 1) then
      allocate (crho(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), stat=istat)
      if (istat /= 0) call error('nciplot', 'could not allocate memory for density cube', faterr)
      allocate (crho_n(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1, 1:nfiles), stat=istat)
      if (istat /= 0) call error('nciplot', 'could not allocate memory for density cube', faterr)
      allocate (cheig(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), stat=istat)
      if (istat /= 0) call error('nciplot', 'could not allocate memory for density cube', faterr)
      allocate (cgrad(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), stat=istat)
      if (istat /= 0) call error('nciplot', 'could not allocate memory for grad', faterr)
      xinc_coarse = xinc ! initialize just in case
      nstep_coarse = nstep ! initialize just in case
   end if
   !===============================================================================!
   ! Allocate memory for coarse grid.
   !===============================================================================!
   if (ind_g .gt. 1) then
      if (allocated(tmp_crho)) deallocate (tmp_crho)
      allocate (tmp_crho(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), stat=istat)
      call move_alloc(tmp_crho, crho)
      if (allocated(tmp_crho_n)) deallocate (tmp_crho_n)
      allocate (tmp_crho_n(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1, 1:nfiles), stat=istat)
      call move_alloc(tmp_crho_n, crho_n)
      if (allocated(tmp_cheigs)) deallocate (tmp_cheigs)
      allocate (tmp_cheigs(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), stat=istat)
      call move_alloc(tmp_cheigs, cheig)
      if (allocated(tmp_cgrad)) deallocate (tmp_cgrad)
      allocate (tmp_cgrad(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), stat=istat)
      call move_alloc(tmp_cgrad, cgrad)
   end if

   !===============================================================================!
   ! Sanity check.
   !===============================================================================!
   if (.not. firstgrid) then
      do i = 1, 3
         if (ubound(rmbox_coarse, i) .ne. (nstep_coarse(i) - 2)) then
            write (*, *) ubound(rmbox_coarse, i)
            write (*, *) nstep_coarse(i) - 2
            call error('nciplot', 'allocations did not work for adaptive grids', faterr)
         end if
      end do
   end if

   !===============================================================================!
   ! Run adaptive grids.
   !===============================================================================!
   if (ispromol) then    ! promolecular densities
      call system_clock(count=c1)
      !$omp parallel do private (x,rho,grad,hess,heigs,hvecs,wk1,wk2,istat,grad2,&
      !$omp dimgrad,intra,rhom,flag,indx,i0,j0,k0) schedule(dynamic)
      do k = 0, nstep(3) - 1
         do j = 0, nstep(2) - 1
            do i = 0, nstep(1) - 1
               x = xinit + (/i, j, k/)*xinc
               if (.not. firstgrid) then
                  flag = .false.
                  do i0 = max(0, i - 1), i
                     do j0 = max(0, j - 1), j
                        do k0 = max(0, k - 1), k
                           ! For each x, look for i, j, k indexes in the previous coarser grid
                           indx = floor(((/i0, j0, k0/)*xinc)/xinc_coarse)
                           indx = (/min(nstep_coarse(1) - 2, indx(1)), min(nstep_coarse(2) - 2, indx(2)), &
                                    min(nstep_coarse(3) - 2, indx(3))/)
                           if ((.not. flag) .and. (.not. (rmbox_coarse(indx(1), indx(2), indx(3))))) then
                              flag = .true.
                              goto 20
                           end if
                        end do
                     end do
                  end do

                  if (.not. flag) then
                     crho(i, j, k) = 100d0
                     cgrad(i, j, k) = 100d0
                     cheig(i, j, k) = 0d0
                     cycle
                  end if
20                continue

               end if

               ! calculate properties at x
               rho_n = 0d0
               call calcprops_pro(x, m, nfiles, rho, rho_n, rhom(1:nfrag), nfrag, autofrag, &
                                  grad, hess, deltag)
               call rs(3, 3, hess, heigs, 0, hvecs, wk1, wk2, istat)
               rho = max(rho, 1d-30)
               grad2 = dot_product(grad, grad)
               dimgrad = sqrt(grad2)/(const*rho**(4.D0/3.D0))
               intra = inter .and. ((any(rhom(1:nfrag) >= sum(rhom(1:nfrag))*rhoparam)) .or. &
                                    (sum(rhom(1:nfrag)) < rhoparam2*rho))
               if (intra) dimgrad = -dimgrad !checks for interatomic, intra is true iff inter and condition hold
               !$omp critical (cubewrite)
               ! write to cube file, dont take into account rho values to low
               if (rho /= 1d-30) then
                  crho(i, j, k) = sign(rho, heigs(2))*100.D0
                  cgrad(i, j, k) = dimgrad
               else
                  crho(i, j, k) = 100d0
                  cgrad(i, j, k) = 100d0
               end if

               do molid = 1, nfiles
                  crho_n(i, j, k, molid) = rhom(molid)
               enddo
               !$omp end critical (cubewrite)

            end do !i=0,nstep(1)-1
         end do ! j=0,nstep(2)-1
      end do  ! k=0,nstep(3)-1
      !$omp end parallel do
      call system_clock(count=c2)
      write (*, "(A, F6.2, A)") ' Time for computing density & RDG = ', real(dble(c2 - c1)/dble(cr), kind=8), ' secs'

   else  ! wavefunction densities
      if (.not. inter) then
         call system_clock(count=c1)
         call calcprops_wfn(xinit, xinc, nstep, m, nfiles, crho, cgrad, cheig)
         !$omp parallel do private (x,rho,grad,hess,heigs,hvecs,wk1,wk2,istat,grad2,&
         !$omp dimgrad,intra,rhom,flag,indx,i0,j0,k0) schedule(dynamic)
         do k = 0, nstep(3) - 1
            do j = 0, nstep(2) - 1
               do i = 0, nstep(1) - 1
                  x = xinit + (/i, j, k/)*xinc
                  if (.not. firstgrid) then
                     ! check if x is used, not removed
                     flag = .false.
                     do i0 = max(0, i - 1), i
                        do j0 = max(0, j - 1), j
                           do k0 = max(0, k - 1), k
                              ! For each x, look for i, j, k indexes in the previous coarser grid
                              indx = floor(((/i0, j0, k0/)*xinc)/xinc_coarse)
                              indx = (/min(nstep_coarse(1) - 2, indx(1)), min(nstep_coarse(2) - 2, indx(2)), &
                                       min(nstep_coarse(3) - 2, indx(3))/)
                              if ((.not. flag) .and. (.not. (rmbox_coarse(indx(1), indx(2), indx(3))))) then
                                 flag = .true.
                                 goto 21
                              end if
                           end do
                        end do
                     end do

                     if (.not. flag) then
                        crho(i, j, k) = 100d0
                        cgrad(i, j, k) = 100d0
                        cheig(i, j, k) = 0d0
                        cycle
                     end if
21                   continue

                  end if
               end do
            end do
         end do
         call system_clock(count=c2)
         write (*, "(A, F6.2, A)") ' Time for computing density & RDG = ', real(dble(c2 - c1)/dble(cr), kind=8), ' secs'
      else !very experimental wfn intermolecular
         call system_clock(count=c1)
         do molid = 1, nfiles
            call calcprops_id_wfn(xinit, xinc, nstep, m, nfiles, molid, crho, cgrad, cheig)
            crho_n(:, :, :, molid) = crho(:, :, :)
         end do
         call calcprops_wfn(xinit, xinc, nstep, m, nfiles, crho, cgrad, cheig)
         do k = 0, nstep(3) - 1
            do j = 0, nstep(2) - 1
               do i = 0, nstep(1) - 1
                  x = xinit + (/i, j, k/)*xinc
                  if (.not. firstgrid) then
                     ! check if x is used, not removed
                     flag = .false.
                     do i0 = max(0, i - 1), i
                        do j0 = max(0, j - 1), j
                           do k0 = max(0, k - 1), k
                              ! For each x, look for i, j, k indexes in the previous coarser grid
                              indx = floor(((/i0, j0, k0/)*xinc)/xinc_coarse)
                              indx = (/min(nstep_coarse(1) - 2, indx(1)), min(nstep_coarse(2) - 2, indx(2)), &
                                       min(nstep_coarse(3) - 2, indx(3))/)
                              if ((.not. flag) .and. (.not. (rmbox_coarse(indx(1), indx(2), indx(3))))) then
                                 flag = .true.
                                 goto 22
                              end if
                           end do
                        end do
                     end do
                     rho = crho(i,j,k)
                     intra = inter .and. (any(abs(crho_n(i, j, k, 1:nfrag))*100d0 >= abs(crho(i, j, k))*rhoparam) .or. &
                                    (sum(abs(crho_n(i, j, k, 1:nfrag))*100d0) < rhoparam2*abs(rho)))
                     if (intra) then !checks for interatomic
                        cgrad(i, j, k) = -cgrad(i, j, k)
                     end if
                     if (.not. flag) then
                        crho(i, j, k) = 100d0
                        cgrad(i, j, k) = 100d0
                        cheig(i, j, k) = 0d0
                        cycle
                     end if
22                   continue

                  end if
               end do
            end do
         end do
         call system_clock(count=c2)
         write (*, "(A, F6.2, A)") ' Time for computing density & RDG = ', real(dble(c2 - c1)/dble(cr), kind=8), ' secs'
      endif !is inter !very experimental wfn intermolecular
   endif !iswfn

   if ((ind_g .le. ng) .or. (ng .eq. 1)) then
      xinc_coarse = xinc ! record increments of the previous coarse grid
      nstep_coarse = nstep
      firstgrid = .false.
      if (allocated(rmbox_coarse)) then
         deallocate (rmbox_coarse)
         allocate (tmp_rmbox(0:nstep_coarse(1) - 2, 0:nstep_coarse(2) - 2, 0:nstep_coarse(3) - 2), stat=istat)
         if (istat /= 0) call error('nciplot', 'could not allocate memory for tmp_rmbox', faterr)
         call build_rmbox_coarse(rhocut, dimcut, ng, ind_g, fginc, tmp_rmbox, crho, cgrad, nstep_coarse)
         call move_alloc(tmp_rmbox, rmbox_coarse)
      else
         allocate (rmbox_coarse(0:nstep_coarse(1) - 2, 0:nstep_coarse(2) - 2, 0:nstep_coarse(3) - 2), stat=istat)
         if (istat /= 0) call error('nciplot', 'could not allocate memory for rmbox_coarse', faterr)
         call build_rmbox_coarse(rhocut, dimcut, ng, ind_g, fginc, rmbox_coarse, crho, cgrad, nstep_coarse)
      end if
      if (allocated(tmp_rmbox)) then
         deallocate (tmp_rmbox)
      end if
   end if

! loop over multi-level grids
   ind_g = ind_g + 1
   if (ind_g .le. ng) then
      goto 12 ! shameful goto to end multilevel grids.
   end if

   !===============================================================================!
   ! Output of .mesh and .sol files.
   !===============================================================================!
   call system_clock(count=c3)
   if (noutput .eq. 4) then
      allocate (vert_use(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1))
      if (noutput .eq. 4) then
         lumesh = 20
         lusol = 21
      else
         lumesh = -1
         lusol = -1
      endif
      if (noutput .eq. 4) then
         open (lumesh, file=trim(oname)//".mesh")
         open (lusol, file=trim(oname)//".sol")
         call write_mesh_file(lumesh, lusol, xinit, xinc, nstep, cgrad, xinc_coarse, rmbox_coarse, &
                              nstep_coarse, vert_use)
         close (lumesh)
         close (lusol)
      endif
      if (allocated(vert_use)) then
         deallocate (vert_use)
      end if
   end if
else ! is a cube file
    ! first we allocate gradient array
    allocate (cgrad(0:nstep(1)-1, 0:nstep(2)-1, 0:nstep(3)-1), stat=istat)
    if (istat /= 0) call error('nciplot', 'could not allocate memory for grad', faterr)
    allocate (crho(0:nstep(1)-1, 0:nstep(2)-1, 0:nstep(3)-1), stat=istat)
    if (istat /= 0) call error('nciplot', 'could not allocate memory for signed density', faterr)
 
    do it1=0,nstep(1)-1 ! first we initialise the reduced gradient array to 100d0
    do it2=0,nstep(2)-1
    do it3=0,nstep(3)-1
    cgrad(it1,it2,it3) = 100d0
    end do
    end do
    end do
    flag_dens_neg = .FALSE.
    DO it1=3,nstep(1)-3 ! since numerical differentiation reaches point +/- 2
      DO it2=3,nstep(2)-3 ! in the following we use the numerical differentiation to compute Hessian on a grid
			DO it3=3,nstep(3)-3
         if (m(1)%cubedens(it1,it2,it3) .LT. 0.0) then ! Density is negative, a messagge appear and values are 0.0
            flag_dens_neg = .TRUE.
            m(1)%cubedens(it1,it2,it3)=abs(m(1)%cubedens(it1,it2,it3))
         end if

			if ((m(1)%cubedens(it1,it2,it3) .LE. rhocut)) THEN!.and. (m(1)%cubedens(it1,it2,it3)  .GT. 1d-6)) THEN ! we restrict computation of gradient & hessian to low density regions
            cgrad(it1,it2,it3) = (SQRT(                                                               &   ! ||grad_p(r)||
            ((m(1)%cubedens(it1+1,it2,it3) - m(1)%cubedens(it1-1,it2,it3)) / (2.*m(1)%xinc0(1)))**2.d0 + &   ! grad_x
   			((m(1)%cubedens(it1,it2+1,it3) - m(1)%cubedens(it1,it2-1,it3)) / (2.*m(1)%xinc0(2)))**2.d0 + &   ! grad_y
   			((m(1)%cubedens(it1,it2,it3+1) - m(1)%cubedens(it1,it2,it3-1)) / (2.*m(1)%xinc0(3)))**2.d0 )) &  ! grad_z
   			/((2.d0*((3.d0*pi**2.d0)**(1.d0/3.d0)))*(m(1)%cubedens(it1,it2,it3)**(4.d0/3.d0)))                   ! (2*((3*pi**2)**1/3) * (p(r)**4/3))
   			
            hess(1,1) = ((m(1)%cubedens(it1+2,it2,it3)) - 2*(m(1)%cubedens(it1+1,it2,it3)) + (m(1)%cubedens(it1,it2,it3)) ) &
                           /(m(1)%xinc0(1) * m(1)%xinc0(1)) 
   			hess(2,2) = ((m(1)%cubedens(it1,it2+2,it3)) - 2*(m(1)%cubedens(it1,it2+1,it3)) + &
                           (m(1)%cubedens(it1,it2,it3)) )/(m(1)%xinc0(2) * m(1)%xinc0(2)) 
   			hess(3,3) = ((m(1)%cubedens(it1,it2,it3+2)) - 2*(m(1)%cubedens(it1,it2,it3+1)) +&
                            (m(1)%cubedens(it1,it2,it3)) )/(m(1)%xinc0(3) * m(1)%xinc0(3)) 
   			hess(1,2) = ((m(1)%cubedens(it1+1,it2+1,it3)) - (m(1)%cubedens(it1+1,it2,it3)) - &
                           (m(1)%cubedens(it1,it2+1,it3)) + (m(1)%cubedens(it1,it2,it3)) )/(m(1)%xinc0(1) * m(1)%xinc0(2)) 
   			hess(1,3) = ((m(1)%cubedens(it1+1,it2,it3+1)) - (m(1)%cubedens(it1+1,it2,it3)) - &
                           (m(1)%cubedens(it1,it2,it3+1)) + (m(1)%cubedens(it1,it2,it3)) )/(m(1)%xinc0(1) * m(1)%xinc0(3)) 
   			hess(2,3) = ((m(1)%cubedens(it1,it2+1,it3+1)) - (m(1)%cubedens(it1,it2+1,it3)) - &
                           (m(1)%cubedens(it1,it2,it3+1)) + (m(1)%cubedens(it1,it2,it3)) )/(m(1)%xinc0(2) * m(1)%xinc0(3)) 
   			hess(2,1) = hess(1,2)
   			hess(3,1) = hess(1,3)
   			hess(3,2) = hess(2,3)
	
   			call RS(3,3,hess,heigs,0,hvecs,wk1,wk2,istat) ! matrix diagonalisation defined in props
   			crho(it1,it2,it3) = sign(real(m(1)%cubedens(it1,it2,it3)),real(heigs(2)))*100.d0
        else
                                                   ! Should be 0.0
            cgrad(it1,it2,it3) = 100d0
            crho(it1,it2,it3) = 0d0
			END IF
			END DO
		END DO
	END DO
        
   call system_clock(count=c3)
   if (flag_dens_neg) then
      write (uout, 137) 
   endif
end if ! isnotcube
   !===============================================================================!
   ! Write .dat file.
   !===============================================================================!
   !REAL :: aux
   aux = 1
   do k = 0, nstep(3) - 1
      do j = 0, nstep(2) - 1
         do i = 0, nstep(1) - 1
            ! fragments for the wfn case
            intra = (cgrad(i, j, k) < 0d0)
            cgrad(i, j, k) = abs(cgrad(i, j, k))
            dimgrad = cgrad(i, j, k)
            rho = crho(i, j, k)/100d0 ! why dividing by 100 ? , true no idea
            ! write the dat file
            if (ludat > 0 .and. .not. intra .and. (abs(rho) < rhocut) .and. (dimgrad < dimcut) .and. &
                abs(rho) > 1d-30) then
               write (ludat, '(1p,E18.10,E18.10)') rho, dimgrad
            endif ! rhocut/dimcut
   
            ! prepare the cube files !!modifJ
!            if (isnotcube) then ! not same cut-off for cube and non-cube inputfiles
             if ((abs(rho) > rhoplot) .or. (dimgrad > dimcut)) then
                 cgrad(i, j, k) = 101d0
             endif !rho cutoff
             if  (intra) then ! intermolecular points also to 100
                  cgrad(i, j, k) = 101d0
             endif
!           else
            if ((abs(rho) == 0.0) .and. (cgrad(i,j,k) .NE. 100.)) then
               cgrad(i,j,k) = 101d0  !NaN values has rho = 0.0, cgrad in that position would be 0
            endif
         end do
      end do
   end do

   !===============================================================================!
   ! Write cube files.
   !===============================================================================!
   if (ludc > 0) call write_cube_body(ludc, nstep, crho)          ! density
   if (lugc > 0) call write_cube_body(lugc, nstep, cgrad)         ! RDG
  
   call system_clock(count=c4)
   write (*, "(A, F6.2, A)") ' Time for writing outputs = ', real(dble(c4 - c3)/dble(cr), kind=8), ' secs'
  
   !===============================================================================!
   ! Integration for promolecular systems. Box removal.
  !===============================================================================!
   if (dointeg) then    ! dointeg
      if (ispromol) then
         !$omp parallel do private (x,rho,grad,hess,heigs,hvecs,wk1,wk2,istat,grad2,&
         !$omp dimgrad,intra,rhom) schedule(dynamic)
         do k = 0, nstep(3) - 2
            do j = 0, nstep(2) - 2
               do i = 0, nstep(1) - 2
                  i1 = (/i, i + 1/)
                  j1 = (/j, j + 1/)
                  k1 = (/k, k + 1/)
                  if (inter) then 
                       do l =1, nfrag     
                          IsInter =((abs(crho_n(i1, j1, k1, l)) .ge. abs(crho(i1, j1, k1))*rhoparam))
                           if (count(IsInter).gt.0) then
                               rmbox_coarse(i, j, k) = .true. !inactive
                           end if
                       enddo    
                  end if
                !  if (((dimgrad > dimcut) .and. .not. rmbox_coarse(i, j, k))) then
                !     rmbox_coarse(i, j, k) = .true. !inactive if dimgrad > dimcut
                !  endif
               end do !i = 0, nstep(3) - 2
            end do !j = 0, nstep(3) - 2
         end do !k = 0, nstep(3) - 2
         !$omp end parallel do
         percent = real(count(rmbox_coarse), kind=8)/(real(size(rmbox_coarse), kind=8))*100d0
         write (*, '(F6.2, A)') percent, '% of small boxes removed for promolecular integration'
      endif  ! ispromol
   endif !dointeg
  
   !===============================================================================!
   ! Integration for non-promolecular systems. Box removal.
  !===============================================================================!
   if (dointeg) then
      if (.not. ispromol) then
         do k = 0, nstep(3) - 2
            do j = 0, nstep(2) - 2
               do i = 0, nstep(1) - 2
                  x = xinit + (/i, j, k/)*xinc
                  rho = abs(crho(i, j, k))/100d0
                  dimgrad = abs(cgrad(i, j, k))
                  if (inter) then
                     if (any(abs(crho_n(i, j, k, 1:nfrag))*100d0 >= abs(crho(i, j, k)*rhoparam))) then
                        rmbox_coarse(i, j, k) = .true.
                     end if
                  end if
                  if (((dimgrad > dimcut) .and. .not. rmbox_coarse(i, j, k))) then
                     rmbox_coarse(i, j, k) = .true. !inactive
                  endif ! rhocut/dimcut
               enddo  !k = 0,nstep(3)-1
            enddo !j = 0,nstep(2)-1
         enddo  !i = 0,nstep(1)-1
         percent = real(count(rmbox_coarse), kind=8)/(real(size(rmbox_coarse), kind=8))*100d0
         write (*, '(F6.2, A)') percent, '% of small boxes removed for density integration'
      endif !not ispromol
   endif !dointeg
  
   !===============================================================================!
   ! Integration and printing, uses the defined rmbox_coarse.
   !===============================================================================!
   ! compute geometry date in the region enclosed by the RDG isosurface:
   ! integral of rho^n (n = 1, 1.5, 2, 2.5, 3, 4/3, 5/3, 0) and rho1*rho2, respectively over the volume and the surface: sum_rhon_vol, sum_rhon_area
   if (dointeg) then
   ! Refine rmbox. Remove boxes with points out of the rhocut range.
   ! Otherwise the total integration is not recover by the integration by ranges.
       do k = 0, nstep(3) - 2
          do j = 0, nstep(2) - 2
             do l = 0, nstep(1) - 2
                 if (.not. (rmbox_coarse(l, j, k))) then  
                      do kk = k,k+1
                         do jj=j,j+1
                             do ll=l,l+1
                                 if ((abs((crho(ll, jj, kk)/100d0)) .gt. rhocut) .and. (cgrad(ll,jj,kk).gt.0)) then
                                    rmbox_coarse(l,j,k)=.true.      
                                 end if    
                              end do 
                           end do
                       end do 
                    end if
                end do 
             end do 
         end do     
 
      call dataGeom(sum_rhon_vol, sum_rhon_area, sum_signrhon_vol, xinc, nstep, crho, crho_n, rmbox_coarse, nfiles)
    !  call dataGeom_points(sum_rhon_vol, sum_signrhon_vol, xinc, nstep, crho, crho_n, &
    !              rmbox_coarse,rmbox_coarse, nfiles)
      write (uout, 117) sum_rhon_vol, sum_signrhon_vol, sum_rhon_area
      call system_clock(count=c5)
      write (*, "(A, F6.2, A)") ' Time for integration by cubes = ', real(dble(c5 - c4)/dble(cr), kind=8), ' secs'
      Total_rho=sum_rhon_vol(1) 
 
   end if
 
   !===============================================================================!
   ! Range integration. Also uses the previously defined rmbox_coarse.
   !===============================================================================!

   if (dorange) then
      allocate (tmp_rmbox(0:nstep(1) - 2, 0:nstep(2) - 2, 0:nstep(3) - 2), stat=istat)
      if (istat /= 0) call error('nciplot', 'could not allocate memory for tmp_rmbox', faterr)
  
      allocate (tmp_rmbox_range_tmp(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), stat=istat)
      if (istat /= 0) call error('nciplot', 'could not allocate memory for tmp_rmbox_range_tmp', faterr) 
  
  
      allocate (tmp_rmbox_range(0:nstep(1) - 2, 0:nstep(2) - 2, 0:nstep(3) - 2,1:nranges), stat=istat)
      if (istat /= 0) call error('nciplot', 'could not allocate memory for tmp_rmbox_range', faterr)
  
      allocate (box_in_range(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1,1:nranges), stat=istat)
      if (istat /= 0) call error('nciplot', 'could not allocate memory for box_in_range', faterr)
  
      allocate (rho_range(1:nranges),stat=istat)    
  
      write (uout, 134)
      !call DoRangeWarning() ! warning: change in active box criterion
  
      tmp_rmbox = .true. 
      tmp_rmbox_range = .true. 
      tmp_rmbox_range_tmp = .true. 
      box_in_range = .true.
  
   
      do i = 1, nranges
         if (srhorange(i, 1) .lt. srhorange(i, 2)) then
            upperbound = srhorange(i, 2)
            lowerbound = srhorange(i, 1)
         else
            upperbound = srhorange(i, 1)
            lowerbound = srhorange(i, 2)
         endif
         if ((upperbound .gt. rhocut) .or. (abs(lowerbound) .gt. rhocut)) then
            call DoRangeWarning2() ! warning: range outside rhocut
         endif
         
    
         do k = 0, nstep(3) - 2
            do j = 0, nstep(2) - 2
               do l = 0, nstep(1) - 2
                   if (.not. (rmbox_coarse(l, j, k))) then  
                       do kk = k,k+1
                          do jj=j,j+1
                             do ll=l,l+1
                                if (((crho(ll, jj, kk)/100d0) .gt. lowerbound) .and. &  
                                         ((crho(ll, jj, kk)/100d0) .lt. upperbound) .and. (cgrad(ll,jj,kk).gt.0)) then 
                                     
                                    box_in_range(ll,jj,kk,i)=.false.    
                                
                              
                               endif            
                             enddo
                          enddo
                       enddo    
                    endif
               enddo
            enddo
          enddo
        enddo
      
    
      ! Integration 
      do i= 1, nranges
         tmp_rmbox_range_tmp(:,:,:)= box_in_range(:,:,:,i)
         call dataGeom_points(sum_rhon_vol, sum_signrhon_vol, xinc, nstep, crho, crho_n, &
                  rmbox_coarse,tmp_rmbox_range_tmp, nfiles)
         write (uout, 135) srhorange(i,1), srhorange(i,2), sum_rhon_vol, sum_signrhon_vol 
         rho_range(i) = sum_rhon_vol(1)
      enddo
      ! Additivy check
      sumrangedensity = sum(rho_range)
      densitydifference = Total_rho - sumrangedensity 
      write (uout, 136) Total_rho, sumrangedensity ,densitydifference
   
   endif ! if range
   call system_clock(count=c6)

   !===============================================================================!
   ! Write .dat and cube files if pprint is set.
   !===============================================================================!
   if (.true.) then
      do k = 0, nstep(3) - 1
         do j = 0, nstep(2) - 1
            do i = 0, nstep(1) - 1
               ! fragments for the wfn case
               intra = (cgrad(i, j, k) < 0d0)
               cgrad(i, j, k) = abs(cgrad(i, j, k))
               dimgrad = cgrad(i, j, k)
               rho = crho(i, j, k)/100d0
               ! write the dat file
               if (ludat > 0 .and. .not. intra .and. (abs(rho) < rhocut) .and. (dimgrad < dimcut) .and. &
                   (abs(rho) > 1d-30) .and. .not. (rmbox_coarse(i, j, k) )) then
                  write (ludat, '(1p,E18.10,E18.10)') rho, dimgrad
               endif ! rhocut/dimcut
 
               ! prepare the cube files
               if ((abs(rho) > rhoplot) .or. (dimgrad > dimcut) .or. (rmbox_coarse(i, j, k) )) then
                  cgrad(i, j, k) = 100d0
               endif !rho cutoff
               if  (intra) then ! intermolecular points also to 100
                  cgrad(i, j, k) = 100d0
               endif
            end do
         end do
      end do
      if (ludc > 0) call write_cube_body(ludc, nstep, crho)          ! density
      if (lugc > 0) call write_cube_body(lugc, nstep, cgrad)         ! RDG
      call system_clock(count=c4)
      write (*, "(A, F6.2, A)") ' Time for writing outputs = ', real(dble(c4 - c3)/dble(cr), kind=8), ' secs'
   end if

   !===============================================================================!
   ! Deallocate grids.
   !===============================================================================!
   if (allocated(tmp_rmbox)) deallocate (tmp_rmbox)
   if (allocated(rmbox_coarse)) deallocate (rmbox_coarse) 
   if (allocated(tmp_rmbox_range)) deallocate (tmp_rmbox_range)
   if (allocated(tmp_rmbox_range_tmp)) deallocate (tmp_rmbox_range_tmp)
   if (allocated(box_in_range)) deallocate (box_in_range)
   if (allocated(rho_range)) deallocate (rho_range)
   if (allocated(rmpoint_coarse)) deallocate (rmpoint_coarse)


   !===============================================================================!
   ! Deallocate grids, close files.
   !===============================================================================!
   if (allocated(crho)) deallocate (crho)
   if (allocated(cheig)) deallocate (cheig)
   if (allocated(cgrad)) deallocate (cgrad)
   if (allocated(crho_n)) deallocate (crho_n)
   if (ludat > 0) close (ludat)

   !===============================================================================!
   ! Write VMD script.
   !===============================================================================!
   if (ligand) then
      nn0 = sum(m(1:udat0 - 1)%n) + 1
      nnf = sum(m(1:udat0 - 1)%n) + m(udat0)%n
   else
      nn0 = 1
      nnf = ntotal
   end if

   if (luvmd > 0) then
      write (luvmd, 114) trim(oname)//"-dens.cube"
      write (luvmd, 115) trim(oname)//"-grad.cube"
      write (luvmd, 116) nn0 - 1, nnf - 1, isordg, 2, 2, 2, -rhoplot*100D0, rhoplot*100D0, 2, 2
      close (luvmd)
   end if

   !===============================================================================!
   ! captar error de no librerias para mencionar por terminal, no matar en caso
   ! NCICLUSTER python script integration
   !===============================================================================!
   if(doclustering) then 
       iounit_p1 = 501
       open(unit=iounit_p1, file='tmp_ncicluster_file', status='replace', action='write')
       ! Write the strings to the file
       write(iounit_p1, '(A)') trim(adjustl(oname))
       ! Close the file
       close(iounit_p1)

      write(command_ncicluster, '(A,A,A,A,F5.2,A,F5.3,A,F5.3)') trim(adjustl(nciplot_home)), "/dev/scripts/ncicluster.py", & 
         " tmp_ncicluster_file", & 
         " --doint True --method kmeans -n 2 --isovalue ", & 
         dimcut, &
         " --outer ", & 
         srhorange(3, 2), &          
         " --inner ", & 
         srhorange(3, 1)    

      py_status = system(trim(adjustl(command_ncicluster)))
      ! Check the return py_status
      select case (py_status)
      case (0)
           print *, "Python script executed successfully."
      case (256)
           print *, "Installation error: Required python library is not installed."
      case (32512)
           print *, "Python script execution failed, remeber to set NCIPLOT_HOME environment variable."
      case default
           print *, "Python script failed with unexpected status code: ", py_status
           print *, "Remember to install all necessary Python libraries"
      end select
      py_status=unlink("tmp_ncicluster_file")
   endif  ! clustering


   !===============================================================================!
   ! Deallocate arrays, call clock end, close output and input files.
   !===============================================================================!
   if (allocated(group)) deallocate (group)
   if (allocated(xmaxat)) deallocate (xmaxat)
   if (allocated(xinitat)) deallocate (xinitat)
   if (allocated(nstepat)) deallocate (nstepat)
   if (allocated(m)) deallocate (m)
   if (allocated(srhorange)) deallocate (srhorange)
   if (allocated(tmp_crho)) deallocate (tmp_crho)
   if (allocated(tmp_crho_n)) deallocate (tmp_crho_n)
   if (allocated(tmp_cheigs)) deallocate (tmp_cheigs)
   if (allocated(tmp_cgrad)) deallocate (tmp_cgrad)
   if (allocated(fginc)) deallocate (fginc)

   call tictac('End')
   if (uin /= stdin) close (uin)
   if (uout /= stdout) close (uout)

   !===============================================================================!
   ! Formats used by NCIPLOT.
   !===============================================================================!
110 format(A, F5.2)

   ! VMD script
114 format('#!/usr/local/bin/vmd', /, &
          '# VMD script written by save_state $Revision: 1.10 $', /, &
          '# VMD version: 1.8.6            ', /, &
          'set viewplist            ', /, &
          'set fixedlist            ', /, &
          '# Display settings            ', /, &
          'display projection   Orthographic            ', /, &
          'display nearclip set 0.000000            ', /, &
          '# load new molecule         ', /, &
          'mol new ', a, ' type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all')
115 format('mol addfile ', a, ' type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all')
116 format('#', /, &
          '# representation of the atoms', /, &
          'mol delrep 0 top', /, &
          'mol representation Lines 1.00000', /, &
          'mol color Name', /, &
          'mol selection {all}', /, &
          'mol material Opaque', /, &
          'mol addrep top', /, &
          'mol representation CPK 1.000000 0.300000 118.000000 131.000000', /, &
          'mol color Name', /, &
          'mol selection {index ', i5, ' to ', i5, ' }', /, &
          'mol material Opaque', /, &
          'mol addrep top', /, &
          '#', /, &
          '# add representation of the surface', /, &
          'mol representation Isosurface ', f8.5, ' 1 0 0 1 1', /, &
          'mol color Volume 0', /, &
          'mol selection {all}', /, &
          'mol material Opaque', /, &
          'mol addrep top', /, &
          'mol selupdate ', i1, ' top 0', /, &
          'mol colupdate ', i1, ' top 0', /, &
          'mol scaleminmax top ', i1, ' ', f8.4, ' ',f7.4, /, &
          'mol smoothrep top ', i1, ' 0', /, &
          'mol drawframes top ', i1, ' {now}', /, &
          'color scale method BGR', /, &
          'set colorcmds { {color Name {C} gray} }', /, &
          '#some more',/)
117 format('                                                     ', / &
          '----------------------------------------------------------------------', /, &
          '                                                    ', / &
          '                     INTEGRATION DATA                        ', /, &
          '                                                     ', / &
          '----------------------------------------------------------------------', /, &
          ' Integration  over the volumes of rho^n                               '/, &
          '----------------------------------------------------------------------', /, &
          ' n=1.0           :', 3X, F15.8, /, &
          ' n=1.5           :', 3X, F15.8, /, &
          ' n=2.0           :', 3X, F15.8, /, &
          ' n=2.5           :', 3X, F15.8, /, &
          ' n=3.0           :', 3X, F15.8, /, &
          ' n=4/3           :', 3X, F15.8, /, &
          ' n=5/3           :', 3X, F15.8, /, &
          ' Volume          :', 3X, F15.8, /, &
          ' rho-sum_i rho_i :', 3X, F15.8, / &
          '                  ', / &
          '---------------------------------------------------------------------', /, &
          ' Integration  over the volumes of sign(lambda2)(rho)^n             '/, &
          '---------------------------------------------------------------------', /, &
          ' n=1.0           :', 3X, F15.8, /, &
          ' n=1.5           :', 3X, F15.8, /, &
          ' n=2.0           :', 3X, F15.8, /, &
          ' n=2.5           :', 3X, F15.8, /, &
          ' n=3.0           :', 3X, F15.8, /, &
          ' n=4/3           :', 3X, F15.8, /, &
          ' n=5/3           :', 3X, F15.8, /, &
          '----------------------------------------------------------------------', /, &
          ' Integration  over the areas of rho^n              '/, &
          '----------------------------------------------------------------------', /, &
          ' n=1.0           :', 3X, F15.8, /, &
          ' n=1.5           :', 3X, F15.8, /, &
          ' n=2.0           :', 3X, F15.8, /, &
          ' n=2.5           :', 3X, F15.8, /, &
          ' n=3.0           :', 3X, F15.8, /, &
          ' n=4/3           :', 3X, F15.8, /, &
          ' n=5/3           :', 3X, F15.8, /, &
          ' Area            :', 3X, F15.8, /, &
          ' rho-sum_i rho_i :', 3X, F15.8, /, &
          '----------------------------------------------------------------------', /, &
          '                   ',/) 

118 format('                                                     ', / &
          '----------------------------------------------------------------------', /, &
          '                                                    ', / &
          '                   INTEGRATION DATA BY POINTS                       ', /, &
          '                                                     ', / &
          '----------------------------------------------------------------------', /, &
          ' Integration  over the volumes of rho^n                               '/, &
          '----------------------------------------------------------------------', /, &
          ' n=1.0           :', 3X, F15.8, /, &
          ' n=1.5           :', 3X, F15.8, /, &
          ' n=2.0           :', 3X, F15.8, /, &
          ' n=2.5           :', 3X, F15.8, /, &
          ' n=3.0           :', 3X, F15.8, /, &
          ' n=4/3           :', 3X, F15.8, /, &
          ' n=5/3           :', 3X, F15.8, /, &
          ' Volume          :', 3X, F15.8, /, &
          ' rho-sum_i rho_i :', 3X, F15.8, / &
          '                  ', / &
          '---------------------------------------------------------------------', /, &
          ' Integration  over the volumes of sign(lambda2)(rho)^n             '/, &
          '---------------------------------------------------------------------', /, &
          ' n=1.0           :', 3X, F15.8, /, &
          ' n=1.5           :', 3X, F15.8, /, &
          ' n=2.0           :', 3X, F15.8, /, &
          ' n=2.5           :', 3X, F15.8, /, &
          ' n=3.0           :', 3X, F15.8, /, &
          ' n=4/3           :', 3X, F15.8, /, &
          ' n=5/3           :', 3X, F15.8, /, &
          '----------------------------------------------------------------------', /, &
          '                   ',/)

120 format(/'-----------------------------------------------------'/ &
           '      Calculation details:'/ &
           '-----------------------------------------------------')
121 format(/, '-----------------------------------------------------'/ &
           '      Operating grid and increments: Grid-', I1/ &
           '-----------------------------------------------------'/ &
           ' x0,y0,z0  = ', f10.4, ' ', f10.4, ' ', f10.4/ &
           ' x1,y1,z1  = ', f10.4, ' ', f10.4, ' ', f10.4/ &
           ' ix,iy,iz  = ', f5.2, '   ', f5.2, '   ', f5.2/ &
           ' nx,ny,nz  = ', i4, '    ', i4, '    ', i4/)

122 format('-----------------------------------------------------'/ &
          '      Writing output in the following units:'/ &
          '-----------------------------------------------------'/)

   ! noutput=1 .or. noutput =3
123 format(' Sign(lambda2)xDensity x Reduced Density Gradient    = ', a,/)

   ! noutput >=2
124 format(' Reduced Density Gradient,RDG      = ', a, / &
          ' Sign(lambda2)xDensity,LS          = ', a, / &
          ' VMD script                        = ', a,/)

130 format('      Using ', a40, ' as LIGAND')
131 format('-----------------------------------------------------'/ &
          '      INPUT INFORMATION:'/ &
          '-----------------------------------------------------')
132 format(/'      MIND YOU'/ &
           '      ONLY ANALYZING INTERMOLECULAR INTERACTIONS     '/)

133 format(/'      MIND YOU'/ &
           '      RUNNING IN PROMOLECULAR MODE     '/)

134 format('                                                                      ', / &
          '----------------------------------------------------------------------', / &
          '                                                                      ', / &
          '               RANGE INTEGRATION DATA                                 ', /, &
          '----------------------------------------------------------------------', /, &
          '                                                                      ')

135 format('----------------------------------------------------------------------', /, &
          ' Interval        :', 2(3X, F15.8) '                     ', /, &
          '                               ', / &
          '----------------------------------------------------------------------', /, &
          ' Integration  over the volumes of rho^n                               '/, &
          '----------------------------------------------------------------------', /, &
          ' n=1.0           :', 3X, F15.8, /, &
          ' n=1.5           :', 3X, F15.8, /, &
          ' n=2.0           :', 3X, F15.8, /, &
          ' n=2.5           :', 3X, F15.8, /, &
          ' n=3.0           :', 3X, F15.8, /, &
          ' n=4/3           :', 3X, F15.8, /, &
          ' n=5/3           :', 3X, F15.8, /, &
          ' Volume          :', 3X, F15.8, /, &
          ' rho-sum_i rho_i :', 3X, F15.8, / &
          '                  ', / &
          '---------------------------------------------------------------------', /, &
          ' Integration  over the volumes of sign(lambda2)(rho)^n             '/, &
          '---------------------------------------------------------------------', /, &
          ' n=1.0           :', 3X, F15.8, /, &
          ' n=1.5           :', 3X, F15.8, /, &
          ' n=2.0           :', 3X, F15.8, /, &
          ' n=2.5           :', 3X, F15.8, /, &
          ' n=3.0           :', 3X, F15.8, /, &
          ' n=4/3           :', 3X, F15.8, /, &
          ' n=5/3           :', 3X, F15.8, /, &
          '---------------------------------------------------------------------', /, &
          '                         ',/) 


136 format(' Additivity Check', /, & 
          '----------------------------------------------------------------------', /, &
          ' Total Density          :', 3X, F15.8 '                     ', /, &
          ' Sum over range         :', 3X, F15.8 '                     ', /, &
          ' Density difference     :', 3X, F15.8 '                     ', /, &
          '----------------------------------------------------------------------', /, &
          '                         ',/) 

137 format('WARNING  -   Negative Density Values Found on Cube File               ', / &
           '             Absolute Values will be assigned                              ',/)


138 format('                                                     ', / &
          '-------------------------------------------------------------------', /, &
          ' CLUSTERING option specified                     ', /, &
          ' NciCluster will be used if all libraries are installed               ', /, &
          ' !!! Remember to set NCIPLOT_HOME variable                     ', /, &
          '-------------------------------------------------------------------', /, &
          '                         ',/) 
contains

   !===============================================================================!
   ! Subroutines for .cube, .mesh and .sol writing.
   !===============================================================================!

   subroutine write_cube_header(lu, l1, l2)

      integer, intent(in) :: lu
      character*(*), intent(in) :: l1, l2

      integer :: i, j

      write (lu, *) trim(l1)
      write (lu, *) trim(l2)
      write (lu, '(I5,3(F12.6))') ntotal, xinit

      write (lu, '(I5,3(F12.6))') nstep(1), xinc(1), 0d0, 0d0
      write (lu, '(I5,3(F12.6))') nstep(2), 0d0, xinc(2), 0d0
      write (lu, '(I5,3(F12.6))') nstep(3), 0d0, 0d0, xinc(3)
      do i = 1, nfiles
         do j = 1, m(i)%n
            write (lu, '(I4,F5.1,F11.6,F11.6,F11.6)') m(i)%z(j), 0d0, m(i)%x(:, j)
         end do
      enddo

   end subroutine write_cube_header

   subroutine write_cube_body(lu, n, c)

      integer, intent(in) :: lu
      integer, intent(in) :: n(3)
      real*8, intent(in) :: c(0:n(1) - 1, 0:n(2) - 1, 0:n(3) - 1)

      integer :: i, j

      do i = 0, n(1) - 1
         do j = 0, n(2) - 1
            write (lu, '(6(1x,e12.5))') (c(i, j, k), k=0, n(3) - 1)
         enddo
      enddo
      close (lu)

   end subroutine write_cube_body

   ! write the .mesh file
   subroutine write_mesh_file(lumesh, lusol, xinit, xinc, nstep, cgrad, xinc_coarse, rmbox_coarse, nstep_coarse, vert_use)
      integer, intent(in) :: lumesh, lusol
      real*8, intent(in) :: xinit(3), xinc(3), xinc_coarse(3)
      integer, intent(in) :: nstep(3), nstep_coarse(3)
      real*8, intent(in) :: cgrad(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1)
      logical, intent(in) :: rmbox_coarse(0:nstep_coarse(1) - 2, 0:nstep_coarse(2) - 2, 0:nstep_coarse(3) - 2)
      integer :: i, j, k, i0, j0, k0, m, count_vert, count_cube, c
      integer :: i1(2), j1(2), k1(2)
      integer :: tetra_cube(1:6, 1:4), ind_cube(1:8), indx(3)
      logical, intent(inout) :: vert_use(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1)
      integer :: ind_vert(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1)
      logical :: flag_use, ohexa, otetra
      ohexa = .false.
      otetra = .true.
      count_vert = 0
      count_cube = 0
      vert_use = .false.
      do i = 0, nstep(1) - 1
         do j = 0, nstep(2) - 1
            do k = 0, nstep(3) - 1
               ! eight neighbor boxes
               flag_use = .false.
               do i0 = max(0, i - 1), i
                  do j0 = max(0, j - 1), j
                     do k0 = max(0, k - 1), k
                        indx = floor(((/i0, j0, k0/)*xinc)/xinc_coarse)
                        indx = (/min(nstep_coarse(1) - 2, indx(1)), min(nstep_coarse(2) - 2, indx(2)), &
                                 min(nstep_coarse(3) - 2, indx(3))/)
                        if (.not. (rmbox_coarse(indx(1), indx(2), indx(3)))) then
                           vert_use(i, j, k) = .true.
                           flag_use = .true.
                           goto 30
                        end if
                     end do
                  end do
               end do
30             continue
               if (flag_use) then
                  count_vert = count_vert + 1
                  ind_vert(i, j, k) = count_vert
               end if
               indx = floor(((/i, j, k/)*xinc)/xinc_coarse)
               indx = (/min(nstep_coarse(1) - 2, indx(1)), min(nstep_coarse(2) - 2, indx(2)), &
                        min(nstep_coarse(3) - 2, indx(3))/)
               if (.not. (rmbox_coarse(indx(1), indx(2), indx(3)))) then
                  count_cube = count_cube + 1
               end if
            end do
         enddo
      enddo

      ! write data
      ! .mesh file
      write (lumesh, "(A)") 'MeshVersionFormatted 2'
      write (lumesh, "(A,/)") 'Dimension 3'
      write (lumesh, "(A)") 'Vertices'
      write (lumesh, '(I20)') count(vert_use)
      ! .sol file
      write (lusol, "(A)") 'MeshVersionFormatted 2'
      write (lusol, "(A,/)") 'Dimension 3'
      write (lusol, "(A)") 'SolAtVertices'
      write (lusol, '(I20)') count(vert_use)
      write (lusol, '(I16, I16)') 1, 1
      do i = 0, nstep(1) - 1
         do j = 0, nstep(2) - 1
            do k = 0, nstep(3) - 1
               if (vert_use(i, j, k)) then
                  write (lumesh, '(1x,3(e16.6), I8)') xinit + (/i, j, k/)*xinc, 1
                  write (lusol, '(1x, 1(f16.6))') cgrad(i, j, k)
               end if
            end do
         enddo
      enddo

      if (ohexa) then
         write (lumesh, "(/,A)") 'Hexaedra '
         write (lumesh, "(I20)") count_cube
         do i = 0, nstep(1) - 2
            do j = 0, nstep(2) - 2
               do k = 0, nstep(3) - 2
                  i1 = (/i, i + 1/)
                  j1 = (/j, j + 1/)
                  k1 = (/k, k + 1/)
                  c = count(vert_use(i1, j1, k1))
                  if (c .eq. 8) then
                     ind_cube = (/ind_vert(i, j, k), ind_vert(i, j, k + 1), ind_vert(i, j + 1, k + 1), ind_vert(i, j + 1, k), &
                        ind_vert(i + 1, j, k), ind_vert(i + 1, j, k + 1), ind_vert(i + 1, j + 1, k + 1), ind_vert(i + 1, j + 1, k)/)
                     write (lumesh, "(1x, 8(I12), I12)") ind_cube, 1
                  end if
               end do
            end do
         end do
      end if

      ! save tetras, by diving each cube into six tetras
      if (otetra) then
         write (lumesh, "(/,A)") 'Tetrahedra'
         write (lumesh, "(I20)") 6*count_cube
         do i = 0, nstep(1) - 2
            do j = 0, nstep(2) - 2
               do k = 0, nstep(3) - 2
                  i1 = (/i, i + 1/)
                  j1 = (/j, j + 1/)
                  k1 = (/k, k + 1/)
                  c = count(vert_use(i1, j1, k1))
                  if (c .eq. 8) then
                     ind_cube = (/ind_vert(i, j, k), ind_vert(i, j, k + 1), ind_vert(i, j + 1, k + 1), ind_vert(i, j + 1, k), &
                        ind_vert(i + 1, j, k), ind_vert(i + 1, j, k + 1), ind_vert(i + 1, j + 1, k + 1), ind_vert(i + 1, j + 1, k)/)
                     call cube2tetra(tetra_cube, ind_cube)
                     do m = 1, 6
                        write (lumesh, "(1x, 4(I16), I16)") tetra_cube(m, :), 1
                     end do
                  end if
               end do
            end do
         end do
      end if

      write (lumesh, '(A)') 'END'
      write (lusol, '(A)') 'END'
   end subroutine write_mesh_file

   !===============================================================================!
   ! Subroutines for multi-level grids.
   !===============================================================================!

   subroutine build_rmbox_coarse(rhocut, dimcut, ng, ind_g, fginc, rmbox_coarse, crho, cgrad, nstep)
      integer, intent(in) :: nstep(3), ng, ind_g
      real*8, intent(in) :: crho(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), &
                            cgrad(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1)
      real*8, intent(in) :: rhocut, dimcut
      real*8, intent(in) :: fginc(ng)
      real*8 :: rhocut0, dimcut0
      integer :: i, j, k
      integer :: i1(2), j1(2), k1(2)
      real*8 :: percent
      logical, intent(out) :: rmbox_coarse(0:nstep(1) - 2, 0:nstep(2) - 2, 0:nstep(3) - 2)
      logical :: flag_grad(2, 2, 2)
      rmbox_coarse = .true.
      rhocut0 = rhocut*fginc(ind_g)
      dimcut0 = dimcut*fginc(ind_g)

      do i = 0, nstep(1) - 2
         do j = 0, nstep(2) - 2
            do k = 0, nstep(3) - 2
               i1 = (/i, i + 1/)
               j1 = (/j, j + 1/)
               k1 = (/k, k + 1/)
               if ((ind_g .eq. ng) .and. (count(crho(i1, j1, k1) .gt. 99) .gt. 0)) then
                  ! small box with removed vertices
                  cycle
               end if

               flag_grad = (((abs(crho(i1, j1, k1))/100) .gt. rhocut0) .or. (cgrad(i1, j1, k1) .gt. dimcut0) .or. &
                            (cgrad(i1, j1, k1) .lt. 0))

               if (count(flag_grad) .lt. 8) then !original criterion
             !  if  (count(flag_grad) .lt.4) then !To be consistent with range integration. Gives irregular s=1 surfaces
                  rmbox_coarse(i, j, k) = .false.  ! active box
               endif

            end do
         end do
      end do
      percent = real(count(rmbox_coarse), kind=8)/(real(size(rmbox_coarse), kind=8))*100d0
      write (*, '(F6.2, A)') percent, '% of small boxes are removed.'
   end subroutine build_rmbox_coarse

   subroutine rmboxtopoint(nstep, rmbox_coarse, rmpoint_coarse)
      integer, intent(in) :: nstep(3)
      integer :: i, j, k, i1, j1, k1
      logical :: checked_point(0:nstep(1)-1, 0: nstep(2)-1, 0:nstep(3)-1)
      logical, intent(in) :: rmbox_coarse(0:nstep(1) - 2, 0:nstep(2) - 2, 0:nstep(3) - 2)
      logical, intent(inout) :: rmpoint_coarse(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1)
      
      !rmpoint_coarse must include the indexes of the points i1,j1,k1, and the corresponding box i,j,k
      !something like rmpoint_coarse(i1,j1,k1,i,j,k)
      !Otherwise shared vertices are include once in the rmpoint_coarse and the number of points and the number 
      !of boxes
      rmpoint_coarse = .true.
      checked_point = .false.
  
      do i = 0, nstep(1) - 2
         do j = 0, nstep(2) - 2
            do k = 0, nstep(3) - 2
               if (.not. rmbox_coarse(i,j,k)) then
                   do i1 = i, i+1
                      do j1 = j, j+1
                         do k1 = k, k+1
                        !    if (.not. checked_point(i1,j1,k1)) then 
                                rmpoint_coarse(i1,j1,k1) = .false.
                                checked_point(i1,j1,k1)  = .true. 
                        !    endif       
                         end do
                      end do  
                   end do
               endif
            end do
         end do
      end do   
      
   end subroutine rmboxtopoint

! write six faces of a cube
   subroutine cube2quad(quad_cube, ind_cube)
      integer, intent(in) :: ind_cube(1:8)
      integer, intent(inout) :: quad_cube(1:6, 1:4)
      integer :: i1, i2, i3, i4, i5, i6, i7, i8
      i1 = ind_cube(1)
      i2 = ind_cube(2)
      i3 = ind_cube(3)
      i4 = ind_cube(4)
      i5 = ind_cube(5)
      i6 = ind_cube(6)
      i7 = ind_cube(7)
      i8 = ind_cube(8)
      quad_cube = reshape((/i1, i2, i3, i4, &
                            i5, i6, i7, i8, &
                            i1, i2, i6, i5, &
                            i4, i3, i7, i8, &
                            i1, i5, i8, i4, &
                            i2, i6, i7, i3/), (/6, 4/), order=(/2, 1/))
   end subroutine cube2quad

   ! divide a cube into six triangles
   subroutine cube2tetra(tetra_cube, ind_cube)
      integer, intent(in) :: ind_cube(1:8)
      integer, intent(inout) :: tetra_cube(1:6, 1:4)
      integer :: i1, i2, i3, i4, i5, i6, i7, i8
      i1 = ind_cube(1)
      i2 = ind_cube(2)
      i3 = ind_cube(3)
      i4 = ind_cube(4)
      i5 = ind_cube(5)
      i6 = ind_cube(6)
      i7 = ind_cube(7)
      i8 = ind_cube(8)
      tetra_cube = reshape((/i1, i2, i3, i7, &
                             i1, i2, i6, i7, &
                             i1, i4, i3, i7, &
                             i1, i4, i8, i7, &
                             i1, i5, i6, i7, &
                             i1, i5, i8, i7/), (/6, 4/), order=(/2, 1/))
   end subroutine cube2tetra

! compute NCI geometry of the region enclosed by RDG isosurface
   subroutine dataGeom(sum_rhon_vol, sum_rhon_area, sum_signrhon_vol, xinc, nstep, crho, crho_n, rmbox_coarse, nfiles)
      real*8, intent(inout) :: sum_rhon_vol(9)
      real*8, intent(inout) :: sum_rhon_area(9)
      real*8, intent(inout) :: sum_signrhon_vol(7)
      real*8, intent(in) :: xinc(3)
      integer, intent(in) :: nstep(3)
      real*8, intent(in) :: crho(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), &
                            crho_n(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1, 1:nfiles)
      logical, intent(in) :: rmbox_coarse(0:nstep(1) - 2, 0:nstep(2) - 2, 0:nstep(3) - 2)
      integer :: i, j, k, n
      integer, intent(in) :: nfiles
      integer :: i1(2), j1(2), k1(2), negative, positive
      real*8 :: sum_signrho, signlambda_2

      sum_rhon_vol = 0
      sum_signrhon_vol = 0
      negative = 0
      positive = 0
      sum_signrho=0
      signlambda_2 = 0
      ! integral of rho^n over the volume of cubes
      do i = 0, nstep(1) - 2
         do j = 0, nstep(2) - 2
            do k = 0, nstep(3) - 2
               if (.not. rmbox_coarse(i, j, k)) then
                  i1 = (/i, i + 1/)
                  j1 = (/j, j + 1/)
                  k1 = (/k, k + 1/)
                  ! Check to avoid passing several times over the same points
                  ! Must be done by points
                  ! if (.not.(Integrated(ii,jj,kk))
                  !
                  ! n = 1, 1.5, 2, 2.5, 3, 4/3, 5/3: sum of rho^n
                  sum_rhon_vol(1) = sum_rhon_vol(1) + sum(abs(crho(i1, j1, k1)/100))*xinc(1)*xinc(2)*xinc(3)/8
                  sum_rhon_vol(2) = sum_rhon_vol(2) + sum(abs(crho(i1, j1, k1)/100)**1.5) &
                                    *xinc(1)*xinc(2)*xinc(3)/8
                  sum_rhon_vol(3) = sum_rhon_vol(3) + sum(abs(crho(i1, j1, k1)/100)**2) &
                                    *xinc(1)*xinc(2)*xinc(3)/8
                  sum_rhon_vol(4) = sum_rhon_vol(4) + sum(abs(crho(i1, j1, k1)/ &
                                                              100)**2.5)*xinc(1)*xinc(2)*xinc(3)/8
                  sum_rhon_vol(5) = sum_rhon_vol(5) + sum(abs(crho(i1, j1, k1)/100)**3) &
                                    *xinc(1)*xinc(2)*xinc(3)/8
                  sum_rhon_vol(6) = sum_rhon_vol(6) + sum(abs(crho(i1, j1, k1)/100)**(1.333)) &
                                    *xinc(1)*xinc(2)*xinc(3)/8
                  sum_rhon_vol(7) = sum_rhon_vol(7) + sum(abs(crho(i1, j1, k1)/100)**(1.666)) &
                                    *xinc(1)*xinc(2)*xinc(3)/8
                  ! n = 0: volume of cubes
                  sum_rhon_vol(8) = sum_rhon_vol(8) + xinc(1)*xinc(2)*xinc(3)
                  ! sum of rho1*rho2
                  do n = 1, nfiles
                     sum_rhon_vol(9) = sum_rhon_vol(9) + sum(crho_n(i1, j1, k1, n)/100* &
                                                      (abs(crho(i1, j1, k1)) - crho_n(i1, j1, k1, n))/100)*xinc(1)*xinc(2)*xinc(3)/8
                  end do

                  !sign_lambda2 x rho
                  sum_signrho = sum(crho(i1, j1, k1))
                  signlambda_2 = sign(1d0, sum_signrho)

                  sum_signrhon_vol(1) = sum_signrhon_vol(1) + signlambda_2*sum(abs(crho(i1, j1, k1)/100))*xinc(1)*xinc(2)*xinc(3)/8
                  sum_signrhon_vol(2) = sum_signrhon_vol(2) + signlambda_2*sum(abs(crho(i1, j1, k1)/100)**1.5) &
                                        *xinc(1)*xinc(2)*xinc(3)/8
                  sum_signrhon_vol(3) = sum_signrhon_vol(3) + signlambda_2*sum(abs(crho(i1, j1, k1)/100)**2) &
                                        *xinc(1)*xinc(2)*xinc(3)/8
                  sum_signrhon_vol(4) = sum_signrhon_vol(4) + signlambda_2*sum(abs(crho(i1, j1, k1)/ &
                                                                                   100)**2.5)*xinc(1)*xinc(2)*xinc(3)/8
                  sum_signrhon_vol(5) = sum_signrhon_vol(5) + signlambda_2*sum(abs(crho(i1, j1, k1)/100)**3) &
                                        *xinc(1)*xinc(2)*xinc(3)/8
                  sum_signrhon_vol(6) = sum_signrhon_vol(6) + signlambda_2*sum(abs(crho(i1, j1, k1)/100)**(1.333)) &
                                        *xinc(1)*xinc(2)*xinc(3)/8
                  sum_signrhon_vol(7) = sum_signrhon_vol(7) + signlambda_2*sum(abs(crho(i1, j1, k1)/100)**(1.666)) &
                                        *xinc(1)*xinc(2)*xinc(3)/8
               end if
            end do
         end do
      end do

      ! integral of rho^n over the surface of cubes
      sum_rhon_area = 0
      do i = 0, nstep(1) - 2
         do j = 0, nstep(2) - 2
            do k = 0, nstep(3) - 2
               if (.not. rmbox_coarse(i, j, k)) then
                  if (i .eq. 0) then
                     i1 = (/i, i/)
                     j1 = (/j, j + 1/)
                     k1 = (/k, k + 1/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(2), xinc(3))
                  elseif (rmbox_coarse(i - 1, j, k)) then
                     i1 = (/i, i/)
                     j1 = (/j, j + 1/)
                     k1 = (/k, k + 1/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(2), xinc(3))
                  elseif ((i .eq. nstep(1) - 2) .or. (rmbox_coarse(i, j, k))) then
                     i1 = (/i + 1, i + 1/)
                     j1 = (/j, j + 1/)
                     k1 = (/k, k + 1/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(2), xinc(3))
                  end if
                  if (j .eq. 0) then
                     i1 = (/i, i + 1/)
                     j1 = (/j, j/)
                     k1 = (/k, k + 1/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(1), xinc(3))
                  elseif (rmbox_coarse(i, j - 1, k)) then
                     i1 = (/i, i + 1/)
                     j1 = (/j, j/)
                     k1 = (/k, k + 1/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(1), xinc(3))
                  elseif ((j .eq. nstep(2) - 2) .or. (rmbox_coarse(i, j, k))) then
                     i1 = (/i, i + 1/)
                     j1 = (/j + 1, j + 1/)
                     k1 = (/k, k + 1/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(1), xinc(3))
                  end if
                  if (k .eq. 0) then
                     i1 = (/i, i + 1/)
                     j1 = (/j, j + 1/)
                     k1 = (/k, k/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(1), xinc(3))
                  elseif (rmbox_coarse(i, j, k - 1)) then
                     i1 = (/i, i + 1/)
                     j1 = (/j, j + 1/)
                     k1 = (/k, k/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(1), xinc(2))
                  elseif ((k .eq. nstep(3) - 2) .or. (rmbox_coarse(i, j, k))) then
                     i1 = (/i, i + 1/)
                     j1 = (/j, j + 1/)
                     k1 = (/k + 1, k + 1/)
                     call compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, xinc(1), xinc(2))
                  end if
               end if
            end do
         end do
      end do

   end subroutine dataGeom 

   subroutine dataGeom_points(sum_rhon_vol, sum_signrhon_vol, xinc, nstep, crho, crho_n, rmbox_coarse,rmpoint_coarse, nfiles)
      real*8, intent(out) :: sum_rhon_vol(9)
      real*8, intent(out) :: sum_signrhon_vol(7)
      real*8, intent(in) :: xinc(3)
      integer, intent(in) :: nstep(3)
      real*8, intent(in) :: crho(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), &
                            crho_n(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1, 1:nfiles)
      logical, intent(in) :: rmbox_coarse(0:nstep(1) - 2, 0:nstep(2) - 2, 0:nstep(3) - 2)
      logical, intent(in) :: rmpoint_coarse(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1)
      integer :: i, j, k, n, i1, j1, k1
      integer, intent(in) :: nfiles
      integer :: negative, positive
      real*8 :: signrho, signlambda_2
      
      sum_rhon_vol = 0
      sum_signrhon_vol = 0
      negative = 0
      positive = 0
      ! integral of rho^n over the points
      do i = 0, nstep(1) - 2
         do j = 0, nstep(2) - 2
            do k = 0, nstep(3) - 2 
               if (.not. rmbox_coarse(i, j, k)) then
                   do i1 = i, i+1
                      do j1 = j, j+1
                         do k1 = k, k+1
                          if (.not. rmpoint_coarse(i1,j1,k1)) then
                            ! n = 1, 1.5, 2, 2.5, 3, 4/3, 5/3: sum of rho^n
                            sum_rhon_vol(1) = sum_rhon_vol(1) + abs(crho(i1, j1, k1)/100)*xinc(1)*xinc(2)*xinc(3)/8
                            sum_rhon_vol(2) = sum_rhon_vol(2) + abs(crho(i1, j1, k1)/100)**1.5 &
                                    *xinc(1)*xinc(2)*xinc(3)/8
                            sum_rhon_vol(3) = sum_rhon_vol(3) + abs(crho(i1, j1, k1)/100)**2 &
                                    *xinc(1)*xinc(2)*xinc(3)/8
                            sum_rhon_vol(4) = sum_rhon_vol(4) + abs(crho(i1, j1, k1)/ &
                                                              100)**2.5*xinc(1)*xinc(2)*xinc(3)/8
                            sum_rhon_vol(5) = sum_rhon_vol(5) + abs(crho(i1, j1, k1)/100)**3 &
                                    *xinc(1)*xinc(2)*xinc(3)/8
                            sum_rhon_vol(6) = sum_rhon_vol(6) + abs(crho(i1, j1, k1)/100)**(1.333) &
                                    *xinc(1)*xinc(2)*xinc(3)/8
                            sum_rhon_vol(7) = sum_rhon_vol(7) + abs(crho(i1, j1, k1)/100)**(1.666) &
                                    *xinc(1)*xinc(2)*xinc(3)/8
                            ! n = 0: volume of cubes
                            sum_rhon_vol(8) = sum_rhon_vol(8) + xinc(1)*xinc(2)*xinc(3)/8

                            !sign_lambda2 x rho
                            signrho = crho(i, j, k)
                            signlambda_2 = sign(1d0, signrho)

                            sum_signrhon_vol(1) = sum_signrhon_vol(1) + (crho(i1, j1, k1)/100)*xinc(1)*xinc(2)*xinc(3)/8
                            sum_signrhon_vol(2) = sum_signrhon_vol(2) + signlambda_2*(abs(crho(i1, j1, k1)/100)**1.5) &
                                        *xinc(1)*xinc(2)*xinc(3)/8
                            sum_signrhon_vol(3) = sum_signrhon_vol(3) + signlambda_2*(abs(crho(i1, j1, k1)/100)**2) &
                                        *xinc(1)*xinc(2)*xinc(3)/8
                            sum_signrhon_vol(4) = sum_signrhon_vol(4) + signlambda_2*(abs(crho(i1, j1, k1)/ &
                                                                                   100)**2.5)*xinc(1)*xinc(2)*xinc(3)/8
                            sum_signrhon_vol(5) = sum_signrhon_vol(5) + signlambda_2*(abs(crho(i1, j1, k1)/100)**3) &
                                        *xinc(1)*xinc(2)*xinc(3)/8
                            sum_signrhon_vol(6) = sum_signrhon_vol(6) + signlambda_2*(abs(crho(i1, j1, k1)/100)**(1.333)) &
                                        *xinc(1)*xinc(2)*xinc(3)/8
                            sum_signrhon_vol(7) = sum_signrhon_vol(7) + signlambda_2*(abs(crho(i1, j1, k1)/100)**(1.666)) &
                                        *xinc(1)*xinc(2)*xinc(3)/8
                         end if    
                         end do 
                     end do 
                   end do  
               end if
            end do
         end do
      end do
   

   end subroutine dataGeom_points

   subroutine compArea(sum_rhon_area, nstep, crho, crho_n, nfiles, i1, j1, k1, a, b)
      ! compute integrals over the boundary surface of the union of cubes
      real*8, intent(inout) :: sum_rhon_area(9)
      real*8, intent(in) :: a, b
      integer, intent(in) :: nstep(3)
      real*8, intent(in) :: crho(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1), &
                            crho_n(0:nstep(1) - 1, 0:nstep(2) - 1, 0:nstep(3) - 1, 1:nfiles)
      integer :: n
      integer, intent(in) :: nfiles
      integer, intent(in) :: i1(2), j1(2), k1(2)
      ! face sides (a, b)
      ! n = 1, 1.5, 2, 2.5, 3, 4/3, 5/3: sum of rho^n
      sum_rhon_area(1) = sum_rhon_area(1) + sum(abs(crho(i1, j1, k1)/100))*a*b/8
      sum_rhon_area(2) = sum_rhon_area(2) + sum(abs(crho(i1, j1, k1)/100)**1.5)*a*b/8
      sum_rhon_area(3) = sum_rhon_area(3) + sum(abs(crho(i1, j1, k1)/100)**2)*a*b/8
      sum_rhon_area(4) = sum_rhon_area(4) + sum(abs(crho(i1, j1, k1)/100)**2.5)*a*b/8
      sum_rhon_area(5) = sum_rhon_area(5) + sum(abs(crho(i1, j1, k1)/100)**3)*a*b/8
      sum_rhon_area(6) = sum_rhon_area(6) + sum(abs(crho(i1, j1, k1)/100)**(1.3333333)) &
                         *a*b/8
      sum_rhon_area(7) = sum_rhon_area(7) + sum(abs(crho(i1, j1, k1)/100)**(1.6666666)) &
                         *a*b/8
      ! area of iso-surface
      sum_rhon_area(8) = sum_rhon_area(8) + a*b
      ! sum of rho1*rho2
      do n = 1, nfiles
         sum_rhon_area(9) = sum_rhon_area(9) + sum(crho_n(i1, j1, k1, n)/100* &
                                                   (abs(crho(i1, j1, k1)) - crho_n(i1, j1, k1, n))/100)*a*b/8
      end do
   end subroutine compArea
 
end program
