c=========================================================================c
      program SAMBVCA
c=========================================================================c
c                                                                         c
c SambVca: Calculation of the Buried Volume of Organometallic Ligands     c
c                                                                         c
c                  RELEASE 2.1 20/07/2019                                 c
c                                                                         c
c Copyright (C) 2008 - 2011 Luigi Cavallo, University of Salerno, Italy   c
c Copyright (C) 2012 - 2019 Luigi Cavallo, KAUST, Saudi Arabia            c
c                                                                         c
c Publications using this tool should cite :                              c
c Poater, A.; Cosenza, B.; Correa, A.; Giudice, S.; Ragone, F.; Scarano,  c
c V.; Cavallo, L. SambVca: A Web Application for the Calculation of the   c
c Buried Volume of c N-Heterocyclic Carbene Ligands                       c
c Eur. J. Inorg. Chem. 2009, 1759                                         c
c                                                                         c
c For any questions/suggestions regarding this code, please contact       c
c Prof. Luigi Cavallo : luigi.cavallo@kaust.edu.sa; lcavallo@unisa.it     c
c                                                                         c
c This program is free software; you can redistribute it and/or modify    c
c it under the terms of the GNU General Public License as published by    c
c the Free Software Foundation; either version 1, or (at your option)     c
c any later version.                                                      c
c                                                                         c
c This program is distributed WITHOUT ANY WARRANTY; without even          c
c the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR     c
c PURPOSE.  See the c GNU General Public License for more details.        c
c                                                                         c
c For the GNU General Public License, see http://www.gnu.org/licenses/    c
c                                                                         c
c Work using this software should cite:                                   c
c Falivene L. et a.. Nat. Chem. 2019, DOI:10.1038/s41557-019-0319-5       c
c                                                                         c
c=========================================================================c

      IMPLICIT NONE

c=========================================================================
c Constant Parameters
c=========================================================================
      Integer, parameter :: NatomMax = 100000   ! max number of atoms
      Integer, parameter :: NpsMax = 10000      ! max number of points for distribution
      Real*8, parameter :: PI = 3.1415926535    ! Pi

c=========================================================================
c Variables to define the system
c=========================================================================
      Integer :: natoms                         ! # of atoms
      Integer :: natoms_new                     ! # of atoms after removing the unnecessary atoms
      Integer :: nTypes                         ! # of atom types

      Real*8, Allocatable :: AtomTypeRadius(:)  ! radius of atom types
      Real*8, Allocatable :: Coord(:)           ! coordinates from input
      Real*8, Allocatable :: Coord_new(:)       ! coordinates after removing unnecessary atoms
      Real*8, Allocatable :: AtomRadius2(:)     ! square of atom radii
      Real*8, Allocatable :: AtomRadius2_new(:) ! square of atom radii after removing unnecessary atoms

      Character*50 :: Title                          ! Frame title, line after natoms line from input
      Character*4, Allocatable :: AtomNames(:)       ! atom names from input
      Character*4, Allocatable :: AtomNames_new(:)   ! record temperary atom names when removing atoms
      Character*4, Allocatable :: AtomTypeName(:)    ! atom names in the data base

c=========================================================================
c Variables to reorient the system
c=========================================================================
      Real*8 :: GC(3)                           ! geometry center
      Real*8 :: x_axis(3)                       ! unit vector for x axis
      Real*8 :: y_axis(3)                       ! unit vector for y axis
      Real*8 :: z_axis(3)                       ! unit vector for z axis
      Real*8 :: z_GC(3)                         ! geometry center of atoms defining z-axis
      Real*8 :: x_temp_GC(3)                    ! geometry center of atoms 
      Real*8 :: bench_vec(3)                    ! align to bench_vec
      Real*8 :: Initial(3)                      ! initial z

c=========================================================================
c Variables to control the calculation
c=========================================================================
      Integer :: RemoveH                        ! if RemoveH = 0 do not remove H atoms. If = 1, remove H atoms
      Integer :: DoMap                          ! if DoMap: =0 do not print surfaces. If = 1 print surfaces
      Integer :: AlignZ                         ! if AlignZ =0, align the molecule along -z; if =1 along +z
      Integer :: IfOverlap                      ! if IfOverlap = 0 voxel is free. If = 1 voxel is buried
      Integer :: IndFrame                       ! index of the molecular frame currently calculated
      Integer :: NumPoints                      ! # of grid points to be scanned: radius / binsize * 2 + 1
      Integer :: NumInside                      ! # of points inside the sphere
      Integer :: nAtomDel                       ! # of atoms to be deleted
      Integer :: nAtomGC                        ! # of atoms to define a geometry center
      Integer :: nAtom_z                        ! # of atoms to define z axis
      Integer :: nAtom_xz                       ! # of atoms to define a second vector

      Integer, allocatable :: IndDel(:)         ! index of the atoms to be deleted
      Integer, allocatable :: IndGC(:)          ! index of the atoms to define a geometry center
      Integer, Allocatable :: IndAxis(:)        ! index of the atoms to define an axis
      Integer, Allocatable :: Indxz(:)          ! index of these atoms

      Integer :: npoints                        ! # of grid points for Vbur distributions
      Integer :: VburDist(NpsMax)               ! distribution of the Vbur among different frames
      Integer :: Vbur4Dist(4, NpsMax)           ! distribution of the Vbur4 among different frames
      Integer :: Vbur8Dist(8, NpsMax)           ! distribution of the Vbur8 among different frames
      Integer :: IfWeb                          ! Hardcoded: If IfWeb = 1 print section headers for the web server 

      Real*8 :: Wgth                            ! weigth of point, 1.0 inside sphere, 0.5 at surface
      Real*8 :: Radius                          ! radius of sphere
      Real*8 :: Radius2                         ! radius of sphere squared
      Real*8 :: Disp                            ! displacement from sphere center
      Real*8 :: BinSize                         ! mesh point length: resolution
      Real*8 :: BinSize2                        ! mesh point length squared
      Real*8 :: BinSize3                        ! mesh point length cubic = volume of voxel
      Real*8 :: Cutoff                          ! cutoff to remove atoms far from the sphere
      Real*8 :: Vburf                           ! free volume
      Real*8 :: Vburo                           ! overlapped volume
      Real*8 :: VburTot                         ! total volume
      Real*8 :: zmin                            ! bottom z surface
      Real*8 :: zmax                            ! top z surface
      Real*8 :: Volume                          ! volume of the sphere
      Real*8 :: space                           ! grid point size 
      Real*8 :: notch                           ! starting point of the plot
      Real*8 :: smin                            ! starting point of one grid
      Real*8 :: smax                            ! ending point of one grid

      Real*8 :: Vbur4f(4)                       ! project free volume into 4 sphere quadrants
      Real*8 :: Vbur8f(8)                       ! project free volume into 8 sphere octants
      Real*8 :: Vbur4o(4)                       ! project occupied volume into 4 sphere quadrants
      Real*8 :: Vbur8o(8)                       ! project occupied volume into 8 sphere octants

c=========================================================================
c Input and output files
c=========================================================================
      Character*4   :: suffix_inp               ! suffix for file : .inp
      Character*4   :: suffix_out               ! suffix for file : .out
      Character*4   :: suffix_dat               ! suffix for file : .dat
      Character*4   :: suffix_xyz               ! suffix for file : .xyz
      Character*100 :: BaseName                 ! the BaseName identifying the job files
      Character*100 :: VburINP                  ! input file for Vbur
      Character*100 :: VburOUT                  ! output file 
      Character*100 :: VburDAT                  ! output file for Vbur values
      Character*100 :: VburDistOUT              ! output file for Vbur_distribution
      Character*100 :: RotatedXYZ               ! output file for rotated molecule
      Character*100 :: Quadrant4                ! project the Vbur into 4 sphere quadrants
      Character*100 :: Quadrant4Dist            ! distribution within 4 sphere quadrants
      Character*100 :: Quadrant8                ! project the Vbur into 8 sphere octants
      Character*100 :: Quadrant8Dist            ! distribution within 8 sphere octants
      Character*100 :: CheckCoord               ! check if the coordinates are successfully manipulated
      Character*100 :: TopSurface               ! top surface 
      Character*100 :: BotSurface               ! bottmon surface 

c=========================================================================
c simple index and temporary variables
c=========================================================================
      Integer :: i, j, k, l, m, n

      Real*8 :: x, y, z, dx, dy, dz
      Real*8 :: dist, dist2
      Real*8 :: temp
      Character*4 :: ATOMNAME
      Character*128 :: WD

c=========================================================================
c setting input and output files
c=========================================================================

      call getarg(1, BaseName)

      if(iargc().eq.0) then
        write(6,*) " "
        write(6,*) "Usage:  "
        write(6,*) "Assuming a myfile.inp input file  "
        write(6,*) "./sambvca21.x myfile "
        stop
      endif

      suffix_inp = ".inp"
      suffix_dat = ".dat"
      suffix_xyz = ".xyz"
      suffix_out = ".out"

c====================================================================================
c Line below  is  a hard coded parameter defining directories and output style.  
c IfWeb=1 triggers options for the webserver. IfWeb=0 to be used in line mode
c====================================================================================
c
      IfWeb=0

      if(IfWeb.eq.0) then
        WD=TRIM(BaseName)
      else
        WD='./temp/'//TRIM(BaseName)//"/"//TRIM(BaseName)
      endif

c====================================================================================
c End defining Web or command line output directory
c====================================================================================

      VburINP        = TRIM(WD)//TRIM(suffix_inp)
      VburOUT        = TRIM(WD)//TRIM(suffix_out)
      VburDAT        = TRIM(WD)//"-VBURvsFRAME"//TRIM(suffix_dat)
      VburDistOUT    = TRIM(WD)//"-VBUR-distrib"//TRIM(suffix_dat)
      RotatedXYZ     = TRIM(WD)//"-rotated"//TRIM(suffix_xyz)
      Quadrant4      = TRIM(WD)//"-QUADvsFRAME"//TRIM(suffix_dat)
      Quadrant4Dist  = TRIM(WD)//"-QUAD-distrib"//TRIM(suffix_dat)
      Quadrant8      = TRIM(WD)//"-OCTvsFRAME"//TRIM(suffix_dat)
      Quadrant8Dist  = TRIM(WD)//"-OCT-distrib"//TRIM(suffix_dat)
      TopSurface     = TRIM(WD)//"-TopSurface"//TRIM(suffix_dat)
      BotSurface     = TRIM(WD)//"-BotSurface"//TRIM(suffix_dat)

      open(101, file=VburINP, status='old')
      open(200, file=VburOUT, status='unknown')
      open(201, file=VburDAT, status='unknown')
      open(202, file=VburDistOUT, status='unknown')
      open(203, file=Quadrant4, status='unknown')
      open(204, file=Quadrant4Dist, status='unknown')
      open(205, file=Quadrant8, status='unknown')
      open(206, file=Quadrant8Dist, status='unknown')
      open(207, file=RotatedXYZ, status='unknown')
      open(208, file=TopSurface, action="write", status="unknown")
      open(209, file=BotSurface, action="write", status="unknown")

      CALL PrintFlag()
c=========================================================================
c Reading the input file
c=========================================================================

c atoms to be deleted 
      read(101,*)nAtomDel
      if(nAtomDel.gt.0)then
        Allocate (IndDel(nAtomDel))
        read(101,*)IndDel(1:nAtomDel)
      endif

c atoms defining the geometry center
      read(101,*)nAtomGC
      Allocate (IndGC(nAtomGC))
      read(101,*)IndGC(1:nAtomGC)

c atoms defining the z-axis
      read(101,*)nAtom_z
      Allocate (IndAxis(nAtom_z))
      read(101,*)IndAxis(1:nAtom_z)

c atoms defining the xz plane
      read(101,*)nAtom_xz
      Allocate (Indxz(nAtom_xz))
      read(101,*)Indxz(1:nAtom_xz)

c radius of the sphere 
      read(101,*)Radius

c displacement of geometry center from the center of sphere
      read(101,*)Disp

c mesh point length
      read(101,*)BinSize

c if remove Hydrogen atoms
      read(101,*)RemoveH

c align molecule to +z or -z
      read(101,*)AlignZ

c if make map
      read(101,*)DoMap

c read atom types and radius
      read(101,*)nTypes

      Allocate (AtomTypeRadius(nTypes))
      Allocate (AtomTypeName(nTypes))

      do i = 1, nTypes
        read(101,*)AtomTypeName(i),AtomTypeRadius(i)
      enddo

      CALL CheckControlParameters(nAtomDel, IndDel, 
     x     nAtomGC, IndGC, nAtom_z, IndAxis, nAtom_xz, Indxz, 
     x     Radius, Disp, BinSize, RemoveH, DoMap, 
     x     AlignZ, nTypes, AtomTypeName, AtomTypeRadius)

c==========================================================================
c Start the big loop over the number of coordinate frames in the input file.
c==========================================================================
      IndFrame = 0
      do while(.true.)

c read coordinates of atoms composing the molecule to be analyzed
        read(101,*,end=901)natoms
        IndFrame = IndFrame+1
        read(101,*,end=901)Title
        Allocate (Coord(natoms*3))
        Allocate (AtomNames(natoms))
        Allocate (AtomRadius2(natoms))
        do i = 1, natoms
          j = (i-1)*3
          read(101,*,end=901)AtomNames(i),
     x    Coord(j+1),Coord(j+2),Coord(j+3)
        enddo
 
c write input coordinates in output:
        write(200,*)
        write(200,'(a39,16x,i5)')
     x  "Checking: No. of atoms of input frame :", natoms

        write(200,*)
        write(200,'(a38)')"Checking: Coordinates of input frame :"
        do i = 1, natoms
          j = (i-1)*3
          write(200,'(a4,1x,3f10.5)')AtomNames(i),
     x    Coord(j+1),Coord(j+2),Coord(j+3)
        enddo

c=========================================================================
c defining the cartesian axis
c=========================================================================
        Call ConstructXYZ(GC, natoms, Coord, IndGC, nAtomGC, z_GC,
     x                    nAtom_z, IndAxis, nAtom_xz, Indxz,
     x                    x_axis, y_axis, z_axis, x_temp_GC)

c align z_axis to 0 0 1
        bench_vec(1) = 0.d0
        bench_vec(2) = 0.d0
        if(AlignZ.eq.0)bench_vec(3) = -1.d0
        if(AlignZ.eq.1)bench_vec(3) = +1.d0
        initial(1) = z_axis(1)
        initial(2) = z_axis(2)
        initial(3) = z_axis(3)

        call ali(Coord, bench_vec, initial,
     x           x_axis, y_axis, z_axis, GC, natoms)

c renew your axis
        Call ConstructXYZ(GC, natoms, Coord, IndGC, nAtomGC, z_GC,
     x                    nAtom_z, IndAxis, nAtom_xz, Indxz,
     x                    x_axis, y_axis, z_axis, x_temp_GC)

c align x_axis to 1 0 0 with fixed z
        bench_vec(1) = 1.d0
        bench_vec(2) = 0.d0
        bench_vec(3) = 0.d0
        initial(1) = x_axis(1)
        initial(2) = x_axis(2)
        initial(3) = 0.d0
        call ali(Coord, bench_vec, initial,
     x           x_axis, y_axis, z_axis, GC, natoms)

c renew your GC
        Call ConstructXYZ(GC, natoms, Coord, IndGC, nAtomGC, z_GC,
     x                    nAtom_z, IndAxis, nAtom_xz, Indxz,
     x                    x_axis, y_axis, z_axis, x_temp_GC)

        Call MoveParticle(Coord, GC, natoms)

c move the particle if Disp > 0.0
        if(Disp.ne.0.0)then
          GC(1) = 0.d0
          GC(2) = 0.d0
          if(AlignZ.eq.0)GC(3) = +Disp
          if(AlignZ.eq.1)GC(3) = -Disp
          Call MoveParticle(Coord, GC, natoms)
        endif

c write XYZ rotated frame in output and -rotated.xyz file
        write(200,*)
        write(200,'(a39)')"Results: Coordinates of rotated frame :"
        call WriteXYZ(200,natoms,Title,Atomnames,Coord)
        call WriteXYZ(207,natoms,Title,Atomnames,Coord)

c remove deleted atoms "AS" and H atoms,  if requested, by working array
        if (nAtomDel.gt.0)then
          do i = 1, nAtomDel
            j = IndDel(i)
            AtomNames(j) = "AS"
          enddo
        endif

        Allocate(Coord_new(natoms*3))
        Allocate(AtomNames_new(natoms))

        ATOMNAME="AS"
        call RemoveAtom(natoms,natoms_new,AtomNames,AtomNames_new,
     x                  Coord,Coord_new,ATOMNAME)

        if(RemoveH.eq.1)then
          ATOMNAME="H"
          call RemoveAtom(natoms,natoms_new,AtomNames,AtomNames_new,
     x                    Coord,Coord_new,ATOMNAME)
        endif

        Deallocate(Coord_new)
        Deallocate(AtomNames_new)

c assign atom radius from the data base
        do i = 1, natoms
          k = 0
          call ToUpperCase(AtomNames(i),ATOMNAME) 
          do j = 1, nTypes
            if (ATOMNAME .eq. AtomTypeName(j))then
              k = 1
              AtomRadius2(i) = AtomTypeRadius(j) * AtomTypeRadius(j)
              goto 800
            endif
          enddo
 800      continue

          if (k.eq.0)then
            write(200,*)"Can't find atom type for atom ", AtomNames(i)
            stop
          endif
        enddo

c remove atoms far away from the sphere
        Cutoff = AtomTypeRadius(1)
        do i = 2, nTypes
          if(Cutoff.lt.AtomTypeRadius(i))Cutoff=AtomTypeRadius(i)
        enddo
        Cutoff = (Cutoff + Radius + 1.0)**2

        Allocate (Coord_new(natoms*3))
        Allocate (AtomRadius2_new(natoms))
        Allocate (AtomNames_new(natoms))
        natoms_new = 0
        do i = 1, natoms
          j = (i-1)*3
          dist = Coord(j+1)**2 + Coord(j+2)**2 + Coord(j+3)**2
          if(dist.le.Cutoff)then
            natoms_new = natoms_new + 1
            k = (natoms_new - 1)*3
            Coord_new(k+1) = Coord(j+1)
            Coord_new(k+2) = Coord(j+2)
            Coord_new(k+3) = Coord(j+3)
            AtomNames_new(natoms_new) = AtomNames(i)
            AtomRadius2_new(natoms_new) = AtomRadius2(i)
          endif
        enddo
        natoms = natoms_new
        do i = 1, natoms
          j = (i-1)*3
          Coord(j+1) = Coord_new(j+1)
          Coord(j+2) = Coord_new(j+2)
          Coord(j+3) = Coord_new(j+3)
          AtomNames(i) = AtomNames_new(i)
          AtomRadius2(i) = AtomRadius2_new(i)
        enddo

        Deallocate(Coord_new)
        Deallocate(AtomNames_new)
        Deallocate(AtomRadius2_new)

c initialize some parameters for the Vbur calculation
c increase marginally square radius of sphere (Radius2) to capture
c grid points that are on the sphere surface.
c 
        Radius2 = Radius * Radius +0.0001*BinSize*BinSize
        BinSize2 = BinSize * BinSize
        BinSize3 = BinSize2 * BinSize
        NumInside = 0
        Vburf = 0.0
        Vburo = 0.0
        VburTot = 0.0

c space is the stepsize for Vbur distribuition over different frames
c for future implementation, currently hardcoded to 1.0
        space = 1.d0
        npoints = int(100.d0/space)
        notch = space / 2.d0

        do i = 1, npoints
          VburDist(i) = 0
        enddo

        do i = 1, 4
          Vbur4f(i) = 0.0
          Vbur4o(i) = 0.0
          do j = 1, npoints
            Vbur4Dist(i,j) = 0
          enddo
        enddo

        do i = 1, 8
          Vbur8f(i) = 0.0
          Vbur8o(i) = 0.0
          do j = 1, npoints
            Vbur8Dist(i,j) = 0
          enddo
        enddo

c=========================================================
c Buried volume calculation starts below by setting number
c of grid points to be sampled on each axis
c check if grid point is inside the sampled sphere
c if grid point is on the surface set volume weight = 0.5
c=========================================================
        NumPoints = int(2.0 * Radius / BinSize + 1.0)

        do i= 1, NumPoints
          x = -Radius + Real(i-1)*BinSize
          do j = 1, NumPoints
            y = -Radius  + Real(j-1)*BinSize
            zmin = Radius * 2.0
            zmax =-Radius * 2.0
            do k = 1, NumPoints
              z = -Radius  + Real(k-1)*BinSize
              IfOverlap = 0
              dist2 = x*x + y*y + z*z
              if(dist2.le.Radius2)then
                NumInside = NumInside + 1
                Wgth = 1.0
                if (abs(dist2-Radius2).lt.0.01*BinSize) Wgth = 0.5
                do l = 1, natoms
                  m = (l-1)*3
                  dx = x - Coord(m+1)
                  dy = y - Coord(m+2)
                  dz = z - Coord(m+3)
                  dist2 = dx*dx + dy*dy + dz*dz
                  if(dist2.lt.AtomRadius2(l)) IfOverlap = 1
                enddo
c
c no overlap of grid point with any atom 
                if(IfOverlap.eq.0)then
                  Vburf = Vburf + BinSize3*Wgth
                  Call Proj4(x, y, z, Vbur4f, BinSize, BinSize3, Wgth)
                  Call Proj8(x, y, z, Vbur8f, BinSize, BinSize3, Wgth)
                endif
c
c grid point overlapping with at least one atom 
                if(IfOverlap.eq.1)then
                  Vburo = Vburo + BinSize3*Wgth
                  Call Proj4(x, y, z, Vbur4o, BinSize, BinSize3, Wgth)
                  Call Proj8(x, y, z, Vbur8o, BinSize, BinSize3, Wgth)
                endif

                VburTot = VburTot + BinSize3*Wgth
              endif

c tune bottom and top surface
              if (IfOverlap.eq.1)then
                if(z.gt.zmax) zmax = z
                if(z.lt.zmin) zmin = z
              endif
            enddo
c make map
            if(DoMap.eq.1)then
              write(208, '(3f8.2)')x,y,zmax
              write(209, '(3f8.2)')x,y,zmin
            endif

          enddo
        enddo

        Volume = Radius * Radius2 * pi * 4.0 / 3.0

        temp = Vburo / VburTot * 100.0
        do j = 1, npoints
          smin = Real(j - 1) * space
          smax = Real(j) * space
          if((temp.gt.smin).and.(temp.lt.smax))
     x          VburDist(j) = VburDist(j) + 1
        enddo

        Call OutputVbur(Volume, NumInside, BinSize3, Vburf, 
     x                  Vburo, VburTot, IndFrame, VburDist, npoints, 
     x                  space, notch, IfWeb)

        write(204,*) '# notch   --        -+        ++        +-   '
        Call OutputQua4(Vbur4f, Vbur4o, IndFrame, IfWeb)
        do j = 1, npoints
          smin = REAL(j - 1) * space
          smax = REAL(j) * space
          do k = 1, 4
            temp = 100.0 * Vbur4o(k) / (Vbur4f(k) + Vbur4o(k))
            if((temp.gt.smin).and.(temp.lt.smax))
     x          Vbur4Dist(k, j) = Vbur4Dist(k, j) + 1.0
          enddo
          write(204,'(f5.1,4f10.1)') smin + notch, 
     x    (REAL(Vbur4Dist(k,j))/REAL(IndFrame), k = 1, 4)
        enddo

        write(206, *)'# notch   ---       -+-      ++-       +--
     x       -++       +++       +-+'
        Call OutputQua8(Vbur8f, Vbur8o, IndFrame, IfWeb)
        do j = 1, npoints
          smin = REAL(j - 1) * space
          smax = REAL(j) * space
          do k = 1, 4
            temp = 100.0 * Vbur8o(k) / (Vbur8f(k) + Vbur8o(k))
            if((temp.gt.smin).and.(temp.lt.smax))
     x         Vbur8Dist(k, j) = Vbur8Dist(k, j) + 1
          enddo
          write(206, '(f5.1,8f10.5)') smin + notch, 
     x    (REAL(Vbur8Dist(k,j))/REAL(IndFrame), k = 1, 8)
        enddo

        Deallocate (Coord)
        Deallocate (AtomNames)
        Deallocate (AtomRadius2)

c Below, end loop over frames
      enddo
 901  continue
 
      write(200,*)
      write(200,*)' No. of frames analyzed : ',IndFrame
      write(200,*)
      write(200,*)'Normal termination'

      if(nAtomDel.gt.0) Deallocate (IndDel)
      Deallocate (IndGC)
      Deallocate (IndAxis)
      Deallocate (Indxz)
      Deallocate (AtomTypeRadius)
      Deallocate (AtomTypeName)

      close(101)
      close(201)
      close(203)
      close(205)
      close(207)

      if(IndFrame.eq.1) then
        close(202,status="delete")
        close(204,status="delete")
        close(206,status="delete")
      else
        close(202)
        close(204)
        close(206)
      endif

      if(DoMap.eq.0) then
        close(208,status="delete")
        close(209,status="delete")
      else
        close(208)
        close(209)
      endif

      end

c=======================================================================
c This subroutine is to output a flag for the program                  c
c=======================================================================
      subroutine PrintFlag()
      Implicit None

      write(200,*)
      write(200,*)'---------------------------------------------------'
      write(200,*)'|                                                 |'
      write(200,*)'|              S A M B V C A  2.1                 |'
      write(200,*)'|                                                 |'
      write(200,*)'|              Release 29-06-2019                 |'
      write(200,*)'|                                                 |'
      write(200,*)'|       Molecular Buried Volume Calculator        |'
      write(200,*)'|                                                 |'
      write(200,*)'|  http://www.molnac.unisa.it/OM-tools/SambVca    |'
      write(200,*)'|                                                 |'
      write(200,*)'|  L. Cavallo & Z. Cao  Email: lcavallo@unisa.it  |'
      write(200,*)'|                                                 |'
      write(200,*)'---------------------------------------------------'
      write(200,*)

      end subroutine PrintFlag

c=======================================================================
c This subroutine is to check if the input file is correctly read      c
c=======================================================================
      subroutine CheckControlParameters(nAtomDel, IndDel, 
     x nAtomGC, IndGC, nAtom_z, IndAxis, nAtom_xz, Indxz, Radius,
     x Disp, BinSize, RemoveH, DoMap, 
     x AlignZ, nTypes, AtomTypeName, AtomTypeRadius)

      Implicit None

      Integer :: i

      Integer :: IfDel                          ! parameter to determine if there is atom to be deleted: > 0, yes
      Integer :: nAtomDel                       ! # of atoms to be deleted
      Integer :: RemoveH                        ! parameter to determine if remove H: =1, remove; =2, not
      Integer :: DoMap                          ! if make map: =1, make; =2, not
      Integer :: AlignZ                         ! =1, align the molecule along z; =2 align the molecule along -z

      Integer :: IndDel(nAtomDel)               ! index of the atoms to be deleted

      Integer :: nAtomGC                        ! # of atoms to define a geometry center
      Integer :: IndGC(nAtomGC)                 ! index of the atoms to define a geometry center

      Integer :: nAtom_z                        ! # of atoms to define x axis
      Integer :: IndAxis(nAtom_z)               ! index of the atoms to define an axis
      Integer :: nAtom_xz                       ! # of atoms to define a second vector
      Integer :: Indxz(nAtom_xz)                ! index of these atoms

      Real*8 :: Radius                          ! radius of the sphere
      Real*8 :: Disp                            ! displacement from the sphere center
      Real*8 :: BinSize                         ! mesh point length: resolution
      Real*8 :: AtomTypeRadius(nTypes)          ! Radius of atoms in database

      Integer :: nTypes                         ! # of atom types in the database

      Character*4 :: AtomTypeName(nTypes)          ! Radius of atoms in database

c End variables definition

      write(200,'(a39,16x,i5)')
     x"Checking: No. of atoms to be removed : ", nAtomDel
      if (nAtomDel.gt.0) then 
         write(200,'(a32)')
     x   "Checking: Atoms to be removed : "
         write(200,'(10i5)')IndDel(1:nAtomDel)
      endif

      write(200,'(a51,4x,i5)')
     x"Checking: No. of atoms defining the sphere center :", nAtomGC 
      write(200,'(a44,11x,i5)')
     x"Checking: Atoms defining the sphere center :"
      write(200,'(10i5)')IndGC(1:nAtomGC)

      write(200,'(a45,10x,i5)')
     x"Checking : No. of atoms defining the Z-axis :", nAtom_z
      write(200,'(10i5)')IndAxis(1:nAtom_z)

      write(200,'(a47,8x,i5)')
     x"Checking : No. of atoms defining the XZ-plane :", nAtom_xz 
      write(200,'(10i5)')Indxz(1:nAtom_xz)

      if (AlignZ.eq.1) write(200,*)
     x   "Checking: Molecule will be aligned along the Z+ direction"
      if (AlignZ.eq.2) write(200,*)
     x   "Checking: Molecule will be aligned along the Z- direction"

      write(200,'(a32,8x,f8.3)')
     x"Checking: Vbur sphere radius : ", Radius

      write(200,'(a38,2x,f8.3)')
     x"Checking: Displacement from origin : " , Disp

      write(200,'(a23,17x,f8.3)')
     x"Checking: Bin width : ", BinSize

      if (RemoveH.eq.1) write(200,*)
     x   "Checking: Removing H atoms"
      if (RemoveH.ne.1) write(200,*)
     x "Checking: Keeping H atoms"

      if (DoMap.eq.1) write(200,*)
     x   "Checking: Outputting steric map"
      if (DoMap.ne.1) write(200,*)
     x   "Checking: No steric map in the output"

      write(200,*)
     x"Checking: No. of atom types in your database :", nTypes
      write(200,*)
     x"Checking: Atom labels and radius in the database :"

      do i = 1,nTypes
        write(200,'(a4,f6.3)')AtomTypeName(i), AtomTypeRadius(i)
      enddo

      end subroutine CheckControlParameters

c=======================================================================
c This subroutine is to construct the original point and the 
c orthogonol axis based on the input information
c=======================================================================
      subroutine ConstructXYZ(GC, natoms, Coord, IndGC, nAtomGC, z_GC,
     x          nAtom_z, IndAxis, nAtom_xz, Indxz, 
     x          x_axis, y_axis, z_axis, x_temp_GC)
      Implicit None

      Integer :: i, j, k
      Integer :: natoms, nAtomGC, nAtom_z, nAtom_xz
      Integer :: IndGC(nAtomGC)
      Integer :: IndAxis(nAtom_z)
      Integer :: Indxz(nAtom_xz)

      Real*8 :: dist
      Real*8 :: Coord(natoms*3)
      Real*8 :: x_axis(3), y_axis(3), z_axis(3)
      Real*8 :: x_temp_GC(3)
      Real*8 :: GC(3)
      Real*8 :: z_GC(3)

c Make geometry center
      do k = 1, 3
        GC(k) = 0.0
      enddo

      do i = 1, nAtomGC
        j = IndGC(i)
        if ((j.le.0).or.(j.gt.natoms))then
          write(200,*)"wrong atoms to define geometry center"
          stop
        endif

        if ((j.gt.0).and.(j.le.natoms))then
          j = (j - 1) * 3
          do k = 1, 3
            GC(k) = GC(k) + Coord(j+k)
          enddo
        endif
      enddo

      do k = 1, 3
        GC(k) = GC(k) / REAL(nAtomGC)
      enddo

c make z-axis
      do k = 1, 3
        z_GC(k) = 0.0
      enddo

      do i = 1, nAtom_z
        j = IndAxis(i)
        if ((j.le.0).or.(j.gt.natoms))then
          write(200,*)"wrong atoms to define z_axis"
          stop
        endif
        if ((j.gt.0).and.(j.le.natoms))then
          j = (j - 1) * 3
          do k = 1, 3
            z_GC(k) = z_GC(k) + Coord(j+k)
          enddo
        endif
      enddo

      do k = 1, 3
        z_GC(k) = z_GC(k) / REAL(nAtom_z)
      enddo
      call UnitVector(z_GC, GC, z_axis, dist)

c second vector in x-z plane
      do k = 1, 3
        x_temp_GC(k) = 0.0
      enddo

      do i = 1, nAtom_xz
        j = Indxz(i)
        if ((j.le.0).or.(j.gt.natoms))then
          write(200,*)"wrong atoms to define x-z plane"
          stop
        endif
        
        if ((j.gt.0).and.(j.le.natoms))then
          j = (j - 1) * 3
          do k = 1, 3
            x_temp_GC(k) = x_temp_GC(k) + Coord(j+k)
          enddo
        endif
      enddo

      do k = 1, 3
        x_temp_GC(k) = x_temp_GC(k) / REAL(nAtom_xz)
      enddo

      Call UnitVector(x_temp_GC, GC, x_axis, dist)

c make y axis: cross product of z and x'
      Call CrossPro(z_axis, x_axis, y_axis)

c make x axis: cross product of y and z
      Call CrossPro(y_axis, z_axis, x_axis)

      end subroutine ConstructXYZ


c=======================================================================
c This subroutine is to get a vector from coordinates of two atoms
c=======================================================================
      subroutine UnitVector(coord1,coord2,vec, dist)
      Implicit None

      Integer :: k
      Real*8 :: dist
      Real*8 :: coord1(3), Coord2(3), vec(3)

      do k = 1, 3
        vec(k) = Coord1(k) - Coord2(k)
      enddo

      dist = sqrt(vec(1)**2+vec(2)**2+vec(3)**2)

      do k = 1, 3
        vec(k) = vec(k) / dist
      enddo

      end subroutine UnitVector

c=======================================================================
c This subroutine is to do the cross product                           c
c=======================================================================
      subroutine CrossPro(vec1, vec2, vec3)
      Implicit None

      Integer :: k
      Real*8 :: dist
      Real*8 :: vec1(3), vec2(3), vec3(3)

      vec3(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
      vec3(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
      vec3(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)

      dist = sqrt (vec3(1)**2+vec3(2)**2+vec3(3)**2)

      do k = 1, 3
        vec3(k) = vec3(k) / dist
      enddo

      end subroutine CrossPro

c=======================================================================
c This subroutine is to align molecule along certain direction         c
c=======================================================================
      subroutine ali(Coord, bench_vec, initial, x_axis, y_axis, z_axis, 
     x          GC, natoms)
      Implicit None

      Integer :: i, j, k
      Integer :: natoms

      Real*8 :: Coord(natoms*3)
      Real*8 :: initial(3)
      Real*8 :: bench_vec(3)
      Real*8 :: x_axis(3)
      Real*8 :: y_axis(3)
      Real*8 :: z_axis(3)
      Real*8 :: GC(3)
      Real*8 :: GC_temp(3)

      Real*8 :: vec1(3), vec2(3), vec3(3)
      Real*8 :: rot(3)
      Real*8 :: Rodrigue(3,3)

      Real*8 :: costheta, sintheta, dl, dist, dist1

c check if the bench_vec and the initial vector are the same, if so, no need to do reorientation
      if((initial(1).eq.bench_vec(1)).and.(initial(2).eq.bench_vec(2))
     x	.and.(initial(3).eq.bench_vec(3)))goto 1001

c get rotation axis
      call CrossPro(initial, bench_vec, rot)

c build Rodrigue's matrix
      costheta = initial(1) * bench_vec(1) + initial(2) * bench_vec(2) 
     x         + initial(3) * bench_vec(3)
      dl = 1.0 - costheta
      sintheta = sqrt(1.0 - costheta * costheta)

      Rodrigue(1,1) = costheta + rot(1) * rot(1) * dl
      Rodrigue(1,2) = rot(1) * rot(2) * dl - rot(3) * sintheta
      Rodrigue(1,3) = rot(2) * sintheta + rot(1) * rot(3) * dl

      Rodrigue(2,1) = rot(3) * sintheta + rot(1) * rot(2) * dl
      Rodrigue(2,2) = costheta + rot(2) * rot(2) * dl
      Rodrigue(2,3) =-rot(1) * sintheta + rot(2) * rot(3) * dl

      Rodrigue(3,1) =-rot(2) * sintheta + rot(1) * rot(3) * dl
      Rodrigue(3,2) = rot(1) * sintheta + rot(2) * rot(3) * dl
      Rodrigue(3,3) = costheta + rot(3) * rot(3) * dl

c set the geometry center to rotate the molecule around
      do k = 1, 3
        GC_temp(k) = Coord(k)
      enddo
      do i = 2, natoms
        j = (i-1) * 3
        GC_temp(1) = Coord(j+1)
        GC_temp(2) = Coord(j+2)
        GC_temp(3) = Coord(j+3)
      enddo
      do k = 1, 3
        GC_temp(k) = GC_temp(k) / REAL(natoms)
      enddo

c rotate the particle
      do i = 1, natoms
        j = (i-1) * 3
        do k = 1, 3
          vec1(k) = Coord(j+k) - GC_temp(k)
        enddo
        dist = sqrt(vec1(1)**2 + vec1(2)**2 + vec1(3)**2)
        do k = 1, 3
          vec1(k) = vec1(k) / dist
        enddo

        call matrixpro(Rodrigue, vec1, vec2, dist1)

        Coord(j+1) = GC_temp(1) + vec2(1) * dist
        Coord(j+2) = GC_temp(2) + vec2(2) * dist
        Coord(j+3) = GC_temp(3) + vec2(3) * dist
      enddo

 1001 continue
      end subroutine ali

c=======================================================================
c This subroutine is to do matrix operation                            c
c=======================================================================
      subroutine matrixPro(m, n, l, dist)
      Implicit None

      Real*8 :: m(3,3)
      Real*8 :: n(3), l(3)

      Real*8 :: dist

      l(1) = m(1,1) * n(1) + m(1,2) * n(2) + m(1,3) * n(3)
      l(2) = m(2,1) * n(1) + m(2,2) * n(2) + m(2,3) * n(3)
      l(3) = m(3,1) * n(1) + m(3,2) * n(2) + m(3,3) * n(3)

      dist = sqrt (l(1)**2 + l(2)**2 + l(3)**2)

      l(1) = l(1) / dist
      l(2) = l(2) / dist
      l(3) = l(3) / dist

      end subroutine matrixPro

c=======================================================================
c This subroutine is to move the geometry center to the origin         c
c=======================================================================
      subroutine MoveParticle(Coord,GC,natoms)
      Implicit None

      Real*8 :: Coord(natoms*3)
      Real*8 :: GC(3)

      Integer :: natoms
      Integer :: i, j

      do i = 1, natoms
        j = (i-1) * 3
        Coord(j+1) = Coord(j+1) - GC(1)
        Coord(j+2) = Coord(j+2) - GC(2)
        Coord(j+3) = Coord(j+3) - GC(3)
      enddo

      GC(1) = 0.0
      GC(2) = 0.0
      GC(3) = 0.0

      end subroutine MoveParticle

c=======================================================================
c This subroutine is to project Vbur into 4 quadrants
c=======================================================================
      subroutine Proj4(x,y,z,Vb4,BS,BS3,W)
      Implicit None
      Real*8 :: x, y, z, BS, BS3, Cut, W
      Real*8 :: Vb4(4)

      Cut = 0.25*BS

      if(abs(x).lt.Cut .and. abs(y).lt.Cut) then
        Vb4(1) = Vb4(1) + W*BS3/4.0
        Vb4(2) = Vb4(2) + W*BS3/4.0
        Vb4(3) = Vb4(3) + W*BS3/4.0
        Vb4(4) = Vb4(4) + W*BS3/4.0
        return
      endif

      if(abs(x).lt.Cut .and. y.lt.0) then
        Vb4(1) = Vb4(1) + W*BS3/2.0
        Vb4(4) = Vb4(4) + W*BS3/2.0
        return
      endif

      if(abs(x).lt.Cut .and. y.gt.0) then
        Vb4(2) = Vb4(2) + W*BS3/2.0
        Vb4(3) = Vb4(3) + W*BS3/2.0
        return
      endif

      if(x.lt.0 .and. abs(y).lt.Cut) then
        Vb4(1) = Vb4(1) + W*BS3/2.0
        Vb4(2) = Vb4(2) + W*BS3/2.0
        return
      endif

      if(x.gt.0 .and. abs(y).lt.Cut) then
        Vb4(3) = Vb4(3) + W*BS3/2.0
        Vb4(4) = Vb4(4) + W*BS3/2.0
        return
      endif

      if((x.lt.0.0).and.(y.lt.0.0))Vb4(1) = Vb4(1) +W*BS3
      if((x.lt.0.0).and.(y.gt.0.0))Vb4(2) = Vb4(2) +W*BS3
      if((x.gt.0.0).and.(y.gt.0.0))Vb4(3) = Vb4(3) +W*BS3
      if((x.gt.0.0).and.(y.lt.0.0))Vb4(4) = Vb4(4) +W*BS3

      end subroutine Proj4

c=======================================================================
c This subroutine is to project Vbur into 8 sphere octants
c=======================================================================
      subroutine Proj8(x,y,z,Vb8,BS,BS3,W)
      Implicit None
      Real*8 :: x, y, z, BS, BS3, Cut, W
      Real*8 :: Vb8(8)

      Cut = 0.25*BS

c check if point is at the origin
      if((abs(x).lt.Cut).and.(abs(y).lt.Cut).and.(abs(z).lt.Cut)) then
        Vb8(1) = Vb8(1) + W*BS3/8.0
        Vb8(2) = Vb8(2) + W*BS3/8.0
        Vb8(3) = Vb8(3) + W*BS3/8.0
        Vb8(4) = Vb8(4) + W*BS3/8.0
        Vb8(5) = Vb8(5) + W*BS3/8.0
        Vb8(6) = Vb8(6) + W*BS3/8.0
        Vb8(7) = Vb8(7) + W*BS3/8.0
        Vb8(8) = Vb8(8) + W*BS3/8.0
        return
      endif

c check if point is on -x axis
      if((x.lt.-Cut).and.(abs(y).lt.Cut).and.(abs(z).lt.Cut)) then    ! -x axis
        Vb8(1) = Vb8(1) + W*BS3/4.0
        Vb8(2) = Vb8(2) + W*BS3/4.0
        Vb8(5) = Vb8(5) + W*BS3/4.0
        Vb8(6) = Vb8(6) + W*BS3/4.0
        return
      endif

c check if point is on +x axis
      if((x.gt.+Cut).and.(abs(y).lt.Cut).and.(abs(z).lt.Cut)) then    ! +x axis
        Vb8(3) = Vb8(3) + W*BS3/4.0
        Vb8(4) = Vb8(4) + W*BS3/4.0
        Vb8(7) = Vb8(7) + W*BS3/4.0
        Vb8(8) = Vb8(8) + W*BS3/4.0
        return
      endif

c check if point is on -y axis
      if((abs(x).lt.Cut).and.(y.lt.-Cut).and.(abs(z).lt.Cut)) then    ! -y axis
        Vb8(1) = Vb8(1) + W*BS3/4.0
        Vb8(4) = Vb8(4) + W*BS3/4.0
        Vb8(5) = Vb8(5) + W*BS3/4.0
        Vb8(8) = Vb8(8) + W*BS3/4.0
        return
      endif

c check if point is on +y axis
      if((abs(x).lt.Cut).and.(y.gt.+Cut).and.(abs(z).lt.Cut)) then    ! +y axis
        Vb8(2) = Vb8(2) + W*BS3/4.0
        Vb8(3) = Vb8(3) + W*BS3/4.0
        Vb8(6) = Vb8(6) + W*BS3/4.0
        Vb8(7) = Vb8(7) + W*BS3/4.0
        return
      endif

c check if point is on -z axis
      if((abs(x).lt.Cut).and.(abs(y).lt.Cut).and.(z.lt.-Cut)) then    ! -z axis
        Vb8(1) = Vb8(1) + W*BS3/4.0
        Vb8(2) = Vb8(2) + W*BS3/4.0
        Vb8(3) = Vb8(3) + W*BS3/4.0
        Vb8(4) = Vb8(4) + W*BS3/4.0
        return
      endif

c check if point is on +z axis
      if((abs(x).lt.Cut).and.(abs(y).lt.Cut).and.(z.gt.+Cut)) then    ! +z axis
        Vb8(5) = Vb8(5) + W*BS3/4.0
        Vb8(6) = Vb8(6) + W*BS3/4.0
        Vb8(7) = Vb8(7) + W*BS3/4.0
        Vb8(8) = Vb8(8) + W*BS3/4.0
        return
      endif

c check if point is on -x-y plane
      if((x.lt.-Cut).and.(y.lt.-Cut).and.(abs(z).lt.Cut)) then    ! -x-y plane
        Vb8(1) = Vb8(1) + W*BS3/2.0
        Vb8(5) = Vb8(5) + W*BS3/2.0
        return
      endif

c check if point is on -x+y plane
      if((x.lt.-Cut).and.(y.gt.+Cut).and.(abs(z).lt.Cut)) then    ! -x+y plane
        Vb8(2) = Vb8(2) + W*BS3/2.0
        Vb8(6) = Vb8(6) + W*BS3/2.0
        return
      endif

c check if point is on +x-y plane
      if((x.gt.+Cut).and.(y.lt.-Cut).and.(abs(z).lt.Cut)) then    ! +x-y plane
        Vb8(4) = Vb8(4) + W*BS3/2.0
        Vb8(8) = Vb8(8) + W*BS3/2.0
        return
      endif

c check if point is on +x+y plane
      if((x.gt.+Cut).and.(y.gt.+Cut).and.(abs(z).lt.Cut)) then    ! +x+y plane
        Vb8(3) = Vb8(3) + W*BS3/2.0
        Vb8(7) = Vb8(7) + W*BS3/2.0
        return
      endif

c check if point is on -x-z plane
      if((x.lt.-Cut).and.(abs(y).lt.Cut).and.(z.lt.-Cut)) then    ! -x-z plane
        Vb8(1) = Vb8(1) + W*BS3/2.0
        Vb8(2) = Vb8(2) + W*BS3/2.0
        return
      endif

c check if point is on -x+z plane
      if((x.lt.-Cut).and.(abs(y).lt.Cut).and.(z.gt.+Cut)) then    ! -x+z plane
        Vb8(5) = Vb8(5) + W*BS3/2.0
        Vb8(6) = Vb8(6) + W*BS3/2.0
        return
      endif

c check if point is on +x-z plane
      if((x.gt.+Cut).and.(abs(y).lt.Cut).and.(z.lt.-Cut)) then    ! +x-z plane
        Vb8(3) = Vb8(3) + W*BS3/2.0
        Vb8(4) = Vb8(4) + W*BS3/2.0
        return
      endif

c check if point is on +x+z plane
      if((x.gt.+Cut).and.(abs(y).lt.Cut).and.(z.gt.+Cut)) then    ! +x+z plane
        Vb8(7) = Vb8(7) + W*BS3/2.0
        Vb8(8) = Vb8(8) + W*BS3/2.0
        return
      endif

c check if point is on -y-z plane
      if((abs(x).lt.Cut).and.(y.lt.-Cut).and.(z.lt.-Cut)) then    ! -y-z plane
        Vb8(1) = Vb8(1) + W*BS3/2.0
        Vb8(4) = Vb8(4) + W*BS3/2.0
        return
      endif

c check if point is on -y+z plane
      if((abs(x).lt.Cut).and.(y.lt.-Cut).and.(z.gt.+Cut)) then    ! -y+z plane
        Vb8(5) = Vb8(5) + W*BS3/2.0
        Vb8(8) = Vb8(8) + W*BS3/2.0
        return
      endif

c check if point is on +y-z plane
      if((abs(x).lt.Cut).and.(y.gt.+Cut).and.(z.lt.-Cut)) then    ! +y-z plane
        Vb8(2) = Vb8(2) + W*BS3/2.0
        Vb8(3) = Vb8(3) + W*BS3/2.0
        return
      endif

c check if point is on +y+z plane
      if((abs(x).lt.Cut).and.(y.gt.+Cut).and.(z.gt.+Cut)) then    ! +y+z plane
        Vb8(6) = Vb8(6) + W*BS3/2.0
        Vb8(7) = Vb8(7) + W*BS3/2.0
        return
      endif

      if((x.lt.0.0).and.(y.le.0.0).and.(z.lt.0.0))Vb8(1) = Vb8(1) +W*BS3 
      if((x.lt.0.0).and.(y.le.0.0).and.(z.ge.0.0))Vb8(5) = Vb8(5) +W*BS3
      if((x.lt.0.0).and.(y.gt.0.0).and.(z.lt.0.0))Vb8(2) = Vb8(2) +W*BS3
      if((x.lt.0.0).and.(y.gt.0.0).and.(z.ge.0.0))Vb8(6) = Vb8(6) +W*BS3
      if((x.ge.0.0).and.(y.ge.0.0).and.(z.lt.0.0))Vb8(3) = Vb8(3) +W*BS3
      if((x.ge.0.0).and.(y.ge.0.0).and.(z.ge.0.0))Vb8(7) = Vb8(7) +W*BS3
      if((x.ge.0.0).and.(y.lt.0.0).and.(z.lt.0.0))Vb8(4) = Vb8(4) +W*BS3
      if((x.ge.0.0).and.(y.lt.0.0).and.(z.ge.0.0))Vb8(8) = Vb8(8) +W*BS3

      end subroutine Proj8

c=======================================================================
c This subroutine is to output general infor. for Vbur calculation     c
c=======================================================================
      subroutine OutputVbur(Volume, NumInside, BinSize3, Vburf,
     x           Vburo, VburTot, IndFrame, VburDist, npoints, 
     x           space, notch, IfWeb)
      Implicit None

      Real*8 :: Volume
      Real*8 :: Vburf
      Real*8 :: Vburo
      Real*8 :: VburTot
      Real*8 :: BinSize3
      Real*8 :: space
      Real*8 :: notch

      Integer :: IfWeb
      Integer :: IndFrame
      Integer :: NumInside
      Integer :: npoints
      Integer :: i

      Integer :: VburDist(npoints)

      write(200,*)
      if (IfWeb.eq.1) write(200,*) '<STATIC_TEXT_1>'
      write(200,'(a27)') 'Results : Volumes in Angs^3'
      write(200,*)
      write(200,'(a30,i15)')' N of voxels examined : ', NumInside
      write(200,'(a30,f15.7)') ' Volume of voxel      : ',BinSize3
      write(200,*)
      if (IfWeb.eq.1) write(200,*) '<STATIC_TEXT_1>'

      write(200,*)
      if (IfWeb.eq.1) write(200,*) '<TABLE_RESULT_1>'
      write(200,*) "  V Free    V Buried   V Total   V Exact"
      write(200,'(4f10.1)') Vburf,Vburo, VburTot, Volume

      write(200,*)
      write(200,'(a20,f10.2,a7)')'Buried volume = ',Vburo,' Angs^3'
      if (IfWeb.eq.1) write(200,*) '<TABLE_RESULT_1>'

      write(200,*)
      if (IfWeb.eq.1) write(200,*) '<TABLE_RESULT_2>'
      write(200,*) "  %V Free   %V Buried  % V Tot/V Ex"
      write(200,'(4f10.1)')
     x  100.0*Vburf/VburTot, 100.0*Vburo/VburTot, 100.0*VburTot/Volume
      if (IfWeb.eq.1) write(200,*) '</TABLE_RESULT_2>'

      write(200,*)
      if (IfWeb.eq.1) write(200,*) '<STATIC_TEXT_2>'
      write(200,'(a35,f8.1)')" The %V Bur of the molecule is: ",
     x  100.0*Vburo / VburTot
      write(201,'(i6,f10.5)')IndFrame, 100.0*Vburo / VburTot
      if (IfWeb.eq.1) write(200,*) '<STATIC_TEXT_2>'
      write(200,*)

      write(202,*) '# notch   n/ntot'
      do i = 1, npoints
        write(202,'(f8.3,f10.6)') REAL(i-1)*space + notch, 
     x        VburDist(i)/Real(IndFrame)
      enddo
      end subroutine OutPutVbur

c=======================================================================
c This subroutine is to output the projected Vbur into 4 quadrants
c=======================================================================
      subroutine OutputQua4(Vbur4f, Vbur4o, IndFrame, IfWeb)
      Implicit None

      Integer :: IndFrame
      Integer :: IfWeb
      Integer :: i
      Real*8 :: Vbur4f(4)
      Real*8 :: Vbur4o(4)

      if (IfWeb.eq.1) write(200,*) '<STATIC_TEXT_3>'
      write(200,*) '                Quadrants analysis'
      if (IfWeb.eq.1) write(200,*) '<STATIC_TEXT_3>'

      if (IfWeb.eq.1) write(200,*) '<TABLE_RESULT_3>'
      write(200,*) 'Quadrant     V f     V b     V t    %V f    %V b'

      write(200,'(a5,4x,5f8.1)')"SW  ", Vbur4f(1), Vbur4o(1),
     x          (Vbur4f(1) + Vbur4o(1)),
     x          100.0 * Vbur4f(1) / (Vbur4f(1) + Vbur4o(1)),
     x          100.0 * Vbur4o(1) / (Vbur4f(1) + Vbur4o(1))
      write(200,'(a5,4x,5f8.1)')"NW  ", Vbur4f(2), Vbur4o(2),
     x          (Vbur4f(2) + Vbur4o(2)),
     x          100.0 * Vbur4f(2) / (Vbur4f(2) + Vbur4o(2)),
     x          100.0 * Vbur4o(2) / (Vbur4f(2) + Vbur4o(2))
      write(200,'(a5,4x,5f8.1)')"NE  ", Vbur4f(3), Vbur4o(3),
     x          (Vbur4f(3) + Vbur4o(3)),
     x          100.0 * Vbur4f(3) / (Vbur4f(3) + Vbur4o(3)),
     x          100.0 * Vbur4o(3) / (Vbur4f(3) + Vbur4o(3))
      write(200,'(a5,4x,5f8.1)')"SE  ", Vbur4f(4), Vbur4o(4),
     x          (Vbur4f(4) + Vbur4o(4)),
     x          100.0 * Vbur4f(4) / (Vbur4f(4) + Vbur4o(4)),
     x          100.0 * Vbur4o(4) / (Vbur4f(4) + Vbur4o(4))

      write(203,*) '# notch   --        -+        ++        +-   '
      write(203,'(i5,4f10.5)')IndFrame, 
     x     (100.0*(Vbur4o(i)/(Vbur4f(i) + Vbur4o(i))), i = 1, 4)
      if (IfWeb.eq.1) write(200,*) '</TABLE_RESULT_3>'

      end subroutine OutputQua4

c=======================================================================
c This subroutine is to output the projected Vbur into 8 octants
c=======================================================================
      subroutine OutputQua8(Vbur8f, Vbur8o, IndFrame, IfWeb)
      Implicit None

      Integer :: IndFrame
      Integer :: IfWeb
      Integer :: i

      Real*8 :: Vbur8f(8)
      Real*8 :: Vbur8o(8)

      write(200,*)
      if(IfWeb.eq.1)write(200,*) '<STATIC_TEXT_4>'
      write(200,*) '                  Octants analysis'
      if(IfWeb.eq.1)write(200,*) '<STATIC_TEXT_4>'

      if(IfWeb.eq.1)write(200,*) '<TABLE_RESULT_4>'
      write(200,*)'Octant       V f      V b    V t    %V f    %V b'

      write(200,'(a5,4x,5f8.1)')" SW-z", Vbur8f(1), Vbur8o(1),
     x          (Vbur8f(1) + Vbur8o(1)),
     x          100.0 * Vbur8f(1) / (Vbur8f(1) + Vbur8o(1)),
     x          100.0 * Vbur8o(1) / (Vbur8f(1) + Vbur8o(1))
      write(200,'(a5,4x,5f8.1)')" NW-z", Vbur8f(2), Vbur8o(2),
     x          (Vbur8f(2) + Vbur8o(2)),
     x          100.0 * Vbur8f(2) / (Vbur8f(2) + Vbur8o(2)),
     x          100.0 * Vbur8o(2) / (Vbur8f(2) + Vbur8o(2))
      write(200,'(a5,4x,5f8.1)')" NE-z", Vbur8f(3), Vbur8o(3), 
     x          (Vbur8f(3) + Vbur8o(3)),
     x          100.0 * Vbur8f(3) / (Vbur8f(3) + Vbur8o(3)),
     x          100.0 * Vbur8o(3) / (Vbur8f(3) + Vbur8o(3))
      write(200,'(a5,4x,5f8.1)')" SE-z", Vbur8f(4), Vbur8o(4),
     x          (Vbur8f(4) + Vbur8o(4)),
     x          100.0 * Vbur8f(4) / (Vbur8f(4) + Vbur8o(4)),
     x          100.0 * Vbur8o(4) / (Vbur8f(4) + Vbur8o(4))
      write(200,'(a5,4x,5f8.1)')" SW+z", Vbur8f(5), Vbur8o(5),
     x          (Vbur8f(5) + Vbur8o(5)),
     x          100.0 * Vbur8f(5) / (Vbur8f(5) + Vbur8o(5)),
     x          100.0 * Vbur8o(5) / (Vbur8f(5) + Vbur8o(5))
      write(200,'(a5,4x,5f8.1)')" NW+z", Vbur8f(6), Vbur8o(6),
     x          (Vbur8f(6) + Vbur8o(6)),
     x          100.0 * Vbur8f(6) / (Vbur8f(6) + Vbur8o(6)),
     x          100.0 * Vbur8o(6) / (Vbur8f(6) + Vbur8o(6))
      write(200,'(a5,4x,5f8.1)')" NE+z", Vbur8f(7), Vbur8o(7),
     x          (Vbur8f(7) + Vbur8o(7)),
     x          100.0 * Vbur8f(7) / (Vbur8f(7) + Vbur8o(7)),
     x          100.0 * Vbur8o(7) / (Vbur8f(7) + Vbur8o(7))
      write(200,'(a5,4x,5f8.1)')" SE+z", Vbur8f(8), Vbur8o(8),
     x          (Vbur8f(8) + Vbur8o(8)),
     x          100.0 * Vbur8f(8) / (Vbur8f(8) + Vbur8o(8)),
     x          100.0 * Vbur8o(8) / (Vbur8f(8) + Vbur8o(8))

      write(205,*) '# notch   ---       -+-      ++-       +--      --+
     x       -++       +++       +-+'
      write(205, '(i5,8f10.5)')IndFrame, 
     x     (100.0*(Vbur8o(i)/(Vbur8f(i) + Vbur8o(i))), i = 1, 8)

      if(IfWeb.eq.1)write(200,*) '</TABLE_RESULT_4>'

      end subroutine OutputQua8

c=======================================================================
c This subroutine converts atom names from lower to upper case
c=======================================================================
      subroutine  ToUpperCase(strIn,strOut)
      implicit none
 
      integer :: i,j
 
      character*4 :: StrIN, StrOut

      do i = 1, len(strIn)
        j = iachar(strIn(i:i))
        if (j>= iachar("a") .and. j<=iachar("z") ) then
          strOut(i:i) = achar(iachar(strIn(i:i))-32)
        else
          strOut(i:i) = strIn(i:i)
        end if
      end do
 
      end subroutine ToUpperCase

c=======================================================================
c This subroutine writes the molecule in file iouint 
c=======================================================================
      subroutine  WriteXYZ(iounit,natoms,Title,AtomNames,Coord)
      implicit none

      Integer :: i,j, iounit
      Integer :: natoms 

      Real*8 :: Coord(natoms*3)

      Character*50 :: Title
      Character*4 :: AtomNames(natoms)

      write(iounit,*)natoms
      write(iounit,'(a39)')Title
      do i = 1, natoms
        j = (i-1)*3
        write(iounit,'(a4,1x,3f10.5)')AtomNames(i),
     x  Coord(j+1),Coord(j+2),Coord(j+3)
      enddo

      end subroutine

c=======================================================================
c This subroutine removes atoms type ATOMNAME from coordinates
c=======================================================================
      subroutine RemoveAtom(natoms,natoms_new,AtomNames,AtomNames_new,
     x           Coord,Coord_new,ATOMNAME)
      implicit none

      Integer :: i,j,k
      Integer :: natoms, natoms_new

      Real*8 :: Coord(natoms*3)
      Real*8 :: Coord_new(natoms*3)

      Character*4 :: ATOMNAME
      Character*4 :: AtomNames(natoms)
      Character*4 :: AtomNames_new(natoms)

      natoms_new = 0
      do i=1, natoms
        if(AtomNames(i).ne.ATOMNAME)then
          j = (i-1)*3
          natoms_new = natoms_new + 1
          k = (natoms_new-1)*3
          Coord_new(k+1) = Coord(j+1)
          Coord_new(k+2) = Coord(j+2)
          Coord_new(k+3) = Coord(j+3)
          AtomNames_new(natoms_new) = AtomNames(i)
        endif
      enddo

      do i = 1, natoms_new
        j = (i-1)*3
        Coord(j+1) = Coord_new(j+1)
        Coord(j+2) = Coord_new(j+2)
        Coord(j+3) = Coord_new(j+3)
        AtomNames(i) = AtomNames_new(i)
      enddo

      natoms = natoms_new

      end subroutine RemoveAtom
