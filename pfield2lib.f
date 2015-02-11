c-----------------------------------------------------------------------
c 2d parallel PIC library for solving field equations
c pfield2lib.f contains procedures to manage guard cells and solve
c              fields equations in fourier space:
c PCGUARD2X copy guard cells in x for 2 component vector array,
c           quadratic interpolation, and distributed data.
c PBGUARD2X copy guard cells in x for 3 component vector array,
c           quadratic interpolation, and distributed data.
c PDGUARD2X copy guard cells in x for scalar array, quadratic
c           interpolation, and distributed data.
c PSCGUARD2 initialize field for 3 component vector array, quadratic
c           interpolation, and distributed data.
c PSCGUARD22 initialize field for 2 component vector array, quadratic
c            interpolation, and distributed data.
c PSGUARD2 initialize field for scalar array, quadratic interpolation,
c          and distributed data.
c PACGUARD2X add guard cells in x for 3 component vector array,
c            quadratic interpolation, and distributed data.
c PACGUARD22X add guard cells in x for 2 component vector array,
c             quadratic interpolation, and distributed data.
c PAGUARD2X add guard cells in x for scalar array, quadratic
c           interpolation, and distributed data.
c PZCGUARD2 zeros out guard cells in extended periodic 3 component
c           vector field, quadratic interpolation, for distributed data.
c PZCGUARD22 zeros out guard cells in extended periodic 2 component
c            vector field, quadratic interpolation, for distributed data
c PZGUARD2 zeros out guard cells in extended periodic scalar field
c          quadratic interpolation, for distributed data.
c PCGUARD2XL copy guard cells in x for 2 component vector array, linear
c            interpolation, for distributed data.
c PBGUARD2XL copy guard cells in x for 3 component vector array, linear
c            interpolation, for distributed data.
c PDGUARD2XL copy guard cells in x for scalar array, linear
c            interpolation, for distributed data.
c PSCGUARD2L initialize field for 3 component vector array, linear
c            interpolation, for distributed data.
c PSCGUARD22L initialize field for 2 component vector array, linear
c             interpolation, for distributed data.
c PSGUARD2L initialize field for scalar array, linear interpolation, for
c           distributed data.
c PACGUARD2XL add guard cells in x for 3 component vector array, linear
c             interpolation, for distributed data.
c PACGUARD22XL add guard cells in x for 2 component vector array, linear
c              interpolation, for distributed data.
c PAGUARD2XL add guard cells in x for scalar array, linear
c            interpolation, for distributed data.
c PZCGUARD2L zeros out guard cells in extended periodic 3 component
c            vector field, linear interpolation, for distributed data.
c PZCGUARD22L zeros out guard cells in extended periodic 2 component
c             vector field, linear interpolation, for distributed data
c PZGUARD2L zeros out guard cells in extended periodic scalar field
c           linear interpolation, for distributed data.
c PPOISP2 solve 2d poisson equation for electric force, potential, or
c         smoothing, for distributed data.
c PPOISP21 solve 1d poisson equation for electric force, potential, or
c          smoothing, for distributed data.
c PPOIS22 solve 2d poisson equation for electric force, for distributed
c         data.
c PPOIS23 solve 2-1/2d poisson equation for electric force, for
c         distributed data.
c PDIVF2 calculates 2d divergence of n component vector in fourier
c        space, for distributed data.
c PGRADF2 calculates 2d gradient of scalar field in fourier space, for
c         distributed data.
c PCURLF2 calculates 2d divergence of 3 component vector in fourier
c         space, for distributed data.
c PCURLF22 calculates 2d divergence of 2 component vector in fourier
c          space, for distributed data.
c PCUPERP2 calculates 2d tranvsere current of 3 component vector in
c          fourier space, for distributed data.
c PCUPERP22 calculates 2d tranvsere current of 2 component vector in
c           fourier space, for distributed data.
c PBPOISP23 solve 2-1/2d vector poisson equation for magnetic force,
c           vector potential, or smoothing, for distributed data.
c PBPOISP22 solve 2d vector poisson equation for magnetic force, vector
c           potential, or smoothing, for distributed data.
c IPBPOISP23 solve 2-1/2d vector poisson equation for magnetic field,
c            for distributed data.
c PMAXWEL2 solve 2d maxwell equation for electric and magnetic fields,
c          for distributed data.
c PEMFIELD2 combines and smooths 2d periodic electric magnetic forces,
c           for distributed data.
c PEMFIELDR2 combines and smooths 2d real electric or magnetic forces
c            for sine-cosine transforms, for distributed data.
c PAVPOT23 calculate 2-1/2d vector potential from magnetic field, for
c          distributed data.
c PAVRPOT32 calculate 2-1/2d radiative part of the vector potential.
c PGTMODES2 extracts selected 2d fourier components from potential
c           array, for distributed data.
c PPTMODES2 places selected 2d fourier components into potential array,
c           for distributed data.
c PGTVMODES2 extracts selected 2d fourier components from vector
c            potential array, for distributed data.
c PPTVMODES2 places selected 2d fourier components into vector potential
c            array, for distributed data.
c PPOYNT2 calculate poynting electromagnetic flux.
c PDPOYNT2 calculate electromagnetic flux in 2-1/2d Darwin field.
c PDPOYNT22 calculate electromagnetic flux in 2d Darwin field.
c PSCFGUARD2 initialize 3 component field with scaled vector array,
c            quadratic interpolation.
c PSCFGUARD22 initialize 2 component field with scaled vector array,
c             quadratic interpolation.
c PSMCGUARD2 initialize field for 4 component tensor array, quadratic
c            interpolation.
c PSMCGUARD22 initialize field for 2 component tensor array, quadratic
c             interpolation.
c PAMCGUARD2X add guard cells in x for n component tensor array,
c             quadratic interpolation, and distributed data.
c PSCFGUARD2L initialize 3 component field with scaled vector array,
c             linear interpolation.
c PSCFGUARD22L initialize 2 component field with scaled vector array,
c              linear interpolation.
c PSMCGUARD2L initialize field for 4 component tensor array, linear
c             interpolation.
c PSMCGUARD22L initialize field for 2 component tensor array, linear
c              interpolation.
c PAMCGUARD2XL add guard cells in x for n component tensor array, linear
c              interpolation, for distributed data.
c PDCUPERP23 calculate 2-1/d transverse derivative of current density
c            from momentum flux.
c PDCUPERP22 calculate 2d transverse derivative of current density from
c            momentum flux.
c PADCUPERP23 calculate 2-1/2d transverse derivative of current density
c             from momentum flux and acceleration density.
c PADCUPERP22 calculate 2d transverse derivative of current density
c             from momentum flux and acceleration density.
c PEPOISP23 solve vector poisson equation for 2-1/2d transverse electric
c           field or force.
c PEPOISP22 solve vector poisson equation for 2d transverse electric
c           field or force.
c PADDQEI2 adds electron and ion densities, for distributed data.
c PADDQEI2X adds electron and ion densities and calculates maximum and
c           minimum plasma frequency, for distributed data.
c PBADDEXT2 adds constant to magnetic field in real space for 2-1/2d
c           code, for distributed data.
c PBADDEXT22 adds constant to magnetic field in real space for 2d code,
c            for distributed data.
c PIMOMENT2 calculates ion momentum for 2-1/2d code from qi*fxy, for
c           distributed data.
c PIMOMENT22 calculates ion momentum for 2d code from qi*fxy, for
c            distributed data.
c PADDVRFIELD2 calculates a = b + c for distributed real vector fields.
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: november 15, 2009
c-----------------------------------------------------------------------
      subroutine PCGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
c replicate extended periodic field
c quadratic interpolation, for distributed data
      implicit none
      real fxy
      integer nyp, nx, nxe, nypmx, nblok
      dimension fxy(2,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, k, l, nyp3
      do 30 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 20 k = 1, nyp3
      do 10 i = 1, 2
      fxy(i,1,k,l) = fxy(i,nx+1,k,l)
      fxy(i,nx+2,k,l) = fxy(i,2,k,l)
      fxy(i,nx+3,k,l) = fxy(i,3,k,l)
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
c replicate extended periodic scalar field
c quadratic interpolation, for distributed data
      implicit none
      real q
      integer nyp, nx, nxe, nypmx, nblok
      dimension q(nxe,nypmx,nblok), nyp(nblok)
c local data
      integer k, l, nyp3
      do 20 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 10 k = 1, nyp3
      q(1,k,l) = q(nx+1,k,l)
      q(nx+2,k,l) = q(2,k,l)
      q(nx+3,k,l) = q(3,k,l)
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBGUARD2X(bxy,nyp,nx,nxe,nypmx,nblok)
c replicate extended periodic vector field
c quadratic interpolation, for distributed data
      implicit none
      real bxy
      integer nyp, nx, nxe, nypmx, nblok
      dimension bxy(3,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, k, l, nyp3
      do 30 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 20 k = 1, nyp3
      do 10 i = 1, 3
      bxy(i,1,k,l) = bxy(i,nx+1,k,l)
      bxy(i,nx+2,k,l) = bxy(i,2,k,l)
      bxy(i,nx+3,k,l) = bxy(i,3,k,l)
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSCGUARD2(cu,nyp,xj0,yj0,zj0,nx,nxe,nypmx,nblok)
c initialize extended periodic field
c quadratic interpolation, for distributed data
      implicit none
      real cu, xj0, yj0, zj0
      integer nyp, nx, nxe, nypmx, nblok
      dimension cu(3,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, j, k, l, nyp3, nx3
c initialize extended field, with zero in the edges
      nx3 = nx + 3
      do 60 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 30 k = 1, nyp(l)
      do 10 j = 1, nx
      cu(1,j+1,k+1,l) = xj0
      cu(2,j+1,k+1,l) = yj0
      cu(3,j+1,k+1,l) = zj0
   10 continue
      do 20 i = 1, 3
      cu(i,1,k+1,l) = 0.
      cu(i,nx+2,k+1,l) = 0.
      cu(i,nx+3,k+1,l) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx3
      do 40 i = 1, 3
      cu(i,j,1,l) = 0.
      cu(i,j,nyp3-1,l) = 0.
      cu(i,j,nyp3,l) = 0.
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSCGUARD22(cu,nyp,xj0,yj0,nx,nxe,nypmx,nblok)
c initialize extended periodic field
c quadratic interpolation, for distributed data
      implicit none
      real cu, xj0, yj0
      integer nyp, nx, nxe, nypmx, nblok
      dimension cu(2,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, j, k, l, nyp3, nx3
c initialize extended field, with zero in the edges
      nx3 = nx + 3
      do 60 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 30 k = 1, nyp(l)
      do 10 j = 1, nx
      cu(1,j+1,k+1,l) = xj0
      cu(2,j+1,k+1,l) = yj0
   10 continue
      do 20 i = 1, 2
      cu(i,1,k+1,l) = 0.
      cu(i,nx+2,k+1,l) = 0.
      cu(i,nx+3,k+1,l) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx3
      do 40 i = 1, 2
      cu(i,j,1,l) = 0.
      cu(i,j,nyp3-1,l) = 0.
      cu(i,j,nyp3,l) = 0.
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSGUARD2(q,nyp,qi0,nx,nxe,nypmx,nblok)
c initialize extended periodic scalar field
c quadratic interpolation, for distributed data
      implicit none
      real q, qi0
      integer nyp, nx, nxe, nypmx, nblok
      dimension q(nxe,nypmx,nblok), nyp(nblok)
c local data
      integer j, k, l, nyp3, nx3
c initialize extended field, with zero in the edges
      nx3 = nx + 3
      do 40 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 20 k = 1, nyp(l)
      do 10 j = 1, nx
      q(j+1,k+1,l) = qi0
   10 continue
      q(1,k+1,l) = 0.
      q(nx+2,k+1,l) = 0.
      q(nx+3,k+1,l) = 0.
   20 continue
      do 30 j = 1, nx3
      q(j,1,l) = 0.
      q(j,nyp3-1,l) = 0.
      q(j,nyp3,l) = 0.
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PACGUARD2X(cu,nyp,nx,nxe,nypmx,nblok)
c accumulate extended periodic vector field
c quadratic interpolation, for distributed data
      implicit none
      real cu
      integer nyp, nx, nxe, nypmx, nblok
      dimension cu(3,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, k, l, nyp3
c accumulate edges of extended field
      do 30 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 20 k = 1, nyp3
      do 10 i = 1, 3
      cu(i,2,k,l) = cu(i,2,k,l) + cu(i,nx+2,k,l)
      cu(i,3,k,l) = cu(i,3,k,l) + cu(i,nx+3,k,l)
      cu(i,nx+1,k,l) = cu(i,nx+1,k,l) + cu(i,1,k,l)
      cu(i,1,k,l) = 0.
      cu(i,nx+2,k,l) = 0.
      cu(i,nx+3,k,l) = 0.
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PACGUARD22X(cu,nyp,nx,nxe,nypmx,nblok)
c accumulate extended periodic vector field
c quadratic interpolation, for distributed data
      implicit none
      real cu
      integer nyp, nx, nxe, nypmx, nblok
      dimension cu(2,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, k, l, nyp3
c accumulate edges of extended field
      do 30 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 20 k = 1, nyp3
      do 10 i = 1, 2
      cu(i,2,k,l) = cu(i,2,k,l) + cu(i,nx+2,k,l)
      cu(i,3,k,l) = cu(i,3,k,l) + cu(i,nx+3,k,l)
      cu(i,nx+1,k,l) = cu(i,nx+1,k,l) + cu(i,1,k,l)
      cu(i,1,k,l) = 0.
      cu(i,nx+2,k,l) = 0.
      cu(i,nx+3,k,l) = 0.
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
c accumulate extended periodic scalar field
c quadratic interpolation, for distributed data
      implicit none
      real q
      integer nyp, nx, nxe, nypmx, nblok
      dimension q(nxe,nypmx,nblok), nyp(nblok)
c local data
      integer k, l, nyp3
c accumulate edges of extended field
      do 20 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 10 k = 1, nyp3
      q(2,k,l) = q(2,k,l) + q(nx+2,k,l)
      q(3,k,l) = q(3,k,l) + q(nx+3,k,l)
      q(nx+1,k,l) = q(nx+1,k,l) + q(1,k,l)
      q(1,k,l) = 0.
      q(nx+2,k,l) = 0.
      q(nx+3,k,l) = 0.
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PZCGUARD2(cu,nyp,nx,nxe,nypmx,nblok)
c zero out guard cells in extended periodic vector field
c quadratic interpolation, for distributed data
      implicit none
      real cu
      integer nyp, nx, nxe, nypmx, nblok
      dimension cu(3,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, j, k, l, nx3
      nx3 = nx + 3
      do 50 l = 1, nblok
c zero out guard cells in x
      do 20 k = 1, nyp(l)
      do 10 i = 1, 3
      cu(i,1,k+1,l) = 0.
      cu(i,nx+2,k+1,l) = 0.
      cu(i,nx+3,k+1,l) = 0.
   10 continue
   20 continue
c zero out guard cells in y
      do 40 j = 1, nx3
      do 30 i = 1, 3
      cu(i,j,1,l) = 0.
      cu(i,j,nyp(l)+2,l) = 0.
      cu(i,j,nyp(l)+3,l) = 0.
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PZCGUARD22(cu,nyp,nx,nxe,nypmx,nblok)
c zero out guard cells in extended periodic vector field
c quadratic interpolation, for distributed data
      implicit none
      real cu
      integer nyp, nx, nxe, nypmx, nblok
      dimension cu(2,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, j, k, l, nx3
      nx3 = nx + 3
      do 50 l = 1, nblok
c zero out guard cells in x
      do 20 k = 1, nyp(l)
      do 10 i = 1, 2
      cu(i,1,k+1,l) = 0.
      cu(i,nx+2,k+1,l) = 0.
      cu(i,nx+3,k+1,l) = 0.
   10 continue
   20 continue
c zero out guard cells in y
      do 40 j = 1, nx3
      do 30 i = 1, 2
      cu(i,j,1,l) = 0.
      cu(i,j,nyp(l)+2,l) = 0.
      cu(i,j,nyp(l)+3,l) = 0.
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PZGUARD2(q,nyp,nx,nxe,nypmx,nblok)
c zero out guard cells in extended periodic scalar field
c quadratic interpolation, for distributed data
      implicit none
      real q
      integer nyp, nx, nxe, nypmx, nblok
      dimension q(nxe,nypmx,nblok), nyp(nblok)
c local data
      integer j, k, l, nx3
      nx3 = nx + 3
      do 30 l = 1, nblok
c zero out guard cells in x
      do 10 k = 1, nyp(l)
      q(1,k+1,l) = 0.
      q(nx+2,k+1,l) = 0.
      q(nx+3,k+1,l) = 0.
   10 continue
c zero out guard cells in y
      do 20 j = 1, nx3
      q(j,1,l) = 0.
      q(j,nyp(l)+2,l) = 0.
      q(j,nyp(l)+3,l) = 0.
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCGUARD2XL(fxy,nyp,nx,nxe,nypmx,nblok)
c replicate extended periodic field
c linear interpolation, for distributed data
      implicit none
      real fxy
      integer nyp, nx, nxe, nypmx, nblok
      dimension fxy(2,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer k, l, nyp1
      do 20 l = 1, nblok
      nyp1 = nyp(l) + 1
      do 10 k = 1, nyp1
      fxy(1,nx+1,k,l) = fxy(1,1,k,l)
      fxy(2,nx+1,k,l) = fxy(2,1,k,l)
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDGUARD2XL(q,nyp,nx,nxe,nypmx,nblok)
c replicate extended periodic scalar field
c linear interpolation, for distributed data
      implicit none
      real q
      integer nyp, nx, nxe, nypmx, nblok
      dimension q(nxe,nypmx,nblok), nyp(nblok)
c local data
      integer k, l, nyp1
      do 20 l = 1, nblok
      nyp1 = nyp(l) + 1
      do 10 k = 1, nyp1
      q(nx+1,k,l) = q(1,k,l)
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBGUARD2XL(bxy,nyp,nx,nxe,nypmx,nblok)
c replicate extended periodic vector field
c linear interpolation, for distributed data
      implicit none
      real bxy
      integer nyp, nx, nxe, nypmx, nblok
      dimension bxy(3,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer k, l, nyp1
      do 20 l = 1, nblok
      nyp1 = nyp(l) + 1
      do 10 k = 1, nyp1
      bxy(1,nx+1,k,l) = bxy(1,1,k,l)
      bxy(2,nx+1,k,l) = bxy(2,1,k,l)
      bxy(3,nx+1,k,l) = bxy(3,1,k,l)
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSCGUARD2L(cu,nyp,xj0,yj0,zj0,nx,nxe,nypmx,nblok)
c initialize extended periodic field
c linear interpolation, for distributed data
      implicit none
      real cu, xj0, yj0, zj0
      integer nyp, nx, nxe, nypmx, nblok
      dimension cu(3,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, j, k, l, nyp1, nx1
c initialize extended field, with zero in the edges
      nx1 = nx + 1
      do 60 l = 1, nblok
      nyp1 = nyp(l) + 1
      do 30 k = 1, nyp(l)
      do 10 j = 1, nx
      cu(1,j,k,l) = xj0
      cu(2,j,k,l) = yj0
      cu(3,j,k,l) = zj0
   10 continue
      do 20 i = 1, 3
      cu(i,nx+1,k,l) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx1
      do 40 i = 1, 3
      cu(i,j,nyp1,l) = 0.
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSCGUARD22L(cu,nyp,xj0,yj0,nx,nxe,nypmx,nblok)
c initialize extended periodic field
c linear interpolation, for distributed data
      implicit none
      real cu, xj0, yj0
      integer nyp, nx, nxe, nypmx, nblok
      dimension cu(2,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, j, k, l, nyp1, nx1
c initialize extended field, with zero in the edges
      nx1 = nx + 1
      do 60 l = 1, nblok
      nyp1 = nyp(l) + 1
      do 30 k = 1, nyp(l)
      do 10 j = 1, nx
      cu(1,j,k,l) = xj0
      cu(2,j,k,l) = yj0
   10 continue
      do 20 i = 1, 2
      cu(i,nx+1,k,l) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx1
      do 40 i = 1, 2
      cu(i,j,nyp1,l) = 0.
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSGUARD2L(q,nyp,qi0,nx,nxe,nypmx,nblok)
c initialize extended periodic scalar field
c linear interpolation, for distributed data
      implicit none
      real q, qi0
      integer nyp, nx, nxe, nypmx, nblok
      dimension q(nxe,nypmx,nblok), nyp(nblok)
c local data
      integer j, k, l, nyp1, nx1
c initialize extended field, with zero in the edges
      nx1 = nx + 1
      do 40 l = 1, nblok
      nyp1 = nyp(l) + 1
      do 20 k = 1, nyp(l)
      do 10 j = 1, nx
      q(j,k,l) = qi0
   10 continue
      q(nx+1,k,l) = 0.
   20 continue
      do 30 j = 1, nx1
      q(j,nyp1,l) = 0.
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PACGUARD2XL(cu,nyp,nx,nxe,nypmx,nblok)
c accumulate extended periodic vector field
c linear interpolation, for distributed data
      implicit none
      real cu
      integer nyp, nx, nxe, nypmx, nblok
      dimension cu(3,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, k, l, nyp1
c accumulate edges of extended field
      do 30 l = 1, nblok
      nyp1 = nyp(l) + 1
      do 20 k = 1, nyp1
      do 10 i = 1, 3
      cu(i,1,k,l) = cu(i,1,k,l) + cu(i,nx+1,k,l)
      cu(i,nx+1,k,l) = 0.
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PACGUARD22XL(cu,nyp,nx,nxe,nypmx,nblok)
c accumulate extended periodic vector field
c linear interpolation, for distributed data
      implicit none
      real cu
      integer nyp, nx, nxe, nypmx, nblok
      dimension cu(2,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, k, l, nyp1
c accumulate edges of extended field
      do 30 l = 1, nblok
      nyp1 = nyp(l) + 1
      do 20 k = 1, nyp1
      do 10 i = 1, 2
      cu(i,1,k,l) = cu(i,1,k,l) + cu(i,nx+1,k,l)
      cu(i,nx+1,k,l) = 0.
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAGUARD2XL(q,nyp,nx,nxe,nypmx,nblok)
c accumulate extended periodic scalar field
c linear interpolation, for distributed data
      implicit none
      real q
      integer nyp, nx, nxe, nypmx, nblok
      dimension q(nxe,nypmx,nblok), nyp(nblok)
c local data
      integer k, l, nyp1
c accumulate edges of extended field
      do 20 l = 1, nblok
      nyp1 = nyp(l) + 1
      do 10 k = 1, nyp1
      q(1,k,l) = q(1,k,l) + q(nx+1,k,l)
      q(nx+1,k,l) = 0.
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PZCGUARD2L(cu,nyp,nx,nxe,nypmx,nblok)
c zero out guard cells in extended periodic vector field
c linear interpolation, for distributed data
      implicit none
      real cu
      integer nyp, nx, nxe, nypmx, nblok
      dimension cu(3,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, j, k, l, nx1
      nx1 = nx + 1
c zero out guard cells in x
      do 50 l = 1, nblok
c zero out guard cells in x
      do 20 k = 1, nyp(l)
      do 10 i = 1, 3
      cu(i,nx+1,k,l) = 0.
   10 continue
   20 continue
c zero out guard cells in y
      do 40 j = 1, nx1
      do 30 i = 1, 3
      cu(i,j,nyp(l)+1,l) = 0.
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PZCGUARD22L(cu,nyp,nx,nxe,nypmx,nblok)
c zero out guard cells in extended periodic vector field
c linear interpolation, for distributed data
      implicit none
      real cu
      integer nyp, nx, nxe, nypmx, nblok
      dimension cu(2,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, j, k, l, nx1
      nx1 = nx + 1
c zero out guard cells in x
      do 50 l = 1, nblok
c zero out guard cells in x
      do 20 k = 1, nyp(l)
      do 10 i = 1, 2
      cu(i,nx+1,k,l) = 0.
   10 continue
   20 continue
c zero out guard cells in y
      do 40 j = 1, nx1
      do 30 i = 1, 2
      cu(i,j,nyp(l)+1,l) = 0.
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PZGUARD2L(q,nyp,nx,nxe,nypmx,nblok)
c zero out guard cells in extended periodic scalar field
c linear interpolation, for distributed data
      implicit none
      real q
      integer nyp, nx, nxe, nypmx, nblok
      dimension q(nxe,nypmx,nblok), nyp(nblok)
c local data
      integer j, k, l, nx1
      nx1 = nx + 1
c zero out guard cells in x
      do 30 l = 1, nblok
c zero out guard cells in x
      do 10 k = 1, nyp(l)
      q(nx+1,k,l) = 0.
   10 continue
c zero out guard cells in y
      do 20 j = 1, nx1
      q(j,nyp(l)+1,l) = 0.
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISP2(q,fx,fy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv
     1,kxp,jblok,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function,
c with periodic boundary conditions for distributed data.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: ffc
c for isign = -1, input: q,ffc,isign,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: fx,fy,we
c approximate flop count is: 33*nxc*nyc + 15*(nxc + nyc)
c for isign = 1, input: q,ffc,isign,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: fx,we
c approximate flop count is: 21*nxc*nyc + 11*(nxc + nyc)
c for isign = 2, input: q,ffc,isign,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: fy
c approximate flop count is: 4*nxc*nyc + 2*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*q(kx,ky)*s(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*q(kx,ky)*s(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
c fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky) = g(kx,ky)*q(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky) = q(kx,ky)*s(kx,ky)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fx(k,j,l) = x component of complex force/charge,
c fy(k,j,l) = y component of complex force/charge,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c jblok = number of data blocks
c kxp = number of data values per block
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c aimag(ffc(k,j,l)) = finite-size particle shape factor s
c real(ffc(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp
      complex q, fx, fy, ffc, zero
      dimension q(nyv,kxp,jblok)
      dimension fx(nyv,kxp,jblok), fy(nyv,kxp,jblok)
      dimension ffc(nyhd,kxp,jblok)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nxh) return
c prepare form factor array
      do 30 l = 1, jblok
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffc(k,j,l) = cmplx(affp,1.)
      else
         ffc(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 100
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nxh) go to 90
      do 80 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 60 j = 1, kxp
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,j,l))*aimag(ffc(k,j,l))
         at2 = dkx*at1
         at3 = dny*float(k - 1)*at1
         fx(k,j,l) = at2*cmplx(aimag(q(k,j,l)),-real(q(k,j,l)))
         fx(k1,j,l) = at2*cmplx(aimag(q(k1,j,l)),-real(q(k1,j,l)))
         fy(k,j,l) = at3*cmplx(aimag(q(k,j,l)),-real(q(k,j,l)))
         fy(k1,j,l) = at3*cmplx(-aimag(q(k1,j,l)),real(q(k1,j,l)))
         wp = wp + at1*(q(k,j,l)*conjg(q(k,j,l)) + q(k1,j,l)*conjg(q(k1,
     1j,l)))
   50    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = real(ffc(1,j,l))*aimag(ffc(1,j,l))
         fx(1,j,l) = dkx*at1*cmplx(aimag(q(1,j,l)),-real(q(1,j,l)))
         fx(k1,j,l) = zero
         fy(1,j,l) = zero
         fy(k1,j,l) = zero
         wp = wp + at1*(q(1,j,l)*conjg(q(1,j,l)))
      endif
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,1,l))*aimag(ffc(k,1,l))
         fx(k,1,l) = zero
         fx(k1,1,l) = zero
         fy(k,1,l) = dny*float(k - 1)*at1*cmplx(aimag(q(k,1,l)),-real(q(
     1k,1,l)))
         fy(k1,1,l) = zero
         wp = wp + at1*(q(k,1,l)*conjg(q(k,1,l)))
   70    continue
         k1 = nyh + 1
         fx(1,1,l) = zero
         fx(k1,1,l) = zero
         fy(1,1,l) = zero
         fy(k1,1,l) = zero
      endif
   80 continue
   90 continue
      we = float(nx)*float(ny)*wp
      return
c calculate potential and sum field energy
  100 if (isign.gt.1) go to 160
      wp = 0.0d0
      if (kstrt.gt.nxh) go to 150
      do 140 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 120 j = 1, kxp
      if ((j+joff).gt.0) then
         do 110 k = 2, nyh
         k1 = ny2 - k
         at2 = real(ffc(k,j,l))
         at1 = at2*aimag(ffc(k,j,l))
         fx(k,j,l) = at2*q(k,j,l)
         fx(k1,j,l) = at2*q(k1,j,l)
         wp = wp + at1*(q(k,j,l)*conjg(q(k,j,l)) + q(k1,j,l)*conjg(q(k1,
     1j,l)))
  110    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = real(ffc(1,j,l))
         at1 = at2*aimag(ffc(1,j,l))
         fx(1,j,l) = at2*q(1,j,l)
         fx(k1,j,l) = zero
         wp = wp + at1*(q(1,j,l)*conjg(q(1,j,l)))
      endif
  120 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 130 k = 2, nyh
         k1 = ny2 - k
         at2 = real(ffc(k,1,l))
         at1 = at2*aimag(ffc(k,1,l))
         fx(k,1,l) = at2*q(k,1,l)
         fx(k1,1,l) = zero
         wp = wp + at1*(q(k,1,l)*conjg(q(k,1,l)))
  130    continue
         k1 = nyh + 1
         fx(1,1,l) = zero
         fx(k1,1,l) = zero
      endif
  140 continue
  150 continue
      we = float(nx)*float(ny)*wp
      return
c calculate smoothing
  160 if (kstrt.gt.nxh) go to 210
      do 200 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 180 j = 1, kxp
      if ((j+joff).gt.0) then
         do 170 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,j,l))
         fy(k,j,l) = at1*q(k,j,l)
         fy(k1,j,l) = at1*q(k1,j,l)
  170    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = aimag(ffc(1,j,l))
         fy(1,j,l) = at1*q(1,j,l)
         fy(k1,j,l) = zero
      endif
  180 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 190 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,1,l))
         fy(k,1,l) = at1*q(k,1,l)
         fy(k1,1,l) = zero
  190    continue
         k1 = nyh + 1
         fy(1,1,l) = cmplx(aimag(ffc(1,1,l))*real(q(1,1,l)),0.)
         fy(k1,1,l) = zero
      endif
  200 continue
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISP21(q,fx,isign,ffc,ax,affp,we,nx,kstrt,kxp,jblok)
c this subroutine solves 1d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function,
c with periodic boundary conditions for distributed data.
c for isign = 0, input: isign,ax,affp,nx,kstrt,kxp,jblok, output: ffc
c for isign = -1, input: q,ffc,isign,nx,kstrt,kxp,jblok, output: fx,we
c approximate flop count is: 15*nxc
c for isign = 1, input: q,ffc,isign,nx,kstrt,kxp,jblok, output: fx,we
c approximate flop count is: 11*nxc
c for isign = 2, input: q,ffc,isign,nx,kstrt,kxp,jblok, output: fx
c approximate flop count is: 2*nxc
c where nxc = (nx/2-1)/nvp, and nvp = number of procs
c if isign < 0, force/charge is calculated using the equations:
c fx(kx) = -sqrt(-1)*kx*g(kx)*q(kx)*s(kx),
c where kx = 2pi*j/nx, and j = fourier mode number,
c g(kx) = (affp/(kx**2)*s(kx,ky),
c s(kx) = exp(-((kx*ax)**2)/2), except for
c fx(kx=pi) = 0, and fx(kx=0) = 0.
c if isign = 1, potential is calculated using the equation:
c fx(kx) = g(kx)*q(kx)
c if isign = 2, smoothing is calculated using the equation:
c fx(kx) = q(kx)*s(kx)
c q(j,l) = complex charge density for fourier mode (jj-1)
c fx(j,l) = x component of complex force/charge,
c for fourier mode (jj-1, where jj = j + kxp*(l - 1)
c jblok = number of data blocks
c kxp = number of data values per block
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c aimag(ffc(j,l)) = finite-size particle shape factor s
c real(ffc(j,l)) = potential green's function g
c for fourier mode (jj-1), where jj = j + kxp*(l - 1)
c ax = half-width of particle in x direction
c affp = normalization constant = nx/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*sum((affp/(kx**2))*|q(kx)*s(kx)|**2)
c nx = system length in x direction
      double precision wp
      complex q, fx, ffc, zero
      dimension q(kxp,jblok), fx(kxp,jblok), ffc(kxp,jblok)
      nxh = nx/2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 30
      if (kstrt.gt.nxh) return
c prepare form factor array
      do 20 l = 1, jblok
      joff = kxp*(l + ks) - 1
      do 10 j = 1, kxp
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = exp(-.5*((dkx*ax)**2))
      if (at1.eq.0.) then
         ffc(j,l) = cmplx(affp,1.)
      else
         ffc(j,l) = cmplx(affp*at2/at1,at2)
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 70
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nxh) go to 60
      do 50 l = 1, jblok
c mode numbers 0 < kx < nx/2
      joff = kxp*(l + ks) - 1
      do 40 j = 1, kxp
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         at1 = real(ffc(j,l))*aimag(ffc(j,l))
         fx(j,l) = dkx*at1*cmplx(aimag(q(j,l)),-real(q(j,l)))
         wp = wp + at1*(q(j,l)*conjg(q(j,l)))
      endif
   40 continue
      if ((l+ks).eq.0) fx(1,l) = zero
   50 continue
   60 continue
      we = float(nx)*wp
      return
c calculate potential and sum field energy
   70 if (isign.gt.1) go to 110
      wp = 0.0d0
      if (kstrt.gt.nxh) go to 100
      do 90 l = 1, jblok
c mode numbers 0 < kx < nx/2
      joff = kxp*(l + ks) - 1
      do 80 j = 1, kxp
      if ((j+joff).gt.0) then
         at2 = real(ffc(j,l))
         at1 = at2*aimag(ffc(j,l))
         fx(j,l) = at2*q(j,l)
         wp = wp + at1*(q(j,l)*conjg(q(j,l)))
      endif
   80 continue
      if ((l+ks).eq.0) fx(1,l) = zero
   90 continue
  100 continue
      we = float(nx)*wp
      return
c calculate smoothing
  110 if (kstrt.gt.nxh) go to 140
      do 130 l = 1, jblok
c mode numbers 0 < kx < nx/2
      joff = kxp*(l + ks) - 1
      do 120 j = 1, kxp
      if ((j+joff).gt.0) then
         at1 = aimag(ffc(j,l))
         fx(j,l) = at1*q(j,l)
      endif
  120 continue
      if ((l+ks).eq.0) then
         fx(1,l) = cmplx(aimag(ffc(1,l))*real(q(1,l)),0.)
      endif
  130 continue
  140 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,k
     1xp,jblok,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions, for distributed data.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: ffc
c for isign /= 0, input: q,ffc,isign,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: fxy,we
c approximate flop count is: 33*nxc*nyc + 15*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
c fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j,l) = x component of complex force/charge,
c fxy(2,k,j,l) = y component of complex force/charge,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c jblok = number of data blocks
c kxp = number of data values per block
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated.
c aimag(ffc(k,j,l)) = finite-size particle shape factor s
c real(ffc(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp
      complex q, fxy, ffc, zero, zt1, zt2
      dimension q(nyv,kxp,jblok), fxy(2,nyv,kxp,jblok)
      dimension ffc(nyhd,kxp,jblok)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nxh) return
c prepare form factor array
      do 30 l = 1, jblok
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffc(k,j,l) = cmplx(affp,1.)
      else
         ffc(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 wp = 0.0d0
      if (kstrt.gt.nxh) go to 90
      do 80 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 60 j = 1, kxp
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,j,l))*aimag(ffc(k,j,l))
         at2 = dkx*at1
         at3 = dny*float(k - 1)*at1
         zt1 = cmplx(aimag(q(k,j,l)),-real(q(k,j,l)))
         zt2 = cmplx(aimag(q(k1,j,l)),-real(q(k1,j,l)))
         fxy(1,k,j,l) = at2*zt1
         fxy(2,k,j,l) = at3*zt1
         fxy(1,k1,j,l) = at2*zt2
         fxy(2,k1,j,l) = -at3*zt2
         wp = wp + at1*(q(k,j,l)*conjg(q(k,j,l)) + q(k1,j,l)*conjg(q(k1,
     1j,l)))
   50    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = real(ffc(1,j,l))*aimag(ffc(1,j,l))
         at3 = dkx*at1
         zt1 = cmplx(aimag(q(1,j,l)),-real(q(1,j,l)))
         fxy(1,1,j,l) = at3*zt1
         fxy(2,1,j,l) = zero
         fxy(1,k1,j,l) = zero
         fxy(2,k1,j,l) = zero
         wp = wp + at1*(q(1,j,l)*conjg(q(1,j,l)))
      endif
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,1,l))*aimag(ffc(k,1,l))
         at2 = dny*float(k - 1)*at1
         zt1 = cmplx(aimag(q(k,1,l)),-real(q(k,1,l)))
         fxy(1,k,1,l) = zero
         fxy(2,k,1,l) = at2*zt1
         fxy(1,k1,1,l) = zero
         fxy(2,k1,1,l) = zero
         wp = wp + at1*(q(k,1,l)*conjg(q(k,1,l)))
   70    continue
         k1 = nyh + 1
         fxy(1,1,1,l) = zero
         fxy(2,1,1,l) = zero
         fxy(1,k1,1,l) = zero
         fxy(2,k1,1,l) = zero
      endif
   80 continue
   90 continue
      we = float(nx)*float(ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,k
     1xp,jblok,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions.  Zeros out z component.
c for distributed data.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: ffc
c for isign /= 0, input: q,ffc,isign,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: fxy,we
c approximate flop count is: 33*nxc*nyc + 15*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
c fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j,l) = x component of complex force/charge,
c fxy(2,k,j,l) = y component of complex force/charge,
c fxy(3,k,j,l) = zero,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c jblok = number of data blocks
c kxp = number of data values per block
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated.
c aimag(ffc(k,j,l)) = finite-size particle shape factor s
c real(ffc(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp
      complex q, fxy, ffc, zero, zt1, zt2
      dimension q(nyv,kxp,jblok), fxy(3,nyv,kxp,jblok)
      dimension ffc(nyhd,kxp,jblok)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nxh) return
c prepare form factor array
      do 30 l = 1, jblok
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffc(k,j,l) = cmplx(affp,1.)
      else
         ffc(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 wp = 0.0d0
      if (kstrt.gt.nxh) go to 90
      do 80 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 60 j = 1, kxp
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,j,l))*aimag(ffc(k,j,l))
         at2 = dkx*at1
         at3 = dny*float(k - 1)*at1
         zt1 = cmplx(aimag(q(k,j,l)),-real(q(k,j,l)))
         zt2 = cmplx(aimag(q(k1,j,l)),-real(q(k1,j,l)))
         fxy(1,k,j,l) = at2*zt1
         fxy(2,k,j,l) = at3*zt1
         fxy(3,k,j,l) = zero
         fxy(1,k1,j,l) = at2*zt2
         fxy(2,k1,j,l) = -at3*zt2
         fxy(3,k1,j,l) = zero
         wp = wp + at1*(q(k,j,l)*conjg(q(k,j,l)) + q(k1,j,l)*conjg(q(k1,
     1j,l)))
   50    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = real(ffc(1,j,l))*aimag(ffc(1,j,l))
         at3 = dkx*at1
         zt1 = cmplx(aimag(q(1,j,l)),-real(q(1,j,l)))
         fxy(1,1,j,l) = at3*zt1
         fxy(2,1,j,l) = zero
         fxy(3,1,j,l) = zero
         fxy(1,k1,j,l) = zero
         fxy(2,k1,j,l) = zero
         fxy(3,k1,j,l) = zero
         wp = wp + at1*(q(1,j,l)*conjg(q(1,j,l)))
      endif
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,1,l))*aimag(ffc(k,1,l))
         at2 = dny*float(k - 1)*at1
         zt1 = cmplx(aimag(q(k,1,l)),-real(q(k,1,l)))
         fxy(1,k,1,l) = zero
         fxy(2,k,1,l) = at2*zt1
         fxy(3,k,1,l) = zero
         fxy(1,k1,1,l) = zero
         fxy(2,k1,1,l) = zero
         fxy(3,k1,1,l) = zero
         wp = wp + at1*(q(k,1,l)*conjg(q(k,1,l)))
   70    continue
         k1 = nyh + 1
         fxy(1,1,1,l) = zero
         fxy(2,1,1,l) = zero
         fxy(3,1,1,l) = zero
         fxy(1,k1,1,l) = zero
         fxy(2,k1,1,l) = zero
         fxy(3,k1,1,l) = zero
      endif
   80 continue
   90 continue
      we = float(nx)*float(ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PDIVF2(f,df,nx,ny,kstrt,ndim,nyv,kxp,jblok)
c this subroutine calculates the divergence in fourier space
c input: all except df, output: df
c approximate flop count is: 16*nxc*nyc + 5*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the divergence is calculated using the equation:
c df(kx,ky) = sqrt(-1)*(kx*fx(kx,ky)+ky*fy(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for df(kx=pi) = df(ky=pi) = df(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c ndim = number of field arrays, must be >= 2
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c kxp = number of data values per block
c jblok = number of data blocks
      complex f, df, zero, zt1
      dimension f(ndim,nyv,kxp,jblok), df(nyv,kxp,jblok)
      if (ndim.lt.2) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate the divergence
      if (kstrt.gt.nxh) return
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         zt1 = dkx*f(1,k,j,l) + dky*f(2,k,j,l)
         df(k,j,l) = cmplx(-aimag(zt1),real(zt1))
         zt1 = dkx*f(1,k1,j,l) - dky*f(2,k1,j,l)
         df(k1,j,l) = cmplx(-aimag(zt1),real(zt1))
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         df(1,j,l) = dkx*cmplx(-aimag(f(1,1,j,l)),real(f(1,1,j,l)))
         df(k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         df(k,1,l) = dky*cmplx(-aimag(f(2,k,1,l)),real(f(2,k,1,l)))
         df(k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         df(1,1,l) = zero
         df(k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGRADF2(df,f,nx,ny,kstrt,ndim,nyv,kxp,jblok)
c this subroutine calculates the gradient in fourier space
c input: all except f, output: f
c approximate flop count is: 12*nxc*nyc + 4*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the gradient is calculated using the equations:
c fx(kx,ky) = sqrt(-1)*kx*df(kx,ky)
c fy(kx,ky) = sqrt(-1)*ky*df(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for fx(kx=pi) = fy(kx=pi) = 0, fx(ky=pi) = fy(ky=pi) = 0,
c and fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c ndim = number of field arrays, must be >= 2
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c kxp = number of data values per block
c jblok = number of data blocks
      complex df, f, zero, zt1
      dimension df(nyv,kxp,jblok), f(ndim,nyv,kxp,jblok)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate the gradient
      if (kstrt.gt.nxh) return
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         zt1 = cmplx(-aimag(df(k,j,l)),real(df(k,j,l)))
         f(1,k,j,l) = dkx*zt1
         f(2,k,j,l) = dky*zt1
         zt1 = cmplx(-aimag(df(k1,j,l)),real(df(k1,j,l)))
         f(1,k1,j,l) = dkx*zt1
         f(2,k1,j,l) = -dky*zt1
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         f(1,1,j,l) = dkx*cmplx(-aimag(df(1,j,l)),real(df(1,j,l)))
         f(2,1,j,l) = zero
         f(1,k1,j,l) = zero
         f(2,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         f(1,k,1,l) = zero
         f(2,k,1,l) = dky*cmplx(-aimag(df(k,1,l)),real(df(k,1,l)))
         f(1,k1,1,l) = zero
         f(2,k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         f(1,1,1,l) = zero
         f(2,1,1,l) = zero
         f(1,k1,1,l) = zero
         f(2,k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCURLF2(f,g,nx,ny,kstrt,nyv,kxp,jblok)
c this subroutine calculates the curl in fourier space
c input: all except g, output: g
c approximate flop count is: 32*nxc*nyc + 10*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the curl is calculated using the equations:
c gx(kx,ky) = sqrt(-1)*ky*fz(kx,ky)
c gy(kx,ky) = -sqrt(-1)*kx*fz(kx,ky)
c gz(kx,ky) = sqrt(-1)*(kx*fy(kx,ky)-ky*fx(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for gx(kx=pi) = gy(kx=pi) = 0, gx(ky=pi) = gy(ky=pi) = 0,
c and gx(kx=0,ky=0) = gy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c kxp = number of data values per block
c jblok = number of data blocks
      complex f, g, zero, zt1, zt2, zt3
      dimension f(3,nyv,kxp,jblok), g(3,nyv,kxp,jblok)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate the curl
      if (kstrt.gt.nxh) return
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         zt1 = cmplx(-aimag(f(3,k,j,l)),real(f(3,k,j,l)))
         zt2 = cmplx(-aimag(f(2,k,j,l)),real(f(2,k,j,l)))
         zt3 = cmplx(-aimag(f(1,k,j,l)),real(f(1,k,j,l)))
         g(1,k,j,l) = dky*zt1
         g(2,k,j,l) = -dkx*zt1
         g(3,k,j,l) = dkx*zt2 - dky*zt3
         zt1 = cmplx(-aimag(f(3,k1,j,l)),real(f(3,k1,j,l)))
         zt2 = cmplx(-aimag(f(2,k1,j,l)),real(f(2,k1,j,l)))
         zt3 = cmplx(-aimag(f(1,k1,j,l)),real(f(1,k1,j,l)))
         g(1,k1,j,l) = -dky*zt1
         g(2,k1,j,l) = -dkx*zt1
         g(3,k1,j,l) = dkx*zt2 + dky*zt3
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         zt1 = cmplx(-aimag(f(3,1,j,l)),real(f(3,1,j,l)))
         zt2 = cmplx(-aimag(f(2,1,j,l)),real(f(2,1,j,l)))
         g(1,1,j,l) = zero
         g(2,1,j,l) = -dkx*zt1
         g(3,1,j,l) = dkx*zt2
         g(1,k1,j,l) = zero
         g(2,k1,j,l) = zero
         g(3,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         zt1 = cmplx(-aimag(f(3,k,1,l)),real(f(3,k,1,l)))
         zt2 = cmplx(-aimag(f(1,k,1,l)),real(f(1,k,1,l)))
         g(1,k,1,l) = dky*zt1
         g(2,k,1,l) = zero
         g(3,k,1,l) = -dky*zt2
         g(1,k1,1,l) = zero
         g(2,k1,1,l) = zero
         g(3,k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         g(1,1,1,l) = zero
         g(2,1,1,l) = zero
         g(3,1,1,l) = zero
         g(1,k1,1,l) = zero
         g(2,k1,1,l) = zero
         g(3,k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCURLF22(f,g,nx,ny,kstrt,nyv,kxp,jblok)
c this subroutine calculates the curl in fourier space
c input: all except g, output: g
c approximate flop count is: 32*nxc*nyc + 10*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the curl is calculated using the equations:
c g(kx,ky) = sqrt(-1)*(kx*fy(kx,ky)-ky*fx(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c kxp = number of data values per block
c jblok = number of data blocks
      complex f, g, zero, zt2, zt3
      dimension f(2,nyv,kxp,jblok), g(nyv,kxp,jblok)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate the curl
      if (kstrt.gt.nxh) return
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         zt2 = cmplx(-aimag(f(2,k,j,l)),real(f(2,k,j,l)))
         zt3 = cmplx(-aimag(f(1,k,j,l)),real(f(1,k,j,l)))
         g(k,j,l) = dkx*zt2 - dky*zt3
         zt2 = cmplx(-aimag(f(2,k1,j,l)),real(f(2,k1,j,l)))
         zt3 = cmplx(-aimag(f(1,k1,j,l)),real(f(1,k1,j,l)))
         g(k1,j,l) = dkx*zt2 + dky*zt3
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         zt2 = cmplx(-aimag(f(2,1,j,l)),real(f(2,1,j,l)))
         g(1,j,l) = dkx*zt2
         g(k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         zt2 = cmplx(-aimag(f(1,k,1,l)),real(f(1,k,1,l)))
         g(k,1,l) = -dky*zt2
         g(k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         g(1,1,l) = zero
         g(k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCUPERP2(cu,nx,ny,kstrt,nyv,kxp,jblok)
c this subroutine calculates the transverse current in fourier space
c input: all, output: cu
c approximate flop count is: 36*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for cux(kx=pi) = cuy(kx=pi) = 0, cux(ky=pi) = cuy(ky=pi) = 0,
c and cux(kx=0,ky=0) = cuy(kx=0,ky=0) = 0.
c cu(i,k,j,l) = i-th component of complex current density and
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c kxp = number of data values per block
c jblok = number of data blocks
      complex cu, zero, zt1
      dimension cu(3,nyv,kxp,jblok)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
      if (kstrt.gt.nxh) return
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         zt1 = at1*(dkx*cu(1,k,j,l) + dky*cu(2,k,j,l))
         cu(1,k,j,l) = cu(1,k,j,l) - dkx*zt1
         cu(2,k,j,l) = cu(2,k,j,l) - dky*zt1
         zt1 = at1*(dkx*cu(1,k1,j,l) - dky*cu(2,k1,j,l))
         cu(1,k1,j,l) = cu(1,k1,j,l) - dkx*zt1
         cu(2,k1,j,l) = cu(2,k1,j,l) + dky*zt1
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         cu(1,1,j,l) = zero
         cu(1,k1,j,l) = zero
         cu(2,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         cu(2,k,1,l) = zero
         cu(1,k1,1,l) = zero
         cu(2,k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         cu(1,1,1,l) = zero
         cu(2,1,1,l) = zero
         cu(1,k1,1,l) = zero
         cu(2,k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCUPERP22(cu,nx,ny,kstrt,nyv,kxp,jblok)
c this subroutine calculates the transverse current in fourier space
c input: all, output: cu
c approximate flop count is: 36*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for cux(kx=pi) = cuy(kx=pi) = 0, cux(ky=pi) = cuy(ky=pi) = 0,
c and cux(kx=0,ky=0) = cuy(kx=0,ky=0) = 0.
c cu(i,k,j,l) = i-th component of complex current density and
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c kxp = number of data values per block
c jblok = number of data blocks
      complex cu, zero, zt1
      dimension cu(2,nyv,kxp,jblok)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
      if (kstrt.gt.nxh) return
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         zt1 = at1*(dkx*cu(1,k,j,l) + dky*cu(2,k,j,l))
         cu(1,k,j,l) = cu(1,k,j,l) - dkx*zt1
         cu(2,k,j,l) = cu(2,k,j,l) - dky*zt1
         zt1 = at1*(dkx*cu(1,k1,j,l) - dky*cu(2,k1,j,l))
         cu(1,k1,j,l) = cu(1,k1,j,l) - dkx*zt1
         cu(2,k1,j,l) = cu(2,k1,j,l) + dky*zt1
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         cu(1,1,j,l) = zero
         cu(1,k1,j,l) = zero
         cu(2,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         cu(2,k,1,l) = zero
         cu(1,k1,1,l) = zero
         cu(2,k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         cu(1,1,1,l) = zero
         cu(2,1,1,l) = zero
         cu(1,k1,1,l) = zero
         cu(2,k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPOISP23(cu,bxy,isign,ffc,ax,ay,affp,ci,wm,nx,ny,kstrt
     1,nyv,kxp,jblok,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with periodic boundary conditions for distributed data.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: ffc
c for isign = -1, input: cu,ffc,isign,ci,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: bxy,wm
c approximate flop count is: 85*nxc*nyc + 36*(nxc + nyc)
c for isign = 1, input: cu,ffc,isign,ci,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: bxy,wm
c approximate flop count is: 63*nxc*nyc + 33*(nxc + nyc)
c for isign = 2, input: cu,ffc,isign,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: bxy
c approximate flop count is: 12*nxc*nyc + 6*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign < 0, magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky)*s(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c             s(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
c bx(ky=pi) = by(ky=pi) = bz(ky=pi) = 0,
c bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
c by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
c bz(kx,ky) = ci*ci*g(kx,ky)*cuz(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky) = cux(kx,ky)*s(kx,ky)
c by(kx,ky) = cuy(kx,ky)*s(kx,ky)
c bz(kx,ky) = cuz(kx,ky)*s(kx,ky)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c jblok = number of data blocks
c kxp = number of data values per block
c kstrt = starting data block number
c aimag(ffc(k,j,l)) = finite-size particle shape factor s
c real(ffc(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, ffc, zero, zt1, zt2, zt3
      dimension cu(3,nyv,kxp,jblok)
      dimension bxy(3,nyv,kxp,jblok)
      dimension ffc(nyhd,kxp,jblok)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 40
      if (kstrt.gt.nxh) return
c prepare form factor array
      do 30 l = 1, jblok
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffc(k,j,l) = cmplx(affp,1.)
      else
         ffc(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 100
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nxh) go to 90
      do 80 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 60 j = 1, kxp
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffc(k,j,l))*aimag(ffc(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         zt1 = cmplx(-aimag(cu(3,k,j,l)),real(cu(3,k,j,l)))
         zt2 = cmplx(-aimag(cu(2,k,j,l)),real(cu(2,k,j,l)))
         zt3 = cmplx(-aimag(cu(1,k,j,l)),real(cu(1,k,j,l)))
         bxy(1,k,j,l) = at2*zt1
         bxy(2,k,j,l) = -at3*zt1
         bxy(3,k,j,l) = at3*zt2 - at2*zt3
         zt1 = cmplx(-aimag(cu(3,k1,j,l)),real(cu(3,k1,j,l)))
         zt2 = cmplx(-aimag(cu(2,k1,j,l)),real(cu(2,k1,j,l)))
         zt3 = cmplx(-aimag(cu(1,k1,j,l)),real(cu(1,k1,j,l)))
         bxy(1,k1,j,l) = -at2*zt1
         bxy(2,k1,j,l) = -at3*zt1
         bxy(3,k1,j,l) = at3*zt2 + at2*zt3
         wp = wp + at1*(cu(1,k,j,l)*conjg(cu(1,k,j,l)) + cu(2,k,j,l)*con
     1jg(cu(2,k,j,l)) + cu(3,k,j,l)*conjg(cu(3,k,j,l)) + cu(1,k1,j,l)*co
     2njg(cu(1,k1,j,l)) + cu(2,k1,j,l)*conjg(cu(2,k1,j,l)) + cu(3,k1,j,l
     3)*conjg(cu(3,k1,j,l)))
   50    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = ci2*real(ffc(1,j,l))*aimag(ffc(1,j,l))
         at2 = dkx*at1
         zt1 = cmplx(-aimag(cu(3,1,j,l)),real(cu(3,1,j,l)))
         zt2 = cmplx(-aimag(cu(2,1,j,l)),real(cu(2,1,j,l)))
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = -at2*zt1
         bxy(3,1,j,l) = at2*zt2
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
         wp = wp + at1*(cu(1,1,j,l)*conjg(cu(1,1,j,l)) + cu(2,1,j,l)*con
     1jg(cu(2,1,j,l)) + cu(3,1,j,l)*conjg(cu(3,1,j,l)))
      endif
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffc(k,1,l))*aimag(ffc(k,1,l))
         at2 = dky*at1
         zt1 = cmplx(-aimag(cu(3,k,1,l)),real(cu(3,k,1,l)))
         zt2 = cmplx(-aimag(cu(1,k,1,l)),real(cu(1,k,1,l)))
         bxy(1,k,1,l) = at2*zt1
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = -at2*zt2
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         wp = wp + at1*(cu(1,k,1,l)*conjg(cu(1,k,1,l)) + cu(2,k,1,l)*con
     1jg(cu(2,k,1,l)) + cu(3,k,1,l)*conjg(cu(3,k,1,l)))
   70    continue
         k1 = nyh + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
      endif
   80 continue
   90 continue
      wm = float(nx)*float(ny)*wp
      return
c calculate vector potential and sum field energy
  100 if (isign.gt.1) go to 160
      wp = 0.0d0
      if (kstrt.gt.nxh) go to 150
      do 140 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 120 j = 1, kxp
      if ((j+joff).gt.0) then
         do 110 k = 2, nyh
         k1 = ny2 - k
         at2 = ci2*real(ffc(k,j,l))
         at1 = at2*aimag(ffc(k,j,l))
         bxy(1,k,j,l) = at2*cu(1,k,j,l)
         bxy(2,k,j,l) = at2*cu(2,k,j,l)
         bxy(3,k,j,l) = at2*cu(3,k,j,l)
         bxy(1,k1,j,l) = at2*cu(1,k1,j,l)
         bxy(2,k1,j,l) = at2*cu(2,k1,j,l)
         bxy(3,k1,j,l) = at2*cu(3,k1,j,l)
         wp = wp + at1*(cu(1,k,j,l)*conjg(cu(1,k,j,l)) + cu(2,k,j,l)*con
     1jg(cu(2,k,j,l)) + cu(3,k,j,l)*conjg(cu(3,k,j,l)) + cu(1,k1,j,l)*co
     2njg(cu(1,k1,j,l)) + cu(2,k1,j,l)*conjg(cu(2,k1,j,l)) + cu(3,k1,j,l
     3)*conjg(cu(3,k1,j,l)))
  110    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = ci2*real(ffc(1,j,l))
         at1 = at2*aimag(ffc(1,j,l))
         bxy(1,1,j,l) = at2*cu(1,1,j,l)
         bxy(2,1,j,l) = at2*cu(2,1,j,l)
         bxy(3,1,j,l) = at2*cu(3,1,j,l)
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
         wp = wp + at1*(cu(1,1,j,l)*conjg(cu(1,1,j,l)) + cu(2,1,j,l)*con
     1jg(cu(2,1,j,l)) + cu(3,1,j,l)*conjg(cu(3,1,j,l)))
      endif
  120 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 130 k = 2, nyh
         k1 = ny2 - k
         at2 = ci2*real(ffc(k,1,l))
         at1 = at2*aimag(ffc(k,1,l))
         bxy(1,k,1,l) = at2*cu(1,k,1,l)
         bxy(2,k,1,l) = at2*cu(2,k,1,l)
         bxy(3,k,1,l) = at2*cu(3,k,1,l)
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         wp = wp + at1*(cu(1,k,1,l)*conjg(cu(1,k,1,l)) + cu(2,k,1,l)*con
     1jg(cu(2,k,1,l)) + cu(3,k,1,l)*conjg(cu(3,k,1,l)))
  130    continue
         k1 = nyh + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
      endif
  140 continue
  150 continue
      wm = float(nx)*float(ny)*wp
      return
c calculate smoothing
  160 if (kstrt.gt.nxh) go to 210
      do 200 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 180 j = 1, kxp
      if ((j+joff).gt.0) then
         do 170 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,j,l))
         bxy(1,k,j,l) = at1*cu(1,k,j,l)
         bxy(2,k,j,l) = at1*cu(2,k,j,l)
         bxy(3,k,j,l) = at1*cu(3,k,j,l)
         bxy(1,k1,j,l) = at1*cu(1,k1,j,l)
         bxy(2,k1,j,l) = at1*cu(2,k1,j,l)
         bxy(3,k1,j,l) = at1*cu(3,k1,j,l)
  170    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = aimag(ffc(1,j,l))
         bxy(1,1,j,l) = at1*cu(1,1,j,l)
         bxy(2,1,j,l) = at1*cu(2,1,j,l)
         bxy(3,1,j,l) = at1*cu(3,1,j,l)
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
      endif
  180 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 190 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,1,l))
         bxy(1,k,1,l) = at1*cu(1,k,1,l)
         bxy(2,k,1,l) = at1*cu(2,k,1,l)
         bxy(3,k,1,l) = at1*cu(3,k,1,l)
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
  190    continue
         k1 = nyh + 1
         at1 = aimag(ffc(1,1,l))
         bxy(1,1,1,l) = cmplx(at1*real(cu(1,1,1,l)),0.)
         bxy(2,1,1,l) = cmplx(at1*real(cu(2,1,1,l)),0.)
         bxy(3,1,1,l) = cmplx(at1*real(cu(3,1,1,l)),0.)
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
      endif
  200 continue
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPOISP22(cu,bxy,bz,isign,ffc,ax,ay,affp,ci,wm,nx,ny,ks
     1trt,nyv,kxp,jblok,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with periodic boundary conditions for distributed data.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: ffc
c for isign = -1, input: cu,ffc,isign,ci,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: bz,wm
c approximate flop count is: 55*nxc*nyc + 24*(nxc + nyc)
c for isign = 1, input: cu,ffc,isign,ci,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: bxy,wm
c approximate flop count is: 45*nxc*nyc + 24*(nxc + nyc)
c for isign = 2, input: cu,ffc,isign,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: bxy
c approximate flop count is: 8*nxc*nyc + 4*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign < 0, magnetic field is calculated using the equations:
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c             s(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bz(kx=pi) = 0, bz(ky=pi) = 0, bz(kx=0,ky=0) = 0.
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
c by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky) = cux(kx,ky)*s(kx,ky)
c by(kx,ky) = cuy(kx,ky)*s(kx,ky)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex vector potential,
c bz(k,j,l) = z component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c jblok = number of data blocks
c kxp = number of data values per block
c kstrt = starting data block number
c aimag(ffc(k,j,l)) = finite-size particle shape factor s
c real(ffc(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, bz, ffc, zero, zt2, zt3
      dimension cu(2,nyv,kxp,jblok)
      dimension bxy(2,nyv,kxp,jblok), bz(nyv,kxp,jblok)
      dimension ffc(nyhd,kxp,jblok)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 40
      if (kstrt.gt.nxh) return
c prepare form factor array
      do 30 l = 1, jblok
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffc(k,j,l) = cmplx(affp,1.)
      else
         ffc(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 100
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nxh) go to 90
      do 80 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 60 j = 1, kxp
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffc(k,j,l))*aimag(ffc(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         zt2 = cmplx(-aimag(cu(2,k,j,l)),real(cu(2,k,j,l)))
         zt3 = cmplx(-aimag(cu(1,k,j,l)),real(cu(1,k,j,l)))
         bz(k,j,l) = at3*zt2 - at2*zt3
         zt2 = cmplx(-aimag(cu(2,k1,j,l)),real(cu(2,k1,j,l)))
         zt3 = cmplx(-aimag(cu(1,k1,j,l)),real(cu(1,k1,j,l)))
         bz(k1,j,l) = at3*zt2 + at2*zt3
         wp = wp + at1*(cu(1,k,j,l)*conjg(cu(1,k,j,l)) + cu(2,k,j,l)*con
     1jg(cu(2,k,j,l)) + cu(1,k1,j,l)*conjg(cu(1,k1,j,l)) + cu(2,k1,j,l)*
     2conjg(cu(2,k1,j,l)))
   50    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = ci2*real(ffc(1,j,l))*aimag(ffc(1,j,l))
         at2 = dkx*at1
         zt2 = cmplx(-aimag(cu(2,1,j,l)),real(cu(2,1,j,l)))
         bz(1,j,l) = at2*zt2
         bz(k1,j,l) = zero
         wp = wp + at1*(cu(1,1,j,l)*conjg(cu(1,1,j,l)) + cu(2,1,j,l)*con
     1jg(cu(2,1,j,l)))
      endif
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffc(k,1,l))*aimag(ffc(k,1,l))
         at2 = dky*at1
         zt2 = cmplx(-aimag(cu(1,k,1,l)),real(cu(1,k,1,l)))
         bz(k,1,l) = -at2*zt2
         bz(k1,1,l) = zero
         wp = wp + at1*(cu(1,k,1,l)*conjg(cu(1,k,1,l)) + cu(2,k,1,l)*con
     1jg(cu(2,k,1,l)))
   70    continue
         k1 = nyh + 1
         bz(1,1,l) = zero
         bz(k1,1,l) = zero
      endif
   80 continue
   90 continue
      wm = float(nx)*float(ny)*wp
      return
c calculate vector potential and sum field energy
  100 if (isign.gt.1) go to 160
      wp = 0.0d0
      if (kstrt.gt.nxh) go to 150
      do 140 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 120 j = 1, kxp
      if ((j+joff).gt.0) then
         do 110 k = 2, nyh
         k1 = ny2 - k
         at2 = ci2*real(ffc(k,j,l))
         at1 = at2*aimag(ffc(k,j,l))
         bxy(1,k,j,l) = at2*cu(1,k,j,l)
         bxy(2,k,j,l) = at2*cu(2,k,j,l)
         bxy(1,k1,j,l) = at2*cu(1,k1,j,l)
         bxy(2,k1,j,l) = at2*cu(2,k1,j,l)
         wp = wp + at1*(cu(1,k,j,l)*conjg(cu(1,k,j,l)) + cu(2,k,j,l)*con
     1jg(cu(2,k,j,l)) + cu(1,k1,j,l)*conjg(cu(1,k1,j,l)) + cu(2,k1,j,l)*
     2conjg(cu(2,k1,j,l)))
  110    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = ci2*real(ffc(1,j,l))
         at1 = at2*aimag(ffc(1,j,l))
         bxy(1,1,j,l) = at2*cu(1,1,j,l)
         bxy(2,1,j,l) = at2*cu(2,1,j,l)
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         wp = wp + at1*(cu(1,1,j,l)*conjg(cu(1,1,j,l)) + cu(2,1,j,l)*con
     1jg(cu(2,1,j,l)))
      endif
  120 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 130 k = 2, nyh
         k1 = ny2 - k
         at2 = ci2*real(ffc(k,1,l))
         at1 = at2*aimag(ffc(k,1,l))
         bxy(1,k,1,l) = at2*cu(1,k,1,l)
         bxy(2,k,1,l) = at2*cu(2,k,1,l)
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         wp = wp + at1*(cu(1,k,1,l)*conjg(cu(1,k,1,l)) + cu(2,k,1,l)*con
     1jg(cu(2,k,1,l)))
  130    continue
         k1 = nyh + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
      endif
  140 continue
  150 continue
      wm = float(nx)*float(ny)*wp
      return
c calculate smoothing
  160 if (kstrt.gt.nxh) go to 210
      do 200 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 180 j = 1, kxp
      if ((j+joff).gt.0) then
         do 170 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,j,l))
         bxy(1,k,j,l) = at1*cu(1,k,j,l)
         bxy(2,k,j,l) = at1*cu(2,k,j,l)
         bxy(1,k1,j,l) = at1*cu(1,k1,j,l)
         bxy(2,k1,j,l) = at1*cu(2,k1,j,l)
  170    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = aimag(ffc(1,j,l))
         bxy(1,1,j,l) = at1*cu(1,1,j,l)
         bxy(2,1,j,l) = at1*cu(2,1,j,l)
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
      endif
  180 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 190 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,1,l))
         bxy(1,k,1,l) = at1*cu(1,k,1,l)
         bxy(2,k,1,l) = at1*cu(2,k,1,l)
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
  190    continue
         k1 = nyh + 1
         at1 = aimag(ffc(1,1,l))
         bxy(1,1,1,l) = cmplx(at1*real(cu(1,1,1,l)),0.)
         bxy(2,1,1,l) = cmplx(at1*real(cu(2,1,1,l)),0.)
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
      endif
  200 continue
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine IPBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,jblok,n
     1yhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field with periodic boundary conditions for distributed data.
c input: cu,ffc,ci,nx,ny,kstrt,nyv,kxp,jblok,nyhd, output: bxy,wm
c approximate flop count is: 85*nxc*nyc + 36*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
c bx(ky=pi) = by(ky=pi) = bz(ky=pi) = 0,
c bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c jblok = number of data blocks
c kxp = number of data values per block
c kstrt = starting data block number
c aimag(ffc(k,j,l)) = finite-size particle shape factor s
c real(ffc(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, ffc, zero, zt1, zt2, zt3
      dimension cu(3,nyv,kxp,jblok)
      dimension bxy(3,nyv,kxp,jblok)
      dimension ffc(nyhd,kxp,jblok)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nxh) go to 50
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffc(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         at1 = at1*aimag(ffc(k,j,l))
         zt1 = cmplx(-aimag(cu(3,k,j,l)),real(cu(3,k,j,l)))
         zt2 = cmplx(-aimag(cu(2,k,j,l)),real(cu(2,k,j,l)))
         zt3 = cmplx(-aimag(cu(1,k,j,l)),real(cu(1,k,j,l)))
         bxy(1,k,j,l) = at2*zt1
         bxy(2,k,j,l) = -at3*zt1
         bxy(3,k,j,l) = at3*zt2 - at2*zt3
         zt1 = cmplx(-aimag(cu(3,k1,j,l)),real(cu(3,k1,j,l)))
         zt2 = cmplx(-aimag(cu(2,k1,j,l)),real(cu(2,k1,j,l)))
         zt3 = cmplx(-aimag(cu(1,k1,j,l)),real(cu(1,k1,j,l)))
         bxy(1,k1,j,l) = -at2*zt1
         bxy(2,k1,j,l) = -at3*zt1
         bxy(3,k1,j,l) = at3*zt2 + at2*zt3
         wp = wp + at1*(cu(1,k,j,l)*conjg(cu(1,k,j,l)) + cu(2,k,j,l)*con
     1jg(cu(2,k,j,l)) + cu(3,k,j,l)*conjg(cu(3,k,j,l)) + cu(1,k1,j,l)*co
     2njg(cu(1,k1,j,l)) + cu(2,k1,j,l)*conjg(cu(2,k1,j,l)) + cu(3,k1,j,l
     3)*conjg(cu(3,k1,j,l)))
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = ci2*real(ffc(1,j,l))
         at2 = dkx*at1
         at1 = at1*aimag(ffc(1,j,l))
         zt1 = cmplx(-aimag(cu(3,1,j,l)),real(cu(3,1,j,l)))
         zt2 = cmplx(-aimag(cu(2,1,j,l)),real(cu(2,1,j,l)))
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = -at2*zt1
         bxy(3,1,j,l) = at2*zt2
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
         wp = wp + at1*(cu(1,1,j,l)*conjg(cu(1,1,j,l)) + cu(2,1,j,l)*con
     1jg(cu(2,1,j,l)) + cu(3,1,j,l)*conjg(cu(3,1,j,l)))
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffc(k,1,l))
         at2 = dky*at1
         at1 = at1*aimag(ffc(k,1,l))
         zt1 = cmplx(-aimag(cu(3,k,1,l)),real(cu(3,k,1,l)))
         zt2 = cmplx(-aimag(cu(1,k,1,l)),real(cu(1,k,1,l)))
         bxy(1,k,1,l) = at2*zt1
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = -at2*zt2
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         wp = wp + at1*(cu(1,k,1,l)*conjg(cu(1,k,1,l)) + cu(2,k,1,l)*con
     1jg(cu(2,k,1,l)) + cu(3,k,1,l)*conjg(cu(3,k,1,l)))
   30    continue
         k1 = nyh + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
      endif
   40 continue
   50 continue
      wm = float(nx)*float(ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PMAXWEL2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,nx,ny,kstrt,ny
     1v,kxp,jblok,nyhd)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with periodic boundary
c conditions.
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 286*nxc*nyc + 84*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the magnetic field is first updated half a step using the equations:
c bx(kx,ky) = bx(kx,ky) - .5*dt*sqrt(-1)*ky*ez(kx,ky)
c by(kx,ky) = by(kx,ky) + .5*dt*sqrt(-1)*kx*ez(kx,ky)
c bz(kx,ky) = bz(kx,ky) - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
c the electric field is then updated a whole step using the equations:
c ex(kx,ky) = ex(kx,ky) + c2*dt*sqrt(-1)*ky*bz(kx,ky)
c                       - affp*dt*cux(kx,ky)*s(kx,ky)
c ey(kx,ky) = ey(kx,ky) - c2*dt*sqrt(-1)*kx*bz(kx,ky)
c                       - affp*dt*cuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = ez(kx,ky) + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
c                       - affp*dt*cuz(kx,ky)*s(kx,ky)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers, except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
c ex(ky=pi) = ey(ky=pi) = ez(ky=pi) = 0,
c ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
c and similarly for bx, by, bz.
c cu(i,k,j,l) = i-th component of complex current density and
c exy(i,k,j,l) = i-th component of complex electric field,
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c aimag(ffc(k,j,l)) = finite-size particle shape factor s
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny*nz**sum((1/affp)*|exyz(kx,ky,kz)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny*nz**sum((c2/affp)*|bxyz(kx,ky,kz)|**2)
c nx/ny = system length in x/y direction
c jblok = number of data blocks
c kxp = number of data values per block
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp, ws
      complex exy, bxy, cu, ffc
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      dimension exy(3,nyv,kxp,jblok), bxy(3,nyv,kxp,jblok)
      dimension cu(3,nyv,kxp,jblok)
      dimension ffc(nyhd,kxp,jblok)
      if (ci.le.0.) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
      if (kstrt.gt.nxh) go to 50
c calculate the electromagnetic fields
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         afdt = adt*aimag(ffc(k,j,l))
c update magnetic field half time step, ky > 0
         zt1 = cmplx(-aimag(exy(3,k,j,l)),real(exy(3,k,j,l)))
         zt2 = cmplx(-aimag(exy(2,k,j,l)),real(exy(2,k,j,l)))
         zt3 = cmplx(-aimag(exy(1,k,j,l)),real(exy(1,k,j,l)))
         zt4 = bxy(1,k,j,l) - dth*(dky*zt1)
         zt5 = bxy(2,k,j,l) + dth*(dkx*zt1)
         zt6 = bxy(3,k,j,l) - dth*(dkx*zt2 - dky*zt3)
c update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt2 = cmplx(-aimag(zt5),real(zt5))
         zt3 = cmplx(-aimag(zt4),real(zt4))
         zt7 = exy(1,k,j,l) + cdt*(dky*zt1) - afdt*cu(1,k,j,l)
         zt8 = exy(2,k,j,l) - cdt*(dkx*zt1) - afdt*cu(2,k,j,l)
         zt9 = exy(3,k,j,l) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,k,j,l)
c update magnetic field half time step and store electric field
         zt1 = cmplx(-aimag(zt9),real(zt9))
         zt2 = cmplx(-aimag(zt8),real(zt8))
         zt3 = cmplx(-aimag(zt7),real(zt7))
         exy(1,k,j,l) = zt7
         exy(2,k,j,l) = zt8
         exy(3,k,j,l) = zt9
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt
     19))
         zt4 = zt4 - dth*(dky*zt1)
         zt5 = zt5 + dth*(dkx*zt1)
         zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
         bxy(1,k,j,l) = zt4
         bxy(2,k,j,l) = zt5
         bxy(3,k,j,l) = zt6
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt
     16))
c update magnetic field half time step, ky < 0
         zt1 = cmplx(-aimag(exy(3,k1,j,l)),real(exy(3,k1,j,l)))
         zt2 = cmplx(-aimag(exy(2,k1,j,l)),real(exy(2,k1,j,l)))
         zt3 = cmplx(-aimag(exy(1,k1,j,l)),real(exy(1,k1,j,l)))
         zt4 = bxy(1,k1,j,l) + dth*(dky*zt1)
         zt5 = bxy(2,k1,j,l) + dth*(dkx*zt1)
         zt6 = bxy(3,k1,j,l) - dth*(dkx*zt2 + dky*zt3)
c update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt2 = cmplx(-aimag(zt5),real(zt5))
         zt3 = cmplx(-aimag(zt4),real(zt4))
         zt7 = exy(1,k1,j,l) - cdt*(dky*zt1) - afdt*cu(1,k1,j,l)
         zt8 = exy(2,k1,j,l) - cdt*(dkx*zt1) - afdt*cu(2,k1,j,l)
         zt9 = exy(3,k1,j,l) + cdt*(dkx*zt2 + dky*zt3) - afdt*cu(3,k1,j,
     1l)
c update magnetic field half time step and store electric field
         zt1 = cmplx(-aimag(zt9),real(zt9))
         zt2 = cmplx(-aimag(zt8),real(zt8))
         zt3 = cmplx(-aimag(zt7),real(zt7))
         exy(1,k1,j,l) = zt7
         exy(2,k1,j,l) = zt8
         exy(3,k1,j,l) = zt9
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt
     19))
         zt4 = zt4 + dth*(dky*zt1)
         zt5 = zt5 + dth*(dkx*zt1)
         zt6 = zt6 - dth*(dkx*zt2 + dky*zt3)
         bxy(1,k1,j,l) = zt4
         bxy(2,k1,j,l) = zt5
         bxy(3,k1,j,l) = zt6
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt
     16))
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         afdt = adt*aimag(ffc(1,j,l))
c update magnetic field half time step
         zt1 = cmplx(-aimag(exy(3,1,j,l)),real(exy(3,1,j,l)))
         zt2 = cmplx(-aimag(exy(2,1,j,l)),real(exy(2,1,j,l)))
         zt5 = bxy(2,1,j,l) + dth*(dkx*zt1)
         zt6 = bxy(3,1,j,l) - dth*(dkx*zt2)
c update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt2 = cmplx(-aimag(zt5),real(zt5))
         zt8 = exy(2,1,j,l) - cdt*(dkx*zt1) - afdt*cu(2,1,j,l)
         zt9 = exy(3,1,j,l) + cdt*(dkx*zt2) - afdt*cu(3,1,j,l)
c update magnetic field half time step and store electric field
         zt1 = cmplx(-aimag(zt9),real(zt9))
         zt2 = cmplx(-aimag(zt8),real(zt8))
         exy(1,1,j,l) = zero
         exy(2,1,j,l) = zt8
         exy(3,1,j,l) = zt9
         ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
         zt5 = zt5 + dth*(dkx*zt1)
         zt6 = zt6 - dth*(dkx*zt2)
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = zt5
         bxy(3,1,j,l) = zt6
         wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
         exy(1,k1,j,l) = zero
         exy(2,k1,j,l) = zero
         exy(3,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         afdt = adt*aimag(ffc(k,1,l))
c update magnetic field half time step
         zt1 = cmplx(-aimag(exy(3,k,1,l)),real(exy(3,k,1,l)))
         zt3 = cmplx(-aimag(exy(1,k,1,l)),real(exy(1,k,1,l)))
         zt4 = bxy(1,k,1,l) - dth*(dky*zt1)
         zt6 = bxy(3,k,1,l) + dth*(dky*zt3)
c update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt3 = cmplx(-aimag(zt4),real(zt4))
         zt7 = exy(1,k,1,l) + cdt*(dky*zt1) - afdt*cu(1,k,1,l)
         zt9 = exy(3,k,1,l) - cdt*(dky*zt3) - afdt*cu(3,k,1,l)
c update magnetic field half time step and store electric field
         zt1 = cmplx(-aimag(zt9),real(zt9))
         zt3 = cmplx(-aimag(zt7),real(zt7))
         exy(1,k,1,l) = zt7
         exy(2,k,1,l) = zero
         exy(3,k,1,l) = zt9
         ws = ws + anorm*(zt7*conjg(zt7) + zt9*conjg(zt9))
         zt4 = zt4 - dth*(dky*zt1)
         zt6 = zt6 + dth*(dky*zt3)
         bxy(1,k,1,l) = zt4
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = zt6
         wp = wp + anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         exy(1,k1,1,l) = zero
         exy(2,k1,1,l) = zero
         exy(3,k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
         exy(1,1,1,l) = zero
         exy(2,1,1,l) = zero
         exy(3,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         exy(1,k1,1,l) = zero
         exy(2,k1,1,l) = zero
         exy(3,k1,1,l) = zero
      endif
   40 continue
   50 continue
      wf = float(nx)*float(ny)*ws
      wm = float(nx)*float(ny)*c2*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PEMFIELD2(fxy,exy,ffc,isign,nx,ny,kstrt,nyv,kxp,jblok,n
     1yhd)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign < 0
c includes additional smoothing
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
      complex fxy, exy, ffc
      dimension fxy(3,nyv,kxp,jblok), exy(3,nyv,kxp,jblok)
      dimension ffc(nyhd,kxp,jblok)
c local data
      integer i, j, k, l, nxh, nyh, ny2, k1
      real at1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      if (kstrt.gt.nxh) return
c add the fields
      if (isign.gt.0) then
         do 50 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
         do 40 j = 1, kxp
         do 20 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,j,l))
         do 10 i = 1, 3
         fxy(i,k,j,l) = fxy(i,k,j,l) + exy(i,k,j,l)*at1
         fxy(i,k1,j,l) = fxy(i,k1,j,l) + exy(i,k1,j,l)*at1
   10    continue
   20    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = aimag(ffc(1,j,l))
         do 30 i = 1, 3
         fxy(i,1,j,l) = fxy(i,1,j,l) + exy(i,1,j,l)*at1
         fxy(i,k1,j,l) = fxy(i,k1,j,l) + exy(i,k1,j,l)*at1
   30    continue
   40    continue
   50    continue
c copy the fields
      else if (isign.lt.0) then
         do 100 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
         do 90 j = 1, kxp
         do 70 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffc(k,j,l))
         do 60 i = 1, 3
         fxy(i,k,j,l) = exy(i,k,j,l)*at1
         fxy(i,k1,j,l) = exy(i,k1,j,l)*at1
   60    continue
   70    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = aimag(ffc(1,j,l))
         do 80 i = 1, 3
         fxy(i,1,j,l) = exy(i,1,j,l)*at1
         fxy(i,k1,j,l) = exy(i,k1,j,l)*at1
   80    continue
   90    continue
  100    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PEMFIELDR2(fxy,exy,ffd,isign,nx,ny,kstrt,nyv,kxp2,j2blo
     1k,nyd)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign < 0
c includes additional smoothing
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyd
      real fxy, exy
      complex ffd
      dimension fxy(3,nyv,kxp2+1,j2blok), exy(3,nyv,kxp2+1,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
c local data
      integer i, j, k, l, ny1
      real at1
      ny1 = ny + 1
      if (kstrt.gt.nx) return
c smooth the em field and add
      if (isign.gt.0) then
         do 70 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
         do 40 j = 1, kxp2
         do 20 k = 1, ny
         at1 = aimag(ffd(k,j,l))
         do 10 i = 1, 3
         fxy(i,k,j,l) = fxy(i,k,j,l) + exy(i,k,j,l)*at1
   10    continue
   20    continue
c mode numbers ky = 0, ny/2
         do 30 i = 1, 3
         fxy(i,ny+1,j,l) = fxy(i,ny+1,j,l) + exy(i,ny+1,j,l)
   30    continue
   40    continue
         do 60 k = 1, ny1
         do 50 i = 1, 3
         fxy(i,k,kxp2+1,l) = fxy(i,k,kxp2+1,l) + exy(i,k,kxp2+1,l)
   50    continue
   60    continue
   70    continue
c copy and smooth the magnetic fields
      else if (isign.lt.0) then
         do 140 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
         do 110 j = 1, kxp2
         do 90 k = 1, ny
         at1 = aimag(ffd(k,j,l))
         do 80 i = 1, 3
         fxy(i,k,j,l) = exy(i,k,j,l)*at1
   80    continue
   90    continue
c mode numbers ky = 0, ny/2
         do 100 i = 1, 3
         fxy(i,ny+1,j,l) = exy(i,ny+1,j,l)
  100    continue
  110    continue
         do 130 k = 1, ny1
         do 120 i = 1, 3
         fxy(i,k,kxp2+1,l) = exy(i,k,kxp2+1,l)
  120    continue
  130    continue
  140    continue
c copy the electric fields
      else if (isign.eq.0) then
         do 210 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
         do 180 j = 1, kxp2
         do 160 k = 1, ny
         do 150 i = 1, 3
         fxy(i,k,j,l) = exy(i,k,j,l)
  150    continue
  160    continue
c mode numbers ky = 0, ny/2
         do 170 i = 1, 3
         fxy(i,ny+1,j,l) = exy(i,ny+1,j,l)
  170    continue
  180    continue
         do 200 k = 1, ny1
         do 190 i = 1, 3
         fxy(i,k,kxp2+1,l) = exy(i,k,kxp2+1,l)
  190    continue
  200    continue
  210    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PAVPOT23(bxy,axy,nx,ny,kstrt,nyv,kxp,jblok)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with periodic boundary conditions
c for distributed data.
c input: bxy,nx,ny,kstrt,nyv,kxp,jblok, output: axy
c approximate flop count is: 38*nxc*nyc + 10*(nxc + nyc)
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the vector potential is calculated using the equations:
c ax(kx,ky) = sqrt(-1)*(ky*bz(kx,ky))/(kx*kx+ky*ky)
c ay(kx,ky) = -sqrt(-1)*(kx*bz(kx,ky))/(kx*kx+ky*ky)
c az(kx,ky) = sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))/(kx*kx+ky*ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
c ax(ky=pi) = ay(ky=pi) = az(ky=pi) = 0,
c ax(kx=0,ky=0) = ay(kx=0,ky=0) = az(kx=0,ky=0) = 0.
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c axy(i,k,j,l) = i-th component of complex vector potential,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c jblok = number of data blocks
c kxp = number of data values per block
c kstrt = starting data block number
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
      complex bxy, axy, zero, zt1, zt2, zt3
      dimension bxy(3,nyv,kxp,jblok), axy(3,nyv,kxp,jblok)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate vector potential
      if (kstrt.gt.nxh) go to 50
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         at2 = dky*at1
         at3 = dkx*at1
         zt1 = cmplx(-aimag(bxy(3,k,j,l)),real(bxy(3,k,j,l)))
         zt2 = cmplx(-aimag(bxy(2,k,j,l)),real(bxy(2,k,j,l)))
         zt3 = cmplx(-aimag(bxy(1,k,j,l)),real(bxy(1,k,j,l)))
         axy(1,k,j,l) = at2*zt1
         axy(2,k,j,l) = -at3*zt1
         axy(3,k,j,l) = at3*zt2 - at2*zt3
         zt1 = cmplx(-aimag(bxy(3,k1,j,l)),real(bxy(3,k1,j,l)))
         zt2 = cmplx(-aimag(bxy(2,k1,j,l)),real(bxy(2,k1,j,l)))
         zt3 = cmplx(-aimag(bxy(1,k1,j,l)),real(bxy(1,k1,j,l)))
         axy(1,k1,j,l) = -at2*zt1
         axy(2,k1,j,l) = -at3*zt1
         axy(3,k1,j,l) = at3*zt2 + at2*zt3
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = 1.0/dkx
         zt1 = cmplx(-aimag(bxy(3,1,j,l)),real(bxy(3,1,j,l)))
         zt2 = cmplx(-aimag(bxy(2,1,j,l)),real(bxy(2,1,j,l)))
         axy(1,1,j,l) = zero
         axy(2,1,j,l) = -at2*zt1
         axy(3,1,j,l) = at2*zt2
         axy(1,k1,j,l) = zero
         axy(2,k1,j,l) = zero
         axy(3,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at2 = 1.0/dky
         zt1 = cmplx(-aimag(bxy(3,k,1,l)),real(bxy(3,k,1,l)))
         zt2 = cmplx(-aimag(bxy(1,k,1,l)),real(bxy(1,k,1,l)))
         axy(1,k,1,l) = at2*zt1
         axy(2,k,1,l) = zero
         axy(3,k,1,l) = -at2*zt2
         axy(1,k1,1,l) = zero
         axy(2,k1,1,l) = zero
         axy(3,k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         axy(1,1,1,l) = zero
         axy(2,1,1,l) = zero
         axy(3,1,1,l) = zero
         axy(1,k1,1,l) = zero
         axy(2,k1,1,l) = zero
         axy(3,k1,1,l) = zero
      endif
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAVRPOT23(axy,bxy,ffc,affp,ci,nx,ny,kstrt,nyv,kxp,jblok
     1,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c the radiative part of the vector potential
c with periodic boundary conditions, for distributed data.
c input: all, output: axy
c approximate flop count is: 68*nxc*nyc + 20*(nxc + nyc)
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the radiative vector potential is updated using the equations:
c ax(kx,ky) = (sqrt(-1)*ky*bz(kx,ky)
c                       - affp*ci2*cux(kx,ky)*s(kx,ky)/(kx*kx+ky*ky)
c ay(kx,ky) = -(sqrt(-1)*kx*bz(kx,ky)
c                       + affp*ci2*cuy(kx,ky)*s(kx,ky))/(kx*kx+ky*ky)
c az(kx,ky) = (sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
c                       - affp*ci2*cuz(kx,ky)*s(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, ci2 = ci*ci
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers, except for
c ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
c ax(ky=pi) = ay(ky=pi) = az(ky=pi) = 0,
c ax(kx=0,ky=0) = ay(kx=0,ky=0) = az(kx=0,ky=0) = 0.
c axy(i,k,j,l) = on entry, i-th component of complex current density cu,
c axy(i,k,j,l) = on exit, i-th component of complex radiative vector
c potential,
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c aimag(ffc(k,j,l)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c kxp = number of data values per block
c jblok = number of data blocks
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, jblok, nyhd
      real affp, ci
      complex axy, bxy, ffc
      dimension axy(3,nyv,kxp,jblok), bxy(3,nyv,kxp,jblok)
      dimension ffc(nyhd,kxp,jblok)
c local data
      integer nxh, nyh, ny2, ks, joff, j, k, k1, l
      real dnx, dny, afc2, dkx, dkx2, dky, at1, at2
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      afc2 = affp*ci*ci
      zero = cmplx(0.,0.)
c calculate the radiative vector potential
      if (kstrt.gt.nxh) go to 50
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = 1.0/(dky*dky + dkx2)
         at2 = afc2*aimag(ffc(k,j,l))
c update radiative vector potential, ky > 0
         zt1 = cmplx(-aimag(bxy(3,k,j,l)),real(bxy(3,k,j,l)))
         zt2 = cmplx(-aimag(bxy(2,k,j,l)),real(bxy(2,k,j,l)))
         zt3 = cmplx(-aimag(bxy(1,k,j,l)),real(bxy(1,k,j,l)))
         axy(1,k,j,l) = at1*(dky*zt1 - at2*axy(1,k,j,l))
         axy(2,k,j,l) = -at1*(dkx*zt1 + at2*axy(2,k,j,l))
         axy(3,k,j,l) = at1*((dkx*zt2 - dky*zt3) - at2*axy(3,k,j,l))
c update radiative vector potential, ky < 0
         zt1 = cmplx(-aimag(bxy(3,k1,j,l)),real(bxy(3,k1,j,l)))
         zt2 = cmplx(-aimag(bxy(2,k1,j,l)),real(bxy(2,k1,j,l)))
         zt3 = cmplx(-aimag(bxy(1,k1,j,l)),real(bxy(1,k1,j,l)))
         axy(1,k1,j,l) = -at1*(dky*zt1 + at2*axy(1,k1,j,l))
         axy(2,k1,j,l) = -at1*(dkx*zt1 + at2*axy(2,k1,j,l))
         axy(3,k1,j,l) = at1*((dkx*zt2 + dky*zt3) - at2*axy(3,k1,j,l))
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = 1.0/dkx2
         at2 = afc2*aimag(ffc(1,j,l))
c update radiative vector potential
         zt1 = cmplx(-aimag(bxy(3,1,j,l)),real(bxy(3,1,j,l)))
         zt2 = cmplx(-aimag(bxy(2,1,j,l)),real(bxy(2,1,j,l)))
         axy(1,1,j,l) = zero
         axy(2,1,j,l) = -at1*(dkx*zt1 + at2*axy(2,1,j,l))
         axy(3,1,j,l) = at1*(dkx*zt2 - at2*axy(3,1,j,l))
         axy(1,k1,j,l) = zero
         axy(2,k1,j,l) = zero
         axy(3,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = 1.0/(dky*dky)
         at2 = afc2*aimag(ffc(k,1,l))
c update radiative vector potential
         zt1 = cmplx(-aimag(bxy(3,k,1,l)),real(bxy(3,k,1,l)))
         zt3 = cmplx(-aimag(bxy(1,k,1,l)),real(bxy(1,k,1,l)))
         axy(1,k,1,l) = at1*(dky*zt1 - at2*axy(1,k,1,l))
         axy(2,k,1,l) = zero
         axy(3,k,1,l) = -at1*(dky*zt3 + at2*axy(3,k,1,l))
         axy(1,k1,1,l) = zero
         axy(2,k1,1,l) = zero
         axy(3,k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         axy(1,1,1,l) = zero
         axy(2,1,1,l) = zero
         axy(3,1,1,l) = zero
         axy(1,k1,1,l) = zero
         axy(2,k1,1,l) = zero
         axy(3,k1,1,l) = zero
      endif
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGTMODES2(pot,pott,nx,ny,it,modesx,modesy,kstrt,nyv,kxp
     1,jblok,nt,modesxpd,modesyd)
c this subroutine extracts lowest order modes from complex array pot
c and stores them into a location in a time history array pott
c modes stored: kx = (kxp*(idproc)+(0,1,...kxp-1)) where idproc=kstrt-1,
c and ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
c except kx = NX/2 is stored at location kxp+1 when idproc=0.
c nx/ny = system length in x/y direction
c it = current time
c modesx/modesy = number of modes to store in x/y direction,
c where modesx <= nx/2+1, modesy <= ny/2+1
c kstrt = starting data block number
c nyv = first dimension of input array pot, nyv >= ny
c kxp = number of data values per block
c jblok = number of data blocks
c nt = first dimension of output array pott, nt >= it
c modesyd = second dimension of array pott, modesyd >= 2*modesy
c modesxpd = third dimension of array pott, modesxpd >= min(modesx,kxp)
c unless modesx = nx/2+1, in which case modesxpd = kxp+1
      implicit none
      integer nx, ny, it, modesx, modesy, kstrt, nyv, kxp, jblok
      integer nt, modesxpd, modesyd
      complex pot, pott
      dimension pot(nyv,kxp,jblok), pott(nt,modesyd,modesxpd,jblok)
c local data
      integer nxh, nyh, jmax, kmax, ny2, j, k, l, j1, k1, ks, joff
      if (it.gt.nt) return
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      kmax = min0(modesy,nyh)
      j1 = kxp + 1
      ks = kstrt - 2
      if (kstrt.gt.nxh) go to 50
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks)
      jmax = modesx - joff
      if (jmax.gt.kxp) then
         jmax = kxp
      else if (jmax.le.0) then
         jmax = 0
      endif
      do 20 j = 1, jmax
      if ((j+joff).gt.1) then
         do 10 k = 2, kmax
         k1 = ny2 - k
         pott(it,2*k-2,j,l) = pot(k,j,l)
         pott(it,2*k-1,j,l) = pot(k1,j,l)
   10    continue
c mode numbers ky = 0, ny/2
         pott(it,1,j,l) = pot(1,j,l)
         if (modesy.gt.nyh) then
            k1 = nyh + 1
            pott(it,ny,j,l) = pot(k1,j,l)
         endif
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, kmax
         k1 = ny2 - k
         pott(it,2*k-2,1,l) = pot(k,1,l)
         pott(it,2*k-1,1,l) = conjg(pot(k,1,l))
         if (modesx.gt.nxh) then
            pott(it,2*k-2,j1,l) = conjg(pot(k1,1,l))
            pott(it,2*k-1,j1,l) = pot(k1,1,l)
         endif
   30    continue
         pott(it,1,1,l) = cmplx(real(pot(1,1,l)),0.0)
         if (modesx.gt.nxh) then
            pott(it,1,j1,l) = cmplx(aimag(pot(1,1,l)),0.0)
         endif
         if (modesy.gt.nyh) then
            k1 = nyh + 1
            pott(it,ny,1,l) = cmplx(real(pot(k1,1,l)),0.0)
            if (modesx.gt.nxh) then
               pott(it,ny,j1,l) = cmplx(aimag(pot(k1,1,l)),0.0)
            endif
         endif
      endif
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPTMODES2(pot,pott,nx,ny,it,modesx,modesy,kstrt,nyv,kxp
     1,jblok,nt,modesxpd,modesyd)
c this subroutine extracts lowest order modes from a location in a time
c history array pott and stores them into complex array pot
c modes stored: kx = (kxp*(idproc)+(0,1,...kxp-1)) where idproc=kstrt-1,
c and ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
c except kx = NX/2 is stored at location kxp+1 when idproc=0.
c nx/ny = system length in x/y direction
c it = current time
c modesx/modesy = number of modes to store in x/y direction,
c where modesx <= nx/2+1, modesy <= ny/2+1
c kstrt = starting data block number
c nyv = first dimension of input array pot, nyv >= ny
c kxp = number of data values per block
c jblok = number of data blocks
c nt = first dimension of output array pott, nt >= it
c modesyd = second dimension of array pott, modesyd  >= 2*modesy
c modesxpd = third dimension of array pott, modesxpd >= min(modesx,kxp)
c unless modesx = nx/2+1, in which case modesxpd = kxp+1
      implicit none
      integer nx, ny, it, modesx, modesy, kstrt, nyv, kxp, jblok
      integer nt, modesxpd, modesyd
      complex pot, pott
      dimension pot(nyv,kxp,jblok), pott(nt,modesyd,modesxpd,jblok)
c local data
      integer nxh, nyh, jmax, kmax, ny2, j, k, l, j1, k1, ks, joff
      complex zero
      if (it.gt.nt) return
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      kmax = min0(modesy,nyh)
      j1 = kxp + 1
      ks = kstrt - 2
      zero = cmplx(0.,0.)
      if (kstrt.gt.nxh) go to 90
      do 80 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks)
      jmax = modesx - joff
      if (jmax.gt.kxp) then
         jmax = kxp
      else if (jmax.le.0) then
         jmax = 0
      endif
      do 30 j = 1, jmax
      if ((j+joff).gt.1) then
         do 10 k = 2, kmax
         k1 = ny2 - k
         pot(k,j,l) = pott(it,2*k-2,j,l)
         pot(k1,j,l) = pott(it,2*k-1,j,l)
   10    continue
         do 20 k = kmax+1, nyh
         k1 = ny2 - k
         pot(k,j,l) = zero
         pot(k1,j,l) = zero
   20    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         pot(1,j,l) = pott(it,1,j,l)
         pot(k1,j,l) = zero
         if (modesy.gt.nyh) then
            pot(k1,j,l) = pott(it,ny,j,l)
         endif
      endif
   30 continue
      do 50 j = jmax+1, kxp
      if ((j+joff).gt.1) then
         do 40 k = 2, nyh
         k1 = ny2 - k
         pot(k,j,l) = zero
         pot(k1,j,l) = zero
   40    continue
         k1 = nyh + 1
c mode numbers ky = 0, ny/2
         pot(1,j,l) = zero
         pot(k1,j,l) = zero
      endif
   50 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 60 k = 2, kmax
         k1 = ny2 - k
         pot(k,1,l) = pott(it,2*k-2,1,l)
         pot(k1,1,l) = zero
         if (modesx.gt.nxh) then
            pot(k1,1,l) = conjg(pott(it,2*k-2,j1,l))
         endif
   60    continue
         do 70 k = kmax+1, nyh
         k1 = ny2 - k
         pot(k,1,l) = zero
         pot(k1,1,l) = zero
   70    continue
         k1 = nyh + 1
         pot(1,1,l) = cmplx(real(pott(it,1,1,l)),0.0)
         pot(k1,1,l) = zero
         if (modesx.gt.nxh) then
            pot(1,1,l) = cmplx(real(pot(1,1,l)),real(pott(it,1,j1,l)))
         endif
         if (modesy.gt.nyh) then
            pot(k1,1,l) = cmplx(real(pott(it,ny,1,l)),0.0)
            if (modesx.gt.nxh) then
               pot(k1,1,l) = cmplx(real(pot(k1,1,l)),real(pott(it,ny,j1,
     1l)))
            endif
         endif
      endif
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGTVMODES2(vpot,vpott,nx,ny,it,modesx,modesy,ndim,kstrt
     1,nyv,kxp,jblok,nt,modesxpd,modesyd)
c this subroutine extracts lowest order modes from complex vector array
c vpot and stores them into a location in a time history array vpott
c modes stored: kx = (kxp*(idproc)+(0,1,...kxp-1)) where idproc=kstrt-1,
c and ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
c except kx = NX/2 is stored at location kxp+1 when idproc=0.
c nx/ny = system length in x/y direction
c it = current time
c modesx/modesy = number of modes to store in x/y direction,
c where modesx <= nx/2+1, modesy <= ny/2+1
c ndim = number of field arrays, must be >= 1
c kstrt = starting data block number
c nyv = second dimension of input array vpot, nyv >= ny
c kxp = number of data values per block
c jblok = number of data blocks
c nt = first dimension of output array vpott, nt >= it
c modesyd = third dimension of array vpott, modesyd >= 2*modesy
c modesxpd = fourth dimension of array vpott,
c modesxpd >= min(modesx,kxp),  unless modesx = nx/2+1,
c in which case modesxpd = kxp+1
      implicit none
      integer nx, ny, it, modesx, modesy, ndim, kstrt, nyv, kxp, jblok
      integer nt, modesxpd, modesyd
      complex vpot, vpott
      dimension vpot(ndim,nyv,kxp,jblok)
      dimension vpott(nt,ndim,modesyd,modesxpd,jblok)
c local data
      integer nxh, nyh, jmax, kmax, ny2, i, j, k, l, j1, k1, ks, joff
      if (it.gt.nt) return
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      kmax = min0(modesy,nyh)
      j1 = kxp + 1
      ks = kstrt - 2
      if (kstrt.gt.nxh) go to 90
      do 80 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks)
      jmax = modesx - joff
      if (jmax.gt.kxp) then
         jmax = kxp
      else if (jmax.le.0) then
         jmax = 0
      endif
      do 40 j = 1, jmax
      if ((j+joff).gt.1) then
         do 20 k = 2, kmax
         k1 = ny2 - k
         do 10 i = 1, ndim
         vpott(it,i,2*k-2,j,l) = vpot(i,k,j,l)
         vpott(it,i,2*k-1,j,l) = vpot(i,k1,j,l)
   10    continue
   20    continue
c mode numbers ky = 0, ny/2
         do 30 i = 1, ndim
         vpott(it,i,1,j,l) = vpot(i,1,j,l)
         if (modesy.gt.nyh) then
            k1 = nyh + 1
            vpott(it,i,ny,j,l) = vpot(i,k1,j,l)
         endif
   30    continue
      endif
   40 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 60 k = 2, kmax
         k1 = ny2 - k
         do 50 i = 1, ndim
         vpott(it,i,2*k-2,1,l) = vpot(i,k,1,l)
         vpott(it,i,2*k-1,1,l) = conjg(vpot(i,k,1,l))
         if (modesx.gt.nxh) then
            vpott(it,i,2*k-2,j1,l) = conjg(vpot(i,k1,1,l))
            vpott(it,i,2*k-1,j1,l) = vpot(i,k1,1,l)
         endif
   50    continue
   60    continue
         do 70 i = 1, ndim
         vpott(it,i,1,1,l) = cmplx(real(vpot(i,1,1,l)),0.0)
         if (modesx.gt.nxh) then
            vpott(it,i,1,j1,l) = cmplx(aimag(vpot(i,1,1,l)),0.0)
         endif
         if (modesy.gt.nyh) then
            k1 = nyh + 1
            vpott(it,i,ny,1,l) = cmplx(real(vpot(i,k1,1,l)),0.0)
            if (modesx.gt.nxh) then
               vpott(it,i,ny,j1,l) = cmplx(aimag(vpot(i,k1,1,l)),0.0)
            endif
         endif
   70    continue
      endif
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPTVMODES2(vpot,vpott,nx,ny,it,modesx,modesy,ndim,kstrt
     1,nyv,kxp,jblok,nt,modesxpd,modesyd)
c this subroutine extracts lowest order modes from a location in a time
c history array vpott and stores them into complex vector array vpot
c modes stored: kx = (kxp*(idproc)+(0,1,...kxp-1)) where idproc=kstrt-1,
c and ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
c except kx = NX/2 is stored at location kxp+1 when idproc=0.
c nx/ny = system length in x/y direction
c it = current time
c modesx/modesy = number of modes to store in x/y direction,
c where modesx <= nx/2+1, modesy <= ny/2+1
c ndim = number of field arrays, must be >= 1
c kstrt = starting data block number
c nyv = second dimension of input array vpot, nyv >= ny
c kxp = number of data values per block
c jblok = number of data blocks
c nt = first dimension of output array vpott, nt >= it
c modesyd = third dimension of array vpott, modesyd  >= 2*modesy
c modesxpd = fourth dimension of array vpott,
c modesxpd >= min(modesx,kxp) unless modesx = nx/2+1,
c in which case modesxpd = kxp+1
      implicit none
      integer nx, ny, it, modesx, modesy, ndim, kstrt, nyv, kxp, jblok
      integer nt, modesxpd, modesyd
      complex vpot, vpott
      dimension vpot(ndim,nyv,kxp,jblok)
      dimension vpott(nt,ndim,modesyd,modesxpd,jblok)
c local data
      integer nxh, nyh, jmax, kmax, ny2, i, j, k, l, j1, k1, ks, joff
      complex zero
      if (it.gt.nt) return
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      kmax = min0(modesy,nyh)
      j1 = kxp + 1
      ks = kstrt - 2
      zero = cmplx(0.,0.)
      if (kstrt.gt.nxh) go to 170
      do 160 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks)
      jmax = modesx - joff
      if (jmax.gt.kxp) then
         jmax = kxp
      else if (jmax.le.0) then
         jmax = 0
      endif
      do 60 j = 1, jmax
      if ((j+joff).gt.1) then
         do 20 k = 2, kmax
         do 10 i = 1, ndim
         k1 = ny2 - k
         vpot(i,k,j,l) = vpott(it,i,2*k-2,j,l)
         vpot(i,k1,j,l) = vpott(it,i,2*k-1,j,l)
   10    continue
   20    continue
         do 40 k = kmax+1, nyh
         k1 = ny2 - k
         do 30 i = 1, ndim
         vpot(i,k,j,l) = zero
         vpot(i,k1,j,l) = zero
   30    continue
   40    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         do 50 i = 1, ndim
         vpot(i,1,j,l) = vpott(it,i,1,j,l)
         vpot(i,k1,j,l) = zero
         if (modesy.gt.nyh) then
            vpot(i,k1,j,l) = vpott(it,i,ny,j,l)
         endif
   50    continue
      endif
   60 continue
      do 100 j = jmax+1, kxp
      if ((j+joff).gt.1) then
         do 80 k = 2, nyh
         k1 = ny2 - k
         do 70 i = 1, ndim
         vpot(i,k,j,l) = zero
         vpot(i,k1,j,l) = zero
   70    continue
   80    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         do 90 i = 1, ndim
         vpot(i,1,j,l) = zero
         vpot(i,k1,j,l) = zero
   90    continue
      endif
  100 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 120 k = 2, kmax
         k1 = ny2 - k
         do 110 i = 1, ndim
         vpot(i,k,1,l) = vpott(it,i,2*k-2,1,l)
         vpot(i,k1,1,l) = zero
         if (modesx.gt.nxh) then
            vpot(i,k1,1,l) = conjg(vpott(it,i,2*k-2,j1,l))
         endif
  110    continue
  120    continue
         do 140 k = kmax+1, nyh
         k1 = ny2 - k
         do 130 i = 1, ndim
         vpot(i,k,1,l) = zero
         vpot(i,k1,1,l) = zero
  130    continue
  140    continue
         k1 = nyh + 1
         do 150 i = 1, ndim
         vpot(i,1,1,l) = cmplx(real(vpott(it,i,1,1,l)),0.0)
         vpot(i,k1,1,l) = zero
         if (modesx.gt.nxh) then
            vpot(i,1,1,l) = cmplx(real(vpot(i,1,1,l)),real(vpott(it,i,1,
     1j1,l)))
         endif
         if (modesy.gt.nyh) then
            vpot(i,k1,1,l) = cmplx(real(vpott(it,i,ny,1,l)),0.0)
            if (modesx.gt.nxh) then
               vpot(i,k1,1,l) = cmplx(real(vpot(i,k1,1,l)),real(vpott(it
     1,i,ny,j1,l)))
            endif
         endif
  150    continue
      endif
  160 continue
  170 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOYNT2(q,exy,bxy,ffc,affp,sx,sy,sz,nx,ny,kstrt,nyv,kxp
     1,jblok,nyhd)
c this subroutine calculates the momentum in the electromagnetic field
c given by the poynting flux.  inputs are the charge density, transverse
c electric field, and magnetic field.  outputs are sx, sy, sz
c equation used is:
c sx = sum((fy(j,k)+exy(2,j,k))*conjg(bxy(3,j,k))-exy(3,j,k)*
c conjg(bxy(2,j,k)))
c sy = sum(exy(3,j,k)*conjg(bxy(1,j,k))-(fx(j,k)+exy(1,j,k))*
c conjg(bxy(3,j,k)))
c sz = sum((fx(j,k)+exy(1,j,k))*conjg(bxy(2,j,k))-(fy(j,k)+exy(2,j,k))*
c conjg(bxy(1,j,k))), where
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*q(kx,ky),
c kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c exy(i,k,j,l) = i-th component of complex transverse electric field,
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c jblok = number of data blocks
c kxp = number of data values per block
c kstrt = starting data block number
c aimag(ffc(k,j,l)) = finite-size particle shape factor s
c real(ffc(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c affp = normalization constant = nx*ny/np, where np=number of particles
c sx/sy/sz = x/y/z components of field momentum
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, jblok, nyhd
      real sx, sy, sz
      complex q, exy, bxy, ffc
      dimension q(nyv,kxp,jblok)
      dimension exy(3,nyv,kxp,jblok), bxy(3,nyv,kxp,jblok)
      dimension ffc(nyhd,kxp,jblok)
c local data
      integer nxh, nyh, ny2, ks, joff, j, k, k1, l
      real dnx, dny, affp, dkx, at1, at2, at3
      complex zt1, zt2
      double precision wx, wy, wz
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
c sum darwin field momentum
      wx = 0.0d0
      wy = 0.0d0
      wz = 0.0d0
      if (kstrt.gt.nxh) go to 50
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,j,l))
         at2 = dkx*at1
         at3 = dny*float(k - 1)*at1
         zt1 = cmplx(aimag(q(k,j,l)),-real(q(k,j,l)))
         zt2 = at3*zt1 + exy(2,k,j,l)
         zt1 = at2*zt1 + exy(1,k,j,l)
         wx = wx + (zt2*conjg(bxy(3,k,j,l))-exy(3,k,j,l)*conjg(bxy(2,k,j
     1,l)))
         wy = wy + (exy(3,k,j,l)*conjg(bxy(1,k,j,l))-zt1*conjg(bxy(3,k,j
     1,l)))
         wz = wz + (zt1*conjg(bxy(2,k,j,l))-zt2*conjg(bxy(1,k,j,l)))
         zt2 = cmplx(aimag(q(k1,j,l)),-real(q(k1,j,l)))
         zt1 = at2*zt2 + exy(1,k1,j,l)
         zt2 = -at3*zt2 + exy(2,k1,j,l)
         wx = wx + (zt2*conjg(bxy(3,k1,j,l))-exy(3,k1,j,l)*conjg(bxy(2,k
     11,j,l)))
         wy = wy + (exy(3,k1,j,l)*conjg(bxy(1,k1,j,l))-zt1*conjg(bxy(3,k
     11,j,l)))
         wz = wz + (zt1*conjg(bxy(2,k1,j,l))-zt2*conjg(bxy(1,k1,j,l)))
   10    continue
c mode numbers ky = 0, ny/2
         at1 = real(ffc(1,j,l))
         at2 = dkx*at1
         zt1 = cmplx(aimag(q(1,j,l)),-real(q(1,j,l)))
         zt2 = exy(2,1,j,l)
         zt1 = at2*zt1 + exy(1,1,j,l)
         wx = wx + (zt2*conjg(bxy(3,1,j,l))-exy(3,1,j,l)*conjg(bxy(2,1,j
     1,l)))
         wy = wy + (exy(3,1,j,l)*conjg(bxy(1,1,j,l))-zt1*conjg(bxy(3,1,j
     1,l)))
         wz = wz + (zt1*conjg(bxy(2,1,j,l))-zt2*conjg(bxy(1,1,j,l)))
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         at1 = real(ffc(k,1,l))
         at3 = dny*float(k - 1)*at1
         zt1 = cmplx(aimag(q(k,1,l)),-real(q(k,1,l)))
         zt2 = at3*zt1 + exy(2,k,1,l)
         zt1 = exy(1,k,1,l)
         wx = wx + (zt2*conjg(bxy(3,k,1,l))-exy(3,k,1,l)*conjg(bxy(2,k,1
     1,l)))
         wy = wy + (exy(3,k,1,l)*conjg(bxy(1,k,1,l))-zt1*conjg(bxy(3,k,1
     1,l)))
         wz = wz + (zt1*conjg(bxy(2,k,1,l))-zt2*conjg(bxy(1,k,1,l)))
   30    continue
      endif
   40 continue
   50 continue
      at1 = 2.0*float(nx)*float(ny)/affp
      sx = at1*wx
      sy = at1*wy
      sz = at1*wz
      return
      end
c-----------------------------------------------------------------------
      subroutine PDPOYNT2(q,cu,ffc,affp,ci,sx,sy,sz,nx,ny,kstrt,nyv,kxp,
     1jblok,nyhd)
c this subroutine calculates the momentum in the darwin field given by
c the poynting flux in 2-1/2d, for distributed data.
c inputs are the charge density and current density.
c outputs are sx, sy, sz
c equation used is:
c sx = sum(fy(j,k)*conjg(bz(j,k))
c sy = sum(-fx(j,k)*conjg(bz(j,k)))
c sz = sum(fx(j,k)*conjg(by(j,k))-fy(j,k)*conjg(bx(j,k))), where
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*q(kx,ky),
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*-kx*cuz(kx,ky)
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))
c kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c cu(i,k,j,l) = i-th component of complex current density,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c jblok = number of data blocks
c kxp = number of data values per block
c kstrt = starting data block number
c aimag(ffc(k,j,l)) = finite-size particle shape factor s
c real(ffc(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c sx/sy/sz = x/y/z components of field momentum
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, jblok, nyhd
      real ci, sx, sy, sz
      complex q, cu, ffc
      dimension q(nyv,kxp,jblok), cu(3,nyv,kxp,jblok)
      dimension ffc(nyhd,kxp,jblok)
c local data
      integer nxh, nyh, ny2, ks, joff, j, k, k1, l
      real dnx, dny, ci2, affp, dkx, at1, at2, at3, at4, at5
      complex zt1, zt2, zt3, zt4, zt5
      double precision wx, wy, wz
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      ci2 = ci*ci
c sum darwin field momentum
      wx = 0.0d0
      wy = 0.0d0
      wz = 0.0d0
      if (kstrt.gt.nxh) go to 50
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,j,l))
         at2 = dkx*at1
         at3 = dny*float(k - 1)*at1
         zt4 = cmplx(aimag(q(k,j,l)),-real(q(k,j,l)))
         zt5 = at3*zt4
         zt4 = at2*zt4
         at4 = ci2*at2
         at5 = ci2*at3
         zt1 = cmplx(-aimag(cu(3,k,j,l)),real(cu(3,k,j,l)))
         zt2 = cmplx(-aimag(cu(2,k,j,l)),real(cu(2,k,j,l)))
         zt3 = cmplx(-aimag(cu(1,k,j,l)),real(cu(1,k,j,l)))
         zt3 = conjg(at4*zt2 - at5*zt3)
         zt2 = conjg(at4*zt1)
         zt1 = conjg(at5*zt1)
         wx = wx + zt5*zt3
         wy = wy - zt4*zt3
         wz = wz - (zt4*zt2+zt5*zt1)
         zt5 = cmplx(aimag(q(k1,j,l)),-real(q(k1,j,l)))
         zt4 = at2*zt5
         zt5 = at3*zt5
         zt1 = cmplx(-aimag(cu(3,k1,j,l)),real(cu(3,k1,j,l)))
         zt2 = cmplx(-aimag(cu(2,k1,j,l)),real(cu(2,k1,j,l)))
         zt3 = cmplx(-aimag(cu(1,k1,j,l)),real(cu(1,k1,j,l)))
         zt3 = conjg(at4*zt2 + at5*zt3)
         zt2 = conjg(at4*zt1)
         zt1 = conjg(at5*zt1)
         wx = wx - zt5*zt3
         wy = wy - zt4*zt3
         wz = wz - (zt4*zt2+zt5*zt1)
   10    continue
c mode numbers ky = 0, ny/2
         at1 = real(ffc(1,j,l))
         at2 = dkx*at1
         zt4 = cmplx(aimag(q(1,j,l)),-real(q(1,j,l)))
         zt4 = at2*zt4
         at2 = ci2*at2
         zt1 = cmplx(-aimag(cu(3,1,j,l)),real(cu(3,1,j,l)))
         zt2 = cmplx(-aimag(cu(2,1,j,l)),real(cu(2,1,j,l)))
         zt3 = conjg(at2*zt2)
         zt2 = conjg(at2*zt1)
         wy = wy - zt4*zt3
         wz = wz - zt4*zt2
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         at1 = real(ffc(k,1,l))
         at3 = dny*float(k - 1)*at1
         zt4 = cmplx(aimag(q(k,1,l)),-real(q(k,1,l)))
         zt5 = at3*zt4
         at3 = ci2*at3
         zt1 = cmplx(-aimag(cu(3,k,1,l)),real(cu(3,k,1,l)))
         zt3 = cmplx(-aimag(cu(1,k,1,l)),real(cu(1,k,1,l)))
         zt1 = conjg(at3*zt1)
         zt3 = conjg(at3*zt3)
         wx = wx - zt5*zt3
         wz = wz - zt5*zt1
   30    continue
      endif
   40 continue
   50 continue
      at1 = 2.0*float(nx)*float(ny)/affp
      sx = at1*wx
      sy = at1*wy
      sz = at1*wz
      return
      end
c-----------------------------------------------------------------------
      subroutine PDPOYNT22(q,cu,ffc,affp,ci,sx,sy,nx,ny,kstrt,nyv,kxp,jb
     1lok,nyhd)
c this subroutine calculates the momentum in the darwin field given by
c the poynting flux in 2d, for distributed data.
c inputs are the charge density and current density.
c outputs are sx, sy
c equation used is:
c sx = sum(fy(j,k)*conjg(bz(j,k))
c sy = sum(-fx(j,k)*conjg(bz(j,k))), where
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*q(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))
c kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c cu(i,k,j,l) = i-th component of complex current density,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c jblok = number of data blocks
c kxp = number of data values per block
c kstrt = starting data block number
c aimag(ffc(k,j,l)) = finite-size particle shape factor s
c real(ffc(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c sx/sy = x/y components of field momentum
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, jblok, nyhd
      real ci, sx, sy
      complex q, cu, ffc
      dimension q(nyv,kxp,jblok), cu(2,nyv,kxp,jblok)
      dimension ffc(nyhd,kxp,jblok)
c local data
      integer nxh, nyh, ny2, ks, joff, j, k, k1, l
      real dnx, dny, ci2, affp, dkx, at1, at2, at3, at4, at5
      complex zt2, zt3, zt4, zt5
      double precision wx, wy
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      ci2 = ci*ci
c sum darwin field momentum
      wx = 0.0d0
      wy = 0.0d0
      if (kstrt.gt.nxh) go to 50
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffc(k,j,l))
         at2 = dkx*at1
         at3 = dny*float(k - 1)*at1
         zt4 = cmplx(aimag(q(k,j,l)),-real(q(k,j,l)))
         zt5 = at3*zt4
         zt4 = at2*zt4
         at4 = ci2*at2
         at5 = ci2*at3
         zt2 = cmplx(-aimag(cu(2,k,j,l)),real(cu(2,k,j,l)))
         zt3 = cmplx(-aimag(cu(1,k,j,l)),real(cu(1,k,j,l)))
         zt3 = conjg(at4*zt2 - at5*zt3)
         wx = wx + zt5*zt3
         wy = wy - zt4*zt3
         zt5 = cmplx(aimag(q(k1,j,l)),-real(q(k1,j,l)))
         zt4 = at2*zt5
         zt5 = at3*zt5
         zt2 = cmplx(-aimag(cu(2,k1,j,l)),real(cu(2,k1,j,l)))
         zt3 = cmplx(-aimag(cu(1,k1,j,l)),real(cu(1,k1,j,l)))
         zt3 = conjg(at4*zt2 + at5*zt3)
         wx = wx - zt5*zt3
         wy = wy - zt4*zt3
   10    continue
c mode numbers ky = 0, ny/2
         at1 = real(ffc(1,j,l))
         at2 = dkx*at1
         zt4 = cmplx(aimag(q(1,j,l)),-real(q(1,j,l)))
         zt4 = at2*zt4
         at2 = ci2*at2
         zt2 = cmplx(-aimag(cu(2,1,j,l)),real(cu(2,1,j,l)))
         zt3 = conjg(at2*zt2)
         wy = wy - zt4*zt3
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         at1 = real(ffc(k,1,l))
         at3 = dny*float(k - 1)*at1
         zt4 = cmplx(aimag(q(k,1,l)),-real(q(k,1,l)))
         zt5 = at3*zt4
         at3 = ci2*at3
         zt3 = cmplx(-aimag(cu(1,k,1,l)),real(cu(1,k,1,l)))
         zt3 = conjg(at3*zt3)
         wx = wx - zt5*zt3
   30    continue
      endif
   40 continue
   50 continue
      at1 = 2.0*float(nx)*float(ny)/affp
      sx = at1*wx
      sy = at1*wy
      return
      end
c-----------------------------------------------------------------------
      subroutine PSCFGUARD2(cus,cu,nyp,q2m0,nx,nxe,nypmx,nblok)
c initialize extended periodic field
c quadratic interpolation, for distributed data
      implicit none
      real cus, cu, q2m0
      integer nyp, nx, nxe, nypmx, nblok
      dimension cus(3,nxe,nypmx,nblok), cu(3,nxe,nypmx,nblok)
      dimension nyp(nblok)
c local data
      integer i, j, k, l, nyp3, nx3
c initialize extended field, with zero in the edges
      nx3 = nx + 3
      do 70 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 40 k = 1, nyp(l)
      do 20 j = 1, nx
      do 10 i = 1, 3
      cus(i,j+1,k+1,l) = -q2m0*cu(i,j+1,k+1,l)
   10 continue
   20 continue
      do 30 i = 1, 3
      cus(i,1,k+1,l) = 0.
      cus(i,nx+2,k+1,l) = 0.
      cus(i,nx+3,k+1,l) = 0.
   30 continue
   40 continue
      do 60 j = 1, nx3
      do 50 i = 1, 3
      cus(i,j,1,l) = 0.
      cus(i,j,nyp3-1,l) = 0.
      cus(i,j,nyp3,l) = 0.
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSCFGUARD22(cus,cu,nyp,q2m0,nx,nxe,nypmx,nblok)
c initialize extended periodic field
c quadratic interpolation, for distributed data
      implicit none
      real cus, cu, q2m0
      integer nyp, nx, nxe, nypmx, nblok
      dimension cus(2,nxe,nypmx,nblok), cu(2,nxe,nypmx,nblok)
      dimension nyp(nblok)
c local data
      integer i, j, k, l, nyp3, nx3
c initialize extended field, with zero in the edges
      nx3 = nx + 3
      do 70 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 40 k = 1, nyp(l)
      do 20 j = 1, nx
      do 10 i = 1, 2
      cus(i,j+1,k+1,l) = -q2m0*cu(i,j+1,k+1,l)
   10 continue
   20 continue
      do 30 i = 1, 2
      cus(i,1,k+1,l) = 0.
      cus(i,nx+2,k+1,l) = 0.
      cus(i,nx+3,k+1,l) = 0.
   30 continue
   40 continue
      do 60 j = 1, nx3
      do 50 i = 1, 2
      cus(i,j,1,l) = 0.
      cus(i,j,nyp3-1,l) = 0.
      cus(i,j,nyp3,l) = 0.
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSMCGUARD2(amu,nyp,x2y2m0,xym0,zxm0,zym0,nx,nxe,nypmx,n
     1blok)
c initialize extended periodic field
c quadratic interpolation, for distributed data
      implicit none
      real amu, x2y2m0, xym0, zxm0, zym0
      integer nyp, nx, nxe, nypmx, nblok
      dimension amu(4,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, j, k, l, nyp3, nx3
c initialize extended field, with zero in the edges
      nx3 = nx + 3
      do 60 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 30 k = 1, nyp(l)
      do 10 j = 1, nx
      amu(1,j+1,k+1,l) = x2y2m0
      amu(2,j+1,k+1,l) = xym0
      amu(3,j+1,k+1,l) = zxm0
      amu(4,j+1,k+1,l) = zym0
   10 continue
      do 20 i = 1, 4
      amu(i,1,k+1,l) = 0.
      amu(i,nx+2,k+1,l) = 0.
      amu(i,nx+3,k+1,l) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx3
      do 40 i = 1, 4
      amu(i,j,1,l) = 0.
      amu(i,j,nyp3-1,l) = 0.
      amu(i,j,nyp3,l) = 0.
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSMCGUARD22(amu,nyp,x2y2m0,xym0,nx,nxe,nypmx,nblok)
c initialize extended periodic field
c quadratic interpolation, for distributed data
      implicit none
      real amu, x2y2m0, xym0
      integer nyp, nx, nxe, nypmx, nblok
      dimension amu(2,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, j, k, l, nyp3, nx3
c initialize extended field, with zero in the edges
      nx3 = nx + 3
      do 60 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 30 k = 1, nyp(l)
      do 10 j = 1, nx
      amu(1,j+1,k+1,l) = x2y2m0
      amu(2,j+1,k+1,l) = xym0
   10 continue
      do 20 i = 1, 2
      amu(i,1,k+1,l) = 0.
      amu(i,nx+2,k+1,l) = 0.
      amu(i,nx+3,k+1,l) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx3
      do 40 i = 1, 2
      amu(i,j,1,l) = 0.
      amu(i,j,nyp3-1,l) = 0.
      amu(i,j,nyp3,l) = 0.
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAMCGUARD2X(amu,nyp,nx,nxe,nypmx,nblok,ndim)
c accumulate extended periodic tensor field
c quadratic interpolation, for distributed data
      implicit none
      real amu
      integer nyp, nx, nxe, nypmx, nblok, ndim
      dimension amu(ndim,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, k, l, nyp3
c accumulate edges of extended field
      do 30 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 20 k = 1, nyp3
      do 10 i = 1, ndim
      amu(i,2,k,l) = amu(i,2,k,l) + amu(i,nx+2,k,l)
      amu(i,3,k,l) = amu(i,3,k,l) + amu(i,nx+3,k,l)
      amu(i,nx+1,k,l) = amu(i,nx+1,k,l) + amu(i,1,k,l)
      amu(i,1,k,l) = 0.
      amu(i,nx+2,k,l) = 0.
      amu(i,nx+3,k,l) = 0.
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSCFGUARD2L(cus,cu,nyp,q2m0,nx,nxe,nypmx,nblok)
c initialize extended periodic field
c linear interpolation, for distributed data
      implicit none
      real cus, cu, q2m0
      integer nyp, nx, nxe, nypmx, nblok
      dimension cus(3,nxe,nypmx,nblok), cu(3,nxe,nypmx,nblok)
      dimension nyp(nblok)
c local data
      integer i, j, k, l, nyp1, nx1
c initialize extended field, with zero in the edges
      nx1 = nx + 1
      do 70 l = 1, nblok
      nyp1 = nyp(l) + 1
      do 40 k = 1, nyp(l)
      do 20 j = 1, nx
      do 10 i = 1, 3
      cus(i,j,k,l) = -q2m0*cu(i,j,k,l)
   10 continue
   20 continue
      do 30 i = 1, 3
      cus(i,nx+1,k,l) = 0.
   30 continue
   40 continue
      do 60 j = 1, nx1
      do 50 i = 1, 3
      cus(i,j,nyp1,l) = 0.
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSCFGUARD22L(cus,cu,nyp,q2m0,nx,nxe,nypmx,nblok)
c initialize extended periodic field
c linear interpolation, for distributed data
      implicit none
      real cus, cu, q2m0
      integer nyp, nx, nxe, nypmx, nblok
      dimension cus(2,nxe,nypmx,nblok), cu(2,nxe,nypmx,nblok)
      dimension nyp(nblok)
c local data
      integer i, j, k, l, nyp1, nx1
c initialize extended field, with zero in the edges
      nx1 = nx + 1
      do 70 l = 1, nblok
      nyp1 = nyp(l) + 1
      do 40 k = 1, nyp(l)
      do 20 j = 1, nx
      do 10 i = 1, 2
      cus(i,j,k,l) = -q2m0*cu(i,j,k,l)
   10 continue
   20 continue
      do 30 i = 1, 2
      cus(i,nx+1,k,l) = 0.
   30 continue
   40 continue
      do 60 j = 1, nx1
      do 50 i = 1, 2
      cus(i,j,nyp1,l) = 0.
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSMCGUARD2L(amu,nyp,x2y2m0,xym0,zxm0,zym0,nx,nxe,nypmx,
     1nblok)
c initialize extended periodic field
c linear interpolation, for distributed data
      implicit none
      real amu, x2y2m0, xym0, zxm0, zym0
      integer nyp, nx, nxe, nypmx, nblok
      dimension amu(4,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, j, k, l, nyp1, nx1
c initialize extended field, with zero in the edges
      nx1 = nx + 1
      do 60 l = 1, nblok
      nyp1 = nyp(l) + 1
      do 30 k = 1, nyp(l)
      do 10 j = 1, nx
      amu(1,j,k,l) = x2y2m0
      amu(2,j,k,l) = xym0
      amu(3,j,k,l) = zxm0
      amu(4,j,k,l) = zym0
   10 continue
      do 20 i = 1, 4
      amu(i,nx+1,k,l) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx1
      do 40 i = 1, 4
      amu(i,j,nyp1,l) = 0.
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSMCGUARD22L(amu,nyp,x2y2m0,xym0,nx,nxe,nypmx,nblok)
c initialize extended periodic field
c linear interpolation, for distributed data
      implicit none
      real amu, x2y2m0, xym0
      integer nyp, nx, nxe, nypmx, nblok
      dimension amu(2,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, j, k, l, nyp1, nx1
c initialize extended field, with zero in the edges
      nx1 = nx + 1
      do 60 l = 1, nblok
      nyp1 = nyp(l) + 1
      do 30 k = 1, nyp(l)
      do 10 j = 1, nx
      amu(1,j,k,l) = x2y2m0
      amu(2,j,k,l) = xym0
   10 continue
      do 20 i = 1, 2
      amu(i,nx+1,k,l) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx1
      do 40 i = 1, 2
      amu(i,j,nyp1,l) = 0.
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAMCGUARD2XL(amu,nyp,nx,nxe,nypmx,nblok,ndim)
c accumulate extended periodic tensor field
c linear interpolation, for distributed data
      implicit none
      real amu
      integer nyp, nx, nxe, nypmx, nblok, ndim
      dimension amu(ndim,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, k, l, nyp1
c accumulate edges of extended field
      do 30 l = 1, nblok
      nyp1 = nyp(l) + 1
      do 20 k = 1, nyp1
      do 10 i = 1, ndim
      amu(i,1,k,l) = amu(i,1,k,l) + amu(i,nx+1,k,l)
      amu(i,nx+1,k,l) = 0.
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp,jblok)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux
c in 2-1/2d with periodic boundary conditions.
c input: all, output: dcu
c approximate flop count is: 45*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the derivative of the current is calculated using the equations:
c dcu(1,kx,ky) = -sqrt(-1)*(kx*vx*vx+ky*vx*vy)
c dcu(2,kx,ky) = -sqrt(-1)*(kx*vx*vy+ky*vy*vy)
c dcu(3,kx,ky) = -sqrt(-1)*(kx*vx*vz+ky*vy*vz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,ky=pi) = dcu(i,kx=0,ky=0) = 0.
c the transverse part is calculated using the equation:
c dcu(1,kx,ky) = dcu(1,kx,ky)-kx*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c dcu(2,kx,ky) = dcu(2,kx,ky)-ky*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c on output:
c dcu(i,k,j,l) = i-th component of transverse part of complex derivative
c of current for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c amu(1,k,j,l) = xx component of complex momentum flux
c amu(2,k,j,l) = xy component of complex momentum flux
c amu(3,k,j,l) = zx component of complex momentum flux
c amu(4,k,j,l) = zy component of complex momentum flux
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = second dimension of field arrays, must be >= ny
c kxp = number of data values per block
c jblok = number of data blocks
      implicit none
      integer nx, ny, kstrt, nyv, kxp, jblok
      complex dcu, amu
      dimension dcu(3,nyv,kxp,jblok), amu(4,nyv,kxp,jblok)
c local data
      integer nxh, nyh, ny2, ks, joff, j, k, k1, l
      real dnx, dny, dkx, dky, dkx2, dky2, dkxy, dkxy2, at1
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
      if (kstrt.gt.nxh) return
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         dky2 = dky*dky
         dkxy = dkx*dky
         dkxy2 = dky2 - dkx2
         at1 = 1.0/(dkx2 + dky2)
         zt1 = cmplx(aimag(amu(1,k,j,l)),-real(amu(1,k,j,l)))
         zt2 = cmplx(aimag(amu(2,k,j,l)),-real(amu(2,k,j,l)))
         zt3 = at1*(dkxy*zt1 + dkxy2*zt2)
         dcu(1,k,j,l) = dky*zt3
         dcu(2,k,j,l) = -dkx*zt3
         zt1 = cmplx(aimag(amu(3,k,j,l)),-real(amu(3,k,j,l)))
         zt2 = cmplx(aimag(amu(4,k,j,l)),-real(amu(4,k,j,l)))
         dcu(3,k,j,l) = dkx*zt1 + dky*zt2
         zt1 = cmplx(aimag(amu(1,k1,j,l)),-real(amu(1,k1,j,l)))
         zt2 = cmplx(aimag(amu(2,k1,j,l)),-real(amu(2,k1,j,l)))
         zt3 = at1*(dkxy*zt1 - dkxy2*zt2)
         dcu(1,k1,j,l) = dky*zt3
         dcu(2,k1,j,l) = dkx*zt3
         zt1 = cmplx(aimag(amu(3,k1,j,l)),-real(amu(3,k1,j,l)))
         zt2 = cmplx(aimag(amu(4,k1,j,l)),-real(amu(4,k1,j,l)))
         dcu(3,k1,j,l) = dkx*zt1 - dky*zt2
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         zt2 = cmplx(aimag(amu(2,1,j,l)),-real(amu(2,1,j,l)))
         dcu(1,1,j,l) = zero
         dcu(2,1,j,l) = dkx*zt2
         zt1 = cmplx(aimag(amu(3,1,j,l)),-real(amu(3,1,j,l)))
         dcu(3,1,j,l) = dkx*zt1
         dcu(1,k1,j,l) = zero
         dcu(2,k1,j,l) = zero
         dcu(3,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         zt2 = cmplx(aimag(amu(2,k,1,l)),-real(amu(2,k,1,l)))
         dcu(1,k,1,l) = dky*zt2
         dcu(2,k,1,l) = zero
         zt2 = cmplx(aimag(amu(4,k,1,l)),-real(amu(4,k,1,l)))
         dcu(3,k,1,l) = dky*zt2
         dcu(1,k1,1,l) = zero
         dcu(2,k1,1,l) = zero
         dcu(3,k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         dcu(1,1,1,l) = zero
         dcu(2,1,1,l) = zero
         dcu(3,1,1,l) = zero
         dcu(1,k1,1,l) = zero
         dcu(2,k1,1,l) = zero
         dcu(3,k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDCUPERP22(dcu,amu,nx,ny,kstrt,nyv,kxp,jblok)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux
c in 2d with periodic boundary conditions.
c input: all, output: dcu
c approximate flop count is: 29*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the derivative of the current is calculated using the equations:
c dcu(1,kx,ky) = -sqrt(-1)*(kx*vx*vx+ky*vx*vy)
c dcu(2,kx,ky) = -sqrt(-1)*(kx*vx*vy+ky*vy*vy)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,ky=pi) = dcu(i,kx=0,ky=0) = 0.
c the transverse part is calculated using the equation:
c dcu(1,kx,ky) = dcu(1,kx,ky)-kx*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c dcu(2,kx,ky) = dcu(2,kx,ky)-ky*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c on output:
c dcu(i,k,j,l) = i-th component of transverse part of complex derivative
c of current for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c amu(1,k,j,l) = xx component of complex momentum flux
c amu(2,k,j,l) = xy component of complex momentum flux
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = second dimension of field arrays, must be >= ny
c kxp = number of data values per block
c jblok = number of data blocks
      implicit none
      integer nx, ny, kstrt, nyv, kxp, jblok
      complex dcu, amu
      dimension dcu(2,nyv,kxp,jblok), amu(2,nyv,kxp,jblok)
c local data
      integer nxh, nyh, ny2, ks, joff, j, k, k1, l
      real dnx, dny, dkx, dky, dkx2, dky2, dkxy, dkxy2, at1
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
      if (kstrt.gt.nxh) return
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         dky2 = dky*dky
         dkxy = dkx*dky
         dkxy2 = dky2 - dkx2
         at1 = 1.0/(dkx2 + dky2)
         zt1 = cmplx(aimag(amu(1,k,j,l)),-real(amu(1,k,j,l)))
         zt2 = cmplx(aimag(amu(2,k,j,l)),-real(amu(2,k,j,l)))
         zt3 = at1*(dkxy*zt1 + dkxy2*zt2)
         dcu(1,k,j,l) = dky*zt3
         dcu(2,k,j,l) = -dkx*zt3
         zt1 = cmplx(aimag(amu(1,k1,j,l)),-real(amu(1,k1,j,l)))
         zt2 = cmplx(aimag(amu(2,k1,j,l)),-real(amu(2,k1,j,l)))
         zt3 = at1*(dkxy*zt1 - dkxy2*zt2)
         dcu(1,k1,j,l) = dky*zt3
         dcu(2,k1,j,l) = dkx*zt3
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         zt2 = cmplx(aimag(amu(2,1,j,l)),-real(amu(2,1,j,l)))
         dcu(1,1,j,l) = zero
         dcu(2,1,j,l) = dkx*zt2
         dcu(1,k1,j,l) = zero
         dcu(2,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         zt2 = cmplx(aimag(amu(2,k,1,l)),-real(amu(2,k,1,l)))
         dcu(1,k,1,l) = dky*zt2
         dcu(2,k,1,l) = zero
         dcu(1,k1,1,l) = zero
         dcu(2,k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         dcu(1,1,1,l) = zero
         dcu(2,1,1,l) = zero
         dcu(1,k1,1,l) = zero
         dcu(2,k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PADCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp,jblok)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux and acceleration density
c in 2-1/2d with periodic boundary conditions.
c input: all, output: dcu
c approximate flop count is: 65*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the derivative of the current is calculated using the equations:
c dcu(1,kx,ky) = dcu(1,kx,ky)-sqrt(-1)*(kx*vx*vx+ky*vx*vy)
c dcu(2,kx,ky) = dcu(2,kx,ky)-sqrt(-1)*(kx*vx*vy+ky*vy*vy)
c dcu(3,kx,ky) = dcu(3,kx,ky)-sqrt(-1)*(kx*vx*vz+ky*vy*vz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,ky=pi) = dcu(i,kx=0,ky=0) = 0.
c the transverse part is calculated using the equation:
c dcu(1,kx,ky) = dcu(1,kx,ky)-kx*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c dcu(2,kx,ky) = dcu(2,kx,ky)-ky*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c on input:
c dcu(i,j,k,l) = complex acceleration density for fourier mode (jj-1,k-1)
c on output:
c dcu(i,k,j,l) = i-th component of transverse part of complex derivative
c of current for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c amu(1,k,j,l) = xx component of complex momentum flux
c amu(2,k,j,l) = xy component of complex momentum flux
c amu(3,k,j,l) = zx component of complex momentum flux
c amu(4,k,j,l) = zy component of complex momentum flux
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = second dimension of field arrays, must be >= ny
c kxp = number of data values per block
c jblok = number of data blocks
      implicit none
      integer nx, ny, kstrt, nyv, kxp, jblok
      complex dcu, amu
      dimension dcu(3,nyv,kxp,jblok), amu(4,nyv,kxp,jblok)
c local data
      integer nxh, nyh, ny2, ks, joff, j, k, k1, l
      real dnx, dny, dkx, dky, dkx2, dky2, dkxy, dkxy2, at1
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
      if (kstrt.gt.nxh) return
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         dky2 = dky*dky
         dkxy = dkx*dky
         dkxy2 = dky2 - dkx2
         at1 = 1.0/(dkx2 + dky2)
         zt1 = cmplx(aimag(amu(1,k,j,l)),-real(amu(1,k,j,l)))
         zt2 = cmplx(aimag(amu(2,k,j,l)),-real(amu(2,k,j,l)))
         zt3 = at1*(dky*dcu(1,k,j,l) - dkx*dcu(2,k,j,l) + dkxy*zt1 + dkx
     1y2*zt2)
         dcu(1,k,j,l) = dky*zt3
         dcu(2,k,j,l) = -dkx*zt3
         zt1 = cmplx(aimag(amu(3,k,j,l)),-real(amu(3,k,j,l)))
         zt2 = cmplx(aimag(amu(4,k,j,l)),-real(amu(4,k,j,l)))
         dcu(3,k,j,l) = dcu(3,k,j,l) + dkx*zt1 + dky*zt2
         zt1 = cmplx(aimag(amu(1,k1,j,l)),-real(amu(1,k1,j,l)))
         zt2 = cmplx(aimag(amu(2,k1,j,l)),-real(amu(2,k1,j,l)))
         zt3 = at1*(dky*dcu(1,k1,j,l) + dkx*dcu(2,k1,j,l) + dkxy*zt1 - d
     1kxy2*zt2)
         dcu(1,k1,j,l) = dky*zt3
         dcu(2,k1,j,l) = dkx*zt3
         zt1 = cmplx(aimag(amu(3,k1,j,l)),-real(amu(3,k1,j,l)))
         zt2 = cmplx(aimag(amu(4,k1,j,l)),-real(amu(4,k1,j,l)))
         dcu(3,k1,j,l) = dcu(3,k1,j,l) + dkx*zt1 - dky*zt2
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         zt2 = cmplx(aimag(amu(2,1,j,l)),-real(amu(2,1,j,l)))
         dcu(1,1,j,l) = zero
         dcu(2,1,j,l) = dcu(2,1,j,l) + dkx*zt2
         zt1 = cmplx(aimag(amu(3,1,j,l)),-real(amu(3,1,j,l)))
         dcu(3,1,j,l) = dcu(3,1,j,l) + dkx*zt1
         dcu(1,k1,j,l) = zero
         dcu(2,k1,j,l) = zero
         dcu(3,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         zt2 = cmplx(aimag(amu(2,k,1,l)),-real(amu(2,k,1,l)))
         dcu(1,k,1,l) = dcu(1,k,1,l) + dky*zt2
         dcu(2,k,1,l) = zero
         zt2 = cmplx(aimag(amu(4,k,1,l)),-real(amu(4,k,1,l)))
         dcu(3,k,1,l) = dcu(3,k,1,l) + dky*zt2
         dcu(1,k1,1,l) = zero
         dcu(2,k1,1,l) = zero
         dcu(3,k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         dcu(1,1,1,l) = zero
         dcu(2,1,1,l) = zero
         dcu(1,k1,1,l) = zero
         dcu(2,k1,1,l) = zero
         dcu(3,k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PADCUPERP22(dcu,amu,nx,ny,kstrt,nyv,kxp,jblok)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux and acceleration density
c in 2d with periodic boundary conditions.
c input: all, output: dcu
c approximate flop count is: 45*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the derivative of the current is calculated using the equations:
c dcu(1,kx,ky) = dcu(1,kx,ky)-sqrt(-1)*(kx*vx*vx+ky*vx*vy)
c dcu(2,kx,ky) = dcu(2,kx,ky)-sqrt(-1)*(kx*vx*vy+ky*vy*vy)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,ky=pi) = dcu(i,kx=0,ky=0) = 0.
c the transverse part is calculated using the equation:
c dcu(1,kx,ky) = dcu(1,kx,ky)-kx*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c dcu(2,kx,ky) = dcu(2,kx,ky)-ky*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c on input:
c dcu(i,j,k,l) = complex acceleration density for fourier mode (jj-1,k-1)
c on output:
c dcu(i,k,j,l) = i-th component of transverse part of complex derivative
c of current for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c amu(1,k,j,l) = xx component of complex momentum flux
c amu(2,k,j,l) = xy component of complex momentum flux
c amu(3,k,j,l) = zx component of complex momentum flux
c amu(4,k,j,l) = zy component of complex momentum flux
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = second dimension of field arrays, must be >= ny
c kxp = number of data values per block
c jblok = number of data blocks
      implicit none
      integer nx, ny, kstrt, nyv, kxp, jblok
      complex dcu, amu
      dimension dcu(2,nyv,kxp,jblok), amu(2,nyv,kxp,jblok)
c local data
      integer nxh, nyh, ny2, ks, joff, j, k, k1, l
      real dnx, dny, dkx, dky, dkx2, dky2, dkxy, dkxy2, at1
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
      if (kstrt.gt.nxh) return
      do 40 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         dky2 = dky*dky
         dkxy = dkx*dky
         dkxy2 = dky2 - dkx2
         at1 = 1.0/(dkx2 + dky2)
         zt1 = cmplx(aimag(amu(1,k,j,l)),-real(amu(1,k,j,l)))
         zt2 = cmplx(aimag(amu(2,k,j,l)),-real(amu(2,k,j,l)))
         zt3 = at1*(dky*dcu(1,k,j,l) - dkx*dcu(2,k,j,l) + dkxy*zt1 + dkx
     1y2*zt2)
         dcu(1,k,j,l) = dky*zt3
         dcu(2,k,j,l) = -dkx*zt3
         zt1 = cmplx(aimag(amu(1,k1,j,l)),-real(amu(1,k1,j,l)))
         zt2 = cmplx(aimag(amu(2,k1,j,l)),-real(amu(2,k1,j,l)))
         zt3 = at1*(dky*dcu(1,k1,j,l) + dkx*dcu(2,k1,j,l) + dkxy*zt1 - d
     1kxy2*zt2)
         dcu(1,k1,j,l) = dky*zt3
         dcu(2,k1,j,l) = dkx*zt3
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         zt2 = cmplx(aimag(amu(2,1,j,l)),-real(amu(2,1,j,l)))
         dcu(1,1,j,l) = zero
         dcu(2,1,j,l) = dcu(2,1,j,l) + dkx*zt2
         dcu(1,k1,j,l) = zero
         dcu(2,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         zt2 = cmplx(aimag(amu(2,k,1,l)),-real(amu(2,k,1,l)))
         dcu(1,k,1,l) = dcu(1,k,1,l) + dky*zt2
         dcu(2,k,1,l) = zero
         dcu(1,k1,1,l) = zero
         dcu(2,k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         dcu(1,1,1,l) = zero
         dcu(2,1,1,l) = zero
         dcu(1,k1,1,l) = zero
         dcu(2,k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PEPOISP23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny,
     1kstrt,nyv,kxp,jblok,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c transverse electric field (or convolution of transverse electric field
c over particle shape), with periodic boundary conditions.
c using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
c A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
c for isign = 0, input: isign,ax,ay,affp,wp0,nx,ny,kstrt,nyv,kxp,jblok,
c nyhd, output: ffe
c for isign /= 0, input: dcu,ffe,isign,affp,ci,nx,ny,kstrt,nyv,kxp,
c jblok,nyhd, output: exy,wf
c approximate flop count is: 59*nxc*nyc + 32*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign = -1, smoothed transverse electric field is calculated
c using the equations:
c ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)*s(kx,ky)
c ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = -ci*ci*g(kx,ky)*dcuz(kx,ky)*s(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
c ex(ky=pi) = ey(ky=pi) = ez(ky=pi) = 0,
c ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
c if isign = 1, unsmoothed transverse electric field is calculated
c using the equations:
c ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)
c ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)
c ez(kx,ky) = -ci*ci*g(kx,ky)*dcuz(kx,ky)
c dcu(i,k,j,l) = i-th component of transverse part of complex derivative
c of current,
c exy(i,k,j,l) = i-th component of complex transverse electric field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c jblok = number of data blocks
c kxp = number of data values per block
c kstrt = starting data block number
c aimag(ffe(k,j,l)) = finite-size particle shape factor s
c real(ffe(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c wp0 = normalized total plasma frequency squared
c ci = reciprical of velocity of light
c transverse electric field energy is also calculated, using
c wf = nx*ny*sum((affp/((kx**2+ky**2)*ci*ci)**2)
c    |dcu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the derivative of current is
c divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
      real ax, ay, affp, wp0, ci, wf
      complex dcu, exy, ffe
      dimension dcu(3,nyv,kxp,jblok)
      dimension exy(3,nyv,kxp,jblok)
      dimension ffe(nyhd,kxp,jblok)
c local data
      integer nxh, nyh, ny2, ks, joff, j, k, k1, l
      real dnx, dny, ci2, wpc, dkx, dky, at1, at2, at3, at4
      complex zero
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 40
      if (kstrt.gt.nxh) return
      wpc = wp0*ci2
c prepare form factor array
      do 30 l = 1, jblok
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffe(k,j,l) = cmplx(affp,1.)
      else
         ffe(k,j,l) = cmplx(affp*at4/(at3 + wpc*at4*at4),at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 100
c calculate smoothed transverse electric field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nxh) go to 90
      do 80 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 60 j = 1, kxp
      if ((j+joff).gt.0) then
         do 50 k = 2, nyh
         k1 = ny2 - k
         at2 = -ci2*real(ffe(k,j,l))
         at1 = at2*aimag(ffe(k,j,l))
         at2 = at2*at2
         exy(1,k,j,l) = at1*dcu(1,k,j,l)
         exy(2,k,j,l) = at1*dcu(2,k,j,l)
         exy(3,k,j,l) = at1*dcu(3,k,j,l)
         exy(1,k1,j,l) = at1*dcu(1,k1,j,l)
         exy(2,k1,j,l) = at1*dcu(2,k1,j,l)
         exy(3,k1,j,l) = at1*dcu(3,k1,j,l)
         wp = wp + at2*(dcu(1,k,j,l)*conjg(dcu(1,k,j,l)) + dcu(2,k,j,l)*
     1conjg(dcu(2,k,j,l)) + dcu(3,k,j,l)*conjg(dcu(3,k,j,l)) + dcu(1,k1,
     2j,l)*conjg(dcu(1,k1,j,l)) + dcu(2,k1,j,l)*conjg(dcu(2,k1,j,l)) + d
     3cu(3,k1,j,l)*conjg(dcu(3,k1,j,l)))
   50    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = -ci2*real(ffe(1,j,l))
         at1 = at2*aimag(ffe(1,j,l))
         at2 = at2*at2
         exy(1,1,j,l) = at1*dcu(1,1,j,l)
         exy(2,1,j,l) = at1*dcu(2,1,j,l)
         exy(3,1,j,l) = at1*dcu(3,1,j,l)
         exy(1,k1,j,l) = zero
         exy(2,k1,j,l) = zero
         exy(3,k1,j,l) = zero
         wp = wp + at2*(dcu(1,1,j,l)*conjg(dcu(1,1,j,l)) + dcu(2,1,j,l)*
     1conjg(dcu(2,1,j,l)) + dcu(3,1,j,l)*conjg(dcu(3,1,j,l)))
      endif
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, nyh
         k1 = ny2 - k
         at2 = -ci2*real(ffe(k,1,l))
         at1 = at2*aimag(ffe(k,1,l))
         at2 = at2*at2
         exy(1,k,1,l) = at1*dcu(1,k,1,l)
         exy(2,k,1,l) = at1*dcu(2,k,1,l)
         exy(3,k,1,l) = at1*dcu(3,k,1,l)
         exy(1,k1,1,l) = zero
         exy(2,k1,1,l) = zero
         exy(3,k1,1,l) = zero
         wp = wp + at2*(dcu(1,k,1,l)*conjg(dcu(1,k,1,l)) + dcu(2,k,1,l)*
     1conjg(dcu(2,k,1,l)) + dcu(3,k,1,l)*conjg(dcu(3,k,1,l)))
   70    continue
         k1 = nyh + 1
         exy(1,1,1,l) = zero
         exy(2,1,1,l) = zero
         exy(3,1,1,l) = zero
         exy(1,k1,1,l) = zero
         exy(2,k1,1,l) = zero
         exy(3,k1,1,l) = zero
      endif
   80 continue
   90 continue
      wf = float(nx)*float(ny)*wp/affp
      return
c calculate unsmoothed transverse electric field and sum field energy
  100 wp = 0.0d0
      if (kstrt.gt.nxh) go to 150
      do 140 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 120 j = 1, kxp
      if ((j+joff).gt.0) then
         do 110 k = 2, nyh
         k1 = ny2 - k
         at2 = -ci2*real(ffe(k,j,l))
         at1 = at2*at2
         exy(1,k,j,l) = at2*dcu(1,k,j,l)
         exy(2,k,j,l) = at2*dcu(2,k,j,l)
         exy(3,k,j,l) = at2*dcu(3,k,j,l)
         exy(1,k1,j,l) = at2*dcu(1,k1,j,l)
         exy(2,k1,j,l) = at2*dcu(2,k1,j,l)
         exy(3,k1,j,l) = at2*dcu(3,k1,j,l)
         wp = wp + at1*(dcu(1,k,j,l)*conjg(dcu(1,k,j,l)) + dcu(2,k,j,l)*
     1conjg(dcu(2,k,j,l)) + dcu(3,k,j,l)*conjg(dcu(3,k,j,l)) + dcu(1,k1,
     2j,l)*conjg(dcu(1,k1,j,l)) + dcu(2,k1,j,l)*conjg(dcu(2,k1,j,l)) + d
     3cu(3,k1,j,l)*conjg(dcu(3,k1,j,l)))
  110    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = -ci2*real(ffe(1,j,l))
         at1 = at2*at2
         exy(1,1,j,l) = at2*dcu(1,1,j,l)
         exy(2,1,j,l) = at2*dcu(2,1,j,l)
         exy(3,1,j,l) = at2*dcu(3,1,j,l)
         exy(1,k1,j,l) = zero
         exy(2,k1,j,l) = zero
         exy(3,k1,j,l) = zero
         wp = wp + at1*(dcu(1,1,j,l)*conjg(dcu(1,1,j,l)) + dcu(2,1,j,l)*
     1conjg(dcu(2,1,j,l)) + dcu(3,1,j,l)*conjg(dcu(3,1,j,l)))
      endif
  120 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 130 k = 2, nyh
         k1 = ny2 - k
         at2 = -ci2*real(ffe(k,1,l))
         at1 = at2*at2
         exy(1,k,1,l) = at2*dcu(1,k,1,l)
         exy(2,k,1,l) = at2*dcu(2,k,1,l)
         exy(3,k,1,l) = at2*dcu(3,k,1,l)
         exy(1,k1,1,l) = zero
         exy(2,k1,1,l) = zero
         exy(3,k1,1,l) = zero
         wp = wp + at1*(dcu(1,k,1,l)*conjg(dcu(1,k,1,l)) + dcu(2,k,1,l)*
     1conjg(dcu(2,k,1,l)) + dcu(3,k,1,l)*conjg(dcu(3,k,1,l)))
  130    continue
         k1 = nyh + 1
         exy(1,1,1,l) = zero
         exy(2,1,1,l) = zero
         exy(3,1,1,l) = zero
         exy(1,k1,1,l) = zero
         exy(2,k1,1,l) = zero
         exy(3,k1,1,l) = zero
      endif
  140 continue
  150 continue
      wf = float(nx)*float(ny)*wp/affp
      return
      end
c-----------------------------------------------------------------------
      subroutine PEPOISP22(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny,
     1kstrt,nyv,kxp,jblok,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c transverse electric field (or convolution of transverse electric field
c over particle shape), with periodic boundary conditions.
c using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
c A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
c for isign = 0, input: isign,ax,ay,affp,wp0,nx,ny,kstrt,nyv,kxp,jblok,
c nyhd, output: ffe
c for isign /= 0, input: dcu,ffe,isign,affp,ci,nx,ny,kstrt,nyv,kxp,
c jblok,nyhd, output: exy,wf
c approximate flop count is: 41*nxc*nyc + 23*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign = -1, smoothed transverse electric field is calculated
c using the equations:
c ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)*s(kx,ky)
c ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)*s(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c ex(kx=pi) = ey(kx=pi) = 0, ex(ky=pi) = ey(ky=pi) = 0, and
c ex(kx=0,ky=0) = ey(kx=0,ky=0) = 0.
c if isign = 1, unsmoothed transverse electric field is calculated
c using the equations:
c ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)
c ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)
c dcu(i,k,j,l) = i-th component of transverse part of complex derivative
c of current,
c exy(i,k,j,l) = i-th component of complex transverse electric field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c jblok = number of data blocks
c kxp = number of data values per block
c kstrt = starting data block number
c aimag(ffe(k,j,l)) = finite-size particle shape factor s
c real(ffe(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c wp0 = normalized total plasma frequency squared
c ci = reciprical of velocity of light
c transverse electric field energy is also calculated, using
c wf = nx*ny*sum((affp/((kx**2+ky**2)*ci*ci)**2)
c    |dcu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the derivative of current is
c divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
      real ax, ay, affp, wp0, ci, wf
      complex dcu, exy, ffe
      dimension dcu(2,nyv,kxp,jblok)
      dimension exy(2,nyv,kxp,jblok)
      dimension ffe(nyhd,kxp,jblok)
c local data
      integer nxh, nyh, ny2, ks, joff, j, k, k1, l
      real dnx, dny, ci2, wpc, dkx, dky, at1, at2, at3, at4
      complex zero
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 40
      if (kstrt.gt.nxh) return
      wpc = wp0*ci2
c prepare form factor array
      do 30 l = 1, jblok
      joff = kxp*(l + ks) - 1
      do 20 j = 1, kxp
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffe(k,j,l) = cmplx(affp,1.)
      else
         ffe(k,j,l) = cmplx(affp*at4/(at3 + wpc*at4*at4),at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 100
c calculate smoothed transverse electric field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nxh) go to 90
      do 80 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 60 j = 1, kxp
      if ((j+joff).gt.0) then
         do 50 k = 2, nyh
         k1 = ny2 - k
         at2 = -ci2*real(ffe(k,j,l))
         at1 = at2*aimag(ffe(k,j,l))
         at2 = at2*at2
         exy(1,k,j,l) = at1*dcu(1,k,j,l)
         exy(2,k,j,l) = at1*dcu(2,k,j,l)
         exy(1,k1,j,l) = at1*dcu(1,k1,j,l)
         exy(2,k1,j,l) = at1*dcu(2,k1,j,l)
         wp = wp + at2*(dcu(1,k,j,l)*conjg(dcu(1,k,j,l)) + dcu(2,k,j,l)*
     1conjg(dcu(2,k,j,l)) + dcu(1,k1,j,l)*conjg(dcu(1,k1,j,l)) + dcu(2,k
     21,j,l)*conjg(dcu(2,k1,j,l)))
   50    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = -ci2*real(ffe(1,j,l))
         at1 = at2*aimag(ffe(1,j,l))
         at2 = at2*at2
         exy(1,1,j,l) = at1*dcu(1,1,j,l)
         exy(2,1,j,l) = at1*dcu(2,1,j,l)
         exy(1,k1,j,l) = zero
         exy(2,k1,j,l) = zero
         wp = wp + at2*(dcu(1,1,j,l)*conjg(dcu(1,1,j,l)) + dcu(2,1,j,l)*
     1conjg(dcu(2,1,j,l)))
      endif
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, nyh
         k1 = ny2 - k
         at2 = -ci2*real(ffe(k,1,l))
         at1 = at2*aimag(ffe(k,1,l))
         at2 = at2*at2
         exy(1,k,1,l) = at1*dcu(1,k,1,l)
         exy(2,k,1,l) = at1*dcu(2,k,1,l)
         exy(1,k1,1,l) = zero
         exy(2,k1,1,l) = zero
         wp = wp + at2*(dcu(1,k,1,l)*conjg(dcu(1,k,1,l)) + dcu(2,k,1,l)*
     1conjg(dcu(2,k,1,l)))
   70    continue
         k1 = nyh + 1
         exy(1,1,1,l) = zero
         exy(2,1,1,l) = zero
         exy(1,k1,1,l) = zero
         exy(2,k1,1,l) = zero
      endif
   80 continue
   90 continue
      wf = float(nx)*float(ny)*wp/affp
      return
c calculate unsmoothed transverse electric field and sum field energy
  100 wp = 0.0d0
      if (kstrt.gt.nxh) go to 150
      do 140 l = 1, jblok
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      joff = kxp*(l + ks) - 1
      do 120 j = 1, kxp
      if ((j+joff).gt.0) then
         do 110 k = 2, nyh
         k1 = ny2 - k
         at2 = -ci2*real(ffe(k,j,l))
         at1 = at2*at2
         exy(1,k,j,l) = at2*dcu(1,k,j,l)
         exy(2,k,j,l) = at2*dcu(2,k,j,l)
         exy(1,k1,j,l) = at2*dcu(1,k1,j,l)
         exy(2,k1,j,l) = at2*dcu(2,k1,j,l)
         wp = wp + at2*(dcu(1,k,j,l)*conjg(dcu(1,k,j,l)) + dcu(2,k,j,l)*
     1conjg(dcu(2,k,j,l)) + dcu(1,k1,j,l)*conjg(dcu(1,k1,j,l)) + dcu(2,k
     21,j,l)*conjg(dcu(2,k1,j,l)))
  110    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = -ci2*real(ffe(1,j,l))
         at1 = at2*at2
         exy(1,1,j,l) = at2*dcu(1,1,j,l)
         exy(2,1,j,l) = at2*dcu(2,1,j,l)
         exy(1,k1,j,l) = zero
         exy(2,k1,j,l) = zero
         wp = wp + at1*(dcu(1,1,j,l)*conjg(dcu(1,1,j,l)) + dcu(2,1,j,l)*
     1conjg(dcu(2,1,j,l)))
      endif
  120 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 130 k = 2, nyh
         k1 = ny2 - k
         at2 = -ci2*real(ffe(k,1,l))
         at1 = at2*at2
         exy(1,k,1,l) = at2*dcu(1,k,1,l)
         exy(2,k,1,l) = at2*dcu(2,k,1,l)
         exy(1,k1,1,l) = zero
         exy(2,k1,1,l) = zero
         wp = wp + at1*(dcu(1,k,1,l)*conjg(dcu(1,k,1,l)) + dcu(2,k,1,l)*
     1conjg(dcu(2,k,1,l)))
  130    continue
         k1 = nyh + 1
         exy(1,1,1,l) = zero
         exy(2,1,1,l) = zero
         exy(1,k1,1,l) = zero
         exy(2,k1,1,l) = zero
      endif
  140 continue
  150 continue
      wf = float(nx)*float(ny)*wp/affp
      return
      end
c-----------------------------------------------------------------------
      subroutine PADDQEI2(qe,qi,nyp,nx,nxe,nypmx,nblok)
c adds electron and ion densities
c assumes guard cells have already been added
c qe/qi = charge density for electrons/ions
c nyp(l) = number of primary gridpoints in particle partition l
c nx = system length in x direction
c nxe = first dimension of field array, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c nblok = number of particle partitions.
      implicit none
      integer nx, nxe, nypmx, nblok
      real qe, qi
      integer nyp
      dimension qe(nxe,nypmx,nblok), qi(nxe,nypmx,nblok)
      dimension nyp(nblok)
c local data
      integer j, k, l
      do 30 l = 1, nblok
      do 20 k = 1, nyp(l)
      do 10 j = 1, nx
      qe(j,k,l) = qe(j,k,l) + qi(j,k,l)
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PADDQEI2X(qe,qi,nyp,qbme,qbmi,wpmax,wpmin,nx,nxe,nypmx,
     1nblok)
c adds electron and ion densities, and calculates maximum and minimum
c plasma frequency.  assumes guard cells have already been added
c qe/qi = charge density for electrons/ions
c nyp(l) = number of primary gridpoints in particle partition l
c qbme/qbmi = charge/mass ratio for electrons/ions
c wpmax/wpmin = maximum/minimum plasma frequency
c nx = system length in x direction
c nxe = first dimension of field array, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c nblok = number of particle partitions.
      implicit none
      integer nx, nxe, nypmx, nblok
      real qbme, qbmi, wpmax, wpmin
      real qe, qi
      integer nyp
      dimension qe(nxe,nypmx,nblok), qi(nxe,nypmx,nblok)
      dimension nyp(nblok)
c local data
      integer j, k, l
      real at1
      double precision sum1
      sum1 = 0.0d0
      wpmax = qbme*qe(1,1,1) + qbmi*qi(1,1,1)
      wpmin = wpmax
      do 30 l = 1, nblok
      do 20 k = 1, nyp(l)
      do 10 j = 1, nx
      at1 = qbme*qe(j,k,l) + qbmi*qi(j,k,l)
      qe(j,k,l) = qe(j,k,l) + qi(j,k,l)
      wpmax = max(wpmax,at1)
      wpmin = min(wpmin,at1)
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBADDEXT2(bxy,nyp,omx,omy,omz,nx,nxe,nypmx,nblok)
c adds constant to magnetic field for 2-1/2d code
c bxy = magnetic field
c nyp(l) = number of primary gridpoints in particle partition l
c omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z
c nx = system length in x direction
c nxe = first dimension of field array, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c nblok = number of particle partitions.
      implicit none
      integer nx, nxe, nypmx, nblok
      real bxy, omx, omy, omz
      integer nyp
      dimension bxy(3,nxe,nypmx,nblok)
      dimension nyp(nblok)
c local data
      integer j, k, l
      do 30 l = 1, nblok
      do 20 k = 1, nyp(l)
      do 10 j = 1, nx
      bxy(1,j,k,l) = bxy(1,j,k,l) + omx
      bxy(2,j,k,l) = bxy(2,j,k,l) + omy
      bxy(3,j,k,l) = bxy(3,j,k,l) + omz
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBADDEXT22(bz,nyp,omz,nx,nxe,nypmx,nblok)
c adds constant to magnetic field for 2d code
c bz = magnetic field
c nyp(l) = number of primary gridpoints in particle partition l
c omz = magnetic field electron cyclotron frequency inz
c nx = system length in x direction
c nxe = first dimension of field array, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c nblok = number of particle partitions.
      implicit none
      integer nx, nxe, nypmx, nblok
      real bz, omz
      integer nyp
      dimension bz(nxe,nypmx,nblok)
      dimension nyp(nblok)
c local data
      integer j, k, l
      do 30 l = 1, nblok
      do 20 k = 1, nyp(l)
      do 10 j = 1, nx
      bz(j,k,l) = bz(j,k,l) + omz
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PIMOMENT2(qi,fxy,nyp,pxi,pyi,pzi,dt,nx,nxe,nypmx,nblok)
c calculate ion momentum for 2-1/2d code from integral of qi*fxy
c assumes guard cells have already been added
      implicit none
      integer nx, nxe, nypmx, nblok
      real pxi, pyi, pzi, dt
      real qi, fxy
      integer nyp
      dimension qi(nxe,nypmx,nblok)
      dimension fxy(3,nxe,nypmx,nblok)
      dimension nyp(nblok)
c local data
      integer j, k, l
      double precision dt1, sum1, sum2, sum3
      sum1 = 0.0d0
      sum2 = 0.0d0
      sum3 = 0.0d0
      do 30 l = 1, nblok
      do 20 k = 1, nyp(l)
      do 10 j = 1, nx
      dt1 = dble(qi(j,k,l))
      sum1 = sum1 + dt1*fxy(1,j,k,l)
      sum2 = sum2 + dt1*fxy(2,j,k,l)
      sum3 = sum3 + dt1*fxy(3,j,k,l)
   10 continue
   20 continue
   30 continue
      pxi = sum1*dt
      pyi = sum2*dt
      pzi = sum3*dt
      return
      end
c-----------------------------------------------------------------------
      subroutine PIMOMENT22(qi,fxy,nyp,pxi,pyi,dt,nx,nxe,nypmx,nblok)
c calculate ion momentum for 2d code from integral of qi*fxy
c assumes guard cells have already been added
      implicit none
      integer nx, nxe, nypmx, nblok
      real pxi, pyi, dt
      real qi, fxy
      integer nyp
      dimension qi(nxe,nypmx,nblok)
      dimension fxy(2,nxe,nypmx,nblok)
      dimension nyp(nblok)
c local data
      integer j, k, l
      double precision dt1, sum1, sum2
      sum1 = 0.0d0
      sum2 = 0.0d0
      do 30 l = 1, nblok
      do 20 k = 1, nyp(l)
      do 10 j = 1, nx
      dt1 = dble(qi(j,k,l))
      sum1 = sum1 + dt1*fxy(1,j,k,l)
      sum2 = sum2 + dt1*fxy(2,j,k,l)
   10 continue
   20 continue
   30 continue
      pxi = sum1*dt
      pyi = sum2*dt
      return
      end
c-----------------------------------------------------------------------
      subroutine PADDVRFIELD2(a,b,c,ndim,nxe,nypmx,nblok)
c this subroutine calculates a = b + c for distributed real vector field
      implicit none
      integer ndim, nxe, nypmx, nblok
      real a, b, c
      dimension a(ndim,nxe,nypmx,nblok)
      dimension b(ndim,nxe,nypmx,nblok), c(ndim,nxe,nypmx,nblok)
c local data
      integer i, j, k, l
      do 40 l = 1, nblok
      do 30 k = 1, nypmx
      do 20 j = 1, nxe
      do 10 i = 1, ndim
      a(i,j,k,l) = b(i,j,k,l) + c(i,j,k,l)
   10 continue
   20 continue
   30 continue
   40 continue
      return
      end
