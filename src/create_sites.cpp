/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "spktype.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "create_sites.h"
#include "app.h"
#include "app_lattice.h"
#include "domain.h"
#include "lattice.h"
#include "region.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

// same as in lattice.cpp

enum{SC_6N,SC_26N,FCC,BCC,DIAMOND,
     FCC_OCTA_TETRA};

enum{BOX,REGION};
enum{DUMMY,IARRAY,DARRAY};

#define DELTALOCAL 1000000
#define DELTABUF 10000
#define EPSILON 0.0001

/* ---------------------------------------------------------------------- */

CreateSites::CreateSites(SPPARKS *spk) : Pointers(spk) {}

/* ---------------------------------------------------------------------- */

void CreateSites::command(int narg, char **arg)
{
  if (app == NULL) 
    error->all(FLERR,"Create_sites command before app_style set");
  if (domain->box_exist == 0) 
    error->all(FLERR,"Create_sites command before simulation box is defined");
  if (app->sites_exist == 1) 
    error->all(FLERR,"Cannot create sites after sites already exist");
  if (domain->lattice == NULL)
    error->all(FLERR,"Cannot create sites with undefined lattice");

  if (narg < 1) error->all(FLERR,"Illegal create_sites command");

  int iarg;
  if (strcmp(arg[0],"box") == 0) {
    style = BOX;
    iarg = 1;
  } else if (strcmp(arg[0],"region") == 0) {
    style = REGION;
    if (narg < 2) error->all(FLERR,"Illegal create_sites command");
    nregion = domain->find_region(arg[1]);
    if (nregion == -1) 
      error->all(FLERR,"Create_sites region ID does not exist");
    iarg = 2;
  } else error->all(FLERR,"Illegal create_sites command");

  // parse optional args

  valueflag = DUMMY;
  nbasis = domain->lattice->nbasis;
  basisflag = new int[nbasis+1];
  basis_ivalue = new int[nbasis+1];
  basis_dvalue = new double[nbasis+1];
  for (int i = 1; i <= nbasis; i++) basisflag[i] = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"value") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal create_sites command");
      valueflag = 1;
      if (strcmp(arg[iarg+1],"site") == 0) {
        valueflag = IARRAY;
        valueindex = 0;
        if (app->iarray == NULL)
          error->all(FLERR,"Creating a quantity application does not support");
      } 
      else if (arg[iarg+1][0] == 'i') {
        valueflag = IARRAY;
        valueindex = atoi(&arg[iarg+1][1]);
        if (valueindex < 1 || valueindex > app->ninteger)
          error->all(FLERR,"Creating a quantity application does not support");
        valueindex--;
      } 
      else if (arg[iarg+1][0] == 'd') {
        valueflag = DARRAY;
        valueindex = atoi(&arg[iarg+1][1]);
        if (valueindex < 1 || valueindex > app->ndouble)
          error->all(FLERR,"Creating a quantity application does not support");
        valueindex--;
      }
      if (valueflag == IARRAY) ivalue = atoi(arg[iarg+2]);
      else dvalue = atof(arg[iarg+2]);
      iarg += 3;
    } 
    else if (strcmp(arg[iarg],"basis") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal create_sites command");
      if (valueflag == DUMMY) 
	      error->all(FLERR,"Must use value option before basis option "
		    "in create_sites command");
      // int ilo,ihi;
      if (nbasis == 0) 
	      error->all(FLERR,"Cannot use create_sites basis with random lattice");
      // int count = 0;
      error->all(FLERR,"Create sites basis not properly implemented");
      // for (int i = ilo; i <= ihi; i++) {
      //   basisflag[i] = 1;
      //   if (valueflag == IARRAY) basis_ivalue[i] = atoi(arg[iarg+2]);
      //   else if (valueflag == DARRAY) basis_dvalue[i] = atof(arg[iarg+2]);
      //   count++;
      // }
      // if (count == 0) error->all(FLERR,"Illegal create_sites command");
      // iarg += 3;
    } else error->all(FLERR,"Illegal create_sites command");
  }

  // error checks
  // full option only allowed for on-lattice, style = BOX, simple lattices
  // create sites, on-lattice

  if (domain->me == 0) {
    if (screen) fprintf(screen,"Creating sites ...\n");
    if (logfile) fprintf(logfile,"Creating sites ...\n");
  }

  app->sites_exist = 1;

  //int dimension = domain->dimension;
  latstyle = domain->lattice->style;

  applattice = (AppLattice *) app;
  latticeflag = 1;
  xlattice = domain->lattice->xlattice;
  ylattice = domain->lattice->ylattice;
  zlattice = domain->lattice->zlattice;
  structured_lattice();
  structured_connectivity();
  ghosts_from_connectivity(applattice,applattice->delpropensity);
  applattice->print_connectivity();
  
  // clean up

  delete [] basisflag;
  delete [] basis_ivalue;
  delete [] basis_dvalue;
}

/* ----------------------------------------------------------------------
   generate sites on structured lattice that fits in simulation box
   loop over entire lattice
   if style = REGION, require a site be in region as well
   each proc keeps sites in its sub-domain
 ------------------------------------------------------------------------- */

void CreateSites::structured_lattice()
{
  // int dimension = domain->dimension;
  int nonperiodic = domain->nonperiodic;
  int xperiodic = domain->xperiodic;
  int yperiodic = domain->yperiodic;
  int zperiodic = domain->zperiodic;

  double boxxlo = domain->boxxlo;
  double boxylo = domain->boxylo;
  double boxzlo = domain->boxzlo;
  double boxxhi = domain->boxxhi;
  double boxyhi = domain->boxyhi;
  double boxzhi = domain->boxzhi;

  double subxlo = domain->subxlo;
  double subylo = domain->subylo;
  double subzlo = domain->subzlo;
  double subxhi = domain->subxhi;
  double subyhi = domain->subyhi;
  double subzhi = domain->subzhi;

  double **basis = domain->lattice->basis;
  int **iarray = app->iarray;
  double **darray = app->darray;

  // in periodic dims:
  // check that simulation box is integer multiple of lattice spacing
  
  nx = static_cast<int> (domain->xprd / xlattice);
  ny = static_cast<int> (domain->yprd / ylattice);
  nz = static_cast<int> (domain->zprd / zlattice);


  if (xperiodic && fabs(nx*xlattice - domain->xprd) > EPSILON*xlattice)
    error->all(FLERR,"Periodic box is not a multiple of lattice spacing");
  if (yperiodic && fabs(ny*ylattice - domain->yprd) > EPSILON*ylattice)
    error->all(FLERR,"Periodic box is not a multiple of lattice spacing");
  if (zperiodic && fabs(nz*zlattice - domain->zprd) > EPSILON*zlattice)
    error->all(FLERR,"Periodic box is not a multiple of lattice spacing");

  // set domain->nx,ny,nz iff style = BOX and system is fully periodic
  // else site IDs may be non-contiguous and/or ordered irregularly

  if (style == BOX && nonperiodic == 0) {
    domain->nx = nx;
    domain->ny = ny;
    domain->nz = nz;
  }

  // simple = 1 if lattice is square or cubic and fills entire box

  int simple = 0;
  if (style == BOX && (latstyle == SC_6N || latstyle == SC_26N)) simple = 1;

  // if dim is periodic:
  //   lattice origin = lower box boundary
  //   loop bounds = 0 to N-1
  // if dim is non-periodic:
  //   lattice origin = 0.0
  //   loop bounds = enough to tile box completely, with all basis atoms
  //   exact boundary checks will be made later
  
  if (xperiodic) {
    xorig = boxxlo;
    xlo = 0;
    xhi = nx-1;
  } 
  else {
    xorig = 0.0;
    xlo = static_cast<int>(std::floor(boxxlo / xlattice));
    xhi = static_cast<int>(std::ceil(boxxhi / xlattice))-1;
//    xlo = static_cast<int> (boxxlo / xlattice);
//    while ((xlo+1)*xlattice > boxxlo) xlo--;
//    xlo++;
    
//    xhi = static_cast<int> (boxxhi / xlattice);
//    while (xhi*xlattice <= boxxhi) xhi++;
//    xhi--;
  }

  if (yperiodic) {
    yorig = boxylo;
    ylo = 0;
    yhi = ny-1;
  } 
  else {
    yorig = 0.0;
    ylo = static_cast<int>(std::floor(boxylo / ylattice));
    yhi = static_cast<int>(std::ceil(boxyhi / ylattice))-1;
//    ylo = static_cast<int> (boxylo / ylattice);
//    while ((ylo+1)*ylattice > boxylo) ylo--;
//    ylo++;
//    yhi = static_cast<int> (boxyhi / ylattice);
//    while (yhi*ylattice <= boxyhi) yhi++;
//    yhi--;
  }
  
  if (zperiodic) {
    zorig = boxzlo;
    zlo = 0;
    zhi = nz-1;
  } 
  else {
    zorig = 0.0;
    zlo = static_cast<int>(std::floor(boxzlo / zlattice));
    zhi = static_cast<int>(std::ceil(boxzhi / zlattice))-1;
//    zlo = static_cast<int> (boxzlo / zlattice);
//    while ((zlo+1)*zlattice > boxzlo) zlo--;
//    zlo++;
//    zhi = static_cast<int> (boxzhi / zlattice);
//    while (zhi*zlattice <= boxzhi) zhi++;
//    zhi--;
  }

  // check for possible overflow of site IDs

  bigint nglobal = (xhi-xlo+1)*nbasis;
  nglobal *= (yhi-ylo+1)*nbasis;
  nglobal *= (zhi-zlo+1)*nbasis;
  if (nglobal > MAXTAGINT) 
    error->all(FLERR,"Site IDs may exceed max ID value");

  // generate xyz coords and store them with site ID
  // tile the simulation box from origin, respecting PBC
  // for non-periodic dims, check if site is within global box
  // for style = REGION, check if site is within region
  // site IDs should be contiguous (1 to N) in any of these cases:
  //   style = BOX and fully periodic
  //   style = BOX and nonperiodic, simple lattice (sq or sc)
  // site IDs may not be contiguous in any of these cases:
  //   style = REGION
  //   box is nonperiodic, lattice is not simple

  int i,j,k,m,nlocal;
  // tagint n,gid,ii,jj,kk;
  tagint n,gid;
  double x,y,z;

  int maxlocal = 0;
  siteijk = NULL;

  int xlo_me = xhi;
  int xhi_me = xlo;
  int ylo_me = yhi;
  int yhi_me = ylo;
  int zlo_me = zhi;
  int zhi_me = zlo;
  
  

  n = 0;
  for (k = zlo; k <= zhi; k++) {
    for (j = ylo; j <= yhi; j++) {
      for (i = xlo; i <= xhi; i++) {
        for (m = 0; m < nbasis; m++) {
          n++;
          gid = n;

          x = (i + basis[m][0])*xlattice + xorig;
          y = (j + basis[m][1])*ylattice + yorig;
          z = (k + basis[m][2])*zlattice + zorig;

          if (nonperiodic) {
            if (!xperiodic && (x < boxxlo || x >= boxxhi)) continue;
            if (!yperiodic && (y < boxylo || y >= boxyhi)) continue;
            if (!zperiodic && (z < boxzlo || z >= boxzhi)) continue;
          }
          if (style == REGION && domain->regions[nregion]->match(x,y,z) == 0) continue;

          if (x < subxlo || x >= subxhi || 
              y < subylo || y >= subyhi || 
              z < subzlo || z >= subzhi) continue;

          applattice->add_site(gid,x,y,z); //also weird that this can grow by delta
          nlocal = app->nlocal; 

          if (nlocal > maxlocal) { // this is weird, isn't the local array size just (subhi-sublo)*(subhi-sublo)*(subhi-sublo)?
            maxlocal += DELTALOCAL; // why are we increasing it by a fixed amount?
            memory->grow(siteijk,maxlocal,4,"create:siteijk");
          }

          siteijk[nlocal-1][0] = i;
          siteijk[nlocal-1][1] = j;
          siteijk[nlocal-1][2] = k;
          siteijk[nlocal-1][3] = m;

          if (valueflag == IARRAY) {
            if (basisflag[m+1])
              iarray[valueindex][nlocal-1] = basis_ivalue[m+1];
            else iarray[valueindex][nlocal-1] = ivalue;
          } 
          else if (valueflag == DARRAY) {
            if (basisflag[m+1]) 
              darray[valueindex][nlocal-1] = basis_dvalue[m+1];
            else darray[valueindex][nlocal-1] = dvalue;
          }

          if (simple) { // I've got no idea what this is doing. Already have subhi/lo?
            xlo_me = MIN(i,xlo_me);
            xhi_me = MAX(i,xhi_me);
            ylo_me = MIN(j,ylo_me);
            yhi_me = MAX(j,yhi_me);
            zlo_me = MIN(k,zlo_me);
            zhi_me = MAX(k,zhi_me);
          }
        }
      }
    }
  }

  // print site count

  tagint nbig = app->nlocal;
  MPI_Allreduce(&nbig,&app->nglobal,1,MPI_SPK_TAGINT,MPI_SUM,world);

  if (domain->me == 0) {
    if (screen)
      fprintf(screen,"  " TAGINT_FORMAT " sites\n",app->nglobal);
    if (logfile)
      fprintf(logfile,"  " TAGINT_FORMAT " sites\n",app->nglobal);
  }

  // for style = BOX and periodic system, check if nglobal is correct

  if (style == BOX && domain->nonperiodic == 0) {
    nbig = nbasis;
    nbig = nbig*nx*ny*nz;
    if (style == BOX && app->nglobal != nbig)
      error->all(FLERR,"Did not create correct number of sites");
  }

  // for simple latice, do further checking and store extents in AppLattice

  if (simple) {

    // convert SPPARKS loop bounds to Stitch integer indices

    int delta;
    if (xperiodic) {
      delta = static_cast<int> (boxxlo);
      xlo_me += delta;
      xhi_me += delta;
    }
    if (yperiodic) {
      delta = static_cast<int> (boxylo);
      ylo_me += delta;
      yhi_me += delta;
    }
    if (zperiodic) {
      delta = static_cast<int> (boxzlo);
      zlo_me += delta;
      zhi_me += delta;
    }

    // check that product of my extents = nlocal

    int flag = 0;
    int mine = (xhi_me-xlo_me+1) * (yhi_me-ylo_me+1) * (zhi_me-zlo_me+1);
    if (mine != app->nlocal) flag = 1;
    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
    if (flagall) error->all(FLERR,"Local simple lattice is not correct");

    // also check that product of global min/max = nglobal
    
    int xminlo = xhi;
    int xmaxhi = xlo;

    MPI_Allreduce(&xlo_me,&xminlo,1,MPI_INT,MPI_MIN,world);
    MPI_Allreduce(&xhi_me,&xmaxhi,1,MPI_INT,MPI_MAX,world);

    int yminlo = xhi;
    int ymaxhi = xlo;
    MPI_Allreduce(&ylo_me,&yminlo,1,MPI_INT,MPI_MIN,world);
    MPI_Allreduce(&yhi_me,&ymaxhi,1,MPI_INT,MPI_MAX,world);

    int zminlo = xhi;
    int zmaxhi = xlo;
    MPI_Allreduce(&zlo_me,&zminlo,1,MPI_INT,MPI_MIN,world);
    MPI_Allreduce(&zhi_me,&zmaxhi,1,MPI_INT,MPI_MAX,world);

    tagint all = (xmaxhi-xminlo+1) * (ymaxhi-yminlo+1) * (zmaxhi-zminlo+1);
    if (all != app->nglobal) 
      error->all(FLERR,"Global simple lattice is not consistent");

    // store lattice bounds info in AppLattice
    // this allows set and dump stitch commands to access the geometry info

    applattice->simple = simple;
    applattice->xlo_simple = xminlo;
    applattice->xhi_simple = xmaxhi;
    applattice->ylo_simple = yminlo;
    applattice->yhi_simple = ymaxhi;
    applattice->zlo_simple = zminlo;
    applattice->zhi_simple = zmaxhi;
    applattice->xlo_me_simple = xlo_me;
    applattice->xhi_me_simple = xhi_me;
    applattice->ylo_me_simple = ylo_me;
    applattice->yhi_me_simple = yhi_me;
    applattice->zlo_me_simple = zlo_me;
    applattice->zhi_me_simple = zhi_me;
  }
}

/* ----------------------------------------------------------------------
   generate site connectivity for on-lattice applications
   respect non-periodic boundaries
   only called for on-lattice models
 ------------------------------------------------------------------------- */

void CreateSites::structured_connectivity()
{
  int i,j,m,max;
  int ineigh,jneigh,kneigh,mneigh;
  tagint gid;
  double xneigh,yneigh,zneigh;

  // int dimension = domain->dimension;
  int nonperiodic = domain->nonperiodic;
  int thermally_insulated = domain->therminsulated;
  int xperiodic = domain->xperiodic;
  int yperiodic = domain->yperiodic;
  int zperiodic = domain->zperiodic;

  double boxxlo = domain->boxxlo;
  double boxylo = domain->boxylo;
  double boxzlo = domain->boxzlo;
  double boxxhi = domain->boxxhi;
  double boxyhi = domain->boxyhi;
  double boxzhi = domain->boxzhi;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  // set maxneigh and allocate idneigh array to store connectivity

  if (latstyle == SC_6N) maxneigh = 6;
  else if (latstyle == SC_26N) maxneigh = 26;
  else if (latstyle == FCC) maxneigh = 12;
  else if (latstyle == BCC) maxneigh = 8;
  else if (latstyle == DIAMOND) maxneigh = 4;
  else if (latstyle == FCC_OCTA_TETRA) maxneigh = 26;
  else error->all(FLERR,"Illegal lattice style");

  memory->create(idneigh,app->nlocal,maxneigh,"create:idneigh");

  // create connectivity offsets

  int nbasis = domain->lattice->nbasis;
  double **basis = domain->lattice->basis;
  memory->create(cmap,nbasis,maxneigh,4,"create:cmap");
  offsets(basis);

  // generate global lattice connectivity for each site
  // for non-periodic dims, site must be in global box and not across boundary
  // for style = REGION, check if site is in region
  // FCC_OCTA_TETRA is special case, # of neighs not same for all sites

  tagint nglobal = app->nglobal;
  int nlocal = app->nlocal;
  tagint *id = app->id;
  uint8_t *numneigh = applattice->numneigh;

  for (i = 0; i < nlocal; i++) {
    numneigh[i] = 0;
    if (latstyle == FCC_OCTA_TETRA) {
      if ((id[i]-1) % 16 < 8) max = maxneigh;
      else max = 14;
    } 
    else max = maxneigh;
    
    for (j = 0; j < max; j++) {

      // ijkm neigh = indices of neighbor site
      // calculated from siteijk and cmap offsets

      m = siteijk[i][3];
      ineigh = siteijk[i][0] + cmap[m][j][0];
      jneigh = siteijk[i][1] + cmap[m][j][1];
      kneigh = siteijk[i][2] + cmap[m][j][2];
      mneigh = cmap[m][j][3];

      // xyz neigh = coords of neighbor site
      // calculated in same manner that structured_lattice() generated coords

      xneigh = (ineigh + basis[mneigh][0])*xlattice + xorig;
      yneigh = (jneigh + basis[mneigh][1])*ylattice + yorig;
      zneigh = (kneigh + basis[mneigh][2])*zlattice + zorig;

      // remap neighbor coords and indices into periodic box via ijk neigh

      if (xperiodic) {
        if (ineigh < 0) {
          xneigh += xprd;
          ineigh += nx;
        }
        if (ineigh >= nx) {
          xneigh -= xprd;
          xneigh = MAX(xneigh,boxxlo);
          ineigh -= nx;
        }
      }
      if (yperiodic) {
        if (jneigh < 0) {
          yneigh += yprd;
          jneigh += ny;
        }
        if (jneigh >= ny) {
          yneigh -= yprd;
          yneigh = MAX(yneigh,boxylo);
          jneigh -= ny;
        }
      }
      if (zperiodic) {
        if (kneigh < 0) {
          zneigh += zprd;
          kneigh += nz;
        }
        if (kneigh >= nz) {
          zneigh -= zprd;
          zneigh = MAX(zneigh,boxzlo);
          kneigh -= nz;
        }
      }
      //Trying to add thermally insulated option, will also affect microstructure simulation
      if (thermally_insulated) {
          if ((xneigh < boxxlo || xneigh >= boxxhi) || (yneigh < boxylo || yneigh >= boxyhi) || (zneigh < boxzlo || zneigh >= boxzhi)) {
              m = siteijk[i][3];
              ineigh = siteijk[i][0];
              jneigh = siteijk[i][1];
              kneigh = siteijk[i][2];
              mneigh = cmap[m][j][3];
              tagint one = 1;   // use this to avoid int overflow in gid calculation
              gid = one * (kneigh - zlo) * (yhi - ylo + 1) * (xhi - xlo + 1) * nbasis +
                  one * (jneigh - ylo) * (xhi - xlo + 1) * nbasis + one * (ineigh - xlo) * nbasis +
                  mneigh + 1;
              idneigh[i][numneigh[i]++] = gid;
              continue;
          }

      // discard neighs that are outside non-periodic box or region

      if (nonperiodic) {
        if (!xperiodic && (xneigh < boxxlo || xneigh >= boxxhi)) continue;
        if (!yperiodic && (yneigh < boxylo || yneigh >= boxyhi)) continue;
        if (!zperiodic && (zneigh < boxzlo || zneigh >= boxzhi)) continue;
      }
      
      }

      if (style == REGION && domain->regions[nregion]->match(xneigh,yneigh,zneigh) == 0) continue;

      // gid = global ID of neighbor
      // calculated in same manner that structured_lattice() generated IDs

      tagint one = 1;   // use this to avoid int overflow in gid calculation

      gid = one*(kneigh-zlo)*(yhi-ylo+1)*(xhi-xlo+1)*nbasis + 
        one*(jneigh-ylo)*(xhi-xlo+1)*nbasis + one*(ineigh-xlo)*nbasis + 
        mneigh + 1;

      if (style == BOX && nonperiodic == 0 && (gid <= 0 || gid > nglobal))
        error->all(FLERR,"Bad neighbor site ID");

      // add gid to neigh list of site I

      idneigh[i][numneigh[i]++] = gid;
    }
  }

  // delete siteijk and connectivity offsets

  memory->destroy(siteijk);
  memory->destroy(cmap);
}

/* ----------------------------------------------------------------------
   set maxneigh and initialize idneigh when lattice created via read_sites
   called from read_sites when it reads in sites and neighbors
 ------------------------------------------------------------------------- */

void CreateSites::read_sites(AppLattice *apl)
{
  int i,j;

  maxneigh = apl->maxneigh;
  memory->create(idneigh,app->nlocal,maxneigh,"create:idneigh");

  uint8_t *numneigh = apl->numneigh;
  int **neighbor = apl->neighbor;
  int nlocal = app->nlocal;

  for (i = 0; i < nlocal; i++)
    for (j = 0; j < numneigh[i]; j++)
      idneigh[i][j] = neighbor[i][j];
}

/* ----------------------------------------------------------------------
   create ghosts sites around local sub-domain
   only called for on-lattice models
   pass applattice as pointer so can call from ReadSites
   numneigh and global neighbor IDs of each owned site are known as input
   add ghost sites for delpropensity layers
   form neigh list for each layer of ghost sites one layer at a time
   when done, delpropensity-1 layers have a full numneigh and neigh list
     last delpropensity layers has a partial numneigh and neigh list
   convert neighbor IDs from global indices to local indices
 ------------------------------------------------------------------------- */

void CreateSites::ghosts_from_connectivity(AppLattice *apl, int delpropensity)
{
  int i,j,k,m,proc,owner_ghost,index_ghost;
  tagint idglobal,idghost,idrecv;
  double x,y,z;
  tagint *id;
  uint8_t *numneigh;
  int **neighbor;
  double **xyz;

  MyHash hash;
  MyIterator loc;

  int me = domain->me;
  int nprocs = domain->nprocs;
  int nlocal = app->nlocal;

  // nchunk = size of one site datum circulated in message

  int nchunk = 7 + maxneigh;

  // setup ring of procs

  int next = me + 1;
  int prev = me -1; 
  if (next == nprocs) next = 0;
  if (prev < 0) prev = nprocs - 1;

  // loop over delpropensity layers to build up layers of ghosts

  int npreviousghost;
  int nghost = 0;

  for (int ilayer = 0; ilayer < delpropensity; ilayer++) {

    // put all sites (owned + current ghosts) in hash
    // key = global ID, value = local index

    id = app->id;

    hash.clear();
    for (i = 0; i < nlocal+nghost; i++)
      hash.insert(std::pair<tagint,int> (id[i],i));

    // make a list of sites I need
    // loop over neighbors of owned + current ghost sites
    // check if site is already an owned or ghost site or already in list
    // if not, add it to new site list and to hash

    double *buf = NULL;
    int nbuf = 0;
    int maxbuf = 0;
    int nsite = 0;

    numneigh = apl->numneigh;

    for (i = 0; i < nlocal+nghost; i++) {
      for (j = 0; j < numneigh[i]; j++) {
	idglobal = idneigh[i][j];
	if (hash.find(idglobal) == hash.end()) {
	  if (nbuf + nchunk >= maxbuf) {
	    maxbuf += DELTABUF;
	    memory->grow(buf,maxbuf,"create:buf");
	  }
	  buf[nbuf] = idglobal;
	  buf[nbuf+1] = -1;
	  nbuf += nchunk;
	  hash.insert(std::pair<tagint,int> (idglobal,nlocal+nghost+nsite));
	  nsite++;
	}
      }
    }

    // maxsize = max buf size on any proc

    int maxsize;
    MPI_Allreduce(&nbuf,&maxsize,1,MPI_INT,MPI_MAX,world);

    memory->grow(buf,maxsize,"create:buf");
    double *bufcopy;
    memory->create(bufcopy,maxsize,"create:bufcopy");

    // cycle site list around ring of procs back to self
    // when receive it, fill in info for any sites I own
    // info = proc, local index, xyz, numneigh, list of global neighbor IDs
    
    MPI_Request request;
    MPI_Status status;

    xyz = app->xyz;

    int size = nbuf;

    for (int loop = 0; loop < nprocs; loop++) {
      if (me != next) {
        MPI_Irecv(bufcopy,maxsize,MPI_DOUBLE,prev,0,world,&request);
        MPI_Send(buf,size,MPI_DOUBLE,next,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&size);
        nsite = size / nchunk;
        memcpy(buf,bufcopy,size*sizeof(double));
      }
      for (int i = 0; i < nsite; i++) {
        m = i * nchunk;
        idrecv = static_cast<tagint> (buf[m++]);
        proc = static_cast<int> (buf[m++]);
        if (proc >= 0) continue;
        loc = hash.find(idrecv);
        if (loc == hash.end() || loc->second >= nlocal) continue;

        j = loc->second;
        buf[m-1] = me;
        buf[m++] = j;
        buf[m++] = xyz[j][0];
        buf[m++] = xyz[j][1];
        buf[m++] = xyz[j][2];
        buf[m++] = numneigh[j];
        for (k = 0; k < numneigh[j]; k++)
          buf[m++] = idneigh[j][k];
      }
    }

    // original site list came back to me around ring
    // realloc idneigh to store neighbor info for these ghost sites
    // extract info for my new layer of ghost sites
    // reset numneigh after each call to add_ghost() in case realloc occurred
    // error if any site is not filled in

    npreviousghost = nghost;
    nghost += nsite;
    memory->grow(idneigh,nlocal+nghost,maxneigh,"create:idneigh");

    for (i = 0; i < nsite; i++) {
      m = i * nchunk;
      idghost = static_cast<tagint> (buf[m++]);
      owner_ghost = static_cast<int> (buf[m++]);
      if (owner_ghost < 0) error->one(FLERR,"Ghost site was not found");
      index_ghost = static_cast<int> (buf[m++]);
      x = buf[m++];
      y = buf[m++];
      z = buf[m++];

      apl->add_ghost(idghost,x,y,z,owner_ghost,index_ghost);
      numneigh = apl->numneigh;

      j = nlocal + npreviousghost + i;
      numneigh[j] = static_cast<int> (buf[m++]);
      for (k = 0; k < numneigh[j]; k++)
	      idneigh[j][k] = static_cast<tagint> (buf[m++]);
    }

    // clean up

    memory->destroy(buf);
    memory->destroy(bufcopy);
  }

  // can now set AppLattice::maxneigh and allocate AppLattice::neighbor

  apl->maxneigh = maxneigh;
  apl->grow(apl->nmax);
  numneigh = apl->numneigh;
  neighbor = apl->neighbor;

  // convert all global neighbors to local indices in AppLattice::neighbor
  // if i is owned or in delpropensity-1 layers, then error if neigh not found
  // if i is ghost in last delpropensity layer, then delete neigh if not found

  for (i = 0; i < nlocal+nghost; i++) {
    j = 0;
    while (j < numneigh[i]) {
      idglobal = idneigh[i][j];
      loc = hash.find(idglobal);
      if (loc != hash.end()) {
        neighbor[i][j] = loc->second;
        j++;
      } 
      else if (i >= nlocal+npreviousghost) {
        numneigh[i]--;
        for (k = j; k < numneigh[i]; k++) idneigh[i][k] = idneigh[i][k+1];
      } 
      else error->one(FLERR,"Ghost connection was not found");
    }
  }

  // no longer need idneigh since AppLattice::neighbor now exists

  memory->destroy(idneigh);
}

/* ---------------------------------------------------------------------- */

void CreateSites::offsets(double **basis)
{
  if (latstyle == SC_6N){
    for (int m = 0; m < nbasis; m++)
      offsets_3d(m,basis,xlattice,xlattice,maxneigh,cmap[m]);
  }
  else if (latstyle == SC_26N){
    fprintf(screen,"sc_26");
    for (int m = 0; m < nbasis; m++)
      offsets_3d(m,basis,xlattice,sqrt(3.0)*xlattice,maxneigh,cmap[m]);
  }
  else if (latstyle == FCC) {
    fprintf(screen,"FCC");
    for (int m = 0; m < nbasis; m++)
      offsets_3d(m,basis,sqrt(2.0)/2.0*xlattice,sqrt(2.0)/2.0*xlattice,
		 maxneigh,cmap[m]);
  }
  else if (latstyle == BCC) {
    for (int m = 0; m < nbasis; m++)
      offsets_3d(m,basis,sqrt(3.0)/2.0*xlattice,sqrt(3.0)/2.0*xlattice,
		 maxneigh,cmap[m]);
  }
  else if (latstyle == DIAMOND) {
    for (int m = 0; m < nbasis; m++)
      offsets_3d(m,basis,sqrt(3.0)/4.0*xlattice,sqrt(3.0)/4.0*xlattice,
		 maxneigh,cmap[m]);
  }
  else if (latstyle == FCC_OCTA_TETRA) {
    fprintf(screen,"FCC_OCTA_TETRA");
    for (int m = 0; m < 4; m++) {
      offsets_3d(m,basis,sqrt(2.0)/2.0*xlattice,sqrt(2.0)/2.0*xlattice,
		 12,&cmap[m][0]);
      offsets_3d(m,basis,0.5*xlattice,0.5*xlattice,6,&cmap[m][12]);
      offsets_3d(m,basis,sqrt(3.0)/4.0*xlattice,sqrt(3.0)/4.0*xlattice,
		 8,&cmap[m][18]);
    }
    for (int m = 4; m < 8; m++) {
      offsets_3d(m,basis,0.5*xlattice,0.5*xlattice,6,&cmap[m][0]);
      offsets_3d(m,basis,sqrt(2.0)/2.0*xlattice,sqrt(2.0)/2.0*xlattice,
		 12,&cmap[m][6]);
      offsets_3d(m,basis,sqrt(3.0)/4.0*xlattice,sqrt(3.0)/4.0*xlattice,
		 8,&cmap[m][18]);
    }
    for (int m = 8; m < nbasis; m++) {
      offsets_3d(m,basis,sqrt(3.0)/4.0*xlattice,sqrt(3.0)/4.0*xlattice,
		 8,&cmap[m][0]);
      offsets_3d(m,basis,0.5*xlattice,0.5*xlattice,6,&cmap[m][8]);
    }
  }
}

/* ---------------------------------------------------------------------- */

void CreateSites::offsets_3d(int ibasis, double **basis, 
                             double cutlo, double cuthi, 
                             int ntarget, int **cmapone)
{
  int i,j,k,m,n;
  double x0,y0,z0,delx,dely,delz,r;

  n = 0;
  x0 = basis[ibasis][0] * xlattice;
  y0 = basis[ibasis][1] * ylattice;
  z0 = basis[ibasis][2] * zlattice;
//  fprintf(screen, "x0=%f y0=%f z0=%f\n",x0,y0,z0);
//  fprintf(screen, "dx=%f dy=%f dz=%f\n",xlattice,ylattice,zlattice);
//  fprintf(screen, "cutlo=%f cuthi=%f\n",cutlo,cuthi);
  for (i = -1; i <= 1; i++) {
    for (j = -1; j <= 1; j++) {
      for (k = -1; k <= 1; k++) {
        for (m = 0; m < nbasis; m++) {
          delx = (i+basis[m][0])*xlattice - x0;
          dely = (j+basis[m][1])*ylattice - y0;
          delz = (k+basis[m][2])*zlattice - z0;
          r = sqrt(delx*delx + dely*dely + delz*delz);
          //fprintf(screen, "r=%f\n",r);
          if (r >= cutlo-EPSILON && r < cuthi+EPSILON) {
            if (n == ntarget) {
              fprintf(screen, "n=%i,ntarget=%i\n",n,ntarget);
              error->all(FLERR,"Incorrect lattice neighbor count: Too many found");
            }
            cmapone[n][0] = i;
            cmapone[n][1] = j;
            cmapone[n][2] = k;
            cmapone[n][3] = m;
            n++;
          }
        }
      }
    }
  }

  if (n != ntarget) {
    fprintf(screen, "n=%i,ntarget=%i\n",n,ntarget);
    error->all(FLERR,"Incorrect lattice neighbor count: Not enough found");
  }
}
