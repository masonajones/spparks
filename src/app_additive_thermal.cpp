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

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <random>
#include <cmath>
#include "string.h"
#include "math.h"
#include "app_additive_thermal.h"
#include "random_fast.h"
#include "error.h"
#include "memory.h"
#include "comm_lattice.h"
#include "solve.h"
#include "lattice.h"
#include "timer.h"
#include "math_const.h"
#include "domain.h"
#include "output.h"
#include <cassert>
using namespace SPPARKS_NS;
using namespace MathConst;
enum{RANDOM};
enum material_state {
  INACTIVE = 0,
  POWDER = 1,
  MOLTEN = 2,
  SOLID = 3
};

// Not used yet
// Intended to make it easier to initiate comms after Comm -> selected is implemented
enum array_names {
  SPIN = 0,
  ACTIVEFLAG = 1,
  MOBILITY = 0,

  // Not sure how to handle these for comms, since we swap the pointers
  TEMPERATURE = 1,
  TEMPERATURE_OLD = 2,

  SOLIDD = 3
};


/* ---------------------------------------------------------------------- */

AppAdditiveThermal::AppAdditiveThermal(SPPARKS *spk, int narg, char **arg) :
  AppPotts(spk,narg,arg)
{
    // only error check for this class, not derived classes
    double x_devIn, y_devIn;

    if (narg != 13 )
        error->all(FLERR,"Illegal app_style command");

    nspins = atoi(arg[1]); //Number of spins/grain IDs
//    velIn = atof(arg[2]); //Velocity of travel (in meters/second)
    time_step = atof(arg[3]); //The size of each FD_step
    path_file_name = arg[4]; //Name of the input file with the proposed scan path  
    x_devIn = atof(arg[5]); //Set the standard deviation for the gaussian source. This will be in meters
    y_devIn = atof(arg[6]); //Std. dev. in y-direction
//    z_devIn = atof(arg[7]); //For a volumetric Gaussian source
//    short_wait_time = atof(arg[8]); //How long to pause during a layer, e.g. while laser is off and moving between rasters
    flux_prefactor = atof(arg[9]); //Total absorbed laser power in Watts
    substrate_height = atoi(arg[10]); //Substrate height in lattice sites
    nsmooth = atoi(arg[11]); //Number of Potts smoothing steps 
    dx = atof(arg[12]); //Lattice spacing (in meters)
 
    time_step_next = time_step;
    time_step_newlayer = time_step;

    //Read in the entire laser path
    path_file();
    done_flag =0;
    new_layer=1;

    path_index = 1;
    ndouble = 4;
    allow_app_update = 1;
    ninteger = 2;
    nlocal_app = 0;
    allow_kmc = 0;
    d_residual = 0;

    x_dev = x_devIn/dx;
    y_dev = y_devIn/dx;

    prefactor_xydev = MY_2PI * x_dev * y_dev * dx * dx * dx;
    recoating_time = 10;
    Kmc = 0.27695; //SPPARKS KMC scaling factor for sq_26 lattice, should not change for different materials

    //Set default values for 304L stainless steel. Can be modified in input file.
    k_solid = 30;
    boundary_temp = 300;
    density = 5706;

    Ko = 0.0000204133;
    Q = 128312;
    Tl = 1723;
    Ts = 1673;
    No = 1e15;
    Tc = 5;
    Tsig = 3;
    sizeNorm = pow(dx,3) * 2;
    sizeSig = pow(dx,3);
    // specific_heat_length = 6;
    solid_front_length = 4;

    //Let's put in default values for array-based parameters
    specific_heat_temps = new double[1];
    specific_heat_vals = new double[1];
    solid_front_coeffs = new double[4];
    
    specific_heat_temps[0] = 300; // not used
    // specific_heat_temps[1] = 600;
    // specific_heat_temps[2] = 1400;
    // specific_heat_temps[3] = 1550;
    // specific_heat_temps[4] = 1650;
    // specific_heat_temps[5] = 1750;
    
    specific_heat_vals[0] = 470.57;
    // specific_heat_vals[1] = 550;
    // specific_heat_vals[2] = 675;
    // specific_heat_vals[3] = 693;
    // specific_heat_vals[4] = 714;
    // specific_heat_vals[5] = 736;
    
    solid_front_coeffs[0] = 1.091e-5;
    solid_front_coeffs[1] = -2.034e-4;
    solid_front_coeffs[2] = 2.74e-3;
    solid_front_coeffs[3] = 1.151e-4;

// DISABLED PARAMETERS - NEED TO BE RE-ENABLED IN CODE AND HERE IF DESIRED
    // T_room = 300;
    //vel = velIn/dx;
    //z_dev = z_devIn/dx;
    //k_powder = 0.3;
    //h = 25;
    //eta = 1;
    //short_wait_time = 2e-7;
    //latent_heat = 285000;

    recreate_arrays();
  
}

/* ----------------------------------------------------------------------
   input script commands unique to this app
------------------------------------------------------------------------- */

void AppAdditiveThermal::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"k_solid") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal solid conductivity command");
    k_solid = atof(arg[0]);
    if (k_solid <= 0) error->all(FLERR,"Illegal solid conductivity");
  } 
  // else if (strcmp(command,"k_powder") == 0) {
  //   if (narg != 1) error->all(FLERR,"Illegal powder conductivity command");
  //   //k_powder = atof(arg[0]);
  //   //if (k_powder <= 0) 
  //   //  error->all(FLERR,"Illegal powder conductivity");
  // }
  else if (strcmp(command,"boundary_temperature") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal boundary temperature command");
    boundary_temp = atof(arg[0]);
    if (boundary_temp <= 0) error->all(FLERR,"Illegal boundary temperature");
  }
  else if (strcmp(command,"density") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal density command");
    density = atof(arg[0]);
    if (density <= 0) error->all(FLERR,"Illegal density");
  }
  // else if (strcmp(command,"t_room") == 0) {
  //   if (narg != 1) error->all(FLERR,"Illegal room temperature command");
  //   T_room = atof(arg[0]);
  //   if (T_room <= 0) error->all(FLERR,"Illegal room temperature");
  // }
  // else if (strcmp(command,"emissivity") == 0) {
  //   if (narg != 1) error->all(FLERR,"Illegal emissivity command");
  //   eta = atof(arg[0]);
  //   if (eta <= 0) error->all(FLERR,"Illegal emissivity");
  // }
  // else if (strcmp(command,"convection_coefficient") == 0) {
  //   if (narg != 1) error->all(FLERR,"Illegal convection coefficient command");
  //   h = atof(arg[0]);
  //   if (h <= 0) error->all(FLERR,"Illegal convection coefficient");
  // }
  // else if (strcmp(command,"Ko") == 0) {
  //   if (narg != 1) error->all(FLERR,"Illegal Arrhenius pre-factor command");
  //   Ko = atof(arg[0]);
  //   if (eta <= 0) error->all(FLERR,"Illegal Arrhenius pre-factor");
  // }
  else if (strcmp(command,"Q") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal Arrhenius exponential-factor command");
    Q = atof(arg[0]);
    if (Q <= 0) error->all(FLERR,"Illegal Q");
  }
  else if (strcmp(command,"liquidus") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal liquidus command");
    Tl = atof(arg[0]);
    if (Tl <= 0) error->all(FLERR,"Illegal liquidus temperature");
  }
  else if (strcmp(command,"solidus") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal solidus command");
    Ts = atof(arg[0]);
    if (Ts <= 0) error->all(FLERR,"Illegal solidus temperature");
  }
  else if (strcmp(command,"nucleation_density") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal nucleation density command");
    No = atof(arg[0]);
    if (No < 0) error->all(FLERR,"Illegal nucleation density");
  }
  else if (strcmp(command,"critical_undercooling") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal critical_undercooling command");
    Tc = atof(arg[0]);
    if (Tc < 0) error->all(FLERR,"Illegal critical undercooling");
  }
  else if (strcmp(command,"undercooling_deviation") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal undercooling_deviation command");
    Tsig = atof(arg[0]);
    if (Tsig < 0) error->all(FLERR,"Illegal undercooling standard deviation");
  }
  else if (strcmp(command,"mean_nuclei_volume") == 0) {
    if (narg !=1) error->all(FLERR,"Illegal mean_nuclei_volume command");
    sizeNorm = atof(arg[0]);
    if (sizeNorm < 0) error->all(FLERR,"Illegal mean nuclei volume");
  }
  else if (strcmp(command,"nuclei_volume_deviation") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal nuclei_volume_deviation command");
    sizeSig = atof(arg[0]);
    if (sizeSig < 0) error->all(FLERR,"Illegal nuclei volume standard deviation");
  }
  else if (strcmp(command,"specific_heat") == 0) {
    delete [] specific_heat_temps;
    delete [] specific_heat_vals;
    if (narg < 1) error->all(FLERR,"Illegal specific_heat command");
    int nHeat;
    nHeat = atoi(arg[0]);
    if (nHeat <= 0) error->all(FLERR,"Illegal specific heat specification");
    specific_heat_temps = new double [nHeat];
    specific_heat_vals = new double [nHeat];
    int j = 0;
    for(int i = 1; i < nHeat * 2 + 1; i = i + 2) {
        specific_heat_temps[j] = atof(arg[i]);
        specific_heat_vals[j] = atof(arg[i + 1]);
        j++;
    }
  }
  // else if (strcmp(command,"latent_heat") == 0) {
  //   if (narg != 1) error->all(FLERR,"Illegal latent_heat command");
  //   latent_heat = atof(arg[0]);
  //   if (latent_heat < 0) 
  //     error->all(FLERR,"Illegal latent heat value");
  // }
	//Need to modify to get solid_front_length
  else if (strcmp(command,"solid_front_vel") == 0) {
    delete [] solid_front_coeffs;
    if (narg < 1) error->all(FLERR,"Illegal solid_front_vel command");
    solid_front_length = atoi(arg[0]);
    if (solid_front_length <= 0) error->all(FLERR,"Illegal solidification front velocity specification");
    solid_front_coeffs = new double [solid_front_length];
    int j = 0;
    for(int i = 1; i < solid_front_length + 1; i++) {
      solid_front_coeffs[j] = atof(arg[i]);
      j++;
    }
  }
  // else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppAdditiveThermal::grow_app()
{
  spin = iarray[0];
  activeFlag = iarray[1]; // try swapping to int8
  MobilityOut = darray[0];
  T = darray[1];
  T_old = darray[2];
  SolidD = darray[3];
//  energy = darray[4];
  
  if (nlocal_app < nlocal) {
    nlocal_app = nlocal;          
  }
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */
void AppAdditiveThermal::init_app()
{
  delete [] sites;
  delete [] unique;
  const double sqrt2 = 1.4142135623731; //std::sqrt(2);//
  const double sqrt3 = 1.7320508075689; //std::sqrt(3);//
  sites = new int[1 + maxneigh];
  unique = new int[1 + maxneigh];
  // Instead of nucleation flags we check the spin # against the proportion of nspins. Cuts down on cache misses.
  // NOTE: this makes visualization weird
  //nucleationFlags = new bool[nspins];
  nucleationTemps = new double[nspins]; // We can make these smaller to match the number of spins that are allowed to nucleate
	nucleationSizes = new double[nspins];
  
  dt_sweep = 1.0/maxneigh;

  // not really sure of the need for this flag/check
  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
      if (spin[i] < 1 || spin[i] > nspins) {
          flag = 1;
      }
      //If we're above the substrate height, randomize the spins
      if(xyz[i][2] > substrate_height) {
          spin[i] = (int) (nspins*ranapp->uniform());
          activeFlag[i] = INACTIVE;
      }
      //If we're less than the value, set activeFlag to the "solid" condition
      else {
          activeFlag[i] = SOLID;
      }
  }
  comm->all();
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
  
	//Initialize the nucleationFlags vector
	if (domain->me==0) {
		nucleation_spins(ranapp);    
	}

	//MPI_Bcast(nucleationFlags,nspins, MPI_C_BOOL,0,world);
  MPI_Bcast(&nucleationCutoff,1, MPI_INT,0,world);
  
//  nucleationTemps = new double[nucleationCutoff];
//	nucleationSizes = new double[nucleationCutoff];

	//Initialize the nucleationTemps and nucleationSizes vectors
	if (domain->me==0) {
			nucleation_init();    
	}

	MPI_Bcast(nucleationTemps,nspins, MPI_DOUBLE,0,world);
	MPI_Bcast(nucleationSizes,nspins, MPI_DOUBLE,0,world);
	
	//Initialize the neighDist array need to fill with good values
	neighDist = new double[26];
	neighDist[0] = sqrt3 * dx;
	neighDist[1] = sqrt2 * dx;
	neighDist[2] = sqrt3 * dx;
	neighDist[3] = sqrt2 * dx;
	neighDist[4] = dx;
	neighDist[5] = sqrt2 * dx;
	neighDist[6] = sqrt3 * dx;
	neighDist[7] = sqrt2 * dx;
	neighDist[8] = sqrt3 * dx;
	neighDist[9] =  sqrt2 * dx;
	neighDist[10] = dx;
	neighDist[11] =  sqrt2 * dx;
	neighDist[12] = sqrt3 * dx;
	neighDist[13] = sqrt3 * dx;
	neighDist[14] =  sqrt2 * dx;
	neighDist[15] = dx;
	neighDist[16] =  sqrt2 * dx;
	neighDist[17] = sqrt3 * dx;
	neighDist[18] = sqrt2 * dx;
	neighDist[19] = sqrt3 * dx;
	neighDist[20] =  sqrt2 * dx;
	neighDist[21] = dx;
	neighDist[22] =  sqrt2 * dx;
	neighDist[23] = sqrt3 * dx;
	neighDist[24] = sqrt2 * dx;
	neighDist[25] = sqrt3 * dx;
 
  // Playing around with methods to encode neighDist data in a more memory friendly way:
  //struct { // 50 bytes -> fits on cache line
  //  const double neighdists [] = {dx, sqrt2*dx, sqrt3*dx};
  //  const uint8_t neighdist_ind [] = {2,1,2,1,0,1,2,1,2,1,0,1,2,2,1,0,1,2,1,2,1,0,1,2,1,2};
  //} neighDist
  
  // Binary encoding of neighdist_ind:
  // Note: not very conducive to real uses
  // 8 bytes (long) + 24 bytes of neighdists = 32 bytes -> half a cache line (not much more useful unless we have something to pack it with)
  // 21210121210122101212101212
  // 0b1001100100011001100100011010010001100110010001100110
	
	//Check that our timestep is short enough to capture solidification behavior
	double max_front_vel = 0;
	int power = solid_front_length -1;
	for(int k = 0; k < solid_front_length; k++) {
		max_front_vel = max_front_vel + solid_front_coeffs[k] * pow(Tl - Ts, power);
		power--;
	}
	if(max_front_vel * time_step > dx) {
		double max_time_step = dx/max_front_vel;
		fprintf(screen,"time_step is too large to capture moving solidification front. Max allowable timestep is %f\n",max_time_step);
		error->all(FLERR,"Illegal time_step");
	}
          
	this->app_update();
}


/* ----------------------------------------------------------------------
   Read in the coordinate file before going into MPI
------------------------------------------------------------------------- */
void AppAdditiveThermal::path_file()
{
   int x_max = 0;
   int y_max = 0;
   int z_max = 0;
   
   std::string line;
   std::ifstream myfile(path_file_name);
   
   //Check if the file exists and is readable.
   if (myfile.good()) {
   	if(domain->me == 0) {
   		fprintf(screen,"Path input file found \n");
   	}
   }
	else {
		error->all(FLERR,"Path input file not found");;
	}

  // new lines will be skipped unless we stop it from happening:    
  myfile.unsetf(std::ios_base::skipws);

  int aNumOfLines = 0;

  std::string aLineStr;
  while (getline(myfile, aLineStr))
  {
      if (!aLineStr.empty())
          line_count++;
  }
	
	//Lets just use a 1D vector with fancy indexing to do our scan_array
	//Each point will be accessed by scan_array[y * sizeX + x], where x is 0-3 and y is 0 - line_count - 1
	//We should also make new arrays for each xyz and value, cause I want to easily access the variables.
	//scan_array = new double[line_count * 5];
  scan_array = new double[5];
	x_scan_array = new double[line_count];
	y_scan_array = new double[line_count];
	z_scan_array = new double[line_count];
	t_scan_array = new double[line_count];
	p_scan_array = new double[line_count];
	
  //Read in the input file. It should always have 5 columns X,Y,Z distance, and pause_flag
  std::ifstream in_file(path_file_name);
  
  for (int row = 0; row < line_count; row++) {
  	getline(in_file, line);
  	if(!in_file.good() )
  		break;
  	
  	std::stringstream iss(line);
  	
  	for (int col = 0; col < 5; col++) {
  		std::string val;
  		getline(iss, val, ',');
  		std::stringstream convertor(val);

  		// convertor >> scan_array[row * 5 + col];  	
      convertor >> scan_array[col];  		
  	}
// Ensure Z height doesn't leave domain, trim path if it does. INTRODUCES BUG
//    if(scan_array[2]>domain->boxzhi) {
//      fprintf(screen,"Laser z height %f higher than domain boundary %f. Truncating scan path.", scan_array[2],domain->boxzhi);
////      x_scan_array[row] = 0.0;
////  	  y_scan_array[row] = 0.0;
////      z_scan_array[row] = 0.0;
////      t_scan_array[row] = scan_array[3];
////      p_scan_array[row] = -1.0;
//      break;
    //}
    x_scan_array[row] = scan_array[0];
  	y_scan_array[row] = scan_array[1];
    z_scan_array[row] = scan_array[2];
    t_scan_array[row] = scan_array[3];
    p_scan_array[row] = scan_array[4];
  }
  delete [] scan_array;
  in_file.close();

// Could scale Power here to cut down on repetitive calculations in main loop.
//  for (int i = 0; i<line_count; i++) {
//    p_scan_array[i] = p_scan_array[i]/prefactor_xydev;
//  }
}



/* ----------------------------------------------------------------------
 perform finite difference on a single site
 ------------------------------------------------------------------------- */

void AppAdditiveThermal::site_event_finitedifference(int i)
{
  static double mult1 = k_solid / (dx * dx);
  static double mult2 = 1/(specific_heat_vals[0] * density);

	double chm_temp;
  if(activeFlag[neighbor[i][13]] != 0) {
    chm_temp = T_old[neighbor[i][13]];
    chm_temp += T_old[neighbor[i][4]];
    chm_temp += T_old[neighbor[i][10]];
    chm_temp += T_old[neighbor[i][12]];
    chm_temp += T_old[neighbor[i][15]];
    chm_temp += T_old[neighbor[i][21]];
    chm_temp += -6 * T_old[i];
    chm_temp *= mult1;
  }
  else {
    chm_temp = T_old[neighbor[i][4]];
    chm_temp += T_old[neighbor[i][10]];
    chm_temp += T_old[neighbor[i][12]];
    chm_temp += T_old[neighbor[i][15]];
    chm_temp += T_old[neighbor[i][21]];
    chm_temp += -5 * T_old[i];
    chm_temp *= mult1;
    double flux_loc = flux_finder(i);
    chm_temp += flux_loc;
  }

  T[i] = T_old[i] + time_step * chm_temp * mult2; // could save 1 more operation if we made mult2 = timestep/(Cp*density) in app update before loop, but probably not worth it
	return;
	
}



/* ----------------------------------------------------------------------
 perform finite difference on the local domain, instead of site by site. WORK IN PROGRESS
 ------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
void AppAdditiveThermal::local_finitedifference()
{
    
  //set values used for convenience
	double Cp = specific_heat_vals[0];
  double flux_loc;
  //double chm_temp = 0; // actually needs to be a simd vector
  // Skip layer 0, nothing happens there
  // First do layer 1
  for(int j=1; j<ny-1; j++) {
      ind = j*nx + 2*nx*ny // something like this
      //_mm_preload(ind) (T and T_old)
      //_mm_preload(ind+nx) (T and T_old)
      //_mm_preload(ind+nx*ny) (T and T_old)
      // below is const
      // assuming j-- already been added
      T[ind] += T_old[ind-1]
      
      for(int i=1; i<nx-1; i+=4) {
        ind = j*nx + 2*nx*ny
        
        //_mm_preload(ind+4) (T and T_old)
        //_mm_preload(ind+nx+4) (T and T_old)
        //_mm_preload(ind+nx*ny+4) (T and T_old)
                
        __mm_mov T[ind]
     	  __mm_mov T_old[ind]
        __mm_add T[ind] T_old[ind] *-6
        
        __mm_mov T_old[ind+1]
        __mm_add T[ind] T_old[ind+1]
        __mm_add T[ind+1] T_old[ind] 
        
        __mm_mov T_old[ind+nx]
        __mm_add T[ind] T_old[ind+nx]
        __mm_add T[ind+nx] T_old[ind]
        
        __mm_mov T_old[ind+nx*ny]
        __mm_add T[ind] T_old[ind+nx*ny]
        __mm_add T[ind+nx*ny] T_old[ind]
        
        __mm_add T[ind] boundary_temp
        
      	//multiply by conduction constants
        __mm_mult T[ind], time_step * k_solid / (dx * dx * (Cp * density)) ;
        __mm_add T[ind] T_old[ind]
      }
    }
    
  for(int k=2; k<nz; k++) {
    for(int j=1; j<ny-1; j++) {
      for(int i=1; i<nx-1; i+=4) {
        flux_loc = 0;
        chm_temp = 0;
       
        if (xyz[i][2] == domain->boxzlo) {
           T[i] = boundary_temp;
           //T_old[i] = boundary_temp;
           return;
        }
      	
      	//First thing, loop through the neighbors of a site and see if any of them are inactive
      	//We need to deal with when a site has more than one inactive value, which isn't that uncommon...
      	//Lets average the contributions from all the orientations
      	for (int j = 0; j < 6; j++) {
      		
      		int s = neighbor[i][good_neigh[j]];
      		
      		//Loop through neighbors, do usual FD if active, and do BCs if not
      		//Treat solid & liquid cases in the same way
      		//Site is Solid (or liquid)
      		
      		//Conduction
      		if(activeFlag[s] != 0) {
          		chm_temp += (T_old[s] - T_old[i]);
      		}
      	  //Do laser flux
      	  //else if(!flux_flag) {
          else {
              //Assuming laser flux only at surface
              //Figure out the flux value for the current heat source.
              flux_loc = flux_finder(i);
              //flux_flag = true;
      	  }
      		
      	}
      	//multiply by conduction constants
        chm_temp *= k_solid / (dx * dx);
      	//Add the laser power/volume
        chm_temp += flux_loc;
       
      	T[i] = T_old[i] + time_step * chm_temp / (Cp * density);
        //T[i] += time_step * chm_temp / (Cp * density);
      }
    }
  }
	return;
	
}

 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
Determine Laser loacation and parameters from path file data based on current time index.
 ------------------------------------------------------------------------- */
void AppAdditiveThermal::position_finder_in() {
  if (p_scan_array[path_index] > 0.0) {
		wait_time = 0.0;
		x_meltspot = x_scan_array[path_index];
		y_meltspot = y_scan_array[path_index];
		z_meltspot = z_scan_array[path_index];
    flux_prefactor = p_scan_array[path_index];
    time_step = time_step_next;
    time_step_next = (t_scan_array[path_index+1] - t_scan_array[path_index]);
    path_index++;
		return;
	}
  else if(p_scan_array[path_index] == 0.0) {
    wait_time = 0.0;
    x_meltspot = x_scan_array[path_index];
		y_meltspot = y_scan_array[path_index];
		z_meltspot = z_scan_array[path_index];
    flux_prefactor = 0.0;
    time_step = time_step_next;
    time_step_next = (t_scan_array[path_index+1] - t_scan_array[path_index]);
    path_index++;
    return;
  }
  else if (p_scan_array[path_index] < 0) {
    if (wait_time > 0) {
      wait_time = (wait_time - time_step);
      if (wait_time <= 0.0) { path_index++; }
      return;
    }
    else if (wait_time <= 0) {
      wait_time = recoating_time;
      wait_time = (wait_time - time_step);
      x_meltspot = 0;
      y_meltspot = 0;
      z_meltspot = z_scan_array[path_index];
      flux_prefactor = 0.0;
      time_step = time_step_next;
      time_step_next = (t_scan_array[path_index+1] - t_scan_array[path_index]);
      if (wait_time <= 0.0) { path_index++; }
      return;
    }
  }
  else {
    //if (domain->me == 0) fprintf(screen, "C5\n"); // function should only be called from rank 0
    fprintf(screen, "C5\n");
    path_index++;
    return;
  }

}

/* ----------------------------------------------------------------------
 Compute the mobility at the specified lattice site. Returns a double
 between 0 and 1 representing the mobility. This is used as normal to compute
 the mobility during solid-state grain growth.
 The MobilityOut array that the values are assigned to is also re-used to control
 smoothing steps immediately after solidification by using negative values
 ------------------------------------------------------------------------- */
double AppAdditiveThermal::compute_mobility(int i) {    
    double mob;
	  mob = exp(-Q/(R*T[i]));
   	return mob; 
}   

/* ----------------------------------------------------------------------
	Simple function to calculate flux at a given lattice site.
	This should represent the flux going into an element the size of of a site
 ------------------------------------------------------------------------- */
double AppAdditiveThermal::flux_finder(int site) {
  if(flux_prefactor==0.0) return 0.0f;
	
	double flux_loc;
  double preFactorTrue = flux_prefactor / (prefactor_xydev);
  static const double x_denom = 1/(x_dev*x_dev);
  static const double y_denom = 1/(y_dev*y_dev);

  //Normalized gaussian. Total flux will stay constant no matter the values of x,y,z_dev
  //flux_loc = preFactorTrue * exp(-0.5 * ((pow((xyz[site][0] - x_meltspot), 2) / (x_dev * x_dev) + pow((xyz[site][1] - y_meltspot), 2) / (y_dev * y_dev))));
  flux_loc = preFactorTrue * exp(-0.5 * ((pow((xyz[site][0] - x_meltspot), 2) * x_denom + pow((xyz[site][1] - y_meltspot), 2) * y_denom)));
  
  // if we split this into two exponents, we might be able to implement a function cache for better performance
  // rectilinear motion will have periods with identical x or y deltas, angular motion will have repetition between subsequent passes (in cases of identical velocity profiles)
  // may also be able to get away with approximations of exp/pow

	return flux_loc;
}



/* ----------------------------------------------------------------------
 iterate through the temperature solver and update phases as needed
 ------------------------------------------------------------------------- */
void AppAdditiveThermal::app_update()
{
    tempMax = 0.0;
    tempMaxAll = 0.0;

    //communicate all sites to make sure it's up-to-date when it starts
    timer->stamp();
    comm->all();
    timer->stamp(TIME_COMM);

    //Calculate the current position
    if(domain->me == 0) {
        position_finder_in();
        //if(z_meltspot>= domain->boxzhi-1) error->all(FLERR, "Layer height out of bounds."); // I guess this should still be checked somewhere, move to path file load at some point...
    }

    timer->stamp(TIME_APP);
    
    constexpr int num_values = 8;
    MPI_Request requests[num_values];
    
    MPI_Ibcast(&wait_time,1,MPI_DOUBLE,0,world, &requests[0]);
    MPI_Ibcast(&done_flag,1,MPI_INT,0,world, &requests[1]);
    MPI_Ibcast(&time_step,1,MPI_DOUBLE,0,world, &requests[2]);
    MPI_Ibcast(&flux_prefactor,1,MPI_DOUBLE,0,world, &requests[3]);
    
    // Probably doesn't need to be transfered, just broadcast done_flag when relevant:
    MPI_Ibcast(&path_index,1,MPI_INT,0,world, &requests[4]); 

    // Don't need transfered if wait time or flux_prefactor<=0:
    MPI_Ibcast(&x_meltspot,1,MPI_DOUBLE,0,world, &requests[5]);
    MPI_Ibcast(&y_meltspot,1,MPI_DOUBLE,0,world, &requests[6]);
    MPI_Ibcast(&z_meltspot,1,MPI_DOUBLE,0,world, &requests[7]);
    
    MPI_Waitall(num_values, requests, MPI_STATUSES_IGNORE);
    
    timer->stamp(TIME_COMM);
    //If we're recoating, run until things get below Ts and then reset all values to room temp
    if(wait_time > 0.01) {
        //See if we're below melting throughout the domain
        tempMax = compute_tempMax();
        MPI_Allreduce(&tempMax,&tempMaxAll,1,MPI_DOUBLE,MPI_MAX,world);
        //If below, set all values to room temp and skip to the end of the wait
        if(tempMaxAll < 1000) {
            if(domain->me == 0) {
                fprintf(screen,"Continuing after waiting %f\n",recoating_time - wait_time);
            }
            path_index++;
            wait_time = -1;
            new_layer = 1;
            nextoutput = output->compute(time,0);
            //MPI_Bcast(&wait_time,1,MPI_DOUBLE,0,world); // probably cheaper to just update these locally
            //MPI_Bcast(&path_index,1,MPI_INT,0,world);
            //MPI_Bcast(&new_layer,1,MPI_INT,0,world);


            __m256d T0_vec = _mm256_set1_pd(boundary_temp);
            //fprintf(screen, "alignment %i, shift by %i elements", (reinterpret_cast<uintptr_t>(&T[i_z2]) % 32), (int)(4-(i_z2 & 0b11)));
            int i = i_z2;
            int check_alignment = 4-(i_z2 & 0b11); // fast modulo with 4
            if(check_alignment) { // Pretty sure this doesn't need an if statement
              for(int j = 0; j < check_alignment; j++) {
                T[i] = boundary_temp;
                //T_old[i] = boundary_temp;
                i++;
              }
            }
            for (; i + 4 < nlocal_zstop; i += 4) { // changed from nlocal
              _mm256_stream_pd(&T[i], T0_vec);
              //_mm256_stream_pd(&T_old[i], T0_vec);
            }
            for (; i < nlocal_zstop; i++) {
               T[i] = boundary_temp;
               //T_old[i] = boundary_temp; // T becomes T_old before next loop, should be way to skip this
            }
            //comm->all(); // Pretty sure this is useless, every process should have reset any temps it can see, nothing else changes here
            return;
        }
    }
    
    // can probably move everything in this loop into the T = boundary temp loop
    if (new_layer==1) {
      bool top_layer = 0;
      for (int i=i_ztop; i<nlocal; i++) { // we can start at nlocal_zstop
          //Stop looping if we're above the meltspot
          if(!top_layer && xyz[i][2]==z_meltspot) {
            i_ztop = i;
            if(domain->me == 0) {
                fprintf(screen,"top start at index %i\n",i_ztop);
            }
            top_layer = 1;
            //activeFlag[i] = 1;
          }
          else if(xyz[i][2] > floor(z_meltspot)) {
            nlocal_zstop = i;
            //fprintf(screen,"Local stop index %i\n",nlocal_zstop);
            break;
          }
          activeFlag[i] = POWDER;
      }
      new_layer=0;
      comm->all(); // only really need to comm activeFlag
    }

    // local_finitedifference();
    
    for (int i=i_z2; i<nlocal_zstop; i++) {
        //If below melt spot, run finite difference
        site_event_finitedifference(i);
        if(T[i] > Tl) {
          tempMax = Tl;
          if(activeFlag[i] != MOLTEN) {   
            //Let's also update the active flag after each FD loop
            //This is also the place to handle nucleation and solidification front impingement
            //Go from solid to molten
            //Or go from powder to molten
            //Also randomize spin and reset cumulative variables
            spin[i] = (int) (nspins * ranapp->uniform());
            SolidD[i] = 0.0;
            MobilityOut[i] = 0;
            activeFlag[i] = MOLTEN;
          }
        }
        //If we're molten, call the mushy_phase function to figure out any phase change
        else if (activeFlag[i] == MOLTEN) {
            mushy_phase(i, ranapp);
            tempMax = MAX(tempMax, T[i]);
        }
        else {
            tempMax = MAX(tempMax, T[i]);
        }  
    }
    
    timer->stamp(TIME_APP);

    //re-sync all the data
    comm->all();
    MPI_Allreduce(&tempMax,&tempMaxAll,1,MPI_DOUBLE,MPI_MAX,world);
    timer->stamp(TIME_COMM);

}
/* ----------------------------------------------------------------------
 We'll often need to restart these simulations. We calculate our spot location by tracking
 an index. Therefore, it'll be easiest to rewrite the path file so that it restarts at our
 final location.
 ------------------------------------------------------------------------- */
void AppAdditiveThermal::path_file_update()
{
	//Output a restart_path file.
	std::string restart_file_str("restart_path.txt");
	FILE* fout = fopen(restart_file_str.c_str(), "w");
	
	double start_distance = t_scan_array[path_index] + d_residual;
	
	//The first index value should be our current location, not array value
	//We might be off by one for the first index
	if(p_scan_array[path_index] == 0) {
		x_scan_array[path_index] = x_meltspot;
		y_scan_array[path_index] = y_meltspot;
		z_scan_array[path_index] = z_meltspot;
	}	
	t_scan_array[path_index] = 0;

	for( int i = path_index; i < line_count; i++) {
		//Update the scan distance
		if(i > path_index) {
			t_scan_array[i] = t_scan_array[i] - start_distance;
		}
		fprintf(fout, "%f,\t%f,\t%f,\t%f,\t%f\n", x_scan_array[i],y_scan_array[i],z_scan_array[i],t_scan_array[i],p_scan_array[i]);
	}
	fclose(fout);
}


/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with no null bin rejection
   flip to random neighbor spin without null bin
   technically this is an incorrect rejection-KMC algorithm
------------------------------------------------------------------------- */
 
  // Possible branchless upgrades:
  // if efinal<=einitial: cmove new state to spin[i] if RNG<=Mobloc
  //    pros: branchless
  //    cons: none (ecxept maybe on intel)
  // change else if to else: cmove new state to spin[i] if RNG<=Mobloc*exp((einitial-efinal)*t_inverse))
  //    pros: branchless
  //    cons: need to be extra sure that logic is correct
  // could probably do something weird with a temp variable, but not sure it would help:
  //    if efinal<=einitial cmove new state to temp, then do RNG<=Mobloc cmov with temp variable
  //    pros: eliminates branches from semi-random efinal check
  //    cons: always does RNG Mobloc check
  //    If efinal is usually less than einitial, then this helps, as the branch predictor would usually do 2nd check anyways
  //    If efinal is usually more than einitial, then this does unnecessary work, though saves branch missprediction time in case of efinal<=einitial
  //    Bottom line: Would need to test. Not sure of temp variable overhed.
  // keep if temp==0, but change to do nothing. Probably no chance of branch miss because always the same. Needed before Check 3.

void AppAdditiveThermal::site_event_rejection(int i, RandomFast *random)
{
  int current_state = spin[i];
  double einitial = site_energy(i,current_state);

  //Assign the local mobility
  double Mobloc = MobilityOut[i];
  
  int j,m,value;
  int nevent = 0;
  bool flip = 0;

  if((Mobloc < 0.0) | (Mobloc > 1.0001)) { // as far as I can tell, a negative mobility can never happen
      MobilityOut[i] = 0;
      return;
  }

  for (j = 0; j < numneigh[i]; j++) {
    int neighid = neighbor[i][j];
    value = spin[neighid];
    // Exclude gas or molten sites from the Potts neighbor tally
    // not sure what value==nspins contributes
    if (value == current_state|| activeFlag[neighid] == 0 || activeFlag[neighid] ==2 || value == nspins) continue;

    for (m = 0; m < nevent; m++) 
      if (value == unique[m]) goto dont_add_neighbor;
    unique[nevent++] = value;
    dont_add_neighbor:;
  }

  if (nevent == 0) return;
  int iran = (int) ((nevent-0.000000001)*random->uniform());
  //if (iran >= nevent) iran = nevent-1;  // this should never happen because int truncates
  int temp_spin = unique[iran];
  double efinal = site_energy(i,temp_spin);

  // accept or reject via Boltzmann criterion
  // Note: temperature here refers to simulation temperature from in.additive, not the local simulated temperature

  if (efinal <= einitial) {
     if (random->uniform() <= Mobloc){
       spin[i] = temp_spin;
       naccept++;
       flip = 1;
     }
  }
  else if (temperature == 0.0) {
    return;
  } 
  else if (random->uniform() <= Mobloc * exp((einitial-efinal)*t_inverse)) { 
    spin[i] = temp_spin;
    naccept++;
    flip = 1;
  }

  // set mask if site could not have changed
  // if site changed, unset mask of sites with affected propensity
  // OK to change mask of ghost sites since never used

  if (Lmask) {
    if (einitial < 0.5*numneigh[i]) mask[i] = 1;
    if (flip)
      for (int j = 0; j < numneigh[i]; j++)
  	    mask[neighbor[i][j]] = 0;
  }
}



/* ----------------------------------------------------------------------
   compute energy of site given state
------------------------------------------------------------------------- */

double AppAdditiveThermal::site_energy(int i, int test_spin)
{
  int isite = test_spin;
  int eng = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (isite != spin[neighbor[i][j]]) eng++;
  return (double) eng;
}

/* ----------------------------------------------------------------------
   Perform evolution for sites in the mushy zone. There are several things that need to happen.
   1. Determine if the site is solid or liquid (from activeFlag)
   2. Determine if the site should nucleate a new grain (from nucleationFlags)
   3. If our current site is liquid & can't nucleate, have it try to switch to a solid neighbor
      (with 4 calculated from the undercooling in someway, not temperature)
   4. If our current site is liquid & can nucleate, check if local undercooling is equal to its 
      critical temp. If so, change the activeFlag value to solid. If not, see if there are any solid sites
      that should capture it.
   5. If our current site is solid, see if it should flip to a neighboring solid value (with
      mobility calculated from undercooling.)
------------------------------------------------------------------------- */
void AppAdditiveThermal::mushy_phase(int i, RandomFast *random){
  	int nevent = 0;
  	int m,value;
    double Tcool = Tl - T[i];
     //For default settings, SolidD[i] =+ (1.091e-5 * pow(Tcool, 3) - 2.034e-4 * pow(Tcool, 2) + 2.74e-3*Tcool + 1.151e-4) * time_step;
       
    //Our site should always be molten and below Tl
    //Check if it's eligible to nucleate
    //if(nucleationFlags[spin[i]]) {
    if(spin[i]<nucleationCutoff && Tcool >= nucleationTemps[spin[i]]) {
        //Can and will nucleate
        activeFlag[i] = SOLID;
        //Don't let nucleated site disappear during smoothing
        SolidD[i] = -nsmooth-2;
        return;
    }
    
    int power = solid_front_length -1;
    for(int k = 0; k < solid_front_length; k++) {
        SolidD[i] += solid_front_coeffs[k] * pow(Tcool, power) * time_step;
        power--;
    }
// Attmept to make this slightly simpler:
// Requires reversing order of coeffs. Introduces a bug.
//    for(int k = 0; k < solid_front_length; k++) {
//        SolidD[i] += solid_front_coeffs[k] * pow(Tcool, k) * time_step;
//    }
   
    if(neighDist[2] <= SolidD[i]) { 
      for (int j = 0; j < numneigh[i]; j++) {
        if(activeFlag[neighbor[i][j]] == 3 ) { //&& neighDist[j] < SolidD[i] // swap to this if not waiting for all neighbors to be close enough
            unique[nevent++] = spin[neighbor[i][j]];	
        }
      }
    }
    
    //If no neighbor is eligible, return before changing anything. Will try next sweep.
    if (nevent == 0) return;
    int iran = (int) ((nevent-0.0000001)*random->uniform()); // Need to double check bounds of RNG to decide if this is necessary
    //if (iran >= nevent) iran = nevent-1;
    spin[i] = unique[iran];
    activeFlag[i] = SOLID;
    SolidD[i] = -1;
    return;

}

/* ----------------------------------------------------------------------
    Only nucleating one site at a time introduce lattice size dependency. Here we will
    use a user-defined nucleation particle size and flip neighboring sites until that size is met
------------------------------------------------------------------------- */
void AppAdditiveThermal::nucleation_particle_flipper(int i, int partRad, RandomFast *random) {
    
    //If one site is big enough to satisfy, skip evertyhing
    if(partRad <= 0) return;
    int nSites = fmin(26,partRad);
    static uint8_t neigh_shells [ ] = {4,10,12,15,21,13, 9,16,3,22,14,11,20,5,1,24,7,18, 0,25,23,2,6,19,17,8};
    int flipped = 0;
    
    //Its hard to go through shells iteravely if the number of neighbors isn't the full 26.
    //Check if this is the case and just loop through the neihgbor list if so
    if(numneigh[i] == 26) {
        //Go through neighbors and nucleate liquid ones
        //Do 1st ones first (this isn't random yet)
        for(int j = 0; j < nSites; j++) {
            int neighid = neighbor[i][neigh_shells[j]];
            if(activeFlag[neighid] == 2) {
                spin[neighid] = spin[i];
                activeFlag[neighid] = 3;
                SolidD[neighid] = -nsmooth -3;
                flipped++;
            }
        }
    }
    else {        
        for(int j = 0; j < numneigh[i]; j++) {
            int neighid = neighbor[i][j];
            if(activeFlag[neighid] == 2) {
                spin[neighid] = spin[i];
                activeFlag[neighid] = 3;
                SolidD[neighid] = -nsmooth -3;
                flipped++;
            }
        }
    }
    //If we didn't fill any sites, or flipped the correct number of sites then return
    if (flipped == 0 || flipped == nSites) {
        return;
    }
    
    int nneigh = 0;
    int possible_neigh[26];
    //If we still haven't satisfied the particle size, pick a neighbor at random and solidify from there.
    //Build a list of same-particle neighbors and pick one randomly
    for(int j =0; j < numneigh[i]; j++) {
        int neighid = neighbor[i][j];
        if(spin[neighid] == spin[i] && activeFlag[neighid] == 3) {
            possible_neigh[nneigh] = j;
            nneigh++;
        }
    }
    //If no possible nieghbors, quit
    if(nneigh == 0) {
        return;
    }
    //Otherwise, randomly pick a possilbe neighbor
    int neighran =  round(((nneigh -1) * random->uniform()));
    nucleation_particle_flipper(neighbor[i][possible_neigh[neighran]],nSites-flipped, random);
    return;    
}


/* ----------------------------------------------------------------------
   Nucleation site initializer. Find the volume of a voxel from dx^3 and Multiply by No.
   This will be the average number of nucleation sites in the voxel. This is also the
   fraction of spins that we want to be able to nucleate new grains. If the value is greater
   than 1 (which we should avoid), allow all spins to nucleate. If not, call a random number
   between zero and one. If the number is less than the fraction, make true. If not, make false.
------------------------------------------------------------------------- */
void AppAdditiveThermal::nucleation_spins(RandomFast *random) {
    double nucleationFraction = dx * dx * dx * No;
    
    //Make all spins nucleation sites. Should avoid this.
    if(nucleationFraction >= 1.0) {
        fprintf(screen,"Nucleation fraction (%f) is greater than 1. Decrease No or increase mesh resolution.\n", nucleationFraction);
//        for (int i = 0; i < nspins; i++) {
//            nucleationFlags[i] = 1;
//        }
    }
    //Do a random number test and allow the spin to nucleate if less than
//    else {
//        fprintf(screen,"Nucleation Fraction is %f \n", nucleationFraction);
//        for (int i = 0; i < nspins; i++) {
//            if(random->uniform() <= nucleationFraction) {
//                nucleationFlags[i] = 1;
//            }
//            else {
//                nucleationFlags[i] = 0;
//            }
//        }
//    }   
    nucleationCutoff = nucleationFraction * nspins;
}

/* ----------------------------------------------------------------------
    The first version of this just initialized critical nucleation temperatures.
    This version will also initialize nucleii size (starting with a normal dist)
------------------------------------------------------------------------- */
void AppAdditiveThermal::nucleation_init() {
    std::normal_distribution<> dist_T{Tc,Tsig};
    std::normal_distribution<> dist_S{sizeNorm,sizeSig};
    std::random_device rd{};
    std::mt19937 gen{rd()};
    
    //Randomly assign a temperature to every spin
    for(int i = 0; i < nspins; i++) {
    //for(int i = 0; i < nucleationCutoff; i++) {
        nucleationTemps[i] = dist_T(gen);
        nucleationSizes[i] = dist_S(gen);
    }
}

/* ----------------------------------------------------------------------
   Evolve simulation time including rejection KMC solver. This function controls
   the relative occurence of solidification/thermal timesteps and solid-state kMC steps
   by comparing the dtMC and time_step variables. time_step is a pre-determined non-constant variable,
   while dtMC depends on lattice size and maximum solid-phase temperature.
 ------------------------------------------------------------------------- 
*/

void AppAdditiveThermal::iterate_rejection(double stoptime)
{
  int i,icolor,nselect,nrange,jset, nMC;
  int *site2i;
  //double tempMax;
  //double tempMaxAll;
  tempMax = 0;
  tempMaxAll = 0;
  double FDElapsed;
  double nucVolume;
  double mobMax = 0;
  double* pointer_swap;
  
  static const double div_dxcubed = 1/pow(dx,3);
  
  // set loop is over:
  // sectors if there are sectors and no colors
  // colors if there are colors and no sectors
  // first nsector sets if there are both sectors and colors

  int nset_loop = nset;
  if (bothflag) nset_loop = nsector;

  int done = 0;
  nMC = 1;
  
  for(int i = 0; i < nlocal; i++) {
    MobilityOut[i] = 0;
  }

  z_meltspot = z_scan_array[0];
  MPI_Bcast(&z_meltspot,1,MPI_DOUBLE,0,world);
  bool top_layer = 0;
  bool z2 = 0;
  for(int i = 0; i < nlocal; i++) {
    if(!z2 && xyz[i][2]!=domain->boxzlo) {
      i_z2 = i;
      if(domain->me == 0) {
          fprintf(screen,"2nd layer start at index %i\n",i_z2);
      }
      z2 = 1;
    }
    else if(!top_layer && xyz[i][2]==z_meltspot) {
      i_ztop = i;
      if(domain->me == 0) {
          fprintf(screen,"top start at index %i\n",i_ztop);
      }
      top_layer = 1;
    }
    else if(xyz[i][2]>z_meltspot) {
      nlocal_zstop = i;
      break;
    }
  }
  
    //Find the highest temperature in the local array and compare with others
  tempMax = boundary_temp;
  //timer->stamp(TIME_APP);
  //MPI_Allreduce(&tempMax,&tempMaxAll,1,MPI_DOUBLE,MPI_MAX,world);
  //timer->stamp(TIME_COMM);
  //Using the global maximum temperature, compute our smallest timestep 
  dtMC = compute_timeMin(tempMaxAll);
  timer->stamp(TIME_APP);
  //This while loop is for the entire simulation run time! (stoptime = "run stoptime")
  while (!done) {
  
//      //Find the highest temperature in the local array and compare with others
//      tempMax = compute_tempMax();
//      timer->stamp(TIME_APP);
//      MPI_Allreduce(&tempMax,&tempMaxAll,1,MPI_DOUBLE,MPI_MAX,world);
//      timer->stamp(TIME_COMM);
//      //Using the global maximum temperature, compute our smallest timestep 
//      dtMC = compute_timeMin(tempMaxAll);
      
      FDElapsed = 0;
      
      for(int i = i_z2; i < nlocal_zstop; i++) {
        MobilityOut[i] = 0;
      }
      mobMax = 0;

      while(dtMC >= FDElapsed) {
        pointer_swap = T_old;
        T_old = T;
        T = pointer_swap;

        //Update temperatures and phases
        app_update();

        //Using the global maximum temperature, compute our smallest timestep
        //Only update dtMC if it is smaller (higher temperature)
        //Basing MC calculation off of the highest temp observed in time_step
        // this can technically be reduced by computing the mob max to test with, then using inverse to caclulate dtMC if true
        if (dtMC > compute_timeMin(tempMaxAll) || mobMax < 1e-8) {
          dtMC = compute_timeMin(tempMaxAll);
          //If new mobMax is larger, use it. This is tested by the enclosing if statement
          mobMax = exp(-Q/(R*tempMaxAll));
        }

        //Add to running mobility values at each site. Multiply mobility by timestep_size
        //which makes the integral of a constant function
        for(int i = i_z2; i < nlocal_zstop; i++) { //nlocal_zstop
          //if(xyz[i][2] > floor(z_meltspot)) break;
          if(activeFlag[i] != SOLID) continue;
          //Also check if we're just solidified and should be "relaxed"
          if(SolidD[i] < 0 && SolidD[i] > -nsmooth -1)    {
            MobilityOut[i] = 1;
            site_event_rejection(i, ranapp);
            SolidD[i]--;
            continue;
          }
          //Check if we just nucleated and need to grow larger
          else if(SolidD[i] == -nsmooth -2) {
            nucVolume = nucleationSizes[spin[i]];
            nucleation_particle_flipper(i, round(nucVolume*div_dxcubed), ranapp);
            SolidD[i] = -nsmooth - 3;
          }
          MobilityOut[i] += time_step * compute_mobility(i);//,ranapp);
        }
        FDElapsed += time_step;
        time += time_step;
        timer->stamp(TIME_SOLVE);
        if (time >= stoptime || path_index > line_count)   done = 1;
        timer->stamp(TIME_OUTPUT);
        if (done) {
          if(domain->me == 0) {
                  fprintf(screen,"done\n");
          }
          break;
        }
      }

      //Compute the normalized mobilities at each site
      //True "mobMax" would be holding the max temp for the entire window
      double mobility_denom = 1/(mobMax * FDElapsed);
      for(int i = i_z2; i < nlocal_zstop; i++) { 
        if(activeFlag[i] == SOLID) MobilityOut[i] = MobilityOut[i] * mobility_denom;
        //MobilityOut[i] = MobilityOut[i]/(mobMax * FDElapsed);
      }  	

    //Do Monte Carlo sweeps
    // Note: this seems more complicated than necessary, I think there is logic for things that will never be used here, but not going to touch for now.
    //for (int j = 0; j<nMC; j++) {
      for (int iset = 0; iset < nset_loop; iset++) {
        if (nprocs > 1) {
          timer->stamp();
          if (sectorflag) comm->sector(iset);
          else comm->all();
          timer->stamp(TIME_COMM);
        }

        if (Lmask) boundary_clear_mask(iset);

        timer->stamp();

        // sectors but no colors (could also be no sectors)
        // random selection of sites in iset
        if (sweepflag == RANDOM) {
          site2i = set[iset].site2i;
          nrange = set[iset].nlocal;
          nselect = set[iset].nselect;
          //std::vector<int> vals = ranapp->m_irandom(nrange,nselect);
          //int *vals = ranapp->m_irandom(nrange,nselect);
          for (i = 0; i < nselect; i++) 
            //sitelist[i] = site2i[vals[i] - 1];
            sitelist[i] = site2i[ranapp->irandom(nrange) - 1];
          (this->*sweep)(nselect,sitelist);
          nattempt += nselect;

          // sectors but no colors, or colors but no sectors
          // ordered sweep over all sites in iset
        } else if (bothflag == 0) {
          //Get rid of nloop because we want this to happen once each time we run it.
          (this->*sweep)(set[iset].nlocal,set[iset].site2i); 
          nattempt += set[iset].nselect;

          // sectors and colors
          // icolor loop is over all colors in a sector
          // jset = set that contains sites of one color in one sector
          // ordered sweep over all sites in jset
        } else {
          for (icolor = 0; icolor < ncolors; icolor++) {
            jset = nsector + iset*ncolors + icolor;
            //Get rid of nloop because we want this to happen once each time we run it.
            (this->*sweep)(set[jset].nlocal,set[jset].site2i);
            nattempt += set[jset].nselect;
          }
        }

        timer->stamp(TIME_SOLVE);

        if (nprocs > 1) {
          if (sectorflag) comm->reverse_sector(iset);
          else comm->all_reverse();
          timer->stamp(TIME_COMM);
        }
      }

      nsweeps++;
    if (dtMC > 1) time += dtMC;
    if (time >= stoptime) done = 1;
    timer->stamp(TIME_OUTPUT);

  }
  nextoutput = output->compute(time+1.0,1);
}

/* ----------------------------------------------------------------------
   Function to determine the highest temperature in a region and compute the corresponding
   timestep. Let's compute timestep after and just return highest temperature
 ------------------------------------------------------------------------- */
double AppAdditiveThermal::compute_tempMax() {
	tempMax = 0;
  for(int i = nlocal_zstop; i >= i_ztop; i--) { 
      if(activeFlag[i] == POWDER) continue;
      //If max temp is above liquidus, just make it liquidus and return
      if(T[i] > Tl) {
            tempMax = Tl;
            return tempMax;
	    }
	    //tempMax = MAX(tempMax, T_old[i]); 
      tempMax = MAX(tempMax, T[i]); 
	}
	return tempMax;
}

/* ----------------------------------------------------------------------
   Given a maximum temperature, compute the minimum MC timestep.
 ------------------------------------------------------------------------- */
double AppAdditiveThermal::compute_timeMin(double tempMax) {
    static const double factor = pow(dx,2) * Kmc/Ko;
    double dtMC;
    dtMC = factor * exp(Q/(R*tempMax));
    return dtMC;
}
