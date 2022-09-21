#include "pluto.h"
#include "ReferenceSystem.h"
#include "Molecules.h"
#include "analysis.h"
#include "boundaries.h"
#include "boundary_fluxes.h"
#include "ReadRestrictionsConfiguration.h"
#include "DustEvolution.h"
#include "PhysicalConstantsCGS.h"
#include "MakemakeTools.h"
#include "Radiation.h"
#include "Molecules.h"
#include "Opacity.h"
#include "Irradiation.h"
#include "AstroTools.h"

#if SELFGRAVITY == YES
#include "SelfGravity.h"
#endif

#if FLUXLIMITEDDIFFUSION == YES
#include "FLD.h"
#endif

#if STELLAREVOLUTION == YES
#include "StellarEvolution.h"
#endif

#if IONIZATION == YES
#include "Ionization.h"
#include "MaterialProperties.h"
#include "DirectRecombination.h"
#endif





/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3, Data *data, Grid *grid, int kl, int jl, int il)
/*
 *
 *
 *
 *********************************************************************** */
{
	static int first_call = 1;
    
    //printf("dx = %e\n", grid[IDIR].dx[2]);
	//printf("dy = %e\n", grid[JDIR].dx[2] * grid[IDIR].x[2]);
	//exit(1);
    
	// ********************
	// * Reference System *
	// ********************
	//
#if RADIATION == YES || IONIZATION == YES || STELLAREVOLUTION == YES
	if(first_call) InitializeReferenceSystem();
#endif
	
    
	// *****************
	// * Density setup *
	// *****************
	//

    us[RHO] = g_inputParam[IC_rho_cgs]/ReferenceDensity;
	
	// ********************
	// * VelocityXY setup *
	// ********************
	//
	us[VX1] = 0.0;
#if COMPONENTS > 1
	us[VX2] = 0.0;
#endif
#if COMPONENTS > 2
	us[VX3] = 0.0;
#endif
    
    double NeutralGasTemperature;
    NeutralGasTemperature = 100.0;
    
    
    // *************
    // * Radiation *
    // *************
    //
#if RADIATION == YES
    data->DustToGasMassRatio[kl][jl][il] = Dust2GasMassRatio;

#if PDVINRT == YES
    data->fun[kl][jl][il]                     = us[PRS];
#endif

#if FLUXLIMITEDDIFFUSION == YES
    data->RadiationEnergyDensity[kl][jl][il]  = a_RadiationConstant * pow(NeutralGasTemperature, 4);
#else
    data->RadiationEnergyDensity[kl][jl][il]  = 0.0;
#endif
    
#if IRRADIATION == YES
    data->IrradiationPowerDensity[kl][jl][il] = 0.0;
#endif

#endif
    
    
    // **************
    // * Ionization *
    // **************
    //
#if IONIZATION == YES
    data->IonizationFraction[kl][jl][il] = 0.0;
    data->NeutralFraction[kl][jl][il]    = 1.0;
    if(urec_min >= 0) data->DirectRecombinationPhotonDensity[kl][jl][il] = urec_min;
    else              data->DirectRecombinationPhotonDensity[kl][jl][il] = 0.0;
    data->IonizationPowerDensity[kl][jl][il] = 0.0;
#endif
	
    
    // ******************
    // * Pressure setup *
    // ******************
    //
    us[PRS] = R_UniversalGasConstant / MolarMass(data,kl,jl,il) * us[RHO] * NeutralGasTemperature;
    
    data->DustTemperature[kl][jl][il] = NeutralGasTemperature;
    data->GasTemperature[kl][jl][il]  = NeutralGasTemperature;

    
    //printf("%e %e\n", us[RHO]*ReferenceDensity, us[PRS]*ReferenceEnergyDensity);
    
  	first_call = 0;
}



/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid
#if STELLAREVOLUTION == YES
               , Star *star
#endif
)
/*
 *
 *
 *********************************************************************** */
{
    static int first_call = 1;
    
    WriteAnalysis(d, grid
#if STELLAREVOLUTION == YES
                  , star
#endif
                  );
    
    first_call = 0;
}



/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*!
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and
 *                    staggered magnetic fields (d->Vs, when used) to
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END,
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
    int kl, jl, il;
    double radius, Rion, x, y, Tgas, dRion;

    Rion = g_inputParam[IBC_ionization_radius_pc] * Parsec;
    dRion = g_inputParam[IBC_delta_ionization_radius_pc] * Parsec;

    if(side == 0){
#if IONIZATION == YES
        if(Rion > 0){
            DOM_LOOP(kl, jl, il){
                radius = grid[IDIR].x[il];
                
                // Step function:
                if(radius <= Rion) x = 1.0;
                else x = 0.0;
                
                // Smooth transition:
//                x = 1.0/M_PI * atan(- 1.0/dRion * (radius - Rion)) + 0.5;
                
                y = 1.0 - x;
                d->IonizationFraction[kl][jl][il] = x;
                d->NeutralFraction[kl][jl][il]    = y;
                Tgas                              = x * IonizedGasTemperature + y * NeutralGasTemperature;
                d->DustTemperature[kl][jl][il]    = Tgas;
                d->GasTemperature[kl][jl][il]     = Tgas;
                d->Vc[PRS][kl][jl][il]            = R_UniversalGasConstant / MolarMass(d,kl,jl,il) * d->Vc[RHO][kl][jl][il] * Tgas;
                
                //                printf("arg = %e\n", - 1.0/dRion * (radius - Rion));
                //                printf("x = %e\n", x);
            }
        }
#endif
    }
    else if(side == X1_BEG){
        BC_Outflow_NoInflow(d, box, side, grid);
#if COMPONENTS > 1
        BC_ZeroGradient(d, box, side, grid, VX2);
#endif
#if COMPONENTS > 2
        BC_ZeroGradient(d, box, side, grid, VX3);
#endif
    }
    else if(side == X1_END){
        BC_Outflow_NoInflow(d, box, side, grid);
#if COMPONENTS > 1
        BC_ZeroGradient(d, box, side, grid, VX2);
#endif
#if COMPONENTS > 2
        BC_ZeroGradient(d, box, side, grid, VX3);
#endif
    }
}



void Test(const Data *data, Grid *grid){
  
 int    kl, jl, il;
 int    i_HII;
 double x, y, r, ustar, urec, R_HII, R_Stroemgren, deviation, devminus, devplus, photon_number;
  
//    R_Stroemgren = StroemgrenRadius(ConstantEUVPhotonRate, Ionization_alpha2, NH * mH_HydrogenMass / mtot * g_inputParam[IC_rho_cgs]/ReferenceDensity);
//    R_Stroemgren = StroemgrenRadius(IrradiationLuminosityPerFrequencyBin[NumberOfFrequencies-1] / (13.6*eV), RecombinationRateIntoAnyButGroundState(IonizedGasTemperature), NH * mH_HydrogenMass / mtot * g_inputParam[IC_rho_cgs]/ReferenceDensity);
//    R_Stroemgren = StroemgrenRadius(IrradiationLuminosityPerFrequencyBin[NumberOfFrequencies-1] / (13.6*eV), RecombinationRateIntoAnyButGroundState(NeutralGasTemperature), NH * mH_HydrogenMass / mtot * g_inputParam[IC_rho_cgs]/ReferenceDensity);
  R_HII        = -1.0;
  i_HII        = -1;
  photon_number        = IrradiationLuminosityPerFrequencyBin[NumberOfFrequencies-1] / (13.6*eV);

 DOM_LOOP(kl,jl,il){
     r     = grid[IDIR].x[il];
     x     = data->IonizationFraction[kl][jl][il];
     y     = data->NeutralFraction[kl][jl][il];
     ustar = data->StellarEUVPhotonDensity[kl][jl][il];
     urec  = data->DirectRecombinationPhotonDensity[kl][jl][il];
         printf("%4d - %e pc: %e cm^-3 %e cm^-3 %e %e\n", il, r/Parsec, ustar*ReferenceNumberDensity, urec*ReferenceNumberDensity, x, y);
     if(x > 0.5){
         i_HII = il;
         R_HII = r;
         //            printf("i_HII = %d\n", i_HII);
         //            printf("R_HII = %e\n", R_HII);
     }
 }
//
//   double ionized_gas_temperature, M_star, R_critical, sound_speed_squared, radius_deviation, trapping_density, expansion_density, limiting_density;
//
//   ionized_gas_temperature = IonizedGasTemperature / ReferenceTemperature;
//   M_star = 20.0 * SolarMass;
//
//   sound_speed_squared = g_gamma * ionized_gas_temperature * (R_UniversalGasConstant / MolarMass(data,kl,jl,il));
//   
//   R_critical = (G_GravityConstant * M_star) / (2.0 * sound_speed_squared);
//
//   deviation = fabs(R_HII                 - R_Stroemgren) / R_Stroemgren;
//   devminus  = fabs(grid[IDIR].x[i_HII-1] - R_Stroemgren) / R_Stroemgren;
//   devplus   = fabs(grid[IDIR].x[i_HII+1] - R_Stroemgren) / R_Stroemgren;
//   
   printf("\n");
//   printf("Test results:\n");
//   printf(" R_Stroemgren    = %e pc\n"    , R_Stroemgren / Parsec);
   printf(" R(%4d)         = %e pc\n", i_HII-1, grid[IDIR].x[i_HII-1] / Parsec);
   printf(" R_HII = R(%4d) = %e pc\n", i_HII  , R_HII                 / Parsec);
   printf(" R(%4d)         = %e pc\n", i_HII+1, grid[IDIR].x[i_HII+1] / Parsec);
    exit(1);
//   printf(" R_critical    = %e pc\n", R_critical / Parsec);
//   printf(" Photon number = %e cgs\n", photon_number);
//   printf(" sound speed = %e cgs \n", sqrt(sound_speed_squared) * ReferenceLength / ReferenceTime);
    
    //    if(deviation > devminus || deviation > devplus){
    //        printf("ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR\n");
    //        printf("ERROR \n");
    //        printf("ERROR  Test failed! \n");
    //        printf("ERROR \n");
    //        printf("ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR\n");
    //        exit(1);
    //    }
    // if(deviation > 1e-2){
    //    printf("ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR\n");
    //    printf("ERROR \n");
    //    printf("ERROR  Test failed! \n");
    //    printf("ERROR \n");
    //    printf("ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR\n");
    //    exit(1);
    // }
    //else{
      //  printf("\n");
      //  printf("Test successful! \n");
      //  printf("\n");
      //  exit(1);
        
    // }

}

