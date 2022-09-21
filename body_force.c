#include "pluto.h"
#include "boundary_fluxes.h"
#include "PhysicalConstantsCGS.h"
#include "ReadRestrictionsConfiguration.h"
#include "ReferenceSystem.h"

#if SELFGRAVITY == YES
#include "SelfGravity.h"
#endif

#if RADIATION == YES
#include "Radiation.h"
#endif

#if IRRADIATION == YES
#include "Irradiation.h"
#endif

#if FLUXLIMITEDDIFFUSION == YES
#include "FLD.h"
#endif

#if RADIATION == YES || IONIZATION == YES
#include "Opacity.h"
#endif

#if IONIZATION == YES
#include "Ionization.h"
#include "DirectRecombination.h"
#include "MaterialProperties.h"
#endif

#if CAKACCELERATION == YES
#include "CAKAcceleration.h"
#endif


#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3, int il, int jl, int kl, const Data *data, Grid *grid, Time_Step *Dts)
/*
 *
 *
 *
 *************************************************************************** */
{
    
    double agrav, Mgrav;
    
    g[IDIR] = 0.0;
    g[JDIR] = 0.0;
    g[KDIR] = 0.0;
    
    
    
    //
    // Central point-mass Gravity:
    //
    agrav    = - G_GravityConstant * g_inputParam[IC_central_mass_Msol] * SolarMass / (x1*x1);
    g[IDIR] += agrav;
    Dts->dt_grav = MIN(Dts->dt_grav, sqrt(- 2.0 * grid[IDIR].dx[il] / agrav));

    
    
    //
    // Irradiation:
    //
#if IRRADIATION == YES
    double aradstar = 0.0;
    if(IrradiationFlag != 0){
        if(StellarRadiativeForceFlag == 0) aradstar = 0.0;
        else if(StellarRadiativeForceFlag == 1) aradstar = data->IrradiationPowerDensity[kl][jl][il] / c_SpeedOfLight / data->Vc[RHO][kl][jl][il];
        else printf("ERROR: StellarRadiativeForceFlag = %d not available.\n", StellarRadiativeForceFlag);
    }
    g[IDIR] += aradstar;
#endif
    
    //
    // Ionization:
    //
#if IONIZATION == YES
    double aionstar = 0.0;
    if(IonizationFlag != 0){
        if(IonizationRadiativeForceFlag == 0) aionstar = 0.0;
        else if(IonizationRadiativeForceFlag == 1) aionstar = data->IonizationPowerDensity[kl][jl][il] / c_SpeedOfLight / data->Vc[RHO][kl][jl][il];
        else printf("ERROR: IonizationRadiativeForceFlag = %d not available.\n", IonizationRadiativeForceFlag);
    }
    g[IDIR] += aionstar;
#endif
    
    
    
    
    
#if FLUXLIMITEDDIFFUSION == YES
    if(RadiationFlag != 0){
        if(ThermalRadiativeForceFlag == 0){
            // No Thermal Radiative Force
            //            aradtherm = 0.0;
        }
        else if(ThermalRadiativeForceFlag == 1){
            
            g[IDIR] += ThermalRadiativeAcceleration(data, grid, kl, jl, il, IDIR);
#if DIMENSIONS > 1
            g[JDIR] += ThermalRadiativeAcceleration(data, grid, kl, jl, il, JDIR);
#endif
#if DIMENSIONS > 2
            g[KDIR] += ThermalRadiativeAcceleration(data, grid, kl, jl, il, KDIR);
#endif
             }
        else if(ThermalRadiativeForceFlag == 2 || ThermalRadiativeForceFlag == 3){
            // thermal radiative force is handled as potential, see below.
            //            aradtherm = 0.0;
        }
        else{
            printf("ERROR: ThermalRadiativeForceFlag = %d not available.\n", ThermalRadiativeForceFlag);
            exit(1);
        }
    }
#endif
    
    
#if IONIZATION == YES
    if(DirectRecombinationFlag == 1 && RecombinationRadiativeForceFlag == 1){

        g[IDIR] += DirectRecombinationAcceleration(data, grid, kl, jl, il, IDIR);
#if DIMENSIONS > 1
        g[JDIR] += DirectRecombinationAcceleration(data, grid, kl, jl, il, JDIR);
#endif
#if DIMENSIONS > 2
        g[KDIR] += DirectRecombinationAcceleration(data, grid, kl, jl, il, KDIR);
#endif
    }
#endif


#if CAKACCELERATION == YES
    double gline[3];

    if(CAKa_T_STAR >= 7500.0 / ReferenceTemperature){

        data->gidir[kl][jl][il] = 0.0;
#if COMPONENTS > 1
        data->gjdir[kl][jl][il] = 0.0;
#endif
#if COMPONENTS > 2
        data->gkdir[kl][jl][il] = 0.0;
#endif

        if((grid[IDIR].lbound != 0 && il>=grid[IDIR].nghost) || grid[IDIR].lbound == 0){
            if(g_inputParam[CAK_ifrc] == 0){
                sobolev(data, grid, kl, jl, il, gline);
            }
            else{
                gcak3d(data, grid, kl, jl, il, gline);
            }
        }

        //
        // if Tgas != Tstar:
        //  Rolf, see e.g. Puls, Springmann, Lennon (2000), Fig. 4 to remember the effect.
        //
        //#if EOS != ISOTHERMAL
        //    double TemperatureCorrectionFactor = 1.0;
        //    RT = UNIT_VELOCITY*UNIT_VELOCITY*(0.599/8.3144621e7);
        //    twind = data->Vc[PRS][kl][jl][il]/data->Vc[RHO][kl][jl][il] * MolarMass(0.0) / R_UniversalGasConstant;
        //    if(twind > tstar/RT){
        //        TemperatureCorrectionFactor = exp(-twind/tstar+1);
        //    }
        //    gline[IDIR] *= TemperatureCorrectionFactor;
        //    gline[JDIR] *= TemperatureCorrectionFactor;
        //    gline[KDIR] *= TemperatureCorrectionFactor;
        //#endif
    }
    else{
        gline[IDIR] = 0.0;
        gline[JDIR] = 0.0;
        gline[KDIR] = 0.0;
    }

    g[IDIR] += gline[IDIR];
    g[JDIR] += gline[JDIR];
    g[KDIR] += gline[KDIR];

    data->gidir[kl][jl][il] = gline[IDIR];
    data->gjdir[kl][jl][il] = gline[JDIR];
    data->gkdir[kl][jl][il] = gline[KDIR];
#endif

    //    if(fabs(aradthermx) > 3.0 * fabs(agrav)){
    //        print("%d %d %d\n", kl, jl, il);
    //        print("agrav      = %+e km/s/dt\n", agrav      * g_dt * ReferenceVelocity * 1e-5);
    //        print("aradstar   = %+e km/s/dt\n", aradstar   * g_dt * ReferenceVelocity * 1e-5);
    //        print("aradthermx = %+e km/s/dt\n", aradthermx * g_dt * ReferenceVelocity * 1e-5);
    //        print(" dtauxl = %+e\n", kapparhoxl*GridCellLength(grid, IDIR, kl,   jl,   il-1));
    //        print(" dtau   = %+e\n", kapparho  *GridCellLength(grid, IDIR, kl,   jl,   il  ));
    //        print(" dtauxr = %+e\n", kapparhoxr*GridCellLength(grid, IDIR, kl,   jl,   il+1));
    //        print(" rhogasxl = %+e\n", rhogasxl);
    //        print(" rhogas   = %+e\n", rhogas  );
    //        print(" rhogasxr = %+e\n", rhogasxr);
    //        // (- Lxl * (ER   - ERxl) / dxl - Lxr * (ERxr - ER  ) / dxr) / rhogas;
    //        print("aradstar   / agrav = %+e\n", aradstar   / agrav);
    //        print("aradthermx / agrav = %+e\n", aradthermx / agrav);
    //    }
    
}
#endif







#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3, int il, int jl, int kl, const Data *data, Grid *grid, Time_Step *Dts)
/*
 *
 *
 *
 *************************************************************************** */
{
    
    double Potential;
    
    Potential = 0.0;
    
    //
    // Self-Gravity of gas in the computational domain:
    //
#if SELFGRAVITY == YES
    //
    // Potential at cell interfaces:
    //
    double dx, fmr, frr, phi, phir;
    
    if(SelfGravityFlag != 0){
        if(g_dir == IDIR){
            dx         = 0.5 * (GridCellLength(grid, IDIR, kl,   jl,   il  ) + GridCellLength(grid, IDIR, kl,   jl,   il+1));
            fmr        = 1.0 - 0.5 * GridCellLength(grid, IDIR, kl,   jl,   il  ) / dx;
            frr        = 1.0 - 0.5 * GridCellLength(grid, IDIR, kl,   jl,   il+1) / dx;
            phi        = data->phi[kl  ][jl  ][il  ];
            phir       = data->phi[kl  ][jl  ][il+1];
        }
        else if(g_dir == JDIR){
            dx         = 0.5 * (GridCellLength(grid, JDIR, kl,   jl,   il  ) + GridCellLength(grid, JDIR, kl,   jl+1, il  ));
            fmr        = 1.0 - 0.5 * GridCellLength(grid, JDIR, kl,   jl,   il  ) / dx;
            frr        = 1.0 - 0.5 * GridCellLength(grid, JDIR, kl,   jl+1, il  ) / dx;
            phi        = data->phi[kl  ][jl  ][il  ];
            phir       = data->phi[kl  ][jl+1][il  ];
        }
        else if(g_dir == KDIR){
            dx         = 0.5 * (GridCellLength(grid, KDIR, kl,   jl,   il  ) + GridCellLength(grid, KDIR, kl+1, jl,   il  ));
            fmr        = 1.0 - 0.5 * GridCellLength(grid, KDIR, kl,   jl,   il  ) / dx;
            frr        = 1.0 - 0.5 * GridCellLength(grid, KDIR, kl+1, jl,   il  ) / dx;
            phi        = data->phi[kl  ][jl  ][il  ];
            phir       = data->phi[kl+1][jl  ][il  ];
        }
        else{
            printf("ERROR:  Value of g_dir unknown. See body_force.c for details.\n");
            exit(1);
        }
        Potential += frr * phir + fmr * phi;
        
#if DEBUGGING > 0
        if(Potential >= 0){
            printf("ERROR: In body_force.c\n");
            printf("ERROR: g_dir=%d; Potential[%d][%d][%d] = %e >= 0\n", g_dir, kl, jl, il, Potential);
            if(g_dir == KDIR){
                printf("ERROR:  phi  = %e >= 0\n", phi);
                printf("ERROR:  phir = %e >= 0\n", phir);
            }
            exit(1);
        }
#endif
        
        //
        // Potential at cell center:
        //
        //if(il == IBEG) printf("%d %e\n", il, x1);
        //Potential += data->phi[kl][jl][il];
    }
#endif
    
    
    
    //
    // Thermal Radiative Force in Eddington approximation:
    //
#if FLUXLIMITEDDIFFUSION == YES
    //
    // Potential at cell interfaces:
    //
    // TODO: specify P_rad etc. at cell interfaces, see above
    
    //
    // Potential at cell center:
    //
    double R, L, rhogas, rhodust, Tdust, Tgas, Trad, kappadust, kappagas, kapparho, ER;
    double ERxl, ERxr, dxl, dxr;
    double ERyl, ERyr, dyl, dyr;
    double ERzl, ERzr, dzl, dzr;
    if(ThermalRadiativeForceFlag == 2){
        Potential  += data->RadiationEnergyDensity[kl][jl][il] / 3.0 / data->Vc[RHO][kl][jl][il];
    }
    else if(ThermalRadiativeForceFlag == 3){
        ER          = data->RadiationEnergyDensity[kl  ][jl  ][il  ];
        ERxl        = data->RadiationEnergyDensity[kl  ][jl  ][il-1];
        ERxr        = data->RadiationEnergyDensity[kl  ][jl  ][il+1];
        rhogas      = data->Vc[RHO][kl  ][jl  ][il  ];
        rhodust     = data->DustToGasMassRatio[kl  ][jl  ][il  ] * rhogas;
        Tdust       = data->DustTemperature[kl  ][jl  ][il  ];
        Tgas        = data->GasTemperature[kl  ][jl  ][il  ];
        Trad        = pow(ER/a_RadiationConstant, 0.25);
        kappadust   = RosselandMeanDustOpacity(Trad   , Tdust  , rhogas  );
        kappagas    = RosselandMeanGasOpacity(Trad    , Tgas   , rhogas   );
        kapparho    = kappadust   * rhodust   + kappagas   * rhogas  ;
        dxl         = 0.5 * (GridCellLength(grid, KDIR, kl,   jl,   il-1) + GridCellLength(grid, KDIR, kl,   jl,   il  ));
        dxr         = 0.5 * (GridCellLength(grid, KDIR, kl,   jl,   il  ) + GridCellLength(grid, KDIR, kl,   jl,   il+1));
        R           = fabs(ERxr - ERxl) / (dxr + dxl);
#if DIMENSIONS > 1
        ERyl        = data->RadiationEnergyDensity[kl  ][jl-1][il  ];
        ERyr        = data->RadiationEnergyDensity[kl  ][jl+1][il  ];
        dyl         = 0.5 * (GridCellLength(grid, KDIR, kl,   jl-1, il  ) + GridCellLength(grid, KDIR, kl,   jl,   il  ));
        dyr         = 0.5 * (GridCellLength(grid, KDIR, kl,   jl,   il  ) + GridCellLength(grid, KDIR, kl,   jl+1, il  ));
        R          += fabs(ERyr - ERyl) / (dyr + dyl);
#endif
#if DIMENSIONS > 2
        ERzl        = data->RadiationEnergyDensity[kl-1][jl  ][il  ];
        ERzr        = data->RadiationEnergyDensity[kl+1][jl  ][il  ];
        dzl         = 0.5 * (GridCellLength(grid, KDIR, kl-1, jl,   il  ) + GridCellLength(grid, KDIR, kl,   jl,   il  ));
        dzr         = 0.5 * (GridCellLength(grid, KDIR, kl,   jl,   il  ) + GridCellLength(grid, KDIR, kl+1, jl,   il  ));
        R          += fabs(ERzr - ERzl) / (dzr + dzl);
#endif
        R          /= kapparho * ER;
        L           = FluxLimiter(R);
        Potential  += L * data->RadiationEnergyDensity[kl][jl][il] / data->Vc[RHO][kl][jl][il];
    }
#endif
    
    return Potential;
}
#endif


