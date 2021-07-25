/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2016 Thomas Henn, Mathias J. Krause,
 *  Marie-Luise Maier
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

/* bifurcation3d.cpp:
 * This example examines a steady particulate flow past a bifurcation. At the inlet,
 * an inflow condition with grad_n u = 0 and rho = 1 is implemented.
 * At both outlets, a Poiseuille profile is imposed on the velocity.
 * After a start time, particles are put into the bifurcation at the
 * inlet and experience a stokes drag force.
 *
 * A publication using the same geometry can be found here:
 * http://link.springer.com/chapter/10.1007/978-3-642-36961-2_5  
 *  *
 */

#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb3D.hh"   // include full template code;
#endif

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
typedef D3Q19<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR> DESCRIPTOR;
#define PARTICLE Particle3D

#define WriteVTK
#define WriteGnuPlot
std::string gnuplotFilename = "gnuplot.dat";

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// set simulation type
// true = resolved particles
// false = subgrid particles
const bool resolved = false;

// Discretization Settings
const int res = 4;                     // resolution of the model
T const charLatticeVelocity = 0.01;     // characteristic lattice velocity

// Time Settings
const T maxPhysT = T( 0.065 );            // max. fluid simulation time in s, SI unit
const T iTwrite = 0.001;                // write out intervall in s

// Domain Settings
T const lengthX = 0.000035*2;               // Cube side length X in m
T const lengthY = 0.000035*2;               // Cube side length Y in m
T const lengthZ = 0.006;                 // Cube side length Z in m

// Fluid Settings
T const physDensity = 1.2;              // Fluid density in kg/(m^3)
T const physViscosity = 1.8E-5;         // Fluid kinematic viscosity in m^2/s

//Particle Settings
T centerX = lengthX*.5;                 // Particle initial postition X
T centerY = lengthY*.5;                 // Particle initial postition Y
T centerZ = lengthZ*.9;                 // Particle initial postition Z
T const sphereDensity = 2500.;          // Particle density in kg/(m^3)
T const sphereDiameter = 0.000035;      // Particle diameter in m

// Characteristic Quantities
T const charPhysLength = sphereDiameter;       // characteristic reference length in m for resolution
T const charPhysVelocity =0.21 ;        // Assumed maximal velocity in m/s        

// Prepare geometry
void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry3D<T>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 2);
  superGeometry.rename(2, 1, 1, 1, 1);

  superGeometry.clean();
  superGeometry.innerClean();

  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter, Dynamics<T, DESCRIPTOR>&
                     designDynamics, SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout, "prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 1, &designDynamics );

  // Material=2 -->bounce back
  sLattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>() );

  clout << "Prepare Lattice ... OK" << std::endl;
  return;
}

//Set Boundary Values for sub-grid particles
void setBoundaryValuesSub( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                        SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout, "setBoundaryValues" );

  if ( iT == 0 ) {
    AnalyticalConst3D<T, T> rho( 1 );
    std::vector<T> velocity( 3, T() );
    AnalyticalConst3D<T, T> uF( velocity );

    sLattice.iniEquilibrium( superGeometry, 1, rho, uF );
    sLattice.iniEquilibrium( superGeometry, 2, rho, uF );

    sLattice.defineRhoU( superGeometry, 1, rho, uF );
    sLattice.defineRhoU( superGeometry, 2, rho, uF );

    // Make the lattice ready for simulation
    sLattice.initialize();
  }
}

//Set Boundary Values for resolved particles
void setBoundaryValuesRes(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                       UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                       SuperGeometry3D<T>& superGeometry)
{
  OstreamManager clout(std::cout, "setBoundaryValues");

  if (iT == 0) {
    AnalyticalConst3D<T, T> zero(0.);
    AnalyticalConst3D<T, T> one(1.);
    sLattice.defineField<descriptors::POROSITY>(superGeometry.getMaterialIndicator({0,1,2}), one);
    // Set initial condition
    AnalyticalConst3D<T, T> ux(0.);
    AnalyticalConst3D<T, T> uy(0.);
    AnalyticalConst3D<T, T> uz(0.);
    AnalyticalConst3D<T, T> rho(1.);
    AnalyticalComposed3D<T, T> u(ux, uy, uz);

    //Initialize all values of distribution functions to their local equilibrium
    sLattice.defineRhoU(superGeometry, 1, rho, u);
    sLattice.iniEquilibrium(superGeometry, 1, rho, u);

    // Make the lattice ready for simulation
    sLattice.initialize();
  }
}

// Outputs simulation results as vtk for sub-grid particles
bool getResultsSub( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, std::size_t iT, T iTwrite,
                 SuperGeometry3D<T>& superGeometry,
                 Timer<double>& timer,
                 SuperParticleSystem3D<T, PARTICLE>& supParticleSystem,
                 T radii, T sphereDensity,
                 SuperParticleSysVtuWriter<T, PARTICLE>& supParticleWriter,
                 PARTICLE<T> &p)
{

  OstreamManager clout( std::cout, "getResults" );
  SuperVTMwriter3D<T> vtmWriter( "settlingSphereSub" );
  SuperVTMwriter3D<T> vtmWriterStartTime( "startingTimeSettlingSphereSub3D" );

  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  vtmWriterStartTime.addFunctor( velocity );
  vtmWriterStartTime.addFunctor( pressure );

  if ( iT == 0 ) {
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
    vtmWriterStartTime.createMasterFile();

    // Print some output of the chosen simulation setup
    clout << "res=" << res <<"; maxTimeSteps="
          << converter.getLatticeTime( maxPhysT ) << "; noOfCuboid="
          << superGeometry.getCuboidGeometry().getNc() << std::endl;
  }

  // Writes the vtk
  if ( iT % converter.getLatticeTime(iTwrite) == 0 ) {
      vtmWriterStartTime.write(iT);
      vtmWriter.write(iT);
  }

  // Writes output on the console
  if (iT% converter.getLatticeTime(iTwrite) == 0 ) {

    // Timer statics
    timer.update( iT );
    timer.printStep();

    // Lattice statistics
    sLattice.getStatistics().print( iT, converter.getPhysTime( iT ) );    
  }

  // Writes output on the console
  if ( (iT% converter.getLatticeTime(iTwrite) == 0 || iT == converter.getLatticeTime( maxPhysT )) ) {

    timer.print( iT );

    // console output number of particles at different material numbers mat
    supParticleSystem.print( {1,2} );
    // console output of escape (E), capture (C) rate for material numbers mat
    supParticleSystem.captureEscapeRate( {2} );
	
    supParticleWriter.write( iT );

    // true as long as certain amount of active particles
    if ( supParticleSystem.globalNumOfActiveParticles() < 1
         && iT > 0.9*converter.getLatticeTime(maxPhysT ) ) {
      return false;
    }
  }
  return true;
}

// Outputs simulation results as vtk and Gnuplot for resolved particles
void getResultsRes(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                SuperGeometry3D<T>& superGeometry, Timer<double>& timer,
                ParticleDynamics3D<T, DESCRIPTOR> particleDynamics,
                SmoothIndicatorF3D<T,T,true> &particle)
{
  OstreamManager clout(std::cout, "getResults");

#ifdef WriteVTK
  SuperVTMwriter3D<T> vtkWriter("sedimentation");
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
  SuperLatticePhysExternalPorosity3D<T, DESCRIPTOR> externalPor(sLattice, converter);
  vtkWriter.addFunctor(velocity);
  vtkWriter.addFunctor(pressure);
  vtkWriter.addFunctor(externalPor);

  if (iT == 0) {
    /// Writes the converter log file
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry(sLattice, superGeometry);
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank3D<T, DESCRIPTOR> rank(sLattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
  }

  if (iT % converter.getLatticeTime(iTwrite) == 0) {
    vtkWriter.write(iT);
  }
#endif


  // Gnuplot constructor (must be static!)
  // for real-time plotting: gplot("name", true) // experimental!
  static Gnuplot<T> gplot( "kinetics" );

  // write pdf at last time step
  if ( iT == int(converter.getLatticeTime( maxPhysT ))-1 ) {
    // writes pdf
    gplot.writePDF();
  }


  /// Writes output on the console
  if (iT % converter.getLatticeTime(iTwrite) == 0) {
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    particleDynamics.print();
    
  //get data to graph
  T physTime = converter.getPhysTime(iT);
  T mass;
  mass = particle.getMass();
  T pos[3], vel[3], acc[3], force[3];
  for (int i=0; i<=2; i++) {
  pos[i] = particle.getPos()[i];
  vel[i] = particle.getVel()[i];
  acc[i] = particle.getAcc()[i];
  force[i] = acc[i] * mass;
  }
  
  clout << "z=" << pos[2];
  clout << "; vel_z=" << vel[2];
  clout << "; mass=" << mass;
  clout << "; acc_z=" << acc[2];
  clout << "; force_z=" << force[2] << endl;
  //rticleDynamics.save("test");
  gplot.setData(physTime, {pos[2], vel[2], acc[2]}, {"z","v_z","a_z"}, "bottom right", {'l','l','l'} );
  gplot.writePNG();
  }
   
    // every (iT%vtkIter) write an png of the plot
    if ( iT%( converter.getLatticeTime( .05 ) ) == 0 ) {
      // writes pngs: input={name of the files (optional), x range for the plot (optional)}
      gplot.writePNG( iT, maxPhysT );
    }
    
    #ifdef WriteGnuPlot
  if (iT % converter.getLatticeTime(iTwrite) == 0) {
    if (singleton::mpi().getRank() == 0) {

      ofstream myfile;
      myfile.open (gnuplotFilename.c_str(), ios::app);
      myfile
          << converter.getPhysTime(iT) << " "
          << std::setprecision(9)
          << particle.getPos()[2] << " "
          << particle.getVel()[2] << " "
          << particle.getMass() << " "
          << particle.getAcc()[2] << " "
          << particle.getAcc()[2] * particle.getMass() << endl;
      myfile.close();
    }
  }
#endif
}

int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===

  olbInit( &argc, &argv );
  if (resolved){
      singleton::directories().setOutputDir( "./tmpRes/" ); // only for resolved particles
  } else {
      singleton::directories().setOutputDir( "./tmpSub/" ); // only for sub-grid particles
  }
  OstreamManager clout( std::cout, "main" );

  UnitConverterFromResolutionAndLatticeVelocity<T,DESCRIPTOR> const converter(
    int {res},                     // resolution
    (T)   charLatticeVelocity,     // charLatticeVelocity
    (T)   charPhysLength,          // charPhysLength: reference length of simulation geometry
    (T)   charPhysVelocity,        // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   physViscosity,           // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   physDensity              // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("settlingSphere3d");

  // === 2nd Step: Prepare Geometry ===

  // Instantiation of a cuboidGeometry
  // values can be changed in the beginning of the file
  std::vector<T> extend(3, T());
  extend[0] = lengthX;
  extend[1] = lengthY;
  extend[2] = lengthZ;
  std::vector<T> origin(3, T());
  IndicatorCuboid3D<T> cuboid(extend, origin);

#ifdef PARALLEL_MODE_MPI
  CuboidGeometry3D<T> cuboidGeometry(cuboid, converter.getConversionFactorLength(), singleton::mpi().getSize());
#else
  CuboidGeometry3D<T> cuboidGeometry(cuboid, converter.getConversionFactorLength(), 7);
#endif
  cuboidGeometry.print();

  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);
  SuperGeometry3D<T> superGeometry(cuboidGeometry, loadBalancer, 2);
  prepareGeometry(converter, superGeometry);

  // === 3rd Step: Prepare Lattice ===

  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

  BGKdynamics<T, DESCRIPTOR> designDynamicsSub( converter.getLatticeRelaxationFrequency(),instances::getBulkMomenta<T, DESCRIPTOR>() );
  PorousParticleBGKdynamics<T, DESCRIPTOR, false> designDynamicsRes(converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>());

  if (resolved){
    prepareLattice(sLattice, converter, designDynamicsRes, superGeometry); // only for resolved particles
  } else  {
    prepareLattice( sLattice, converter, designDynamicsSub, superGeometry ); // only for sub-grid particles
  }

  // === 3.1 Step: Particles ===
  clout << "Prepare Particles ..." << std::endl;

  // SuperParticleSystems3D
  SuperParticleSystem3D<T, PARTICLE> supParticleSystem( superGeometry );
  // define which properties are to be written in output data
  SuperParticleSysVtuWriter<T, PARTICLE> supParticleWriter( supParticleSystem,
      "particles", SuperParticleSysVtuWriter<T, PARTICLE>::particleProperties::
      velocity
      | SuperParticleSysVtuWriter<T, PARTICLE>::particleProperties::mass
      | SuperParticleSysVtuWriter<T, PARTICLE>::particleProperties::radius
      | SuperParticleSysVtuWriter<T, PARTICLE>::particleProperties::active );

  SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR> getVel( sLattice, converter );
  
  // material numbers where particles should be reflected
  std::set<int> boundMaterial = {2};
  auto materialBoundary = make_shared
                          < MaterialBoundary3D<T, PARTICLE>
                          > ( superGeometry, boundMaterial );

  //Forces acting on the particle (Drag,Gravity and bouyancy)
  auto stokesDragForce = make_shared
                         < StokesDragForce3D<T, PARTICLE, DESCRIPTOR>
                         > ( getVel, converter );


  auto WeightForce = make_shared 
                        < WeightForce3D< T, PARTICLE >
                        > (std::vector<T>({0,0,-1}),9.81);
 
  auto buoyancyForce = make_shared
  								< BuoyancyForce3D< T, PARTICLE, DESCRIPTOR >
  								> (converter, std::vector<T>({0,0,-1}), 9.81);

  // Adding forces and boundary to Particle System
  supParticleSystem.addForce( WeightForce );
  supParticleSystem.addForce( stokesDragForce );
  supParticleSystem.addForce( buoyancyForce );

  supParticleSystem.addBoundary( materialBoundary );
  supParticleSystem.setOverlap( 2. * converter.getConversionFactorLength() );

  //add 1 particle with indicator 1 and fixed location
  PARTICLE<T> p({centerX,centerY,centerZ},4. / 3. * M_PI * std::pow( (sphereDiameter/2), 3 ) * sphereDensity, (sphereDiameter/2),1);
  supParticleSystem.addParticle(p);

  clout << "Prepare Particles ... OK" << std::endl;
  
  // === 4th Step: Main Loop with Timer ===

  Timer<double> timer( converter.getLatticeTime( maxPhysT ),
                            superGeometry.getStatistics().getNvoxel() );

  timer.start();

  std::size_t iT = 0;

  if (!resolved){ // only for sub-grid particles
    sLattice.communicate();
    
    for ( ; iT <= converter.getLatticeTime( maxPhysT ); ++iT ) {
      setBoundaryValuesSub( sLattice, converter, iT, superGeometry );
      //sLattice.collideAndStream();
      supParticleSystem.simulate( converter.getConversionFactorTime() );
    
    if (iT % converter.getLatticeTime(iTwrite) == 0) {
      supParticleSystem.getOutput("pos_output", iT, converter.getConversionFactorTime(),1);
      supParticleSystem.getOutput("vel_output", iT, converter.getConversionFactorTime(),2);
      supParticleSystem.getOutput("force_output", iT, converter.getConversionFactorTime(),16);
    }
      
      if ( !getResultsSub( sLattice, converter, iT,
                        iTwrite, superGeometry, timer,
                        supParticleSystem, sphereDiameter, sphereDensity,
                        supParticleWriter, p) ) {
        break;
      }
    }
  }

  // Create Particle Dynamics
  ParticleDynamics3D<T, DESCRIPTOR> particleDynamics(sLattice, converter, superGeometry, lengthX, lengthY, lengthZ, {.0, .0, -9.81 * (1. - physDensity / sphereDensity)});

  // Create Sphere Indicator
  //T epsilon = 0.5*converter.getConversionFactorLength();
  T epsilon = 0;
  SmoothIndicatorSphere3D<T, T, true> particleIndicator ({centerX,centerY,centerZ}, 0.5*sphereDiameter, epsilon, sphereDensity, {0.,0.,0.});

  SuperField3D<T,DESCRIPTOR,POROSITY> superExtPorosity(superGeometry, sLattice, sLattice.getOverlap());
  SuperField3D<T,DESCRIPTOR,VELOCITY_NUMERATOR> superExtNumerator(superGeometry, sLattice, sLattice.getOverlap());
  SuperField3D<T,DESCRIPTOR,VELOCITY_DENOMINATOR> superExtDenominator(superGeometry, sLattice, sLattice.getOverlap());
  particleDynamics.addParticle( particleIndicator );

  /// === 5th Step: Definition of Initial and Boundary Conditions ===
  if (resolved){ // only for resolved particles
    particleDynamics.print();

    setBoundaryValuesRes(sLattice, converter, 0, superGeometry);

    clout << "MaxIT: " << converter.getLatticeTime(maxPhysT) << std::endl;
    for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT)+10; ++iT) {
      particleDynamics.simulateTimestep("verlet");
      getResultsRes(sLattice, converter, iT, superGeometry, timer, particleDynamics,particleIndicator);
      sLattice.collideAndStream();
      superExtPorosity.communicate();
      superExtNumerator.communicate();
      superExtDenominator.communicate();
    }
  }


  timer.stop();
  timer.printSummary();
}
