/*    Copyright (c) 2010-2017, Delft University of Technology
*    All rigths reserved
*
*    This file is part of the Tudat. Redistribution and use in source and
*    binary forms, with or without modification, are permitted exclusively
*    under the terms of the Modified BSD license. You should have received
*    a copy of the license with this file. If not, please or visit:
*    http://tudat.tudelft.nl/LICENSE.
*/

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>

#include <spaceProjectAssignment2/inputOutputDirectories.h>

namespace tudat
{


///                                                                 ////
/// Barebones class for aerodynamic guidance of the space shuttle   ////
///                                                                 ////
///
class SpaceShuttleAerodynamicGuidance: public aerodynamics::AerodynamicGuidance
{
public:

   //! Constructor
   SpaceShuttleAerodynamicGuidance( const simulation_setup::NamedBodyMap& bodyMap )
   {
       vehicleFlightConditions_ = bodyMap.at( "STS" )->getFlightConditions( );
       coefficientInterface_ = bodyMap.at("STS")-> getAerodynamicCoefficientInterface( );
       ////                                                                                ////
       ////    Retrieve relevant objects from environment and set as member variables      ////
       ////                                                                                ////
    }
       double airspeed    = vehicleFlightConditions_->getCurrentAirspeed( );
       double density     = vehicleFlightConditions_->getCurrentDensity( ) ;
       double radius      = vehicleFlightConditions_->getCurrentAltitude( ) + 6378137.0;
       double vehiclemass = 165000 * 0.45;
       double gravitydown = 9.81*(6378137.0/radius)*(6378137.0/radius);
       double gamma       = vehicleFlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::flight_path_angle );
       double latitude    = vehicleFlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::latitude_angle );
       double heading     = vehicleFlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( reference_frames::heading_angle );
       const double rotationEarth = ( 2.0 * mathematical_constants::PI ) / 86164.0916;
       const double refereceArea = 2690.0 * 0.3048 * 0.3048;
       double liftforce;

       // call objects from environment like angles and altitude here!
       // Have a look at aerodynamic guidance and vehicleFlightconditions_ header files, not cpp file!

       // Cl and S for liftforce
       // gravity downward component

   //! The aerodynamic angles are to be computed here
   void updateGuidance( const double time )
   {
       // Define input to aerodynamic coefficients: take care of order of input (this depends on how the coefficients are created)!
       std::vector< double > currentAerodynamicCoefficientsInput_;
       currentAerodynamicCoefficientsInput_.push_back( currentAngleOfAttack_ );
       currentAerodynamicCoefficientsInput_.push_back( currentAngleOfSideslip_ );
       currentAerodynamicCoefficientsInput_.push_back( vehicleFlightConditions_->getCurrentMachNumber( ) );

       // Update and retrieve current aerodynamic coefficients
       coefficientInterface_->updateCurrentCoefficients( currentAerodynamicCoefficientsInput_ );
       Eigen::Vector3d currentAerodynamicCoefficients = coefficientInterface_->getCurrentForceCoefficients( );

       liftforce   = 0.5 * density * refereceArea * airspeed * airspeed * currentAerodynamicCoefficients(2);

       // functions for angle of attack, angle of sideslip and current bank angle
       if ( vehicleFlightConditions_->getCurrentMachNumber( ) < 6.0 ) {
           currentAngleOfAttack_ = 10.0 * mathematical_constants::PI / 180.0;
       }

       else if (vehicleFlightConditions_->getCurrentMachNumber( ) > 12.0 || vehicleFlightConditions_->getCurrentAltitude( ) > 150.0E3 ) {
           currentAngleOfAttack_ = 40.0 * mathematical_constants::PI / 180.0;
       }

       else {
           currentAngleOfAttack_ = ( 5.0 * vehicleFlightConditions_->getCurrentMachNumber( ) - 20.0 ) * mathematical_constants::PI / 180.0;
       }

       currentAngleOfSideslip_ = 0.0;

// define terms of the RHS of the change in flight path angle equation
        double RHS1 = liftforce / vehiclemass;
        double RHS2 = - ( gravitydown + pow(airspeed,2) / radius ) * cos( gamma ) +  2 * rotationEarth * airspeed * cos( latitude ) * cos( heading ) +
                pow(rotationEarth,2) * radius * cos( latitude ) * ( cos( latitude ) * cos (gamma ) + sin( gamma ) * sin( latitude ) * cos ( heading) );

        if (fabs(RHS1) >= fabs(RHS2))
        { currentBankAngle_ = acos(-RHS2/RHS1);}
        else if (RHS2 < 0) {
        currentBankAngle_ = 0.0;
        }
        else {
            currentBankAngle_ = mathematical_constants::PI; // or 90 deg? Because it doesn't make sense to have a bank angle bigger than 90 deg.
        }

    }




private:

   boost::shared_ptr< aerodynamics::FlightConditions > vehicleFlightConditions_;
   boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > coefficientInterface_;
   boost::shared_ptr< simulation_setup::Body > Body_;

};

}

//! Execute propagation of orbits of STS during entry.
int main( )
{
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


   using namespace tudat::ephemerides;
   using namespace tudat::interpolators;
   using namespace tudat::numerical_integrators;
   using namespace tudat::spice_interface;
   using namespace tudat::simulation_setup;
   using namespace tudat::basic_astrodynamics;
   using namespace tudat::orbital_element_conversions;
   using namespace tudat::propagators;
   using namespace tudat::aerodynamics;
   using namespace tudat::basic_mathematics;
   using namespace tudat::input_output;
   using namespace tudat;

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   // Load Spice kernels.
   spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
   spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );

   // Set simulation start epoch.
   const double simulationStartEpoch = 0.0;
   const double simulationEndEpoch = 5000.0;

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   // Define simulation body settings.
   std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
           getDefaultBodySettings( { "Earth" } );
   bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< simulation_setup::ConstantEphemerisSettings >(
               Eigen::Vector6d::Zero( ), "SSB", "J2000" );
   bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
   bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< AtmosphereSettings >( nrlmsise00 );

   // Create Earth object
   simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

   // Create vehicle objects.
   bodyMap[ "STS" ] = boost::make_shared< simulation_setup::Body >( );

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


   // Create vehicle  coefficients
   std::map< int, std::string > forceCoefficientFiles;
   forceCoefficientFiles[ 0 ] =
            tudat_applications::getStsInputPath( ) + "STS_CD.dat";
   forceCoefficientFiles[ 2 ] =
            tudat_applications::getStsInputPath( ) + "STS_CL.dat";

   boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
           simulation_setup::readTabulatedAerodynamicCoefficientsFromFiles(
               forceCoefficientFiles, 2690.0 * 0.3048 * 0.3048,
               boost::assign::list_of( aerodynamics::angle_of_attack_dependent )(  aerodynamics::mach_number_dependent ),
               true, true );

   bodyMap[ "STS" ]->setAerodynamicCoefficientInterface(
               createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "STS" ) );
   bodyMap[ "STS" ]->setConstantBodyMass( 165000 * 0.45 );

   // Finalize body creation.
   setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   // Define propagator settings variables.
   SelectedAccelerationMap accelerationMap;
   std::vector< std::string > bodiesToPropagate;
   std::vector< std::string > centralBodies;

   // Define propagation settings.
   std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSTS;

   accelerationsOfSTS[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                    basic_astrodynamics::central_gravity ) );

   accelerationMap[  "STS" ] = accelerationsOfSTS;
   bodiesToPropagate.push_back( "STS" );
   centralBodies.push_back( "SSB" );

   // set guidance
   boost::shared_ptr< aerodynamics::AerodynamicGuidance > aerodynamicGuidance =
           boost::make_shared< SpaceShuttleAerodynamicGuidance>( bodyMap); //Create user-defined guidance object here
   setGuidanceAnglesFunctions( aerodynamicGuidance, bodyMap);

   basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
               bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

   ///                                                                         ///
   /// Create acceleration models and set guidance (in that order!!) here      ///
   ///                                                                         ///

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////             DEFINE INITIAL STATE                   ////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   // Set spherical elements for STS.
   Eigen::Vector6d stsSphericalEntryState;
   stsSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
           spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
   stsSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) = 0.3;
   stsSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) = 1.2;
   stsSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 7.45E3;
   stsSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
           -1.2 * mathematical_constants::PI / 180.0;
   stsSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = 0.6;

   // Convert sts state from spherical elements to Cartesian elements.
   Eigen::Vector6d systemInitialState = transformStateToGlobalFrame(
               convertSphericalOrbitalToCartesianState(
                               stsSphericalEntryState ),
               simulationStartEpoch, bodyMap.at( "Earth" )->getRotationalEphemeris( ) );

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////     CREATE INTEGRATION/PROPAGATION SETTINGS; PROPAGATE ORBIT; SAVE OUTPUT      ////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;

   dependentVariablesList.push_back(
                   boost::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "STS", "Earth" ) );
//   dependentVariablesList.push_back(
//                   boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >( "STS", 0 ) ); // latitude
//   dependentVariablesList.push_back(
//                   boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >( "STS", 1 ) ); // longitude
//   dependentVariablesList.push_back(
//                   boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >( "STS", 6 ) ); // bank angle
//   dependentVariablesList.push_back(
//                   boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >( "STS", 3 ) ); // flight path angle
//   dependentVariablesList.push_back(
//                   boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >( "STS", 4 ) ); // angle of attack
   dependentVariablesList.push_back(
                   boost::make_shared< SingleDependentVariableSaveSettings >( mach_number_dependent_variable, "STS", "Earth" ) );
   dependentVariablesList.push_back(
                   boost::make_shared< SingleDependentVariableSaveSettings >( airspeed_dependent_variable, "STS", "Earth" ) );
   dependentVariablesList.push_back(
                   boost::make_shared< SingleDependentVariableSaveSettings >( aerodynamic_moment_coefficients_dependent_variable, "STS" ) );


 // Create object with list of dependent variables
 //  boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
 //          boost::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

   boost::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
           boost::make_shared< SingleDependentVariableSaveSettings >(
               altitude_dependent_variable, "STS", "Earth" );

   boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
           boost::make_shared< PropagationDependentVariableTerminationSettings >(
               terminationDependentVariable, 10.0E3, true );

   boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings, cowell,
              boost::make_shared< DependentVariableSaveSettings >( dependentVariablesList ) );

   // create prop and integrator settings
//   boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings =
//           boost::make_shared< TranslationalStatePropagatorSettings< double > >
//           ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch);

   const double fixedStepSize = 30.0;
   boost::shared_ptr< IntegratorSettings< > > integratorSettings =
           boost::make_shared< IntegratorSettings< > >
           ( rungeKutta4, simulationStartEpoch, fixedStepSize );

   // integrate and propagate


   // Create simulation object and propagate dynamics.
   SingleArcDynamicsSimulator< > dynamicsSimulator(
               bodyMap, integratorSettings, propagatorSettings, true, false, false );
   std::map< double, Eigen::VectorXd > cartesianIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
   std::map< double, Eigen::VectorXd > keplerianIntegrationResult;

   // Compute map of Kepler elements
   Eigen::Vector6d currentCartesianState;
   for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = cartesianIntegrationResult.begin( );
        stateIterator != cartesianIntegrationResult.end( ); stateIterator++ )
   {
       // Retrieve current Cartesian state (convert to Moon-centered frame if needed)
       currentCartesianState = stateIterator->second;

       keplerianIntegrationResult[ stateIterator->first ] =
               convertCartesianToKeplerianElements(
                   currentCartesianState, bodyMap[ "Earth" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
   }

   // Save data

   input_output::writeDataMapToTextFile( cartesianIntegrationResult,
                                         "stsOrbitCart.dat",
                                         tudat_applications::getStsOutputPath( ),
                                         "",
                                         std::numeric_limits< double >::digits10,
                                         std::numeric_limits< double >::digits10,
                                         "," );



   ///                                                                                          ///
   ///    Create integration and propagation settings here, then propagate and save dynamics    ///
   ///                                                                                          ///

   std::string outputDirectory = tudat_applications::getStsOutputPath( );

   // Final statement.
   // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
   return EXIT_SUCCESS;
}
