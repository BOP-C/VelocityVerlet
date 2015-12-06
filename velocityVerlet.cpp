/* 	
	Oscillatory-ballistic motion regularities of a gravitational pendulum
	Submitted to Nonlinear Dynamics

	Authors: Sebastian Micluta-Campeanu and Tiberius O. Cheche
	Faculty of Physics, University of Bucharest, Bucharest 077125, Romania
	Corresponding author: Tiberius O. Cheche; e-mail: cheche@gate.sinica.edu.tw
	Date: December 4, 2015

	Description:
	The program simulates the trajectory and the number of 'collisions' (the event when the pendulum string 
	is suddenly tensed during its motion) with the following initial conditions: 
	the gravitational pendulum is launched from the equilibrium vertical position 
	with velocity perpendicular to the string and initial angular oscillation amplitude between
	90deg and 180deg(characterized by gamma = omega0^2 * l / g,
	where omega0 is the initial angular amplitude and l is the length of the string
	and g is the gravitational acceleration; gamma takes values between 2 and 5).
	We consider that the radial component of velocity vanishes at each collision  event.

	Input:
	1. without command line arguments
		gamma (from stdin), 
		mass, length of the wire,
		simulation time, time step, number of steps between file writes
		(in S.I.) (in a file named "input.dat")

	2. with command line arguments
	If the program is with command line arguments, they should respect the following rules:
		If the program is run with 2 arguments, the first argument will be a real number 
			between 2 and 5 and the second argument will be a pozitive integer greater than 0.
			The first argument represents the gamma parameter defined in the description above.
			The second argument represents the maximum number of collisions.
		If the program is run with 3 arguments, the first 2 arguments will respect the rules above 
			and the third argument will represent a file name for the phase space portrait.
	
	Output:
	If the program is run without arguments, the output mode is selected
		at runtime by the user and the output will be stored in a file named "trajectory.dat".
	If the program is run with 2 arguments, the program switches to output mode 0(nothing in "trajectory.dat") and 
		the simulation results will be appended in a file named "results.dat". 
		The results will have the following format:
		timeStep gamma noCollisions remainingEnergy AngleAtTheLastCollision simulationTime
	If the program is run with 3 arguments, in addition to the above-mentioned behaviour,
		the program will switch to output mode 2(phase space portrait) and the result will be 
		storred in a file with the name given by the third argument.

	For detailed information at each collision switch the verbose flag to true.
	The logging information will be stored in a file named "log.dat".

*/ 

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "point.h"
#include "common.h"

// constants
const long double PI = 3.14159265358979323846264338327950;
const long double g = 1;

// files
std::ifstream in("input.dat");
std::ofstream out("trajectory.dat");
std::ofstream logFile("log.dat");
std::ofstream logResults("results.dat", std::ofstream::out | std::ofstream::app);

// flags
bool verbose = false;

int main(int argc, char *argv[])
{
	long double angle, omega0, gamma, mass, l, dt, dtheta, simTime, t;
	int fileWrite, outputMode, maxNumberOfCollisions, signChange;
	char* fileName;

	out.precision(std::numeric_limits<long double>::digits10);
	std::cout.precision(std::numeric_limits<long double>::digits10);
	logFile.precision(std::numeric_limits<long double>::digits10);
	logResults.precision(std::numeric_limits<long double>::digits10);
	
	if (argc >= 3)	// if the input is not given as an argument, use std::cin
	{
		outputMode = 0;
		verbose = false;
		gamma = atof(argv[1]);
		maxNumberOfCollisions = atoi(argv[2]);
		if (gamma <= 2 || gamma >= 5 || maxNumberOfCollisions <=0 )
		{
			std::cout << "Arguments:\nvelocityVerlet [gamma] [maxNumberOfCollisions] {filename}\ngamma must be between 2 and 5 and maxNumberOfCollisions must be >=0";
			return 1;
		}
		if (argc == 4)
		{
			outputMode = 2;
			fileName = argv[3];
			out.close();
			out.open(fileName, std::ofstream::out);
		}
		angle = 1.5 * PI;
	}
	else
	{
		try {
			std::cout << "Possible output modes:\n";
			std::cout << "\
0. no output\n\
1. trajectory equation\n\
2. phase space\n\
3. E(t)\n\
4. x(t)\n\
5. y(t)\n\
6. theta(t)\n\
7. dtheta / dt (t)\n\
8. E(nrCollisions)\n";
			std::cout << "Output mode(0 - 17): ";
			std::cin >> outputMode;
			if (outputMode < 0 || outputMode > 8 || !std::cin.good())
				throw outputMode;
		}
		catch (int choice)
		{
			std::cerr << "Choice " << choice << " is unavailable. Try again.\n";
			return 1;
		}
		angle = 0;
		// transform angle to radians & adjust to theta being measured from Ox
		angle = (angle * PI) / 180 + 1.5 * PI;
		try {
			std::cout << "Gamma(between 2 and 5): ";
			std::cin >> gamma;
			if (gamma <= 2 || gamma >= 5 || !std::cin.good())
				throw gamma;
			std::cout << "Maximum number of collisions: ";
			std::cin >> maxNumberOfCollisions;
			if (maxNumberOfCollisions <= 0 || !std::cin.good())
				throw maxNumberOfCollisions;
		}
		catch (long double invalidGamma) {
			std::cout << invalidGamma << " is not a valid value for gamma. Try again.\n";
			return 1;
		}
		catch (int invalidMax) {
			std::cout << invalidMax << " is not a valid value for maximum number of collisions. Try again.\n";
			return 1;
		}
	}
	
	// mass, length of the wire, 
	// simulation time, time step, number of steps between file writes
	in >> mass >> l >> simTime >> dt >> fileWrite;
	omega0 = sqrt(gamma * g / l);

	if (argc >= 3)
	{
		logResults << dt << '\t' << gamma << '\t';
	}

	signChange = 0;
	dtheta = 0;
	t = 0;
	point *p = new point(l, angle, omega0, mass, dt, dtheta, t, outputMode);

	for (long long i = 0; t < simTime && p->numberOfCollisions < maxNumberOfCollisions; i++)
	{
		long double old = dtheta;
		t += p->dt;
		try
		{
			p->move(dtheta, angle, t);

			// file size optimization
			if (i % fileWrite == 0)
				p->write(dtheta, p->dt, t);
		}
		catch (const char* message)
		{
			std::cout << message << "at t = " << t << '\n';
			if (argc >= 3)
				logResults << "err:\t" << message << '\t' << p->numberOfCollisions << '\t' << p->E() << '\t' << atan2(p->lastXc, -p->lastYc) + 1.5 * PI << '\t' << t << '\n';
			delete p;
			return 0;
		}
		catch (...) {
			std::cerr << "Exception occurred";
		}

		if (p->E() < mass * g * l)	// stop simulation if there isn't enough energy for another collision
		{
			if (old * dtheta < 0)
				signChange++;
			if (signChange == 3)
			{
				//make sure that the ellipse in phase space is complete
				break;
			}
		}
	}

	if (argc >= 3)
		logResults << p->numberOfCollisions << '\t' << p->E() << '\t' << atan2(p->lastXc, -p->lastYc) + 1.5 * PI << '\t' << t << '\n';
	else
		std::cout << "Number of collisions: " << p->numberOfCollisions << '\n' << "Energy: " << p->E() << '\n';

	delete p;
	return 0;
}
