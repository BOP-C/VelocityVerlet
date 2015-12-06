/*
	Oscillatory-ballistic motion regularities of a gravitational pendulum
	Submitted to Nonlinear Dynamics

	Authors: Sebastian Micluta-Campeanu and Tiberius O. Cheche
	Faculty of Physics, University of Bucharest, Bucharest 077125, Romania
	Corresponding author: Tiberius O. Cheche; e-mail: cheche@gate.sinica.edu.tw
	Date: December 4, 2015

	Description: Point class implementation.
*/

#include <iostream>

#include "point.h"
#include "common.h"

// constructor
point::point(long double l, long double angle, long double omega, long double m, long double Dt, long double &dtheta, long double &t, int output)
	:
	mass(m),
	T(0),
	dt(Dt),
	collision(false),
	lastXc(0),
	lastYc(0),
	currentOutputMode(output),
	defaultOutputMode(output),
	numberOfCollisions(0)
{
	r = l;
	theta = angle;
	omega0 = omega;
	v.x = -omega0 * l * sin(angle);
	v.y = -omega0 * l * cos(angle);

	// Initialization for Verlet method

	updateCartesian();
	computeAcceleration(dtheta, angle, t);

	if (currentOutputMode != 2)	// dtheta can not be written at this point because it was not computed
		write(0, dt, 0);
}

// current total energy
long double point::E() const
{
	return mass * g * (y + r) + mass * v.mod() * v.mod() / 2;
}

void point::write(long double dtheta, long double dt, long double t) const
{
	switch (currentOutputMode)
	{
	case 0:
		// no output
		break;
	case 1:
		// trajectory equation
		out << x << ' ' << y << '\n';
		if (collision)
			out << "\n\n\n";
		break;
	case 2:
		// phase space
		out << theta << ' ' << dtheta / dt << '\n';
		break;
	case 3:
		// E(t) (using real velocity)
		out << t << ' ' << E() << '\n';
		break;
	case 4:
		// x(t)
		out << t << ' ' << x << '\n';
		break;
	case 5:
		// y(t)
		out << t << ' ' << y << '\n';
		break;
	case 6:
		// theta(t)
		out << t << ' ' << theta << '\n';
		break;
	case 7:
		// dtheta / dt (t)
		out << t << ' ' << dtheta / dt << '\n';
		break;
	case 8:
		// E(nrCollisions)
		out << numberOfCollisions << ' ' << E() << '\n';
		break;
	default:
		throw "Invalid output mode!\n";
		break;
	}
}

void point::computeAcceleration(long double &dtheta, long double &alpha, long double &t)
{
	long double oldT = T;
	long double delta = v.mod() * dt;					// distance since last step
	long double epsilon = r - sqrt(x * x + y * y);		// distance from the point to the circle of radius r ( r - sqrt(x ^ 2 + y ^ 2) )
	long double oldTheta = theta;
	a[0] = a[1];
	updateTheta();

	T = mass * g * (-3 * sin(theta) + 2 * sin(alpha) + r / g * omega0 * omega0);

	if (T < 0)
	{
		T = 0;
		if (currentOutputMode == 2 || currentOutputMode == 4 || currentOutputMode == 10 || currentOutputMode == 11)	// interrupt output when the wire is not stretched
			currentOutputMode = 0;
	}

	if (epsilon > delta && oldT == 0 && t)
	{
		T = 0;
	}

	if (oldT && T == 0)
	{
		// log event
		if (verbose) {
			logFile << "The wire is no longer stretched at t = " << t << '\n' << "Coordinates:\nx: " << x << "\ny: " << y << "\nEnergy: " << E() << "\n\n";
		}
		if (currentOutputMode == 1)
			out << "\n\n\n";
	}

	// the radial component of the velocity vanishes when the wire is stretched again
	if (oldT == 0 && T > 0 && t)
	{
		// update the number of collisions
		numberOfCollisions++;
		collision = true;
		// log event
		if (verbose) {
			logFile << "\nCollision number: " << numberOfCollisions << "\n\n" << "The wire is stretched again at t = " << t << ".\nLength error: " << epsilon << ".\n";
		}
		v.getTangentialComponent(theta);
		
		// recalculate thread tension
		reInitialize(dtheta, alpha, t);

		lastXc = x;
		lastYc = y;
		if (currentOutputMode == 0)
			currentOutputMode = defaultOutputMode;
	}

	dtheta = theta - oldTheta;
	// correction for avoiding anomalies in phase space
	if (dtheta > 2 * PI || dtheta + 1e-3 > 2 * PI)
		dtheta -= 2 * PI;
	if (-dtheta > 2 * PI || -dtheta + 1e-3 > 2 * PI)
		dtheta += 2 * PI;

	if (sqrt(x * x + y * y) - r > r / 10)
		throw "The thread was broken ";

	a[1].x = (-T * cos(theta)) / mass;
	a[1].y = -T * sin(theta) / mass - g;
}

// Move the point using velocity Verlet 
void point::move(long double &dtheta, long double &angle, long double &t)
{
	update(a, v, dt);
	computeAcceleration(dtheta, angle, t);
	if (!collision)
		v.update(a, dt);
	else
		collision = false;
}

// Re-initialize velocity Verlet
void point::reInitialize(long double &dtheta, long double &alpha, long double &t)
{
	// update omega0 according to the current velocity & alpha to the current angle 
	omega0 = v.mod() / r; 
	alpha = theta;
	
	// Recalculate the thread tension because omega0 was modified
	T = mass * g * (-3 * sin(theta) + 2 * sin(alpha) + r / g * omega0 * omega0);

	write(dtheta, dt, t);
	if (verbose) {
		logFile << "Coordinates: \nx: " << x << "\ny: " << y << '\n' << "Energy after movement: " << E() << "\n\n";
	}
}
