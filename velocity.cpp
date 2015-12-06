/*
	Oscillatory-ballistic motion regularities of a gravitational pendulum
	Submitted to Nonlinear Dynamics

	Authors: Sebastian Micluta-Campeanu and Tiberius O. Cheche
	Faculty of Physics, University of Bucharest, Bucharest 077125, Romania
	Corresponding author: Tiberius O. Cheche; e-mail: cheche@gate.sinica.edu.tw
	Date: December 4, 2015

	Description: Velocity class implementation.
	At the collision point the velocity vector can have multiple 
	orientations. In order to get the tangential component
	one must know this orientation.
	Those orientations can be organized using the values of theta and beta,
	where theta is the current angle and betta is the angle between the
	string and the velocity vector.
	We use A-D to denote in which quadrant is theta and 1-8 for the different 
	intervals in which beta takes values.
*/

#include <iostream>

#include "velocity.h"
#include "common.h"

void velocity::update(acceleration a[2], long double dt)
{
	x += (a[0].x + a[1].x) * dt / 2;
	y += (a[0].y + a[1].y) * dt / 2;
}

void velocity::getTangentialComponent(long double theta)
{
	long double beta, v_tg;

	beta = atan2(y, x);
	if (beta < 0)
		beta += 2 * PI;
	if (verbose) {
		logFile << "beta: " << beta << '\n' << "theta: " << theta << '\n';
	}
																	  // the first quadrant
	if (theta >= 0 && theta < PI / 2)
	{
		if (verbose)
		{	
			logFile << "first: ";
		}
		if ((beta - 1.5 * PI) > theta)	// case A8
		{
			v_tg = -mod() * sin(beta - theta);
			x = -v_tg * sin(theta);
			y = -v_tg * cos(theta);
			if (verbose)
			{
				logFile << "A8\n";
			}
			return;
		}
	}
	// the second quadrant
	if (theta >= PI / 2 && theta < PI)
	{
		if (verbose)
		{
			logFile << "second: ";
		}								// only B5 is possible
		if (beta > theta && (beta - PI / 2) < theta)	// cases B4 & B5
		{
			v_tg = mod() * sin(beta - theta);
			x = -v_tg * sin(theta);
			y = v_tg * cos(theta);
			if (verbose)
			{
				logFile << "B4 & B5\n";
			}
			return;
		}
	}
	// the third quadrant
	if (theta >= PI && theta < 1.5 * PI)
	{
		if (verbose)
		{
			logFile << "third: ";
		}
		if (beta > theta && (beta - PI / 2) < theta)	// cases C6 & C7
		{
			v_tg = mod() * sin(beta - theta);
			x = -v_tg * sin(theta);
			y = v_tg * cos(theta);
			if (verbose)
			{
				logFile << "C6 & C7\n";
			}
			return;
		}
		if (beta < theta)	// case C5
		{
			v_tg = -mod() * sin(beta - theta);
			x = -v_tg * sin(theta);
			y = v_tg * cos(theta);
			if (verbose)
			{
				logFile << "C5\n";
			}
			return;
		}
	}
	// the fourth quadrant
	if (theta >= 1.5 * PI && theta <= 2 * PI)
	{
		if (verbose)
		{
			logFile << "fourth: ";
		}
		if (beta > (theta - PI / 2) && beta < theta)	// cases D6 & D7
		{
			v_tg = -mod() * sin(beta - theta);
			x = v_tg * sin(theta);
			y = -v_tg * cos(theta);
			if (verbose)
			{
				logFile << "D6 & D7\n";
			}
			return;
		}
	}

	std::cerr << "Precision limit reached for the given parameters\n";
	throw "Precision limit reached for the given parameters ";
}