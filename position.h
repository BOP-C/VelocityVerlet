/*
	Oscillatory-ballistic motion regularities of a gravitational pendulum
	Submitted to Nonlinear Dynamics

	Authors: Sebastian Micluta-Campeanu and Tiberius O. Cheche
	Faculty of Physics, University of Bucharest, Bucharest 077125, Romania
	Corresponding author: Tiberius O. Cheche; e-mail: cheche@gate.sinica.edu.tw
	Date: December 4, 2015

	Description: Position class header. Provides coordinates and methods 
	for modyfinig them.
*/

#pragma once

#include "acceleration.h"
#include "velocity.h"

class position
{
protected:
	long double x, y;												// Carthesian coordinates
	long double r, theta;											// length of the string and angle(measured from Ox)
	long double omega0;												// Initial angular velocity
	void update(acceleration a[2], velocity v, long double dt);		// velocity Verlet
	void updateCartesian();
	void updateTheta();
public:
	position() { x = y = r = theta = 0; };
};