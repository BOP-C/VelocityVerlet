/*
	Oscillatory-ballistic motion regularities of a gravitational pendulum
	Submitted to Nonlinear Dynamics

	Authors: Sebastian Micluta-Campeanu and Tiberius O. Cheche
	Faculty of Physics, University of Bucharest, Bucharest 077125, Romania
	Corresponding author: Tiberius O. Cheche; e-mail: cheche@gate.sinica.edu.tw
	Date: December 4, 2015

	Description: Velocity class header.
*/

#pragma once

#include <cmath>
#include "acceleration.h"

class velocity
{
public:
	velocity() { x = y = 0; };
	void update(acceleration a[2], long double dt);
	void getTangentialComponent(long double theta);				// this method is used to obtain the tangential component of velocity at the collision points
	long double mod() const { return sqrt(x * x + y * y); }
	long double x, y;
};