/*
	Oscillatory-ballistic motion regularities of a gravitational pendulum
	Submitted to Nonlinear Dynamics

	Authors: Sebastian Micluta-Campeanu and Tiberius O. Cheche
	Faculty of Physics, University of Bucharest, Bucharest 077125, Romania
	Corresponding author: Tiberius O. Cheche; e-mail: cheche@gate.sinica.edu.tw
	Date: December 4, 2015

	Description: Position class implementation.
*/

#include <iostream>

#include "position.h"
#include "common.h"

// velocity Verlet
void position::update(acceleration a[2], velocity v, long double dt)
{
	x += v.x * dt + a[1].x * dt * dt / 2;
	y += v.y * dt + a[1].y * dt * dt / 2;
}

void position::updateCartesian()
{
	if (theta < 0)
		theta += 2 * PI;
	x = r * cos(theta);
	y = r * sin(theta);
	return;
}

void position::updateTheta()
{
	theta = atan2(x, -y) + 1.5 * PI;
	if (theta > 2 * PI)
		theta -= 2 * PI;
}