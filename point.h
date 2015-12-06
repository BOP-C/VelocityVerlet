/*
	Oscillatory-ballistic motion regularities of a gravitational pendulum
	Submitted to Nonlinear Dynamics

	Authors: Sebastian Micluta-Campeanu and Tiberius O. Cheche
	Faculty of Physics, University of Bucharest, Bucharest 077125, Romania
	Corresponding author: Tiberius O. Cheche; e-mail: cheche@gate.sinica.edu.tw
	Date: December 4, 2015

	Description: Point class header. Point inhehits coordinates and the 
	methods for modifying them from position class.
*/

#pragma once

#include "position.h"
#include "velocity.h"
#include "acceleration.h"

class point : public position
{
private:
	velocity v;
	acceleration a[2];
	long double T;
	long double mass;
	bool collision;
	int currentOutputMode;
	const int defaultOutputMode;
	void reInitialize(long double &dtheta, long double &alpha, long double &t);
	void computeAcceleration(long double &dtheta, long double &alpha, long double &t);
public:
	void move(long double &dtheta, long double &angle, long double &t);
	long double dt;
	int numberOfCollisions;
	long double lastXc;			// last collision's abscissa
	long double lastYc;			// last collision's ordinate
	long double E() const;
	void write(long double dtheta, long double dt, long double t) const;
	point(long double l, long double angle, long double omega, long double m, long double Dt, long double &dtheta, long double &t, int output);	// constructor
};