/*
	Oscillatory-ballistic motion regularities of a gravitational pendulum
	Submitted to Nonlinear Dynamics

	Authors: Sebastian Micluta-Campeanu and Tiberius O. Cheche
	Faculty of Physics, University of Bucharest, Bucharest 077125, Romania
	Corresponding author: Tiberius O. Cheche; e-mail: cheche@gate.sinica.edu.tw
	Date: December 4, 2015

	Description: Common variables.
*/

#pragma once
#include <fstream>

// constants
extern const long double PI;
extern const long double g;

// files
extern std::ifstream in;
extern std::ofstream out;
extern std::ofstream logFile;
extern std::ofstream logResults;

// flags
extern bool verbose;