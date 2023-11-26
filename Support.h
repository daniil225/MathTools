#pragma once
#include <vector>
#include<iostream>

namespace SolverSLAU {
	
	// Binary operation
	std::vector<double> operator+(const std::vector<double> &vec1, const std::vector<double> &vec2);
	std::vector<double> operator*(const double lamba, const std::vector<double>& vec);
	std::vector<double> operator*(const std::vector<double>& vec, const double lambda);
	// —кал€рное произведение 
	double operator*(const std::vector<double>& vec1, const std::vector<double> &vec2);
	std::vector<double> operator-(const std::vector<double>& vec1, const std::vector<double>& vec2);

	// Unar opiration
	std::vector<double>& operator+=(std::vector<double>& vec1, const std::vector<double>& vec2);
	std::vector<double>& operator*=(std::vector<double>& vec1, const double lambda);
	std::vector<double>& operator-=(std::vector<double>& vec1, const std::vector<double>& vec2);

	std::ostream& operator<<(std::ostream& cnt, const std::vector<double>& vec);
	
};