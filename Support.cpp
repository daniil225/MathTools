#include "Support.h"

namespace SolverSLAU {

	// Binary operation
	std::vector<double> operator+(const std::vector<double>& vec1, const std::vector<double>& vec2)
	{
		std::vector<double> res;
		int size = vec1.size();
		res.resize(size);

		for (int i = 0; i < vec1.size(); i++)
			res[i] = vec1[i] + vec2[i];

		return res;
	}

	std::vector<double> operator*(const double lamba, const std::vector<double>& vec)
	{
		std::vector<double> res;
		int size = vec.size();
		res.resize(size);

		for (int i = 0; i < size; i++)
			res[i] = vec[i] * lamba;

		return res;
	}
	
	std::vector<double> operator*(const std::vector<double>& vec, const double lambda)
	{
		return operator*(lambda, vec);
	}

	// —кал€рное произведение 
	double operator*(const std::vector<double>& vec1, const std::vector<double>& vec2)
	{
		double res = 0.0;
		int size = vec1.size();

		for (int i = 0; i < size; i++)
			res += vec1[i] * vec2[i];

		return res;
	}

	std::vector<double> operator-(const std::vector<double>& vec1, const std::vector<double>& vec2)
	{
		std::vector<double> res;
		int size = vec1.size();
		res.resize(size);
		
		for (int i = 0; i < size; i++)
			res[i] = vec1[i] - vec2[i];

		return res;
	}

	// Unar opiration
	std::vector<double>& operator+=(std::vector<double>& vec1, const std::vector<double>& vec2)
	{
		int size = vec1.size();

		for (int i = 0; i < size; i++)
			vec1[i] += vec2[i];

		return vec1;
	}
	std::vector<double>& operator*=(std::vector<double>& vec1, const double lambda)
	{
		int size = vec1.size();

		for (int i = 0; i < size; i++)
			vec1[i] *= lambda;

		return vec1;
	}

	std::vector<double>& operator-=(std::vector<double>& vec1, const std::vector<double>& vec2)
	{
		int size = vec1.size();

		for (int i = 0; i < size; i++)
			vec1[i] -= vec2[i];

		return vec1;
	}

	std::ostream& operator<<(std::ostream& cnt, const std::vector<double>& vec)
	{
		int size = vec.size();
		cnt << "[";
		for (int i = 0; i < size; i++)
			cnt << vec[i] << " ";
		cnt << "]\n";

		return cnt;
	}


};