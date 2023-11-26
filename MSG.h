#pragma once

#include "DataStruct.h"
#include "ISolverSLAU.h"
#include "Factorization.h"



namespace SolverSLAU {

	//using namespace Factorization;

	//  ����������� ���������� ��� ������� ����  
	//	NONE - ��� ����������� ����������� �������� 
	//	DIAGONALFACTORIZATION - ������������ ������������������ 
	//	LLTFACTORIZATION - �������� ������������ ��������� 
	//	LUSQFACTORIZATION - �������� ������������ � ����������� ���������� ���������� ������� LU(sq)
	enum class MSGResolverModification
	{
		NONE,
		DIAGONALFACTORIZATION,
		LLTFACTORIZATION,
		//LUSQFACTORIZATION
	};


	/****************************************************************************************************************************/
	/******************************************************* MSG *********************************************************************/

	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	class MSG : public ISolverSLAU<OutData, Matrix, SLAU>
	{
	private:

		MSGResolverModification modification;

		// ��������� ��� ��������������� �������� ��� ���������� ������� 
		//MSGResolverModification::NONE
		std::vector<double> rk;
		std::vector<double> Azk;
		std::vector<double> zk;

		/* ��� ������������ LLT */
		Factorization::SPARSE_ROW_COL_SYMETRIC_FACT FactMatr;

		// xk - ��� ������ � �������� ������ 
		double f_norm = 0.0; // ����� ������� f
		double rk_norm = 0.0; // ����� ������� �������
		double non_repan = 0.0; // ������� 

		// ������������� �������� ����� ������������ ��������������� ����� �������� ������ 
		void MSG_Classic(OutData& x);

		// ������������ ������������������ ��� ������� ������� 
		void MSG_Diagonal(OutData& x);

		void MSG_LLTFact(OutData& x);


	public:

		MSG(SLAU<Matrix>& data, MSGResolverModification modification_ = MSGResolverModification::NONE) :
			ISolverSLAU<OutData, Matrix, SLAU>(data), modification(modification_) {}


		// ���������� ���������� �������� ������� ������� ���� � ��������� ����������� 
		void solve(OutData& out);

		~MSG() {}

	};


	/* Private section  */
	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	void  MSG<OutData, Matrix, SLAU>::MSG_Classic(OutData& outData)
	{
		// ��������� ������������ �������� ������ ����������� ���������� 
		int N = this->Data.f.size();
		int maxiter = this->Data.maxiter;
		double eps = this->Data.eps;
		// ����� ������ �� ��� �������� ������ 
		std::vector<double>& f = this->Data.f;
		std::vector<double>& xk = outData.x;
		f_norm = std::sqrt(f * f); // ����� ������� f

		if (this->Data.x.size() == 0)
		{
			// ���� ���������� ����������� ��� �� 
			xk.resize(N);
		}
		else
		{
			//  ���� ��������� ����������� ��� �� ����, ��
			xk.resize(N);
			std::copy(this->Data.x.begin(), this->Data.x.end(), xk.begin());
		}

		rk.resize(N);
		Azk.resize(N);
		zk.resize(N);

		// ���� ���������� ����������� ��� 
		if (this->Data.x.size() == 0)
		{
			std::copy(f.begin(), f.end(), rk.begin()); // r0 = f0
		}
		else
		{
			Azk = this->MultA(xk); // A*x0
			rk = f - Azk; // r0 = f - A*x0
		}

		std::copy(rk.begin(), rk.end(), zk.begin()); // z0 = r0
		non_repan = std::sqrt(rk * rk) / f_norm; // ������������� ������� 

		for (int k = 1; k < maxiter; k++)
		{
			Azk = this->MultA(zk);

			double dot_rk = rk * rk;
			double alphak = dot_rk / (Azk * zk);
			xk = xk + alphak * zk;
			rk = rk - alphak * Azk;

			double bettak = rk * rk / dot_rk;

			zk = rk + bettak * zk;

			non_repan = std::sqrt(rk * rk) / f_norm;
			std::cout << "Iteration: " << k << " Non-repan: " << non_repan << "\n";

			if (non_repan < eps) break;
		}

		// ������������� �������� ������� 
		rk = f - this->MultA(xk);
		non_repan = std::sqrt(rk * rk) / f_norm;

		std::cout << "�������� ������� r = f - Ax: " << non_repan << "\n";

	}


	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	void  MSG<OutData, Matrix, SLAU>::MSG_Diagonal(OutData& outData)
	{
		// ��������� �������� ������ ����������� ���������� � ������������ ������������������� 
		int N = this->Data.f.size();
		int maxiter = this->Data.maxiter;
		double eps = this->Data.eps;
		// ����� ������ �� ��� �������� ������ 
		std::vector<double>& di = this->Data.matr.di;
		std::vector<double>& f = this->Data.f;
		std::vector<double>& xk = outData.x;
		f_norm = std::sqrt(f * f); // ����� ������� f

		if (this->Data.x.size() == 0)
		{
			// ���� ���������� ����������� ��� �� 
			xk.resize(N);
		}
		else
		{
			//  ���� ��������� ����������� ��� �� ����, ��
			xk.resize(N);
			std::copy(this->Data.x.begin(), this->Data.x.end(), xk.begin());
		}

		rk.resize(N);
		Azk.resize(N);
		zk.resize(N);
		// ���� ���������� ����������� ��� 
		if (this->Data.x.size() == 0)
		{
			std::copy(f.begin(), f.end(), rk.begin()); // r0 = f0
		}
		else
		{
			Azk = this->MultA(xk); // A*x0
			rk = f - Azk; // r0 = f - A*x0
		}

		// ���� ���������� ������ ������� ��� �������������� ��������� ������� �� ������ 
		auto MultD = [](const std::vector<double> &di,const std::vector<double>& vec)
		{
			int size = vec.size();
			std::vector<double> res(size);

			for (int i = 0; i < size; i++)
				res[i] = vec[i] / di[i];
			return res;
		};

		zk = MultD(di,rk);

		for (int k = 1; k < maxiter; k++)
		{
			Azk = this->MultA(zk);
			double dot_rk = MultD(di, rk) * rk;
			double alphak = dot_rk /( Azk * zk);

			xk = xk + alphak * zk;
			rk = rk - alphak * Azk;

			std::vector<double> Mrk = MultD(di, rk);
			double bettak = (Mrk * rk) / dot_rk;
			zk = Mrk + bettak * zk;


			non_repan = std::sqrt(rk * rk) / f_norm;
			std::cout << "Iteration: " << k << " Non-repan: " << non_repan << "\n";

			if (non_repan < eps) break;
		}

		// ������������� �������� ������� 
		rk = f - this->MultA(xk);
		non_repan = std::sqrt(rk * rk) / f_norm;

		std::cout << "�������� ������� r = f - Ax: " << non_repan << "\n";
		
	}

	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	void  MSG<OutData, Matrix, SLAU>::MSG_LLTFact(OutData& outData)
	{
		// ��������� �������� ������ ����������� ���������� � ������������ ������������������� 
		int N = this->Data.f.size();
		int maxiter = this->Data.maxiter;
		double eps = this->Data.eps;
		// ����� ������ �� ��� �������� ������ 
		std::vector<double>& di = this->Data.matr.di;
		std::vector<double>& f = this->Data.f;
		std::vector<double>& xk = outData.x;
		f_norm = std::sqrt(f * f); // ����� ������� f

		if (this->Data.x.size() == 0)
		{
			// ���� ���������� ����������� ��� �� 
			xk.resize(N);
		}
		else
		{
			//  ���� ��������� ����������� ��� �� ����, ��
			xk.resize(N);
			std::copy(this->Data.x.begin(), this->Data.x.end(), xk.begin());
		}

		rk.resize(N);
		Azk.resize(N);
		zk.resize(N);
		// ���� ���������� ����������� ��� 
		if (this->Data.x.size() == 0)
		{
			std::copy(f.begin(), f.end(), rk.begin()); // r0 = f0
		}
		else
		{
			Azk = this->MultA(xk); // A*x0
			rk = f - Azk; // r0 = f - A*x0
		}

		
		// �������� ���������� 
		FactMatr = Factorization::LLTFact(this->Data.matr);
		

		zk = this->M_1x(FactMatr, rk);


		for (int k = 1; k < maxiter; k++)
		{
			Azk = this->MultA(zk);
			double dot_rk = this->M_1x(FactMatr, rk) * rk;
			double alphak = dot_rk / (Azk * zk);

			xk = xk + alphak * zk;
			rk = rk - alphak * Azk;

			std::vector<double> tmp = this->M_1x(FactMatr, rk);

			double bettak = (tmp * rk) / dot_rk;

			zk = tmp + bettak * zk;


			non_repan = std::sqrt(rk * rk) / f_norm;
			std::cout << "Iteration: " << k << " Non-repan: " << non_repan << "\n";

			if (non_repan < eps) break;
		}


		// ������������� �������� ������� 
		rk = f - this->MultA(xk);
		non_repan = std::sqrt(rk * rk) / f_norm;

		std::cout << "�������� ������� r = f - Ax: " << non_repan << "\n";

	}

	/* Public section */
	// ���������� ���������� �������� ������� ������� ���� � ��������� ����������� 
	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	void MSG<OutData, Matrix, SLAU>::solve(OutData& out)
	{
		switch (modification)
		{
		case MSGResolverModification::NONE:
			MSG_Classic(out);
			break;
		case MSGResolverModification::DIAGONALFACTORIZATION:
			MSG_Diagonal(out);
			break;
		case MSGResolverModification::LLTFACTORIZATION:
			MSG_LLTFact(out);
			break;
		default:
			std::cout << "�������������� ��� �����������\n";
			break;
		}
	}
};