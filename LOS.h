#pragma once
#include"ISolverSLAU.h"
#include "Factorization.h"
#include "DataStruct.h"

namespace SolverSLAU {
	/****************************************************************************************************************************/
	/******************************************************* LOS ***************************************************************/
	using namespace Factorization;

	enum class LOSResolverModification
	{
		NONE,
		DIAGONALFACTORIZATION,
		LLTFACTORIZATION,
		//LUSQFACTORIZATION
	};


	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	class LOS : public ISolverSLAU<OutData, Matrix, SLAU>
	{
	private:

		/* ��������������� ��������� ������ */
		std::vector<double> rk;
		std::vector<double> Ark;
		std::vector<double> zk;
		std::vector<double> pk;

		/* ��������� ��� ������������ ������������ */
		std::vector<double> D;
		std::vector<double> SAQrk;
		/* ��� ������������������ �������� ����������� ��������� */
		SPARSE_ROW_COL_SYMETRIC_FACT FactMatr;
		std::vector<double> Qrk; 


		// ��������� ��� ��������������� �������� ��� ���������� ������� 
		void LOS_Classic(OutData& out);

		void LOS_Diagonal(OutData& out);

		void LOS_LLTFact(OutData& out);

		LOSResolverModification modification;
	public:

		LOS(SLAU<Matrix>& data, LOSResolverModification modification_ = LOSResolverModification::NONE) : 
			ISolverSLAU<OutData, Matrix, SLAU>(data), modification(modification_) {}

		// ���������� ���������� �������� ������� ������� ���� � ��������� ����������� 
		void solve(OutData& out);

	};


	/* Private section  */
	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	void LOS<OutData, Matrix, SLAU>::LOS_Classic(OutData& out)
	{
		// ��������� ������������ �������� ������ ����������� ���������� 
		int N = this->Data.f.size();
		int maxiter = this->Data.maxiter;
		double eps = this->Data.eps;
		// ����� ������ �� ��� �������� ������ 
		std::vector<double>& f = this->Data.f;
		std::vector<double>& xk = out.x;

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
		Ark.resize(N);
		zk.resize(N);
		pk.resize(N);


		// ���� ���������� ����������� ��� 
		if (this->Data.x.size() == 0)
		{
			std::copy(f.begin(), f.end(), rk.begin()); // r0 = f0
		}
		else
		{
			Ark = this->MultA(xk); // A*x0
			rk = f - Ark; // r0 = f - A*x0
		}

		std::copy(rk.begin(), rk.end(), zk.begin()); // z0 = r0

		pk = this->MultA(rk); // p0 = Ar0
		std::copy(pk.begin(), pk.end(), Ark.begin());

		double dot_rk = rk * rk; //  ������� ����� ��������� ������������ 
		for (int k = 1; k < maxiter; k++)
		{
			double dot_pk = pk * pk;
			double alphak = (pk * rk) / dot_pk;
			// xk
			xk = xk + alphak * zk;
			// rk
			rk = rk - alphak * pk;
			// Bettak
			Ark = this->MultA(rk);
			double bettak = -1.0 * (pk * Ark) / dot_pk;
			zk = rk + bettak * zk;
			pk = Ark + bettak * pk;

			// (rk; rk)
			dot_rk -= alphak * alphak * dot_pk;
			std::cout << "Iteration: " << k << " (rk; rk) = " << dot_rk << "\n";
			if (dot_rk < eps) break;
		}
		// ������������� �������� ������� 
		rk = f - this->MultA(xk);
		dot_rk = rk * rk;

		std::cout << "�������� ������� r = f - Ax: " << dot_rk << "\n";


	}
	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	void LOS<OutData, Matrix, SLAU>::LOS_Diagonal(OutData& out)
	{
		// ��������� ������������ �������� ������ ����������� ���������� 
		int N = this->Data.f.size();
		int maxiter = this->Data.maxiter;
		double r0_norm2 = 0.0;
		double eps = this->Data.eps* this->Data.eps;
		// ����� ������ �� ��� �������� ������ 
		std::vector<double>& f = this->Data.f;
		std::vector<double>& xk = out.x;
		std::vector<double>& di = this->Data.matr.di;

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
		Ark.resize(N);
		zk.resize(N);
		pk.resize(N);

		// ������� ������������� ������������������ 
		D = DiagonalFact(di);
		// ���� ���������� ������ ������� ��� �������������� ��������� ������� �� ������ 
		auto MultD = [](const std::vector<double>& di, const std::vector<double>& vec)
		{
			int size = vec.size();
			std::vector<double> res(size);

			for (int i = 0; i < size; i++)
				res[i] = vec[i] / di[i];
			return res;
		};

		// ���� ���������� ����������� ��� 
		if (this->Data.x.size() == 0)
		{
			std::copy(f.begin(), f.end(), rk.begin()); // r0 = f0
		}
		else
		{
			Ark = this->MultA(xk); // A*x0
			rk = f - Ark; // r0 = f - A*x0
		}

		rk = MultD(D, rk); // r0 � ���������
		r0_norm2 = rk * rk;

		zk = MultD(D, rk);

		pk = MultD(D,(this->MultA(zk)));
		SAQrk.resize(N);
		
		double dot_rk = rk * rk; //  ������� ����� ��������� ������������ 
		for (int k = 1; k < maxiter; k++)
		{
			double dot_pk = pk * pk;
			double alphak = (pk * rk) / dot_pk;
			// xk
			xk = xk + alphak * zk;
			// rk
			rk = rk - alphak * pk;
			// Bettak
			SAQrk = MultD(D, this->MultA( MultD(D, rk) ) );
			double bettak = -1.0 * (pk*SAQrk) / dot_pk;
			zk = MultD(D, rk) + bettak * zk;
			pk = SAQrk + bettak * pk;

			// (rk; rk)
			dot_rk -= alphak * alphak * dot_pk;
			std::cout << "Iteration: " << k << " (rk; rk) = " << dot_rk/r0_norm2 << "\n";
			if (dot_rk/r0_norm2 < eps) break;
		}
		
		// ������������� �������� ������� 
		rk = f - this->MultA(xk);
		dot_rk = rk * rk;

		std::cout << "�������� ������� r = f - Ax: " << dot_rk/r0_norm2 << "\n";

	}

	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	void LOS<OutData, Matrix, SLAU>::LOS_LLTFact(OutData& out)
	{
		// ��������� ������������ �������� ������ ����������� ���������� 
		int N = this->Data.f.size();
		int maxiter = this->Data.maxiter;
		double r0_norm2 = 0.0;
		double eps = this->Data.eps * this->Data.eps;
		// ����� ������ �� ��� �������� ������ 
		std::vector<double>& f = this->Data.f;
		std::vector<double>& xk = out.x;

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
		Ark.resize(N);
		zk.resize(N);
		pk.resize(N);

		// ���� ���������� ����������� ��� 
		if (this->Data.x.size() == 0)
		{
			std::copy(f.begin(), f.end(), rk.begin()); // r0 = f0
		}
		else
		{
			Ark = this->MultA(xk); // A*x0
			rk = f - Ark; // r0 = f - A*x0
		}

		// ������������ �������� ����������� ���������
		FactMatr = LLTFact(this->Data.matr);

		/* ���������� ����� ����������� ��� ���� ������� */
		rk = this->normal(FactMatr, rk);

		zk = this->reverse(FactMatr, rk);

		pk = this->normal(FactMatr,this->MultA(zk));

		r0_norm2 = rk * rk;
		double dot_rk = rk * rk; //  ������� ����� ��������� ������������ 
		SAQrk.resize(N);
		Qrk.resize(N);

		for (int k = 1; k < maxiter; k++)
		{
			double dot_pk = pk * pk;
			double alphak = (pk * rk) / dot_pk;
			// xk
			xk = xk + alphak * zk;
			// rk
			rk = rk - alphak * pk;
			// Bettak
			Qrk = this->reverse(FactMatr, rk);
			SAQrk = this->normal(FactMatr,this->MultA(Qrk));

			double bettak = -1.0 * (pk * SAQrk) / dot_pk;
			zk = Qrk + bettak * zk;
			pk = SAQrk + bettak * pk;

			// (rk; rk)
			dot_rk -= alphak * alphak * dot_pk;
			std::cout << "Iteration: " << k << " (rk; rk) = " << dot_rk/r0_norm2 << " < " << eps << "\n";
			if (dot_rk / r0_norm2 < eps) break;
		}

		// ������������� �������� ������� 
		rk = f - this->MultA(xk);
		dot_rk = rk * rk;

		std::cout << "�������� ������� r = f - Ax: " << dot_rk / r0_norm2 << "\n";
	}

	/* Public section */

	template<typename OutData, typename Matrix, template<typename Matrix> class SLAU>
	void LOS<OutData, Matrix, SLAU>::solve(OutData& out)
	{
		switch (modification)
		{
		case LOSResolverModification::NONE:
			LOS_Classic(out);
			break;
		case LOSResolverModification::DIAGONALFACTORIZATION:
			LOS_Diagonal(out);
			break;
		case LOSResolverModification::LLTFACTORIZATION:
			LOS_LLTFact(out);
			break;
		default:
			std::cout << "�������������� ��� �����������\n";
			break;
		}
	}

};