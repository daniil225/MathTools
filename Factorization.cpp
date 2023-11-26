#include "Factorization.h"
#include <cmath>
#include <iostream>


namespace Factorization
{
	std::vector<double> DiagonalFact(const std::vector<double>& di)
	{
		int size = di.size();
		std::vector<double> fact_di(size);

		for (int i = 0; i < size; i++)
			fact_di[i] = 1 / std::sqrt(di[i]);
		return fact_di;
	}


	SPARSE_ROW_COL_SYMETRIC_FACT LLTFact(SPARSE_ROW_COL_SYMETRIC_MATRIX& matr)
	{
        int N = matr.di.size();
        int size = matr.gg.size();
        std::vector<double>& di = matr.di;
        std::vector<double>& gg = matr.gg;
        std::vector<int>& ig = matr.ig;
        std::vector<int>& jg = matr.jg;

       
        SPARSE_ROW_COL_SYMETRIC_FACT res;
       
        res.di.resize(N);
        res.gg.resize(size);
        
        std::copy(gg.begin(), gg.end(), res.gg.begin());
        std::copy(di.begin(), di.end(), res.di.begin());

        

        std::vector<double>& di_l = res.di;
        std::vector<double>& gg_l = res.gg;

        
        for (size_t i = 0; i < N; i++) {
            double sum_diag = 0;

            for (size_t j = ig[i]; j < ig[i + 1]; j++) {
                double sum = 0;
                int jk = ig[jg[j]];
                int ik = ig[i];
                while ((ik < j) && (jk < ig[jg[j] + 1]))
                {
                    int l = jg[jk] - jg[ik];
                    if (l == 0) {
                        sum += gg_l[jk] * gg_l[ik];
                        ik++; jk++;
                    }
                    jk += (l < 0);
                    ik += (l > 0);
                }
                gg_l[j] -= sum;
                gg_l[j] /= di_l[jg[j]];
                sum_diag += gg_l[j] * gg_l[j];
            }
            di_l[i] -= sum_diag;
            di_l[i] = sqrt(abs(di_l[i]));
        }

        return res;
	}


	//SPARSE_ROW_COL_FACT LUsqFact(SPARSE_ROW_COL& matr)
	//{

	//}


	//SPARSE_ROW_COL_FACT LUFact(SPARSE_ROW_COL& matr)
	//{

	//}
};
