#include "MSG.h"
#include "LOS.h"
#include <iostream>


int main()
{

	setlocale(LC_ALL, "RU_ru");
	InitData<SPARSE_ROW_COL_SYMETRIC> inData;
	inData.MatrInitFile.di = "di.txt";
	inData.MatrInitFile.ig = "ig.txt";
	inData.MatrInitFile.jg = "jg.txt";
	inData.MatrInitFile.gg = "ggl.txt";
	inData.kuslau = "kuslau.txt";
	inData.f = "pr.txt";


	SLAU<SPARSE_ROW_COL_SYMETRIC_MATRIX> slau;
	slau.Load(inData);
	
	OutData out;
	SolverSLAU::LOS<OutData, SPARSE_ROW_COL_SYMETRIC_MATRIX,SLAU> solver(slau, SolverSLAU::LOSResolverModification::LLTFACTORIZATION);
	solver.solve(out);

	out.file_to_save = "res.txt";
	out.save();

	return 0;
}