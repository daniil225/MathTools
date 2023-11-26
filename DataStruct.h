#pragma once
#include "config.h"
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
// Общие структуры для хранения данных их инициализации 

enum class MatrixType {
	SPARSE_ROW_COL_SYMETRIC_MATRIX,
	SPARSE_ROW_COL_MATRIX,
	PROFILE_ROW_COL_SYMETRIC_MATRIX,
	PROFILE_ROW_COL_MATRIX
};



// MatrixType::SPARSE_ROW_COL_SYMETRIC
struct SPARSE_ROW_COL_SYMETRIC
{
	std::string gg = "";
	std::string ig = "";
	std::string jg = "";
	std::string di = "";
};

struct SPARSE_ROW_COL_SYMETRIC_MATRIX
{
	MatrixType type = MatrixType::SPARSE_ROW_COL_SYMETRIC_MATRIX;
	std::vector<double> gg;
	std::vector<int> ig;
	std::vector<int> jg;
	std::vector<double> di;
};

////////////////////////////////////////////////////////////////////////


struct SPARSE_ROW_COL
{
	std::string ggu = "";
	std::string ggl = "";
	std::string ig = "";
	std::string jg = "";
	std::string di = "";
};

struct SPARSE_ROW_COL_MATRIX
{
	MatrixType type = MatrixType::SPARSE_ROW_COL_MATRIX;
	std::vector<double> ggu;
	std::vector<double> ggl;
	std::vector<int> ig;
	std::vector<int> jg;
	std::vector<double> di;
};


////////////////////////////////////////////////////////////////////////

template<typename MatrixInitStruct>
struct InitData
{
	// std::string FileFormat; - Формат файла txt - обычный текстовый или файлы прямого доступа 
	// MatrixType type = MatrixType::SPARSE_ROW_COL_SYMETRIC; // По умолчанию разряженный строчно симметричный вормат 
	int maxiter; // Максимальное число итераций 
	double eps; // Норма вектра невязки 

	std::string x = ""; // Вектор начального приближения (Не обязательный параметр) 
	std::string f = ""; //  Вектор правой части  
	std::string kuslau = ""; // Размерность, невязка, число итераций 
	MatrixInitStruct MatrInitFile; // Структура для инициализации матрицы 
};


// Структура для инициализации решателя СЛАУ итерационными методами
// Самостоятельный класс, который предоставляет метод для своей инициализаии 
template<typename Matrix>
struct SLAU
{
private:
	int N; // Размерность СЛАУ
	int size; // Размер массивов gg, ggu, jg, ggl (в зависимости от формата) 
	
	// Профильный и разряженный формат отличаются лишь наличием массива JG или его отсутсвия все Остальное точно такое же
	// Наличие полей для загрузки указывается стандартом и является резервом 
	// Создадим методы для загрузки каждой компоненты отдельно 
	
	void Load_commonPart(const std::string &kuslau, const std::string &f, const std::string &x);
	void Load_di(const std::string& filename);
	void Load_ig(const std::string& filename);

	#ifdef _SPARSE_
	
	#ifdef _SPARSE_ROW_COL_
	void Load_ggl(const std::string& filename);
	void Load_ggu(const std::string& filename);
	#endif
		
	#ifdef _SPARSE_ROW_COL_SYMETRIC_
	void Load_gg(const std::string& filename);
	#endif
	
	void Load_jg(const std::string& filename);
	
	#endif
	

public:

	int maxiter; // Максимальное число итераций 
	double eps;
	std::vector<double> f; // Вектор правой части 
	std::vector<double> x; // Начальное приближение 

	Matrix matr; // Матрица решаемой слау

	// Конструкторы не нужны интерфейс открытый => если например происходит сборка матрицы, то можно просто собирать сразу структуру СЛАУ 

	// Дружественная Функция Загрузки данных для СЛАУ из файлов (  не является компонентой стуктуры матрицы ) 
	template<typename MatrixInitStructType>
	void Load(InitData<MatrixInitStructType> &);

	// Распечатка всей информации о СЛАУ 
	template<typename T>
	friend std::ostream& operator<<(std::ostream& cnt, const SLAU<T>& slau);
};

template<typename Matrix>
std::ostream& operator<<(std::ostream& cnt, const SLAU<Matrix>& slau)
{
	if (slau.N > 0)
	{
		cnt << "СЛАУ: \n";
		cnt << "Размерность = " <<  slau.N << "  eps = " << slau.eps << "  maxiter = " << slau.maxiter << "\n";
		cnt << "di = [";
		for (int i = 0; i < slau.N; i++)
			cnt << slau.matr.di[i] << " ";
		cnt << "]\n";

		cnt << "f = [";
		for (int i = 0; i < slau.N; i++)
			cnt << slau.f[i] << " ";
		cnt << "]\n";

		cnt << "ig = [";
		for (int i = 0; i < slau.size+1; i++)
			cnt << slau.matr.ig[i] << " ";
		cnt << "]\n";

		#ifdef _SPARSE_
		cnt << "jg = [";
		for (int i = 0; i < slau.size; i++)
			cnt << slau.matr.jg[i] << " ";
		cnt << "]\n";

		#ifdef _SPARSE_ROW_COL_SYMETRIC_
		cnt << "gg = [";
		for (int i = 0; i < slau.size; i++)
			cnt << slau.matr.gg[i] << " ";
		cnt << "]\n";
		#endif

		#ifdef _SPARSE_ROW_COL_
		cnt << "ggl = [";
		for (int i = 0; i < slau.size; i++)
			cnt << slau.matr.ggl[i] << " ";
		cnt << "]\n";

		cnt << "ggu = [";
		for (int i = 0; i < slau.size; i++)
			cnt << slau.matr.ggu[i] << " ";
		cnt << "]\n";
		#endif

		#endif
	}
	else
	{
		cnt << "СЛАУ пока не инициализировалась\n";
	}

	return cnt;
}

/* Private section */

template<typename Matrix>
void SLAU<Matrix>::Load_commonPart(const std::string& kuslau, const std::string& f, const std::string& x)
{
	std::ifstream in;

	// Общая часть для всех типов матриц 
	if (kuslau != "")
	{
		in.open(kuslau);
		in >> N >> maxiter >> eps;
		in.close();
		std::cout << "kuslau Load\n";
	}

	// Загрузка правой части 
	if (f != "")
	{
		in.open(f);
		this->f.resize(N);
		for (int i = 0; i < N; i++)
			in >> this->f[i];
		in.close();
		std::cout << "f Load\n";
	}

	// Загружаем начальное приближение если оно есть 
	if (x != "")
	{
		in.open(x);
		this->x.resize(N);
		for (int i = 0; i < N; i++)
			in >> this->x[i];
		in.close();
		std::cout << "x Load\n";
	}
}

#ifdef _SPARSE_ROW_COL_SYMETRIC_
template<typename Matrix>
void SLAU<Matrix>::Load_gg(const std::string& filename)
{
	std::ifstream in;
	in.open(filename);
	matr.gg.resize(size);

	for (int i = 0; i < size; i++)
		in >> matr.gg[i];

	in.close();
	std::cout << "gg Load\n";
}
#endif

template<typename Matrix>
void SLAU<Matrix>::Load_ig(const std::string& filename)
{
	std::ifstream in;
	in.open(filename);
	bool is_C_format = false; // true == C++ Формат, false == fortran format 
	matr.ig.resize(N + 1);
	// Проверка на формат, если начало с 1 то это формат Фортрана если нет то С++
	in >> matr.ig[0];
	if (matr.ig[0] == 0) is_C_format = true;
	else matr.ig[0]--;

	for (int i = 1; i < N + 1; i++)
		in >> matr.ig[i];
	if (!is_C_format)
		for (int i = 1; i < N + 1; i++)
			matr.ig[i]--;

	in.close();

	size = matr.ig[matr.ig.size() - 1];
	std::cout << "ig Load\n";

}

#ifdef _SPARSE_
template<typename Matrix>
void SLAU<Matrix>::Load_jg(const std::string& filename)
{
	std::ifstream in;
	bool is_C_format = false; // true == start from 0 false overwise
	in.open(filename);
	matr.jg.resize(size);
	in >> matr.jg[0];
	if (matr.jg[0] == 0) is_C_format = true;
	else matr.jg[0]--;

	for (int i = 1; i < size; i++)
		in >> matr.jg[i];

	if (!is_C_format)
		for (int i = 1; i < size; i++)
			matr.jg[i]--;

	in.close();

	std::cout << "jg Load\n";
}
#endif

#ifdef _SPARSE_ROW_COL_
template<typename Matrix>
void SLAU<Matrix>::Load_ggl(const std::string& filename)
{
	std::ifstream in;
	in.open(filename);
	matr.ggl.resize(size);

	for (int i = 0; i < size; i++)
		in >> matr.ggl[i];

	in.close();
	std::cout << "ggl Load\n";
}

template<typename Matrix>
void SLAU<Matrix>::Load_ggu(const std::string& filename)
{
	std::ifstream in;
	in.open(filename);
	matr.ggu.resize(size);

	for (int i = 0; i < size; i++)
		in >> matr.ggu[i];

	in.close();
	std::cout << "ggu Load\n";
}
#endif

template<typename Matrix>
void SLAU<Matrix>::Load_di(const std::string& filename)
{
	std::ifstream in;
	in.open(filename);
	matr.di.resize(N);
	for (int i = 0; i < N; i++)
		in >> matr.di[i];
	in.close();

	std::cout << "di Load\n";
}


template<typename Matrix>
template<typename MatrixInitStructType>
void SLAU<Matrix>::Load(InitData<MatrixInitStructType>& inData)
{
	
	Load_commonPart(inData.kuslau, inData.f, inData.x);

	switch (matr.type)
	{
	case MatrixType::SPARSE_ROW_COL_SYMETRIC_MATRIX:
		Load_ig(inData.MatrInitFile.ig);
		Load_di(inData.MatrInitFile.di);

		#ifdef _SPARSE_ROW_COL_SYMETRIC_
		Load_gg(inData.MatrInitFile.gg);
		#endif

		#ifdef _SPARSE_
		Load_jg(inData.MatrInitFile.jg);
		#endif
		break;
	case MatrixType::SPARSE_ROW_COL_MATRIX:
		Load_ig(inData.MatrInitFile.ig);
		Load_di(inData.MatrInitFile.di);
		#ifdef _SPARSE_ROW_COL_
		Load_ggu(inData.MatrInitFile.gg);
		Load_ggl(inData.MatrInitFile.gg);
		#endif

		#ifdef _SPARSE_
		Load_jg(inData.MatrInitFile.jg);
		#endif

		break;
	default:
		std::cout << "Undefined Matrix type\n";
		break;
	}
}


//// Выходные данные решения СЛАУ 

struct OutData
{
	std::vector<double> x;
	std::string file_to_save = ""; // Файл для сохранения данных 

	// Сохранить данные в файл 
	void save(const std::string &delimeter = " ") const;

	friend std::ostream& operator<<(std::ostream& cnt, const OutData& out);
};