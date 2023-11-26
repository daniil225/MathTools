#include "DataStruct.h"


void OutData::save(const std::string& delimeter) const
{
	if (file_to_save != "")
	{
		std::ofstream out;
		out.open(file_to_save);
		for (int i = 0; i < x.size(); i++)
			out << x[i] << delimeter;
		out.close();
		std::cout << "x сохранен в файл " << file_to_save << "\n";
	}
	else
	{
		std::cout << "Файл для сохранения не задан\n";
	}
}

std::ostream& operator<<(std::ostream& cnt, const OutData& out)
{
	cnt << "x = [";
	for (int i = 0; i < out.x.size(); i++)
		cnt << out.x[i] << " ";
	cnt << "]\n";

	return cnt;
}