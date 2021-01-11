#include <itp/core>
using namespace std;

vector<int> getACF(const vector<vector<int>>& data)
{
	
	int nframe = data.size();
	int nmol = data[0].size();
	vector<int> acf(nframe, 0);

	for (int m = 0; m < nmol; m++)
	{
		for (int i = 0; i < nframe; i++)
		{
			if (data[i][m] == 1)
			{
				for (int j = 0; j < nframe - i; j++)
				{
					if (data[i + j][m] == 1)
					{
						acf[j]++;
					}
					else
					{
						break;
					}
				}
			}
		}
	}
	return acf;
}

bool contain(const vector<int>& vec, int value)
{
	return std::find(vec.begin(), vec.end(), value) != vec.end();
}

int main()
{
	vector<int> data = { 1, 2, 3, 5 };

	fmt::print("{}\n", contain(data, 4));
	

}