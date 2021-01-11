#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <map>
#include <vector>
#include <sstream>
using namespace std;

void getIndexFromSel(const string& NdxFnm, int* Ngx, int** Index, char** Grpname, const vector<string>& selection)
{
	std::ifstream fid(NdxFnm);
	if (!fid)
	{
		printf("Onfly Error: Connot open ndx file: %s, Exit!", NdxFnm.c_str());
		std::exit(-2);
	}
	std::string line, key;
	std::map<std::string, std::vector<int>> index;
	size_t n1, n2;
	std::stringstream ss;
	int num;

	while (std::getline(fid, line))
	{
		if (line.empty())
			continue;
		n1 = line.find('[');
		if (n1 != std::string::npos)
		{
			n2 = line.find(']');
			key = line.substr(n1 + 1, n2 - n1 - 1);
			key = " " + key + " ";
			key.erase(0, key.find_first_not_of(" "));
			key.erase(key.find_last_not_of(" ") + 1);
			while (true)
			{
				if (index.find(key) != index.end())
					key += "_repeat";
				else
					break;
			}
			continue;
		}
		ss.str(line);
		while (ss >> num)
		{
			index[key].push_back(num - 1);
		}
		ss.clear();
	}

	vector<int> tmp;
	for (int i = 0; i < selection.size(); i++)
	{
		tmp = index[selection[i]];
		Index[i] = new int[tmp.size()];
		std::copy_n(tmp.data(), tmp.size(), Index[i]);
		Ngx[i] = tmp.size();
		Grpname[i] = new char[tmp.size()+1];
		std::copy_n(selection[i].begin(), selection[i].size(), Grpname[i]);
		Grpname[i][selection[i].size()] = '\0';
	}
}


int main()
{
	string fnm = "index.ndx";
	int ngrps = 1;
	int* ngx;
	int** index;
	char** grpname;
	vector<string> sel = {"System", "NA"};

	ngx = new int[ngrps];
	index = new int* [ngrps];
	grpname = new char* [ngrps];

	getIndexFromSel(fnm, ngx, index, grpname, sel);


}