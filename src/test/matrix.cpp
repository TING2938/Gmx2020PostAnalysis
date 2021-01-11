/*
#include <iostream>
using namespace std;

template <typename T>
class CMatrix {
	T** data;
	int row, col;

public:
	CMatrix() {
		row = 0; 
		col = 0; 
		data = NULL;
	}

	CMatrix(int i, int j) {
		data = new T*[i];
		for (int m = 0; m < i; m++) {
			data[m] = new T[j];
		}
		row = i;
		col = j;
	}

	~CMatrix() {
		if (data != NULL)
		{
			for (int m = 0; m < row; m++) {
				delete data[m];
			}
			delete data;
		}
	}

	CMatrix<T>& operator= (CMatrix<T>& a) {
		row = a.row;
		col = a.col;
		if (data != NULL)
		{
			for (int m = 0; m < row; m++) {
				delete data[m];
			}
			delete data;
		}
		data = new T * [row];
		for (int m = 0; m < row; m++) {
			data[m] = new T[col];
		}
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				data[i][j] = a.data[i][j];
			}
		}
		return *this;
	}

	CMatrix<T>& operator+=(CMatrix<T>& a) {
		if (row != a.row || col != a.col) {
			cout << "Error" << endl;
		}
		else {
			for (int i = 0; i < row; i++) {
				for (int j = 0; j < col; j++) {
					data[i][j] += a.data[i][j];
				}
			}
		}
		return *this;
	}

	template <typename S>
	friend ostream& operator<< (ostream& os, CMatrix<S>& a);


	T& operator()(int a, int b) {
		return data[a][b];
	}

	T* operator[](int a) {
		return data[a];
	}
	
};


template <typename S>
ostream& operator<< (ostream& os, CMatrix<S>& a) {
	for (int i = 0; i < a.row; i++) {
		for (int j = 0; j < a.col; j++) {
			os << a.data[i][j] << " ";
		}
		os << endl;
	}
	//os << "value";
	return os;
}

int main() {
	CMatrix<double> A(3, 4);
	CMatrix<double> B(3, 4);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 4; j++) {
			A(i, j) = i + j;
		}
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 4; j++) {
			B[i][j] = i - j;
		}
	}
	cout << "A:" << endl << A;
	cout << "B:" << endl << B;
	CMatrix<double> C;
	C = A;
	cout << "C=A:" << endl << C;
	A += B;
	cout << "A+=B:" << endl << A;
	cout << "A(2,3):" << A(2, 3) << endl;
	cout << "B[2][3]:" << B[2][3] << endl;
	return 0;
	
}
*/

/*
#include <iostream>
#include <string>
using namespace std;

class KaoTi
{
protected:
	double score;
	string question;
	string answer;
public:
	virtual void input();
	virtual void show();
	virtual double pt();
};

class tiankong :public KaoTi {
public:
	virtual void input();
	virtual void show();
	virtual double pt();
};

class jianda :public KaoTi {
public:
	virtual void input();
	virtual void show();
	virtual double pt();
};

class chengxu :public KaoTi {
public:
	virtual void input();
	virtual void show();
	virtual double pt();
};

int main() {
	KaoTi* A[3];
	A[0] = new tiankong;
	A[1] = new jianda;
	A[2] = new chengxu;
	double sore = 0;
	for (int i = 0; i < 3; i++) {
		A[i]->input();
		sore += A[i]->pt();
		A[i]->show();
	}
	return 0;
}
*/

#include <iostream>

void func(int (*arr)[4][5][6])
{
	arr[0][0][0][0] = 12;
}

int main()
{
	int array[3][4][5][6];
	func(array);
	std::cout << array[0][0][0][0] << std::endl;
}