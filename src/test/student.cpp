/*
#include <iostream>
#include <string>
using namespace std;


class person
{
protected:
	string NO;
	string name;
public:
	person(){}
	virtual void input() = 0;
	virtual void show() = 0;
};

class student: public person
{
	string CNO;
	double score;
public:
	virtual void input() {
		cout << "\nplease input student information:" << endl;
		cout << "please input NO:" << endl;
		cin >> NO;
		cout << "please input name:" << endl;
		cin >> name;
		cout << "please input class NO:" << endl;
		cin >> CNO;
		cout << "please input score:" << endl;
		cin >> score;
	}
	virtual void show() {
		cout << "\nstudent information:" << endl;
		cout << "NO:" << NO << endl;
		cout << "name:" << name << endl;
		cout << "class NO:" << CNO << endl;
		cout << "score:" << score << endl;
	}
};
class teacher : public person
{
	string title;
	string depart;
public:
	virtual void input() {
		cout << "\nplease input teacher information:" << endl;
		cout << "please input NO:" << endl;
		cin >> NO;
		cout << "please input name:" << endl;
		cin >> name;
		cout << "please input title:" << endl;
		cin >> title;
		cout << "please input department:" << endl;
		cin >> depart;
	}
	virtual void show() {
		cout << "\nteacher information:" << endl;
		cout << "NO:" << NO << endl;
		cout << "name:" << name << endl;
		cout << "title:" << title << endl;
		cout << "department:" << depart << endl;
	}
};

void showInfo(person& per)
{
	per.show();
}

void inputInfo(person& per)
{
	per.input();
}

int main() {
	student a;
	teacher b;

	inputInfo(a);
	inputInfo(b);

	showInfo(a);
	showInfo(b);

	return 0;
}
*/

/*
#include <iostream>
using namespace std;

class A {
public:
	virtual void Print(int a)
	{
		cout << "A::P1 " << "a = " << a << endl;
	}
	virtual void Print(float a, double b = 2.8)
	{
		cout << "A::P2 " << "a = " << a << " , b = " << b << endl;
	}
};
class B : public A {
public:
	virtual void Print(int a)
	{
		cout << "B::P1 " << "a = " << a << endl;
	}
	virtual void Print(float a, float b)
	{
		cout << "B::P2 " << "a = " << a << " , b = " << b << endl;
	}
};
void Show(A* p)
{
	p->Print(2);
	p->Print(2, 1.9);
}

int main()
{
	A* pa = new A;
	B* pb = new B;
	Show(pa);
	Show(pb);
	delete pa;
	delete pb;
}
*/

/*
#include <iostream>
#include <string>
using namespace std;

class OBJ1
{
public:
	OBJ1() { cout << "1-OBJ1类构造" << endl; }
};

class OBJ2
{
public:
	OBJ2() { cout << "2-OBJ2类构造" << endl; }
};

class Base1
{
public:
	Base1() { cout << "3-Base1类构造" << endl; }
};
class Base2 : virtual public Base1
{
public:
	Base2() { cout << "4-Base2类构造" << endl; }
};

class Base3
{
public:
	Base3() { cout << "5-Base3类构造" << endl; }
};

class Base4
{
public:
	Base4() { cout << "6-Base4类构造" << endl; }
};

class DerivedA : public Base1, virtual public Base2,
	public Base3, virtual public Base4
{
public:
	DerivedA() :Base4(), Base3(), Base2(), Base1(), obj2(), obj1()
	{
		cout << "7-派生类构造成功" << endl;
	}
protected:
	OBJ1 obj1;
	OBJ2 obj2;
	Base1 base1;
	Base2 base2;
};

int main(int argc, char* argv[])
{
	DerivedA aa;

	cout << "8-end" << endl;
	return 0;
}
*/

#include <iostream>
#include <cmath>
using namespace std;
class Cshape {
public:
	virtual void input() = 0;
	virtual void area() = 0;
	virtual void show() = 0;
};

class Ellipse : public Cshape {
	double a=4, b, s;
public:
	virtual void input() {
		cout << "\nplease input Ellipse information:" << endl;
		cout << "please input a:" << endl;
		cin >> a;
		cout << "please input b:" << endl;
		cin >> b;
	}
	virtual void area() {
		s = 3.1415926 * a * b;
	}
	virtual void show() {
		cout << "\nEllipse information:" << endl;
		cout << "a=" << a << endl;
		cout << "b=" << b << endl;
		cout << "area=" << s << endl;
	}
};

class Rectangle : public Cshape {
	double a, b, s;
public:
	virtual void input() {
		cout << "\nplease input Rectangle information:" << endl;
		cout << "please input a:" << endl;
		cin >> a;
		cout << "please input b:" << endl;
		cin >> b;
	}
	virtual void area() {
		s = a * b;
	}
	virtual void show() {
		cout << "\nRectangle information:" << endl;
		cout << "a=" << a << endl;
		cout << "b=" << b << endl;
		cout << "area=" << s << endl;
	}
};

class Square : public Cshape {
	double a, s;
public:
	virtual void input() {
		cout << "\nplease input Square information:" << endl;
		cout << "please input a:" << endl;
		cin >> a;
		
	}
	virtual void area() {
		s = a * a;
	}
	virtual void show() {
		cout << "\nSquare information:" << endl;
		cout << "a=" << a << endl;
		cout << "area=" << s << endl;
	}
};

class Triangle : public Cshape {
	double a, b, c, s;
public:
	virtual void input() {
		cout << "\nplease input Triangle information:" << endl;
		cout << "please input a:" << endl;
		cin >> a;
		cout << "please input b:" << endl;
		cin >> b;
		cout << "please input c:" << endl;
		cin >> c;
	}
	virtual void area() {
		double p = (a + b + c) / 2;
		s = sqrt(p * (p - a) * (p - b) * (p - c));
	}
	virtual void show() {
		cout << "\nTriangle information:" << endl;
		cout << "a=" << a << endl;
		cout << "b=" << b << endl;
		cout << "c=" << c << endl;
		cout << "area=" << s << endl;
	}
};

/*void inputinfo(Cshape* a) {
	a->input();
}
void calcarea(Cshape* a) {
	a->area();
}
void showinfo(Cshape* a) {
	a->show();
}*/


int main() {
	Cshape* a[4];
	a[0] = new Ellipse;
	a[1] = new Rectangle;
	a[2] = new Square;
	a[3] = new Triangle;
	for (int i = 0; i < 4; i++) {
		a[i]->input();
		a[i]->area();
		a[i]->show();
	}
	return 0;
}