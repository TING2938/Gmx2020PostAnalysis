#include <cstring>
#include <iostream>
using namespace std;
class String
{
	char* str;
	int len;
public:
	String() { len = 0; str = nullptr; }
	String(const char* a) {
		len = strlen(a);
		str = new char[len];
		for (int i = 0; i < len; i++) {
			str[i] = a[i];
		}
	}
	String& operator= (const String& n) {
		if (len < n.len) {
			delete[] str;
			str = new char[n.len];
		}

		for (int i = 0; i < n.len; i++) {
			str[i] = n.str[i];
		}
		return *this;
	}
	String(const String& n) {
		if (len < n.len) {
			delete[] str;
			str = new char[n.len];
		}
		len = n.len;
		for (int i = 0; i < n.len; i++) {
			str[i] = n.str[i];
		}
	}
	~String() {
		delete[] str;
	}

	void print() const {
		for (int i = 0; i < len; i++) {
			cout << str[i];
		}
		cout << endl;
	}

	friend String operator + (const String& a, const String& b);
};

String operator + (const String& a, const String& b) {
	String n;
	n.str = new char[a.len + b.len];
	n.len = a.len + b.len;
	for (int i = 0; i < a.len; i++) {
		n.str[i] = a.str[i];
	}
	
	for (int i = 0; i < b.len; i++) {
		n.str[a.len + i] = b.str[i];
	}
	
	return n;
}

int main() {
	String a("aaa"),b("bbb");
	String c = a + b;
	c = a;
	a.print();
	b.print();
	c.print();
	return 0;
}