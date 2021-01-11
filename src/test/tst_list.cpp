/* learn chain table */

#include <iostream>

class Node;
class List;

class Node
{
public:
	int elem;
	Node* link;
};


Node* newNode()
{
	return new Node;
}

void deleteNode(Node* p)
{
	delete p;
}

class List
{
public:
	using iterator = Node * ;

	List()
	{
		_first = nullptr;
		_size = 0;
	}

	~List()
	{
		auto p = _first;
		while (p) {
			p = _first->link;
			deleteNode(_first);
			_first = p;
		}
	}

	void append(int value)
	{
		auto push = newNode();
		push->elem = value;
		push->link = nullptr;

		if (_first != nullptr) { 
			auto b = end();
			b->link = push;
		} else {
			_first = push;
		}
		++_size; 
	}

	iterator begin()
	{
		return _first;
	}

	iterator end()
	{
		
		Node* sol = _first; 
		for (size_t i = 1; i < _size; ++i) {
			sol = _first->link;
		}
		return sol;
	}

	void print()
	{
		auto p = _first;
		while (p) {
			std::cout << p->elem << ' ';
			p = p->link;
		}
		std::cout << std::endl;
	}
private:
	Node* _first;
	size_t _size;
};

int main()
{
	List l;
	l.append(2);
	l.append(4);
	l.print();
}
