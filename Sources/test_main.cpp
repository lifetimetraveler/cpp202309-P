#include"user.h"
using namespace std;



int main() {
	string name;
	string sequence;
	User user1;
	cout << "�̸�";
	cin >> name;
	cout << "����";
	cin >> sequence;
	user1.SetName(name);
	user1.SetSequence(sequence);

	cout << "�̸�:" << user1.GetName() << endl;
	cout << "����: " << user1.GetSequence() << endl;
	return 0;
}