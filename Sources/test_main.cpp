#include"user.h"
using namespace std;



int main() {
	string name;
	string sequence;
	User user1;
	cout << "이름";
	cin >> name;
	cout << "서열";
	cin >> sequence;
	user1.SetName(name);
	user1.SetSequence(sequence);

	cout << "이름:" << user1.GetName() << endl;
	cout << "서열: " << user1.GetSequence() << endl;
	return 0;
}