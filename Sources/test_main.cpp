#include"user.h"
using namespace std;



int main() {
	string name;
	string sequence;
	User user1;
	cout << "�̸�";
	name = "����ȣȣ";
	//cin >> name;
	sequence = "abcabc";
	//cin >> sequence;
	user1.SetName(name);
	user1.SetSequence(sequence);

	cout << "�̸�:" << user1.GetName() << endl;
	cout << "����: " << user1.GetSequence() << endl;

	user1.FrameSetting();
cout << "dna2: " << user1.GetDna2() << endl;
cout << "dna3: " << user1.GetDna3() << endl;
	cout << "dna4:" << user1.GetRDna1() << endl;
	cout << "dna5:" << user1.GetRDna2()<<endl;
	cout << "dna6:" << user1.GetRDna3() << endl;
	

	cout << "dna1:" << user1.GetDna1() << endl;
	return 0;
}