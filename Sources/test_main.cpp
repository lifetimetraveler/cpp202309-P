#include"user.h"
using namespace std;



int main() {
	string name;
	string sequence;
	int length;
	int length1;
	User user1;
	Orf user1o;
	cout << "�̸�";
	name = "����ȣȣ";
	cout << "DNA����:";
	//cin >> name;
	sequence = "abcabcabcabcabcabc";
	//cin >> sequence;
	cout << "ã�����ϴ� ORF�� ũ�� : ";
	//cin>>length>>length1;
	length = 100;
	length1 = 1000;

	user1.SetName(name);
	user1.SetSequence(sequence);

	cout << "�̸�:" << user1.GetName() << endl;
	cout << "����: " << user1.GetSequence() << endl;

	user1.FrameSetting();
	cout << "dna2: " << user1.GetDna2() << endl;
	cout << "dna3: " << user1.GetDna3() << endl;
	cout << "dna4:" << user1.GetRDna1() << endl;
	cout << "dna5:" << user1.GetRDna2() << endl;
	cout << "dna6:" << user1.GetRDna3() << endl;
	cout << "dna1:" << user1.GetDna1() << endl;
	user1o.TransferSeq(user1);
	cout << "����" << endl;
	for (int i = 0; i < user1o.orf1.size(); i++)
	{
		cout<<user1o.orf1[i]<<" ";
	}
	return 0;
}