#include"user.h"
using namespace std;



int main() {
	string name;
	string sequence;
	int length;
	int length1;
	User user1;
	Orf user1o;
	cout << "이름";
	name = "하하호호";
	cout << "DNA서열:";
	//cin >> name;
	sequence = "abcabcabcabcabcabc";
	//cin >> sequence;
	cout << "찾고자하는 ORF의 크기 : ";
	//cin>>length>>length1;
	length = 100;
	length1 = 1000;

	user1.SetName(name);
	user1.SetSequence(sequence);

	cout << "이름:" << user1.GetName() << endl;
	cout << "서열: " << user1.GetSequence() << endl;

	user1.FrameSetting();
	cout << "dna2: " << user1.GetDna2() << endl;
	cout << "dna3: " << user1.GetDna3() << endl;
	cout << "dna4:" << user1.GetRDna1() << endl;
	cout << "dna5:" << user1.GetRDna2() << endl;
	cout << "dna6:" << user1.GetRDna3() << endl;
	cout << "dna1:" << user1.GetDna1() << endl;
	user1o.TransferSeq(user1);
	cout << "벡터" << endl;
	for (int i = 0; i < user1o.orf1.size(); i++)
	{
		cout<<user1o.orf1[i]<<" ";
	}
	return 0;
}