#include"user.h"
using namespace std;



int main() {
	string name;//������ �̸�
	string sequence;//������ ����
	
	User user;//User�� �̸�, ����, ORF����, �ٸ� reading frame������
	          //�����ϴ� Ŭ����
	
	cout << "�̸�";
	//name = "����ȣȣ";
	cout << "DNA����:";
	cin >> name;
	//sequence = "AGCCAAAAAATGGGGAAAGGGAAACCCAAAGGGTGATTTGGGAAA";
	cin >> sequence;
	cout << "ã�����ϴ� ORF�� ũ�� : ";
	cin >> user.length1 >> user.length2;//ORF�� ������ �Է�
	//user.length1 = 1;
	//user.length2 = 100;

	user.SetName(name);//��ü���� �̸��ʱ�ȭ
	user.SetSequence(sequence);//��ü���� ���� �ʱ�ȭ

	cout << "�̸�:" << user.GetName() << endl;
	cout << "����: " << user.GetSequence() << endl;

	user.FrameSetting();
	cout << "dna2: " << user.GetDna2() << endl;
	cout << "dna3: " << user.GetDna3() << endl;
	cout << "dna4:" << user.GetRDna1() << endl;
	cout << "dna5:" << user.GetRDna2() << endl;
	cout << "dna6:" << user.GetRDna3() << endl;
	cout << "dna1:" << user.GetDna1() << endl;

	Orf user_seq1(user.GetDna1());//Orf��ü�� �м��� ���� ����
	user_seq1.TransferSeq();//������ ��Ʈ�� ����ȭ

	
	cout << "����" << endl;
	for (int i = 0; i < user_seq1.orf1.size(); i++)
	{
		cout<<user_seq1.orf1[i]<<" ";
	}
	cout << endl;
	user_seq1.IndexFinder();//�ε��� ã��
	user_seq1.OrfFinder(user);//orf ã��


	
	//�Ϸ��� orf ���
	for (int i = 0; i < user_seq1.complete_orf[0].size(); i++)
	{
		cout << user_seq1.complete_orf[0][i] << " ";
	}



	return 0;
}