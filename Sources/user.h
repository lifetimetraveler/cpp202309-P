#include<iostream>
#include<string>
#include<vector>
using namespace std;

class User {
	string user_name;//������ �̸�
	string user_dna_sequence;//������ DNA������ �޴� ����

	//dna1~6�� DNA ������ �յڷ� ���� ���(2)*DNA�� �� �� ������ ���������ϳ� (3)
	// DNA�� codon�� �� DNA������ ���ΰ��̱⶧����
	//�������� ���� �� �ִ�. 
	
public:
	string dna1;
	string dna2;
	string dna3;
	string reverse_dna1;
	string reverse_dna2;
	string reverse_dna3;
	void SetName(string name);//�̸��� �ʱ�ȭ�ϴ� �Լ�
	void SetSequence(string sequence);//������ �ʱ�ȭ�ϴ� �Լ�
	void FrameSetting();//codon�� frame�� �����ϴ� �Լ�
	string GetName();//�̸��� ��ȯ�ϴ� �Լ�
	string GetSequence();//������ ��ȯ�ϴ� �Լ�
};

class Orf
{

public:
	vector<string> orf1;
	vector<string> orf2;
	vector<string> orf3;
	vector<string> reverse_orf1;
	vector<string> reverse_orf2;
	vector<string> reverse_orf3;
	void TransferSeq(User user);

};