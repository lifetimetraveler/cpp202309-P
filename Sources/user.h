#include<iostream>
#include<string>
using namespace std;

class User {
	string user_name;//������ �̸�
	string user_dna_sequence;//������ DNA������ �޴� ����

	//dna1~6�� DNA ������ �յڷ� ���� ���(2)*DNA�� �� �� ������ ���������ϳ� (3)
	// DNA�� codon�� �� DNA������ ���ΰ��̱⶧����
	//�������� ���� �� �ִ�. 
	string dna1;
	string dna2;
	string dna3;
	string dna4;
	string dna5;
	string dna6;
public:
	void SetName(string name);//�̸��� �ʱ�ȭ�ϴ� �Լ�
	void SetSequence(string sequence);//������ �ʱ�ȭ�ϴ� �Լ�
	void FrameSetting();//codon�� frame�� �����ϴ� �Լ�
	string GetName();//�̸��� ��ȯ�ϴ� �Լ�
	string GetSequence();//������ ��ȯ�ϴ� �Լ�
};