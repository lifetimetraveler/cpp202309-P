#include<iostream>
#include<string>
#include<vector>
using namespace std;

class User {
	string user_name;//������ �̸�
	string user_dna_sequence;//������ DNA������ �޴� ����
    string dna1;
	string dna2;
	string dna3;
	string reverse_dna1;
	string reverse_dna2;
	string reverse_dna3;
	//dna1~6�� DNA ������ �յڷ� ���� ���(2)*DNA�� �� �� ������ ���������ϳ� (3)
	// DNA�� codon�� �� DNA������ ���ΰ��̱⶧����
	//�������� ���� �� �ִ�. 
	
public:
	int length1=0;
	int length2=0;


	void SetName(string name);//�̸��� �ʱ�ȭ�ϴ� �Լ�
	void SetSequence(string sequence);//������ �ʱ�ȭ�ϴ� �Լ�
	void FrameSetting();//codon�� frame�� �����ϴ� �Լ�
	string GetName();//�̸��� ��ȯ�ϴ� �Լ�
	string GetSequence();//������ ��ȯ�ϴ� �Լ�
	string GetDna1();
	string GetDna2();
	string GetDna3();
	string GetRDna1();
	string GetRDna2();
	string GetRDna3();
};

class Orf
{

public:
	
	vector<vector<string>> complete_orf;
	vector<string> orf1;
	vector<int> atg_index;
	vector<int> tga_index;
	vector<int> taa_index;
	vector<int> tag_index;

	vector<string> orf2;
	vector<string> orf3;
	vector<string> reverse_orf1;
	vector<string> reverse_orf2;
	vector<string> reverse_orf3;


	void TransferSeq(User user);
	void IndexFinder();
	void OrfFinder(User user);
	
};