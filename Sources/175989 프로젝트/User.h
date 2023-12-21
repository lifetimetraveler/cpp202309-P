#include<iostream>
#include<string>
#include<vector>
using namespace std;

//������ ORf�� �ε����� �����ϱ� ����
class SavedIndex
{public:
	int case_num = -1;//���� �ڵ��� ��ž �ڵ��� ����� ��. 0,1,2 �� ���� �߿� �ϳ��� ����
	//�ؼ� ����
    int start_index=-1;// ���� ���� �ε���
	int stop_index = -1;// ���� ���� �ε���
	


	SavedIndex(int case_num, int start_index, int stop_index)
	{
		this->case_num = case_num;
		this->start_index = start_index;
		this->stop_index = stop_index;
	}
};
//�߸� orf���� ��Ʈ���� orf���� ����� ��ġ�� �����ϱ� ���� 
class IntronIndex
{
	int intron_start;
	int intron_stop;
};

//������ �̸�, ������ ����, ã���� �ϴ� orf�� ���̹���, reading frame�� �ٲ� ������
//  ������ ���� Ŭ����
class User {
	string user_name;//������ �̸�
	string user_dna_sequence;//������ DNA������ �޴� ����
    string dna1;//������, ���� reading frame
	string dna2;//������, ���� reading frame+1
	string dna3;//������, ���� reading frame+2
	string reverse_dna1;//������, reading frame
	string reverse_dna2;//������, reading frame+1
	string reverse_dna3;//������, reading frame+2
	//dna1~6�� DNA ������ �յڷ� ���� ���(2)*DNA�� �� �� ������ ���������ϳ� (3)
	// DNA�� codon�� �� DNA������ ���ΰ��̱⶧����
	// �������� ���� �� �ִ�. ���⿡ �������� ��츦 ��� ���ϸ� 6���� ���
	
public:
	int length1=0;//ã���� �ϴ� orf�� ������ �ּ� ��
	int length2=0;//ã���� �ϴ� orf�� ������ �ִ� ��


	void SetName(string name);//�̸��� �ʱ�ȭ�ϴ� �Լ�
	void SetSequence(string sequence);//������ �ʱ�ȭ�ϴ� �Լ�
	void FrameSetting();//���� ������ reading frame�� �޸��Ͽ� �����ϴ� �Լ�
	string GetName();//�̸��� ��ȯ�ϴ� �Լ�
	string GetSequence();//������ ��ȯ�ϴ� �Լ�
	string GetDna1();//dna1�� ��ȯ�ϴ� �Լ�
	string GetDna2();//dna2�� ��ȯ�ϴ� �Լ�
	string GetDna3();//dna3�� ��ȯ�ϴ� �Լ�
	string GetRDna1();//reaverse_dna1�� ��ȯ�ϴ� �Լ�
	string GetRDna2();//reaverse_dna2�� ��ȯ�ϴ� �Լ�
	string GetRDna3();//reaverse_dna3�� ��ȯ�ϴ� �Լ�
};

//�ϳ��� reading frame���� orf�� ã�� �װ��� �����ϴ� Ŭ����
class Orf
{
	string original_seq;//�޾ƿ� �м��� ������ ����
    
public:
	vector<SavedIndex> complete_index;//������ orf�� �ε��� ���� ����//�׽�Ʈ�� ���� cout���� ������ �ۺ�
	vector<int> complete_score;//���� ������ �����ϱ� ���� �迭
	Orf(string original_seq) { this->original_seq = original_seq; }//������
	vector<string> orf1;//ó�������� �Ų������ϱ����� ��Ʈ������ ������ ����� ����
	vector<vector<string>> complete_orf;//orf�� ������ ���͸� �����ϱ����� ���� ����
	vector<vector<string>> intron_removed;//��Ʈ���� ���ŵ� 
	vector<vector<string>> protein;//amino acid���·� ��ȯ�� ORF�� �����ϱ� ����

	//orf���� atg(start codon)�� ����
	vector<int> atg_index;

    //orf���� stop codon�� �ִ� �迭�� ����
	vector<int> tga_index;
	vector<int> taa_index;
	vector<int> tag_index;

	string GetOriSeq();



	//������ ����ȭ�ϴ� �Լ�
	void TransferSeq();
	//����ȭ�� �������� stop, start codon�� ã�� �Լ�
	void IndexFinder();
	//ã�Ƴ� �ε����� ���� ������ �����ϴ� �Լ�
	void OrfFinder(User user);
	//��Ʈ���� ã�Ƴ��� �Լ�, ���ŵ�
	void IntronFinder();
	//amino acid���·� ��ȯ�ϱ� ���� �Լ�
	void CodonDecipher();
	//kozak score�� ����ϴ� �Լ�. score�� ��ȯ.
	void KozakCalculator();
	

	
	
};

