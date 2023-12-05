#include"user.h"
#include <fstream>
using namespace std;



int main() {
	string name;//������ �̸�
	string sequence;//������ ����
	ifstream sequence_file("sequence.txt");//������ ���� ����
	User user;//User�� �̸�, ����, ORF����, �ٸ� reading frame������
	//�����ϴ� Ŭ����

//�̸��Է�, �׽�Ʈ�� ���� ������ �̸��� ����, ������ ����߽��ϴ�.
	cout << "�̸�: " << endl;
	//cin >> name;
	name = "����ȣȣ";

	//���� �Է�
	cout << "DNA����:" << endl;
	if (!sequence_file)//���� �߻� �� �޼���//�̻��� ������ ���� �Է�
	{
		cerr << "������ �ҷ����� ���߽��ϴ�.";
	}
	else
	{
		sequence_file >> sequence;
	}

	//sequence = "AGTTTTAAAGGGCCCTTTAAAGGGCCCGTAAAAAGCCAAAAAATGGGGAAAGGGAAACCCAAAGGGTGATAATAGAAAAGATAATAGTGGGTTTCCCGGGGGGAAAGGGGTAAATGGGGAAAGGGAAACCCAAAGGGTGAAAATGGGGAAAGGGAAACCCAAAGGGTGAAGTTTTAAAGGGCCCTTTAAAGGGCCCGTAAAGTTTTAAAGGGCCCTTTAAAGGGCCCGTAAA";

	//ORF ���� �Է�
	cout << "ã�����ϴ� ORF�� ũ�� : " << endl;
	//cin >> user.length1 >> user.length2;//ORF�� ������ �Է�
	user.length1 = 1;
	user.length2 = 100;
	//length1�� �� ū ���� �Է��� ���, ���� �ٲ��ش�.
	if (user.length2 < user.length1)
	{
		int temp_length = user.length1;
		user.length1 = user.length2;
		user.length2 = temp_length;
	}

	user.SetName(name);//��ü���� �̸��ʱ�ȭ
	user.SetSequence(sequence);//��ü���� ���� �ʱ�ȭ
	//������ ������ ���
	cout << "�̸� :" << user.GetName() << endl;
	cout << "���� : " << user.GetSequence() << endl;
	cout << "���� : " << user.length1 << "," << user.length2 << endl;
	//6������ ������ readingframe���� ������ �����ϴ� �Լ�
	user.FrameSetting();
	//������ �Ǿ����� Ȯ���ϰ� ���� ����ϴ� �ڵ�
	cout << "dna2: " << user.GetDna2() << endl;
	cout << "dna3: " << user.GetDna3() << endl;
	cout << "dna4:" << user.GetRDna1() << endl;
	cout << "dna5:" << user.GetRDna2() << endl;
	cout << "dna6:" << user.GetRDna3() << endl;
	cout << "dna1:" << user.GetDna1() << endl;


	//�� 6���� ������ ���� ���� ������� �м�. 6�� �ݺ�
	Orf user_seq1(user.GetDna1());//Orf��ü�� �м��� ���� ���޹� Orf ��ü ����
	user_seq1.TransferSeq();//������ ��Ʈ�� ����ȭ
	cout << "����ȭdna1: ";
	for (int i = 0; i < user_seq1.orf1.size(); i++)//��Ʈ�� ����ȭ�� ������ ���.
	{
		cout << user_seq1.orf1[i] << " ";
	}
	cout << endl;

	user_seq1.IndexFinder();//start codon�� stop codon�� �ش��ϴ� �ε��� ã��
	user_seq1.OrfFinder(user);//orf ã��

	//�Ϸ��� orf ���
	if (user_seq1.complete_orf.empty())//���ͺ�������� ���� �ȳ� ���
	{
		cout << "ORF�� Ž������ �ʾҽ��ϴ�.";
	}
	else//���Ͱ� ���ִٸ� ORF�� ���
	{//��Ʈ�� ���ߺ��Ϳ� �����ϴ� ������ ���.
		cout << "ã�� ORF: ";
		for (int k = 0; k < user_seq1.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq1.complete_orf[k].size(); i++)
			{

				cout << user_seq1.complete_orf[k][i] << " ";
			}cout << endl << "          ";
		}
	}cout << endl << endl;

	Orf user_seq2(user.GetDna2());//Orf��ü�� �м��� ���� ����
	user_seq2.TransferSeq();//������ ��Ʈ�� ����ȭ
	cout << "����ȭdna2: ";
	for (int i = 0; i < user_seq2.orf1.size(); i++)//��Ʈ�� ����ȭ�� ������ ���.
	{
		cout << user_seq2.orf1[i] << " ";
	}
	cout << endl;

	user_seq2.IndexFinder();//start codon�� stop codon�� �ش��ϴ� �ε��� ã��
	user_seq2.OrfFinder(user);//orf ã��

	if (user_seq2.complete_orf.empty())
	{
		cout << "ORF�� Ž������ �ʾҽ��ϴ�.";
	}
	else//���Ͱ� ���ִٸ� ORF�� ���
	{
		cout << "ã�� ORF: ";
		for (int k = 0; k < user_seq2.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq2.complete_orf[k].size(); i++)
			{

				cout << user_seq2.complete_orf[k][i] << " ";
			}cout << endl << "          ";
		}
	}cout << endl << endl;


	Orf user_seq3(user.GetDna3());//Orf��ü�� �м��� ���� ����
	user_seq3.TransferSeq();//������ ��Ʈ�� ����ȭ
	cout << "����ȭdna3: ";
	for (int i = 0; i < user_seq3.orf1.size(); i++)//��Ʈ�� ����ȭ�� ������ ���.
	{
		cout << user_seq3.orf1[i] << " ";
	}
	cout << endl;

	user_seq3.IndexFinder();//start codon�� stop codon�� �ش��ϴ� �ε��� ã��
	user_seq3.OrfFinder(user);//orf ã��

	if (user_seq3.complete_orf.empty())
	{
		cout << "ORF�� Ž������ �ʾҽ��ϴ�.";
	}
	else//���Ͱ� ���ִٸ� ORF�� ���
	{
		cout << "ã�� ORF: ";
		for (int k = 0; k < user_seq3.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq3.complete_orf[k].size(); i++)
			{

				cout << user_seq3.complete_orf[k][i] << " ";
			}cout << endl << "          ";
		}
	}cout << endl << endl;


	Orf user_seq4(user.GetRDna1());//Orf��ü�� �м��� ���� ����
	user_seq4.TransferSeq();//������ ��Ʈ�� ����ȭ
	cout << "����ȭdna4: ";
	for (int i = 0; i < user_seq4.orf1.size(); i++)//��Ʈ�� ����ȭ�� ������ ���.
	{
		cout << user_seq4.orf1[i] << " ";
	}
	cout << endl;
	user_seq4.IndexFinder();//start codon�� stop codon�� �ش��ϴ� �ε��� ã��
	user_seq4.OrfFinder(user);//orf ã��

	if (user_seq4.complete_orf.empty())
	{
		cout << "ORF�� Ž������ �ʾҽ��ϴ�.";
	}
	else//���Ͱ� ���ִٸ� ORF�� ���
	{
		cout << "ã�� ORF: ";
		for (int k = 0; k < user_seq4.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq4.complete_orf[k].size(); i++)
			{

				cout << user_seq4.complete_orf[k][i] << " ";
			}cout << endl << "          ";
		}
	}cout << endl << endl;


	Orf user_seq5(user.GetRDna2());//Orf��ü�� �м��� ���� ����
	user_seq5.TransferSeq();//������ ��Ʈ�� ����ȭ
	cout << "����ȭdna5: ";
	for (int i = 0; i < user_seq5.orf1.size(); i++)//��Ʈ�� ����ȭ�� ������ ���.
	{
		cout << user_seq5.orf1[i] << " ";
	}
	cout << endl;
	user_seq5.IndexFinder();//start codon�� stop codon�� �ش��ϴ� �ε��� ã��
	user_seq5.OrfFinder(user);//orf ã��

	if (user_seq5.complete_orf.empty())
	{
		cout << "ORF�� Ž������ �ʾҽ��ϴ�.";
	}
	else//���Ͱ� ���ִٸ� ORF�� ���
	{
		cout << "ã�� ORF: ";
		for (int k = 0; k < user_seq5.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq5.complete_orf[k].size(); i++)
			{

				cout << user_seq5.complete_orf[k][i] << " ";
			}cout << endl << "          ";
		}
	}cout << endl << endl;


	Orf user_seq6(user.GetRDna3());//Orf��ü�� �м��� ���� ����
	user_seq6.TransferSeq();//������ ��Ʈ�� ����ȭ
	cout << "����ȭdna6: ";
	for (int i = 0; i < user_seq6.orf1.size(); i++)//��Ʈ�� ����ȭ�� ������ ���.
	{
		cout << user_seq6.orf1[i] << " ";
	}
	cout << endl;
	user_seq6.IndexFinder();//start codon�� stop codon�� �ش��ϴ� �ε��� ã��
	user_seq6.OrfFinder(user);//orf ã��

	if (user_seq6.complete_orf.empty())
	{
		cout << "ORF�� Ž������ �ʾҽ��ϴ�.";
	}
	else//���Ͱ� ���ִٸ� ORF�� ���
	{
		cout << "ã�� ORF: ";
		for (int k = 0; k < user_seq6.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq6.complete_orf[k].size(); i++)
			{

				cout << user_seq6.complete_orf[k][i] << " ";
			}cout << endl << "          ";
		}
	}cout << endl << endl;


	//���˿��� �ε��� �׽�Ʈ
	for(int i=0;i<user_seq1.complete_index.size();i++)
	{
		cout << user_seq1.complete_index[i].case_num << " ";
		cout << user_seq1.complete_index[i].start_index << " ";
		cout << user_seq1.complete_index[i].stop_index<< " ";
	}

	//��Ʈ�� ���� �׽�Ʈ
	user_seq1.IntronFinder();
	if (user_seq1.intron_removed.empty())//���ͺ�������� ���� �ȳ� ���
	{
		cout << "ORF�� Ž������ �ʾҽ��ϴ�.";
	}
	else//���Ͱ� ���ִٸ� ORF�� ���
	{//��Ʈ�� ���ߺ��Ϳ� �����ϴ� ������ ���.
		cout << "ã�� ORF: ";
		for (int k = 0; k < user_seq1.intron_removed.size(); k++)
		{
			for (int i = 0; i < user_seq1.intron_removed[k].size(); i++)
			{

				cout << user_seq1.intron_removed[k][i] << " ";
			}cout << endl << "          ";
		}
	}cout << endl << endl;

	return 0;
}