#include"user.h"
using namespace std;



void User::SetName(string name)//����� �̸��ʱ�ȭ
{
	user_name = name;
}
void User::SetSequence(string sequence)//����� DNA�����ʱ�ȭ
{
	user_dna_sequence = sequence;
}
void User::FrameSetting()//reading frame�� �ٲ㼭 ������ �����ϴ� �Լ�
{
	string temp_dna=user_dna_sequence;//DNA������ frame �� ������ ���� �ӽú���;
	

	dna1 = user_dna_sequence;//������ ������, reading frame DNA ����

	//�� �� �ϳ��� ������ frameshift
	dna2 = temp_dna.erase(0,1);//������ ������, +1 reding frame DNA ����

	//�� �� ���� ������ frameshift
	dna3 = temp_dna.erase(0,1);//������ ������, +2 reding frame DNA ����

	for (int i=user_dna_sequence.length()-1;i>=0; i--)
	{
		reverse_dna1 += user_dna_sequence[i];//������ ������, reading frame DNA ����
	}

	temp_dna = reverse_dna1;//�������� ���� ����

	 //�� �� �ϳ��� ������ frameshift
	reverse_dna2 = temp_dna.erase(0, 1);//������ ������, +1 reding frame DNA ����

	 //�� �� ���� ������ frameshift
	reverse_dna3 = temp_dna.erase(0, 1);//������ ������, +2 reding frame DNA ����


}
string User::GetName() {//����� �̸���ȯ
	return user_name;
}
string User::GetSequence() {//����� DNA������ȯ
	return user_dna_sequence;
}
//����� reading frame ��ȯ�� ������ ��ȯ�ϴ� �Լ���
string User::GetDna1() {
	return dna1;
}
string User::GetDna2() {
	return dna2;
}
string User::GetDna3() {
	return dna3;
}
string User::GetRDna1() {
	return reverse_dna1;
}
string User::GetRDna2() {
	return reverse_dna2;
}
string User::GetRDna3() {
	return reverse_dna3;
}


//��Ʈ�������� ������ ������ ������ ��Ʈ�� ���Ϳ� �����Ͽ�
//���Ŀ� ������ �����Ѵ�.
void Orf::TransferSeq()//��, ��Ʈ�� ����ȭ�ϴ� �Լ�
{
	for (int i=0; i < original_seq.length();i+=3)
	{
		string triplet = original_seq.substr(i,3);//�Ʒ��� �ڵ带 ����ص� �ȴ�.
		//to_string(original_seq[i])+to_string(original_seq[i+1])+to_string(original_seq[i+2])
		//������ �̷��� ����Ϸ��� string�� ���̺��� ū �ε��� ���� ����ϰ� �ȴ�.
		//�� ��쿡�� 3�� ����� �ǵ��� ������ �ٵ��� �ؼ� ���ǻ� string�� �����ִ�
		//��Ʈ�� Ŭ���� ����Լ��� ����Ͽ���.
		//substr�� ������ ����� ���̻� ��ȯ���� �ʴ´�.
		orf1.push_back(triplet);
	}
}
//������ ������ start codon�� stop codon�� �ε����� �����ϴ� �Լ�
void Orf::IndexFinder()
{//ATG(start codon)�� ��Ʈ���� ã�Ƽ� �׶��� �ε��� ���� ���Ϳ� ����
	for (int i = 0; i < orf1.size(); i++)
	{
		if (orf1[i] == "ATG")
			atg_index.push_back(i);
	}
	//TGA(stop codon)�� ��Ʈ���� ã�Ƽ� �׶��� �ε��� ���� ���Ϳ� ����
	for (int i = 0; i < orf1.size(); i++)
	{
		if (orf1[i] == "TGA")
			tga_index.push_back(i);
	}
	//TAA(stop codon)�� ��Ʈ���� ã�Ƽ� �׶��� �ε��� ���� ���Ϳ� ����
	for (int i = 0; i < orf1.size(); i++)
	{
		if (orf1[i] == "TAA")
			taa_index.push_back(i);
	}
	//TAG(stop codon)�� ��Ʈ���� ã�Ƽ� �׶��� �ε��� ���� ���Ϳ� ����
	for (int i = 0; i < orf1.size(); i++)
	{
		if (orf1[i] == "TAG")
			tag_index.push_back(i);
	}
}

//����� �ε����� ���� ���� ����(ATG)���� ��(TAA,TAG,TGA)������ ������� ������
//���� �� ������ ��Ʈ�� ���ߺ��Ϳ� �����Ѵ�.
void Orf::OrfFinder(User user)
{
	int length1 = user.length1;//����ڰ� �Է��� ������ �����´�.
	int length2 = user.length2;
	for (int i = 0; i < atg_index.size(); i++)//����� �����ڵ��� �ִ� �ε����� ����ϱ�����
	{//start codon�� �����ϱ� ������ ������ stop codon�� ����. ������ stop codon�� ������ for������ Ž��
		for (int j = 0; j < tga_index.size(); j++)//stop codon�� ����� �ε����� ����ϱ�����
		{//stop codon�� TGA�� ��� 
			if (atg_index[i] < tga_index[j])//���� �ڵ��� stop codon���� �տ� ���� ��
			{//���ۺ��� �������� ������ ũ�Ⱑ ����ڰ� �Է��� ������ ��ġ�ϴ� ���
				if (tga_index[j] - atg_index[i]+1 >= length1 && tga_index[j] - atg_index[i] +1<= length2)
				{//�� ���ۺ��� ������(ORF)�� ��Ʈ�� ���� ���Ϳ� ����(push_back())�Ѵ�.
					vector<string> temp_com_orf(orf1.begin() + atg_index[i], orf1.begin() + tga_index[j]+1);
					complete_orf.push_back(temp_com_orf);
				}
			}
		}
		for (int j = 0; j < taa_index.size(); j++)//stop codon�� �ε����� ����ϱ� ����
		{//���� for���� ������ ����, ������ stop codon�� TAA�� ���
			if (atg_index[i] < taa_index[j])
			{
				if (taa_index[j] - atg_index[i] + 1 >= length1 && taa_index[j] - atg_index[i] + 1 <= length2)
				{
					vector<string> temp_com_orf(orf1.begin() + atg_index[i], orf1.begin() + taa_index[j] + 1);
					complete_orf.push_back(temp_com_orf);
				}
			}
		}
		for (int j = 0; j < tag_index.size(); j++)//stop codon�� �ε����� ����ϱ� ����
		{//���� for���� ������ ����, ������ stop codon�� TAG�� ���
			if (atg_index[i] < tag_index[j])
			{
				if (tag_index[j] - atg_index[i] + 1 >= length1 && tag_index[j] - atg_index[i] + 1 <= length2)
				{
					vector<string> temp_com_orf(orf1.begin() + atg_index[i], orf1.begin() + tag_index[j] + 1);
					complete_orf.push_back(temp_com_orf);
				}
			}
		}
	}
}

