#include"user.h"
using namespace std;



void User::SetName(string name)
{
	user_name = name;
}
void User::SetSequence(string sequence)
{
	user_dna_sequence = sequence;
}
void User::FrameSetting()
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
string User::GetName() {
	return user_name;
}
string User::GetSequence() {
	return user_dna_sequence;
}

