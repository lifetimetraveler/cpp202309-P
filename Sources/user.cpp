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



void Orf::TransferSeq(User user)
{
	string temp_orf = user.GetDna1();
	for (int i=0; i < temp_orf.length();i+=3)
	{
		string triplet = temp_orf.substr(i,3);
		orf1.push_back(triplet);
	}
}
void Orf::OrfFinder()
{
	string 
}

