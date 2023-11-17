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
	string temp_dna=user_dna_sequence;//DNA서열의 frame 별 저장을 위한 임시변수;
	

	dna1 = user_dna_sequence;//유저의 정방향, reading frame DNA 저장

	//맨 앞 하나를 지워서 frameshift
	dna2 = temp_dna.erase(0,1);//유저의 정방향, +1 reding frame DNA 저장

	//맨 앞 둘을 지워서 frameshift
	dna3 = temp_dna.erase(0,1);//유저의 정방향, +2 reding frame DNA 저장

	for (int i=user_dna_sequence.length()-1;i>=0; i--)
	{
		reverse_dna1 += user_dna_sequence[i];//유저의 역방향, reading frame DNA 저장
	}

	temp_dna = reverse_dna1;//역방향을 위한 저장

	 //맨 앞 하나를 지워서 frameshift
	reverse_dna2 = temp_dna.erase(0, 1);//유저의 역방향, +1 reding frame DNA 저장

	 //맨 앞 둘을 지워서 frameshift
	reverse_dna3 = temp_dna.erase(0, 1);//유저의 역방향, +2 reding frame DNA 저장


}
string User::GetName() {
	return user_name;
}
string User::GetSequence() {
	return user_dna_sequence;
}

