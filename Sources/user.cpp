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


//스트링형태의 서열을 세개씩 나눠서 스트링 벡터에 저장하여
//추후에 가공을 쉽게한다.
void Orf::TransferSeq()//
{
	for (int i=0; i < original_seq.length();i+=3)
	{
		string triplet = original_seq.substr(i,3);
		orf1.push_back(triplet);
	}
}
//가공한 서열의 start codon과 stop codon의 인덱스를 저장하는 함수
void Orf::IndexFinder()
{
	for (int i = 0; i < orf1.size(); i++)
	{
		if (orf1[i] == "ATG")
			atg_index.push_back(i);
	}
	for (int i = 0; i < orf1.size(); i++)
	{
		if (orf1[i] == "TGA")
			tga_index.push_back(i);
	}
	for (int i = 0; i < orf1.size(); i++)
	{
		if (orf1[i] == "TAA")
			taa_index.push_back(i);
	}
	for (int i = 0; i < orf1.size(); i++)
	{
		if (orf1[i] == "TAG")
			tag_index.push_back(i);
	}
}
void Orf::OrfFinder(User user)
{
	int length1 = user.length1;
	int length2 = user.length2;
	for (int i = 0; i < atg_index.size(); i++)
	{
		for (int j = 0; j < tga_index.size(); j++)
		{
			if (atg_index[i] < tga_index[j])
			{
				if (tga_index[j] - atg_index[i]+1 >= length1 && tga_index[j] - atg_index[i] +1<= length2)
				{
					vector<string> temp_com_orf(orf1.begin() + atg_index[i], orf1.begin() + tga_index[j]+1);
					complete_orf.push_back(temp_com_orf);
				}
			}
		}
		for (int j = 0; j < taa_index.size(); i++)
		{
			if (atg_index[i] < taa_index[j])
			{
				if (taa_index[j] - atg_index[i] + 1 >= length1 && taa_index[j] - atg_index[i] + 1 <= length2)
				{
					vector<string> temp_com_orf(orf1.begin() + atg_index[i], orf1.begin() + taa_index[j] + 1);
					complete_orf.push_back(temp_com_orf);
				}
			}
		}
		for (int j = 0; j < tag_index.size(); i++)
		{
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

