#include"user.h"
using namespace std;



void User::SetName(string name)//사용자 이름초기화
{
	user_name = name;
}
void User::SetSequence(string sequence)//사용자 DNA서열초기화
{
	user_dna_sequence = sequence;
}
void User::FrameSetting()//reading frame을 바꿔서 서열을 저장하는 함수
{
	string temp_dna=user_dna_sequence;//DNA서열의 frame 별 저장을 위한 임시변수;
	

	dna1 = user_dna_sequence;//유저의 정방향, reading frame DNA 저장

	//맨 앞 하나를 지워서 frameshift
	dna2 = temp_dna.erase(0,1);//유저의 정방향, +1 reding frame DNA 저장

	//맨 앞 둘을 지워서 frameshift
	dna3 = temp_dna.erase(0,1);//유저의 정방향, +2 reding frame DNA 저장

	for (int i=(int)user_dna_sequence.length()-1;i>=0; i--)
	{
		reverse_dna1 += user_dna_sequence[i];//유저의 역방향, reading frame DNA 저장
	}

	temp_dna = reverse_dna1;//역방향을 위한 저장

	 //맨 앞 하나를 지워서 frameshift
	reverse_dna2 = temp_dna.erase(0, 1);//유저의 역방향, +1 reding frame DNA 저장

	 //맨 앞 둘을 지워서 frameshift
	reverse_dna3 = temp_dna.erase(0, 1);//유저의 역방향, +2 reding frame DNA 저장
}
string User::GetName() {//사용자 이름반환
	return user_name;
}
string User::GetSequence() {//사용자 DNA서열반환
	return user_dna_sequence;
}
//저장된 reading frame 변환된 서열을 반환하는 함수들
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
void Orf::TransferSeq()//즉, 스트링 벡터화하는 함수
{
	for (int i=0; i < original_seq.length();i+=3)
	{
		string triplet = original_seq.substr(i,3);//아래의 코드를 사용해도 된다.
		//to_string(original_seq[i])+to_string(original_seq[i+1])+to_string(original_seq[i+2])
		//하지만 이렇게 사용하려면 string의 길이보다 큰 인덱스 값을 사용하게 된다.
		//그 경우에는 3의 배수가 되도록 서열을 다듬어야 해서 편의상 string에 원래있는
		//스트링 클래스 멤버함수를 사용하였다.
		//substr은 범위를 벗어나면 더이상 반환하지 않는다.
		orf1.push_back(triplet);
	}
}
//가공한 서열의 start codon과 stop codon의 인덱스를 저장하는 함수
void Orf::IndexFinder()
{//ATG(start codon)인 스트링을 찾아서 그때의 인덱스 값을 벡터에 저장
	for (int i = 0; i < orf1.size(); i++)
	{
		if (orf1[i] == "ATG")
			atg_index.push_back(i);
	}
	//TGA(stop codon)인 스트링을 찾아서 그때의 인덱스 값을 벡터에 저장
	for (int i = 0; i < orf1.size(); i++)
	{
		if (orf1[i] == "TGA")
			tga_index.push_back(i);
	}
	//TAA(stop codon)인 스트링을 찾아서 그때의 인덱스 값을 벡터에 저장
	for (int i = 0; i < orf1.size(); i++)
	{
		if (orf1[i] == "TAA")
			taa_index.push_back(i);
	}
	//TAG(stop codon)인 스트링을 찾아서 그때의 인덱스 값을 벡터에 저장
	for (int i = 0; i < orf1.size(); i++)
	{
		if (orf1[i] == "TAG")
			tag_index.push_back(i);
	}
}

//저장된 인덱스의 값에 따라 시작(ATG)부터 끝(TAA,TAG,TGA)까지가 사용자의 범위에
//들어가면 그 서열을 스트링 이중벡터에 저장한다.
void Orf::OrfFinder(User user)
{
	int length1 = user.length1;//사용자가 입력한 범위를 가져온다.
	int length2 = user.length2;
	for (int i = 0; i < atg_index.size(); i++)//저장된 시작코돈이 있는 인덱스를 사용하기위함
	{//start codon은 동일하기 때문에 세개의 stop codon이 공유. 하위의 stop codon을 세개의 for문으로 탐색
		for (int j = 0; j < tga_index.size(); j++)//stop codon이 저장된 인덱스를 사용하기위함
		{//stop codon이 TGA인 경우 
			if (atg_index[i] < tga_index[j])//시작 코돈이 stop codon보다 앞에 있을 때
			{//시작부터 끝까지의 서열의 크기가 사용자가 입력한 범위와 일치하는 경우
				if (tga_index[j] - atg_index[i]+1 >= length1 && tga_index[j] - atg_index[i] +1<= length2)
				{//그 시작부터 끝까지(ORF)를 스트링 이중 벡터에 저장(push_back())한다.
					vector<string> temp_com_orf(orf1.begin() + atg_index[i], orf1.begin() + tga_index[j]+1);
					complete_orf.push_back(temp_com_orf);
					
					//ORF가 추출 될 때마다의 인덱스 값을 저장
					SavedIndex temp_saved(0, atg_index[i], tga_index[j]);
					complete_index.push_back(temp_saved);
						
				}
			}
		}
		for (int j = 0; j < taa_index.size(); j++)//stop codon의 인덱스를 사용하기 위함
		{//위의 for문과 동일한 원리, 하지만 stop codon이 TAA인 경우
			if (atg_index[i] < taa_index[j])
			{
				if (taa_index[j] - atg_index[i] + 1 >= length1 && taa_index[j] - atg_index[i] + 1 <= length2)
				{
					vector<string> temp_com_orf(orf1.begin() + atg_index[i], orf1.begin() + taa_index[j] + 1);
					complete_orf.push_back(temp_com_orf);

					SavedIndex temp_saved(1, atg_index[i], taa_index[j]);
					complete_index.push_back(temp_saved);
				}
			}
		}
		for (int j = 0; j < tag_index.size(); j++)//stop codon의 인덱스를 사용하기 위함
		{//위의 for문과 동일한 원리, 하지만 stop codon이 TAG인 경우
			if (atg_index[i] < tag_index[j])
			{
				if (tag_index[j] - atg_index[i] + 1 >= length1 && tag_index[j] - atg_index[i] + 1 <= length2)
				{
					vector<string> temp_com_orf(orf1.begin() + atg_index[i], orf1.begin() + tag_index[j] + 1);
					complete_orf.push_back(temp_com_orf);

					SavedIndex temp_saved(2, atg_index[i], tag_index[j]);
					complete_index.push_back(temp_saved);
				}
			}
		}
	}
}
//인트론을 제거하는 함수
void Orf::IntronFinder()
{
	intron_removed = complete_orf;//찾은 orf를 인트론이 잘린 변수가 담길 곳에 복사

	for (int a = 0; a < intron_removed.size(); a++)//찾은 ORF의 수만큼 반복
	{

		for(int b=0;b<intron_removed[a].size()-5;b++)//끝까지 못가게 하여 없는 인덱스를 넣지 않도록
		{//전체 서열을 탐색
            //탐색하다가 인트론의 시작 부분을 만나면
			if (intron_removed[a][b] == "GTA" && intron_removed[a][b+1] == "AGT")
			{
				for (int c = b + 1; c < intron_removed[a].size(); c++)//그 이후의 서열을 탐색
				{   //인트론의 끝부분에 해당하는 부분을 찾으면
					if(intron_removed[a][c-3]=="TTT"&&intron_removed[a][c-2]=="TTT")
					{
						if (intron_removed[a][c - 1][0] == 'T' && intron_removed[a][c - 1][1] == 'T')
						{
							if (intron_removed[a][c] == "CAG")
							{
								for (int i = b; i <= c; i++)//찾아낸 인트론의 처음과 끝을,
								{
									intron_removed[a][i] = "intron";//해당하는 범위를 모두 인트론으로 바꿈
								}
								break;//완료한 후 for문을 탈출 하여 가장 처음으로 만난 인트론의 끝부분까지만
								//인트론으로 바꿈.
							}
						}
					}
				}
			}
		}
	}
	
}

void Orf::CodonDecipher()
{
	protein = intron_removed;
	for (int a = 0; a < protein.size(); a++)
	{
		for(int b=0;b<protein[a].size();b++)//각각 저장된 orf를 탐색
		{
			string temp_codon = protein[a][b];
			if (temp_codon == "UUU" || "UUC")
				protein[a][b] = " F ";
			else if(temp_codon == "UUA" || "UUG"||"CUU"||"CUC"||"CUA"||"CUG")
				protein[a][b] = " L ";
			else if (temp_codon == "AUU" || "AUC" || "AUA")
				protein[a][b] = " L ";
			else if (temp_codon == "AUG")
				protein[a][b] = " M ";
			else if (temp_codon == "GUU" || "GUC" || "GUA" || "GUG")
				protein[a][b] = " V ";
			else if (temp_codon == "UCU" || "UCC" || "UCA" || "UCG" || "AGU" || "AGC")
				protein[a][b] = " S ";
			else if (temp_codon == "CCU" || temp_codon == "CCC" || temp_codon == "CCA" || temp_codon == "CCG")
				protein[a][b] = " P ";
			else if (temp_codon == "ACU" || temp_codon == "ACC" || temp_codon == "ACA" || temp_codon == "ACG")
				protein[a][b] = " T ";
			else if (temp_codon == "GCU" || temp_codon == "GCC" || temp_codon == "GCA" || temp_codon == "GCG")
				protein[a][b] = " A ";
			else if (temp_codon == "UAU" || temp_codon == "UAC")
				protein[a][b] = " Y ";
			else if (temp_codon == "CAU" || temp_codon == "CAC")
				protein[a][b] = " H ";
			else if (temp_codon == "CAA" || temp_codon == "CAG")
				protein[a][b] = " Q ";
			else if (temp_codon == "AAU" || temp_codon == "AAC")
				protein[a][b] = " N ";
			else if (temp_codon == "AAA" || temp_codon == "AAG")
				protein[a][b] = " K ";
			else if (temp_codon == "GAU" || temp_codon == "GAC")
				protein[a][b] = " D ";
			else if (temp_codon == "GAA" || temp_codon == "GAG")
				protein[a][b] = " E ";
			else if (temp_codon == "UGU" || temp_codon == "UGC")
				protein[a][b] = " C ";
			else if (temp_codon == "UGG")
				protein[a][b] = " W ";
			else if (temp_codon == "CGU" || temp_codon == "CGC" || temp_codon == "CGA" || temp_codon == "CGG" || temp_codon == "AGA" || temp_codon == "AGG")
				protein[a][b] = " R ";
			else if (temp_codon == "GGU" || temp_codon == "GGC" || temp_codon == "GGA" || temp_codon == "GGG")
				protein[a][b] = " G ";
			else if (temp_codon == "UAA" || temp_codon == "UAG" || temp_codon == "UGA")
				protein[a][b] = "stop";

		}
	}
}
