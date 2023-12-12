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

	for (int i=(int)user_dna_sequence.length()-1;i>=0; i--)
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
					
					//ORF�� ���� �� �������� �ε��� ���� ����
					SavedIndex temp_saved(0, atg_index[i], tga_index[j]);
					complete_index.push_back(temp_saved);
						
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

					SavedIndex temp_saved(1, atg_index[i], taa_index[j]);
					complete_index.push_back(temp_saved);
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

					SavedIndex temp_saved(2, atg_index[i], tag_index[j]);
					complete_index.push_back(temp_saved);
				}
			}
		}
	}
}
//��Ʈ���� �����ϴ� �Լ�
void Orf::IntronFinder()
{
	intron_removed = complete_orf;//ã�� orf�� ��Ʈ���� �߸� ������ ��� ���� ����

	for (int a = 0; a < intron_removed.size(); a++)//ã�� ORF�� ����ŭ �ݺ�
	{

		for(int b=0;b<intron_removed[a].size()-5;b++)//������ ������ �Ͽ� ���� �ε����� ���� �ʵ���
		{//��ü ������ Ž��
            //Ž���ϴٰ� ��Ʈ���� ���� �κ��� ������
			if (intron_removed[a][b] == "GTA" && intron_removed[a][b+1] == "AGT")
			{
				for (int c = b + 1; c < intron_removed[a].size(); c++)//�� ������ ������ Ž��
				{   //��Ʈ���� ���κп� �ش��ϴ� �κ��� ã����
					if(intron_removed[a][c-3]=="TTT"&&intron_removed[a][c-2]=="TTT")
					{
						if (intron_removed[a][c - 1][0] == 'T' && intron_removed[a][c - 1][1] == 'T')
						{
							if (intron_removed[a][c] == "CAG")
							{
								for (int i = b; i <= c; i++)//ã�Ƴ� ��Ʈ���� ó���� ����,
								{
									intron_removed[a][i] = "intron";//�ش��ϴ� ������ ��� ��Ʈ������ �ٲ�
								}
								break;//�Ϸ��� �� for���� Ż�� �Ͽ� ���� ó������ ���� ��Ʈ���� ���κб�����
								//��Ʈ������ �ٲ�.
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
		for(int b=0;b<protein[a].size();b++)//���� ����� orf�� Ž��
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
