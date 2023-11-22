#include<iostream>
#include<string>
#include<vector>
using namespace std;

class User {
	string user_name;//유저의 이름
	string user_dna_sequence;//유저의 DNA서열을 받는 변수
    string dna1;
	string dna2;
	string dna3;
	string reverse_dna1;
	string reverse_dna2;
	string reverse_dna3;
	//dna1~6은 DNA 서열을 앞뒤로 읽은 경우(2)*DNA의 맨 앞 서열을 무엇으로하냐 (3)
	// DNA의 codon은 세 DNA서열이 묶인것이기때문에
	//세가지로 읽을 수 있다. 
	
public:
	int length1=0;
	int length2=0;


	void SetName(string name);//이름을 초기화하는 함수
	void SetSequence(string sequence);//서열을 초기화하는 함수
	void FrameSetting();//codon의 frame을 설정하는 함수
	string GetName();//이름을 반환하는 함수
	string GetSequence();//서열을 반환하는 함수
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