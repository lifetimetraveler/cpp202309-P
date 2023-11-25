#include<iostream>
#include<string>
#include<vector>
using namespace std;


//유저의 이름, 유저의 서열, 찾고자 하는 orf의 길이범위, reading frame을 바꾼 서열의
//  저장을 위한 클래스
class User {
	string user_name;//유저의 이름
	string user_dna_sequence;//유저의 DNA서열을 받는 변수
    string dna1;//정방향, 원래 reading frame
	string dna2;//정방향, 원래 reading frame+1
	string dna3;//정방향, 원래 reading frame+2
	string reverse_dna1;//역방향, reading frame
	string reverse_dna2;//역방향, reading frame+1
	string reverse_dna3;//역방향, reading frame+2
	//dna1~6은 DNA 서열을 앞뒤로 읽은 경우(2)*DNA의 맨 앞 서열을 무엇으로하냐 (3)
	// DNA의 codon은 세 DNA서열이 묶인것이기때문에
	// 세가지로 읽을 수 있다. 여기에 역방향일 경우를 모두 합하면 6가지 경우
	
public:
	int length1=0;//찾고자 하는 orf의 범위의 최소 값
	int length2=0;//찾고자 하는 orf의 범위의 최대 값


	void SetName(string name);//이름을 초기화하는 함수
	void SetSequence(string sequence);//서열을 초기화하는 함수
	void FrameSetting();//본래 서열의 reading frame을 달리하여 저장하는 함수
	string GetName();//이름을 반환하는 함수
	string GetSequence();//서열을 반환하는 함수
	string GetDna1();//dna1을 반환하는 함수
	string GetDna2();//dna2을 반환하는 함수
	string GetDna3();//dna3을 반환하는 함수
	string GetRDna1();//reaverse_dna1을 반환하는 함수
	string GetRDna2();//reaverse_dna2을 반환하는 함수
	string GetRDna3();//reaverse_dna3을 반환하는 함수
};

//하나의 reading frame에서 orf를 찾고 그것을 저장하는 클래스
class Orf
{
	

public:
string original_seq;//받아온 분석할 서열을 저장
	Orf(string original_seq) { this->original_seq = original_seq; }//생성자
	vector<string> orf1;//처리과정을 매끄럽게하기위해 스트링으로 세개씩 나누어서 저장
	vector<vector<string>> complete_orf;//orf를 추출한 벡터를 저장하기위한 이중 벡터

	//orf에서 atg(start codon)을 저장
	vector<int> atg_index;

    //orf에서 stop codon이 있는 배열을 저장
	vector<int> tga_index;
	vector<int> taa_index;
	vector<int> tag_index;

	//서열을 벡터화하는 함수
	void TransferSeq();
	//벡터화한 서열에서 stop, start codon을 찾는 함수
	void IndexFinder();
	//찾아낸 인덱스에 맞춰 서열을 저장하는 함수
	void OrfFinder(User user);
	
};