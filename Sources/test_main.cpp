#include"user.h"
using namespace std;



int main() {
	string name;//유저의 이름
	string sequence;//유저의 서열
	
	User user;//User의 이름, 서열, ORF범위, 다른 reading frame서열을
	          //저장하는 클래스
	
	cout << "이름";
	//name = "하하호호";
	cout << "DNA서열:";
	cin >> name;
	//sequence = "AGCCAAAAAATGGGGAAAGGGAAACCCAAAGGGTGATTTGGGAAA";
	cin >> sequence;
	cout << "찾고자하는 ORF의 크기 : ";
	cin >> user.length1 >> user.length2;//ORF의 범위를 입력
	//user.length1 = 1;
	//user.length2 = 100;

	user.SetName(name);//객체내의 이름초기화
	user.SetSequence(sequence);//객체내의 서열 초기화

	cout << "이름:" << user.GetName() << endl;
	cout << "서열: " << user.GetSequence() << endl;

	user.FrameSetting();
	cout << "dna2: " << user.GetDna2() << endl;
	cout << "dna3: " << user.GetDna3() << endl;
	cout << "dna4:" << user.GetRDna1() << endl;
	cout << "dna5:" << user.GetRDna2() << endl;
	cout << "dna6:" << user.GetRDna3() << endl;
	cout << "dna1:" << user.GetDna1() << endl;

	Orf user_seq1(user.GetDna1());//Orf객체에 분석할 서열 전달
	user_seq1.TransferSeq();//서열의 스트링 벡터화

	
	cout << "벡터" << endl;
	for (int i = 0; i < user_seq1.orf1.size(); i++)
	{
		cout<<user_seq1.orf1[i]<<" ";
	}
	cout << endl;
	user_seq1.IndexFinder();//인덱스 찾기
	user_seq1.OrfFinder(user);//orf 찾기


	
	//완료한 orf 출력
	for (int i = 0; i < user_seq1.complete_orf[0].size(); i++)
	{
		cout << user_seq1.complete_orf[0][i] << " ";
	}



	return 0;
}