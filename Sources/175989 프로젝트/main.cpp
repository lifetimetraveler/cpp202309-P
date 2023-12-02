#include"user.h"

using namespace std;



int main() {
	string name;//유저의 이름
	string sequence;//유저의 서열

	User user;//User의 이름, 서열, ORF범위, 다른 reading frame서열을
	//저장하는 클래스

//이름입력, 테스트를 위해 임의의 이름과 서열, 범위를 사용했습니다.
	cout << "이름: " << endl;
	//cin >> name;
	name = "하하호호";

	//서열 입력
	cout << "DNA서열:" << endl;
	//cin >> sequence;
	sequence = "AGTTTTAAAGGGCCCTTTAAAGGGCCCGTAAAAAGCCAAAAAATGGGGAAAGGGAAACCCAAAGGGTGATAATAGAAAAGATAATAGTGGGTTTCCCGGGGGGAAAGGGGTAAATGGGGAAAGGGAAACCCAAAGGGTGAAAATGGGGAAAGGGAAACCCAAAGGGTGAAGTTTTAAAGGGCCCTTTAAAGGGCCCGTAAAGTTTTAAAGGGCCCTTTAAAGGGCCCGTAAA";

	//ORF 범위 입력
	cout << "찾고자하는 ORF의 크기 : " << endl;
	//cin >> user.length1 >> user.length2;//ORF의 범위를 입력
	user.length1 = 1;
	user.length2 = 100;
	//length1에 더 큰 수를 입력한 경우, 둘을 바꿔준다.
	if (user.length2 < user.length1)
	{
		int temp_length = user.length1;
		user.length1 = user.length2;
		user.length2 = temp_length;
	}

	user.SetName(name);//객체내의 이름초기화
	user.SetSequence(sequence);//객체내의 서열 초기화
	//저장한 값들을 출력
	cout << "이름 :" << user.GetName() << endl;
	cout << "서열 : " << user.GetSequence() << endl;
	cout << "범위 : " << user.length1 << "," << user.length2 << endl;
	//6가지의 가능한 readingframe으로 서열을 저장하는 함수
	user.FrameSetting();
	//저장이 되었는지 확인하가 위해 출력하는 코드
	cout << "dna2: " << user.GetDna2() << endl;
	cout << "dna3: " << user.GetDna3() << endl;
	cout << "dna4:" << user.GetRDna1() << endl;
	cout << "dna5:" << user.GetRDna2() << endl;
	cout << "dna6:" << user.GetRDna3() << endl;
	cout << "dna1:" << user.GetDna1() << endl;


	//각 6개의 서열에 대해 같은 방법으로 분석. 6번 반복
	Orf user_seq1(user.GetDna1());//Orf객체에 분석할 서열 전달및 Orf 객체 생성
	user_seq1.TransferSeq();//서열의 스트링 벡터화
	cout << "벡터화dna1: ";
	for (int i = 0; i < user_seq1.orf1.size(); i++)//스트링 벡터화된 서열을 출력.
	{
		cout << user_seq1.orf1[i] << " ";
	}
	cout << endl;

	user_seq1.IndexFinder();//start codon과 stop codon에 해당하는 인덱스 찾기
	user_seq1.OrfFinder(user);//orf 찾기

	//완료한 orf 출력
	if (user_seq1.complete_orf.empty())//벡터비어있으면 관련 안내 출력
	{
		cout << "ORF가 탐색되지 않았습니다.";
	}
	else//벡터가 차있다면 ORF를 출력
	{//스트링 이중벡터에 존재하는 모든것을 출력.
		cout << "찾은 ORF: ";
		for (int k = 0; k < user_seq1.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq1.complete_orf[k].size(); i++)
			{

				cout << user_seq1.complete_orf[k][i] << " ";
			}cout << endl << "          ";
		}
	}cout << endl << endl;

	Orf user_seq2(user.GetDna2());//Orf객체에 분석할 서열 전달
	user_seq2.TransferSeq();//서열의 스트링 벡터화
	cout << "벡터화dna2: ";
	for (int i = 0; i < user_seq2.orf1.size(); i++)//스트링 벡터화된 서열을 출력.
	{
		cout << user_seq2.orf1[i] << " ";
	}
	cout << endl;

	user_seq2.IndexFinder();//start codon과 stop codon에 해당하는 인덱스 찾기
	user_seq2.OrfFinder(user);//orf 찾기

	if (user_seq2.complete_orf.empty())
	{
		cout << "ORF가 탐색되지 않았습니다.";
	}
	else//벡터가 차있다면 ORF를 출력
	{
		cout << "찾은 ORF: ";
		for (int k = 0; k < user_seq2.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq2.complete_orf[k].size(); i++)
			{

				cout << user_seq2.complete_orf[k][i] << " ";
			}cout << endl << "          ";
		}
	}cout << endl << endl;


	Orf user_seq3(user.GetDna3());//Orf객체에 분석할 서열 전달
	user_seq3.TransferSeq();//서열의 스트링 벡터화
	cout << "벡터화dna3: ";
	for (int i = 0; i < user_seq3.orf1.size(); i++)//스트링 벡터화된 서열을 출력.
	{
		cout << user_seq3.orf1[i] << " ";
	}
	cout << endl;

	user_seq3.IndexFinder();//start codon과 stop codon에 해당하는 인덱스 찾기
	user_seq3.OrfFinder(user);//orf 찾기

	if (user_seq3.complete_orf.empty())
	{
		cout << "ORF가 탐색되지 않았습니다.";
	}
	else//벡터가 차있다면 ORF를 출력
	{
		cout << "찾은 ORF: ";
		for (int k = 0; k < user_seq3.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq3.complete_orf[k].size(); i++)
			{

				cout << user_seq3.complete_orf[k][i] << " ";
			}cout << endl << "          ";
		}
	}cout << endl << endl;


	Orf user_seq4(user.GetRDna1());//Orf객체에 분석할 서열 전달
	user_seq4.TransferSeq();//서열의 스트링 벡터화
	cout << "벡터화dna4: ";
	for (int i = 0; i < user_seq4.orf1.size(); i++)//스트링 벡터화된 서열을 출력.
	{
		cout << user_seq4.orf1[i] << " ";
	}
	cout << endl;
	user_seq4.IndexFinder();//start codon과 stop codon에 해당하는 인덱스 찾기
	user_seq4.OrfFinder(user);//orf 찾기

	if (user_seq4.complete_orf.empty())
	{
		cout << "ORF가 탐색되지 않았습니다.";
	}
	else//벡터가 차있다면 ORF를 출력
	{
		cout << "찾은 ORF: ";
		for (int k = 0; k < user_seq4.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq4.complete_orf[k].size(); i++)
			{

				cout << user_seq4.complete_orf[k][i] << " ";
			}cout << endl << "          ";
		}
	}cout << endl << endl;


	Orf user_seq5(user.GetRDna2());//Orf객체에 분석할 서열 전달
	user_seq5.TransferSeq();//서열의 스트링 벡터화
	cout << "벡터화dna5: ";
	for (int i = 0; i < user_seq5.orf1.size(); i++)//스트링 벡터화된 서열을 출력.
	{
		cout << user_seq5.orf1[i] << " ";
	}
	cout << endl;
	user_seq5.IndexFinder();//start codon과 stop codon에 해당하는 인덱스 찾기
	user_seq5.OrfFinder(user);//orf 찾기

	if (user_seq5.complete_orf.empty())
	{
		cout << "ORF가 탐색되지 않았습니다.";
	}
	else//벡터가 차있다면 ORF를 출력
	{
		cout << "찾은 ORF: ";
		for (int k = 0; k < user_seq5.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq5.complete_orf[k].size(); i++)
			{

				cout << user_seq5.complete_orf[k][i] << " ";
			}cout << endl << "          ";
		}
	}cout << endl << endl;


	Orf user_seq6(user.GetRDna3());//Orf객체에 분석할 서열 전달
	user_seq6.TransferSeq();//서열의 스트링 벡터화
	cout << "벡터화dna6: ";
	for (int i = 0; i < user_seq6.orf1.size(); i++)//스트링 벡터화된 서열을 출력.
	{
		cout << user_seq6.orf1[i] << " ";
	}
	cout << endl;
	user_seq6.IndexFinder();//start codon과 stop codon에 해당하는 인덱스 찾기
	user_seq6.OrfFinder(user);//orf 찾기

	if (user_seq6.complete_orf.empty())
	{
		cout << "ORF가 탐색되지 않았습니다.";
	}
	else//벡터가 차있다면 ORF를 출력
	{
		cout << "찾은 ORF: ";
		for (int k = 0; k < user_seq6.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq6.complete_orf[k].size(); i++)
			{

				cout << user_seq6.complete_orf[k][i] << " ";
			}cout << endl << "          ";
		}
	}cout << endl << endl;





	return 0;
}