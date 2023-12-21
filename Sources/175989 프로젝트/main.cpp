#include"user.h"
#include <fstream>
using namespace std;



int main() {
	string name;//유저의 이름
	string sequence;//유저의 서열
	ifstream sequence_file("sequence.txt");//유저의 서열 파일
	ofstream result("output.txt", ios::out | ios::trunc);//분석결과를 출력하기 위함.
	User user;//User의 이름, 서열, ORF범위, 다른 reading frame서열을
	//저장하는 클래스

//이름입력, 테스트를 위해 임의의 이름과 서열, 범위를 사용했습니다.
	cout << "이름: " << endl;
	//cin >> name;
	name = "하하호호";

	//서열 입력
	cout << "DNA서열:" << endl;
	if (!sequence_file)//오류 발생 시 메세지//이상이 없으면 서열 입력
	{
		cerr << "파일을 불러오지 못했습니다.";
	}
	else
	{
		sequence_file >> sequence;
	}
	sequence_file.close();

	//sequence = "AGTTTTAAAGGGCCCTTTAAAGGGCCCGTAAAAAGCCAAAAAATGGGGAAAGGGAAACCCAAAGGGTGATAATAGAAAAGATAATAGTGGGTTTCCCGGGGGGAAAGGGGTAAATGGGGAAAGGGAAACCCAAAGGGTGAAAATGGGGAAAGGGAAACCCAAAGGGTGAAGTTTTAAAGGGCCCTTTAAAGGGCCCGTAAAGTTTTAAAGGGCCCTTTAAAGGGCCCGTAAA";

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
	cout << "dna1:" << user.GetDna1() << endl << endl;


	//각 6개의 서열에 대해 같은 방법으로 분석. 6번 반복
	Orf user_seq1(user.GetDna1());//Orf객체에 분석할 서열 전달및 Orf 객체 생성
	user_seq1.TransferSeq();//서열의 스트링 벡터화
	cout << "벡터화dna1: " << endl;
	for (int i = 0; i < user_seq1.orf1.size(); i++)//스트링 벡터화된 서열을 출력.
	{
		cout << user_seq1.orf1[i] << " ";
	}
	cout << endl;

	user_seq1.IndexFinder();//start codon과 stop codon에 해당하는 인덱스 찾기
	user_seq1.OrfFinder(user);//orf 찾기

	//완료한 orf 출력
	cout << "찾은 ORF: " << endl;
	if (user_seq1.complete_orf.empty())//벡터비어있으면 관련 안내 출력
	{
		cout << "ORF가 탐색되지 않았습니다." << endl;
	}
	else//벡터가 차있다면 ORF를 출력
	{//스트링 이중벡터에 존재하는 모든것을 출력.

		for (int k = 0; k < user_seq1.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq1.complete_orf[k].size(); i++)
			{

				cout << user_seq1.complete_orf[k][i] << " ";
			}cout << endl;
		}
	}cout << endl << endl;

	//인트론 제거
	cout << "인트론 가공 후: " << endl;
	user_seq1.IntronFinder();
	if (user_seq1.intron_removed.empty())//벡터비어있으면 관련 안내 출력
	{
		cout << "ORF가 존재하지 않습니다." << endl;
	}
	else//벡터가 차있다면 ORF를 출력
	//스트링 이중벡터에 존재하는 모든것을 출력.
	{

		//찾은 ORF 서열의 수만큼 반복
		for (int k = 0; k < user_seq1.intron_removed.size(); k++)
		{//각 ORF에서 string 벡터를 출력
			for (int i = 0; i < user_seq1.intron_removed[k].size(); i++)//
			{

				cout << user_seq1.intron_removed[k][i] << " ";
			}cout << endl;
		}
	}

	//codon 변환 및 확인
	user_seq1.CodonDecipher();//인트론 제거된 ORF를 단백질화 시키는 함수
	cout << "codon 해독 후: " << endl;
	if (user_seq1.protein.empty())//벡터비어있으면 관련 안내 출력
	{
		cout << "단백질 서열이 존재하지 않습니다.";
	}
	else//벡터가 차있다면 ORF를 출력
	//스트링 이중벡터에 존재하는 모든것을 출력.
	{

		//찾은 ORF 서열의 수만큼 반복
		for (int k = 0; k < user_seq1.protein.size(); k++)
		{//각 ORF에서 string 벡터를 출력
			for (int i = 0; i < user_seq1.protein[k].size(); i++)//
			{

				cout << user_seq1.protein[k][i] << " ";
			}cout << endl;
		}
	}
	cout << endl;
	user_seq1.KozakCalculator();
	cout << "Kozak score : ";
	for (int i = 0; i < user_seq1.complete_index.size(); i++)
	{
		if (user_seq1.complete_score[i] == 7777)
		{
			cout << "계산불가" << " ";
		}
		else
		{
			cout << user_seq1.complete_score[i] << " ";
		}
	}
	cout << endl << endl;


	Orf user_seq2(user.GetDna2());//Orf객체에 분석할 서열 전달
	user_seq2.TransferSeq();//서열의 스트링 벡터화
	cout << "벡터화dna2: " << endl;
	for (int i = 0; i < user_seq2.orf1.size(); i++)//스트링 벡터화된 서열을 출력.
	{
		cout << user_seq2.orf1[i] << " ";
	}
	cout << endl;

	user_seq2.IndexFinder();//start codon과 stop codon에 해당하는 인덱스 찾기
	user_seq2.OrfFinder(user);//orf 찾기
	cout << "찾은 ORF: " << endl;
	if (user_seq2.complete_orf.empty())
	{
		cout << "ORF가 탐색되지 않았습니다." << endl;
	}
	else//벡터가 차있다면 ORF를 출력
	{

		for (int k = 0; k < user_seq2.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq2.complete_orf[k].size(); i++)
			{

				cout << user_seq2.complete_orf[k][i] << " ";
			}cout << endl;
		}
	}cout << endl << endl;
	//인트론 제거
	cout << "인트론 가공 후: " << endl;
	user_seq2.IntronFinder();
	if (user_seq2.intron_removed.empty())//벡터비어있으면 관련 안내 출력
	{
		cout << "ORF가 존재하지 않습니다." << endl;
	}
	else//벡터가 차있다면 ORF를 출력
	{//스트링 이중벡터에 존재하는 모든것을 출력.

		//찾은 ORF 서열의 수만큼 반복
		for (int k = 0; k < user_seq2.intron_removed.size(); k++)
		{//각 ORF에서 string 벡터를 출력
			for (int i = 0; i < user_seq2.intron_removed[k].size(); i++)//
			{

				cout << user_seq2.intron_removed[k][i] << " ";
			}cout << endl;
		}
	}
	//codon 변환 및 확인
	user_seq2.CodonDecipher();//인트론 제거된 ORF를 단백질화 시키는 함수
	cout << "codon 해독 후: " << endl;
	if (user_seq2.protein.empty())//벡터비어있으면 관련 안내 출력
	{
		cout << "단백질 서열이 존재하지 않습니다.";
	}
	else//벡터가 차있다면 ORF를 출력
	//스트링 이중벡터에 존재하는 모든것을 출력.
	{

		//찾은 ORF 서열의 수만큼 반복
		for (int k = 0; k < user_seq2.protein.size(); k++)
		{//각 ORF에서 string 벡터를 출력
			for (int i = 0; i < user_seq2.protein[k].size(); i++)//
			{

				cout << user_seq2.protein[k][i] << " ";
			}cout << endl;
		}
	}
	cout << endl;
	user_seq2.KozakCalculator();
	cout << "Kozak score : ";
	for (int i = 0; i < user_seq2.complete_index.size(); i++)
	{
		if (user_seq2.complete_score[i] == 7777)
		{
			cout << "계산불가" << " ";
		}
		else
		{
			cout << user_seq2.complete_score[i] << " ";
		}
	}
	cout << endl << endl;



	Orf user_seq3(user.GetDna3());//Orf객체에 분석할 서열 전달
	user_seq3.TransferSeq();//서열의 스트링 벡터화
	cout << "벡터화dna3: " << endl;
	for (int i = 0; i < user_seq3.orf1.size(); i++)//스트링 벡터화된 서열을 출력.
	{
		cout << user_seq3.orf1[i] << " ";
	}
	cout << endl;

	user_seq3.IndexFinder();//start codon과 stop codon에 해당하는 인덱스 찾기
	user_seq3.OrfFinder(user);//orf 찾기
	cout << "찾은 ORF: " << endl;
	if (user_seq3.complete_orf.empty())
	{
		cout << "ORF가 탐색되지 않았습니다." << endl;
	}
	else//벡터가 차있다면 ORF를 출력
	{

		for (int k = 0; k < user_seq3.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq3.complete_orf[k].size(); i++)
			{

				cout << user_seq3.complete_orf[k][i] << " ";
			}cout << endl;
		}
	}cout << endl << endl;
	//인트론 제거
	cout << "인트론 가공 후: " << endl;
	user_seq3.IntronFinder();
	if (user_seq3.intron_removed.empty())//벡터비어있으면 관련 안내 출력
	{
		cout << "ORF가 존재하지 않습니다." << endl;
	}
	else//벡터가 차있다면 ORF를 출력
	{//스트링 이중벡터에 존재하는 모든것을 출력.

		//찾은 ORF 서열의 수만큼 반복
		for (int k = 0; k < user_seq3.intron_removed.size(); k++)
		{//각 ORF에서 string 벡터를 출력
			for (int i = 0; i < user_seq3.intron_removed[k].size(); i++)
			{

				cout << user_seq3.intron_removed[k][i] << " ";
			}cout << endl;
		}
	}
	//codon 변환 및 확인
	user_seq3.CodonDecipher();//인트론 제거된 ORF를 단백질화 시키는 함수
	cout << "codon 해독 후: " << endl;
	if (user_seq3.protein.empty())//벡터비어있으면 관련 안내 출력
	{
		cout << "단백질 서열이 존재하지 않습니다.";
	}
	else//벡터가 차있다면 ORF를 출력
	//스트링 이중벡터에 존재하는 모든것을 출력.
	{

		//찾은 ORF 서열의 수만큼 반복
		for (int k = 0; k < user_seq3.protein.size(); k++)
		{//각 ORF에서 string 벡터를 출력
			for (int i = 0; i < user_seq3.protein[k].size(); i++)//
			{

				cout << user_seq3.protein[k][i] << " ";
			}cout << endl;
		}
	}
	cout << endl;
	user_seq3.KozakCalculator();
	cout << "Kozak score : ";
	for (int i = 0; i < user_seq3.complete_index.size(); i++)
	{
		if (user_seq3.complete_score[i] == 7777)
		{
			cout << "계산불가" << " ";
		}
		else
		{
			cout << user_seq3.complete_score[i] << " ";
		}
	}
	cout << endl << endl;


	Orf user_seq4(user.GetRDna1());//Orf객체에 분석할 서열 전달
	user_seq4.TransferSeq();//서열의 스트링 벡터화
	cout << "벡터화dna4: " << endl;
	for (int i = 0; i < user_seq4.orf1.size(); i++)//스트링 벡터화된 서열을 출력.
	{
		cout << user_seq4.orf1[i] << " ";
	}
	cout << endl;
	user_seq4.IndexFinder();//start codon과 stop codon에 해당하는 인덱스 찾기
	user_seq4.OrfFinder(user);//orf 찾기
	cout << "찾은 ORF: " << endl;
	if (user_seq4.complete_orf.empty())
	{
		cout << "ORF가 탐색되지 않았습니다." << endl;
	}
	else//벡터가 차있다면 ORF를 출력
	{

		for (int k = 0; k < user_seq4.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq4.complete_orf[k].size(); i++)
			{

				cout << user_seq4.complete_orf[k][i] << " ";
			}cout << endl;
		}
	}cout << endl << endl;
	//인트론 제거
	cout << "인트론 가공 후: " << endl;
	user_seq4.IntronFinder();
	if (user_seq4.intron_removed.empty())//벡터비어있으면 관련 안내 출력
	{
		cout << "ORF가 존재하지 않습니다." << endl;
	}
	else//벡터가 차있다면 ORF를 출력
	{//스트링 이중벡터에 존재하는 모든것을 출력.

		//찾은 ORF 서열의 수만큼 반복
		for (int k = 0; k < user_seq4.intron_removed.size(); k++)
		{//각 ORF에서 string 벡터를 출력
			for (int i = 0; i < user_seq4.intron_removed[k].size(); i++)//
			{

				cout << user_seq4.intron_removed[k][i] << " ";
			}cout << endl;
		}
	}
	//codon 변환 및 확인
	user_seq4.CodonDecipher();//인트론 제거된 ORF를 단백질화 시키는 함수
	cout << "codon 해독 후: " << endl;
	if (user_seq4.protein.empty())//벡터비어있으면 관련 안내 출력
	{
		cout << "단백질 서열이 존재하지 않습니다.";
	}
	else//벡터가 차있다면 ORF를 출력
	//스트링 이중벡터에 존재하는 모든것을 출력.
	{

		//찾은 ORF 서열의 수만큼 반복
		for (int k = 0; k < user_seq4.protein.size(); k++)
		{//각 ORF에서 string 벡터를 출력
			for (int i = 0; i < user_seq4.protein[k].size(); i++)//
			{

				cout << user_seq4.protein[k][i] << " ";
			}cout << endl;
		}
	}
	cout << endl;
	user_seq4.KozakCalculator();
	cout << "Kozak score : ";
	for (int i = 0; i < user_seq4.complete_index.size(); i++)
	{
		if (user_seq4.complete_score[i] == 7777)
		{
			cout << "계산불가" << " ";
		}
		else
		{
			cout << user_seq4.complete_score[i] << " ";
		}
	}
	cout << endl << endl;



	Orf user_seq5(user.GetRDna2());//Orf객체에 분석할 서열 전달
	user_seq5.TransferSeq();//서열의 스트링 벡터화
	cout << "벡터화dna5: " << endl;
	for (int i = 0; i < user_seq5.orf1.size(); i++)//스트링 벡터화된 서열을 출력.
	{
		cout << user_seq5.orf1[i] << " ";
	}
	cout << endl;
	user_seq5.IndexFinder();//start codon과 stop codon에 해당하는 인덱스 찾기
	user_seq5.OrfFinder(user);//orf 찾기
	cout << "찾은 ORF: " << endl;
	if (user_seq5.complete_orf.empty())
	{
		cout << "ORF가 탐색되지 않았습니다." << endl;
	}
	else//벡터가 차있다면 ORF를 출력
	{

		for (int k = 0; k < user_seq5.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq5.complete_orf[k].size(); i++)
			{

				cout << user_seq5.complete_orf[k][i] << " ";
			}cout << endl;
		}
	}cout << endl << endl;
	//인트론 제거
	cout << "인트론 가공 후: " << endl;
	user_seq5.IntronFinder();
	if (user_seq5.intron_removed.empty())//벡터비어있으면 관련 안내 출력
	{
		cout << "ORF가 존재하지 않습니다." << endl;
	}
	else//벡터가 차있다면 ORF를 출력
	{//스트링 이중벡터에 존재하는 모든것을 출력.

		//찾은 ORF 서열의 수만큼 반복
		for (int k = 0; k < user_seq5.intron_removed.size(); k++)
		{//각 ORF에서 string 벡터를 출력
			for (int i = 0; i < user_seq5.intron_removed[k].size(); i++)//
			{

				cout << user_seq5.intron_removed[k][i] << " ";
			}cout << endl;
		}
	}
	//codon 변환 및 확인
	user_seq5.CodonDecipher();//인트론 제거된 ORF를 단백질화 시키는 함수
	cout << "codon 해독 후: " << endl;
	if (user_seq5.protein.empty())//벡터비어있으면 관련 안내 출력
	{
		cout << "단백질 서열이 존재하지 않습니다.";
	}
	else//벡터가 차있다면 ORF를 출력
	//스트링 이중벡터에 존재하는 모든것을 출력.
	{

		//찾은 ORF 서열의 수만큼 반복
		for (int k = 0; k < user_seq5.protein.size(); k++)
		{//각 ORF에서 string 벡터를 출력
			for (int i = 0; i < user_seq5.protein[k].size(); i++)//
			{

				cout << user_seq5.protein[k][i] << " ";
			}cout << endl;
		}
	}
	cout << endl;
	user_seq5.KozakCalculator();
	cout << "Kozak score : ";
	for (int i = 0; i < user_seq5.complete_index.size(); i++)
	{
		if (user_seq5.complete_score[i] == 7777)
		{
			cout << "계산불가" << " ";
		}
		else
		{
			cout << user_seq5.complete_score[i] << " ";
		}
	}
	cout << endl << endl;




	Orf user_seq6(user.GetRDna3());//Orf객체에 분석할 서열 전달
	user_seq6.TransferSeq();//서열의 스트링 벡터화
	cout << "벡터화dna6: " << endl;
	for (int i = 0; i < user_seq6.orf1.size(); i++)//스트링 벡터화된 서열을 출력.
	{
		cout << user_seq6.orf1[i] << " ";
	}
	cout << endl;
	user_seq6.IndexFinder();//start codon과 stop codon에 해당하는 인덱스 찾기
	user_seq6.OrfFinder(user);//orf 찾기
	cout << "찾은 ORF: " << endl;
	if (user_seq6.complete_orf.empty())
	{
		cout << "ORF가 탐색되지 않았습니다." << endl;
	}
	else//벡터가 차있다면 ORF를 출력
	{

		for (int k = 0; k < user_seq6.complete_orf.size(); k++)
		{
			for (int i = 0; i < user_seq6.complete_orf[k].size(); i++)
			{

				cout << user_seq6.complete_orf[k][i] << " ";
			}cout << endl;
		}
	}cout << endl << endl;

	//인트론 제거
	cout << "인트론 가공 후: " << endl;
	user_seq6.IntronFinder();
	if (user_seq6.intron_removed.empty())//벡터비어있으면 관련 안내 출력
	{
		cout << "ORF가 존재하지 않습니다." << endl;
	}
	else//벡터가 차있다면 ORF를 출력
	{//스트링 이중벡터에 존재하는 모든것을 출력.

		//찾은 ORF 서열의 수만큼 반복
		for (int k = 0; k < user_seq6.intron_removed.size(); k++)
		{//각 ORF에서 string 벡터를 출력
			for (int i = 0; i < user_seq6.intron_removed[k].size(); i++)//
			{

				cout << user_seq6.intron_removed[k][i] << " ";
			}cout << endl;
		}
	}
	//codon 변환 및 확인
	user_seq6.CodonDecipher();//인트론 제거된 ORF를 단백질화 시키는 함수
	cout << "codon 해독 후: " << endl;
	if (user_seq6.protein.empty())//벡터비어있으면 관련 안내 출력
	{
		cout << "단백질 서열이 존재하지 않습니다.";
	}
	else//벡터가 차있다면 ORF를 출력
	//스트링 이중벡터에 존재하는 모든것을 출력.
	{

		//찾은 ORF 서열의 수만큼 반복
		for (int k = 0; k < user_seq6.protein.size(); k++)
		{//각 ORF에서 string 벡터를 출력
			for (int i = 0; i < user_seq6.protein[k].size(); i++)//
			{

				cout << user_seq6.protein[k][i] << " ";
			}cout << endl;
		}
	}
	cout << endl;
	user_seq6.KozakCalculator();
	cout << "Kozak score : ";
	for (int i = 0; i < user_seq6.complete_index.size(); i++)
	{
		if (user_seq6.complete_score[i] == 7777)
		{
			cout << "계산불가" << " ";
		}
		else
		{
			cout << user_seq6.complete_score[i] << " ";
		}
	}
	cout << endl << endl;

	//오알에프 인덱스 테스트
	for (int i = 0; i < user_seq1.complete_index.size(); i++)
	{
		cout << user_seq1.complete_index[i].case_num << " ";
		cout << user_seq1.complete_index[i].start_index << " ";
		cout << user_seq1.complete_index[i].stop_index << " ";
	}



	//결과 출력
	string out_sequence = user.GetSequence();
	for (char ori : out_sequence)//유저에게 받은 서열 출력
	{
		result << ori;
	}
	result << endl;
	//찾은 cds 파일로 출력
	for (int a = 0; a < user_seq1.protein.size(); a++)
	{
		for (int i = 0; i < user_seq1.complete_index[a].start_index; i++)
				{
					result << "   ";
				}
		for (int b = 0; b < user_seq1.protein[a].size(); b++)
		{
			result << user_seq1.protein[a][b];
		}result << endl;
	}
	for (int a = 0; a < user_seq2.protein.size(); a++)
	{
		result << " ";
		for (int i = 0; i < user_seq2.complete_index[a].start_index; i++)
		{
			result << "   ";
		}
		for (int b = 0; b < user_seq2.protein[a].size(); b++)
		{
			result << user_seq2.protein[a][b];
		}result << endl;
	}
	for (int a = 0; a < user_seq3.protein.size(); a++)
	{
		result << "  ";
		for (int i = 0; i < user_seq3.complete_index[a].start_index; i++)
		{
			result << "   ";
		}
		for (int b = 0; b < user_seq3.protein[a].size(); b++)
		{
			result << user_seq3.protein[a][b];
		}result << endl;
	}
	for (int a = 0; a < user_seq4.protein.size(); a++)
	{
		for (int i = 0; i < out_sequence.length() - 3 * (user_seq4.complete_index[a].stop_index+1); i++)
		{
			result << " ";
		}
		for (int b = user_seq4.protein[a].size()-1; b>=0; b--)
		{
			result << user_seq4.protein[a][b];
		}result << endl;
	}
	for (int a = 0; a < user_seq5.protein.size(); a++)
	{
		for (int i = 0; i < out_sequence.length() - 3 * (user_seq5.complete_index[a].stop_index+1) - 1; i++)
		{
			result << " ";
		}
		for (int b = user_seq5.protein[a].size() - 1; b >= 0; b--)
		{
			result << user_seq5.protein[a][b];
		}result << endl;
	}
	for (int a = 0; a < user_seq6.protein.size(); a++)
	{
		for (int i = 0; i < out_sequence.length() - 3 * (user_seq6.complete_index[a].stop_index+1) - 2; i++)
		{
			result << " ";
		}
		for (int b = user_seq6.protein[a].size() - 1; b >= 0; b--)
		{
			result << user_seq6.protein[a][b];
		}result << endl;
	}












	return 0;
}