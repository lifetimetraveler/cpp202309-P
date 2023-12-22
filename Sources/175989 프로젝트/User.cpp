#include "user.h"
using namespace std;

void User::SetName(string name)  // 사용자 이름초기화
{
  user_name = name;
}
void User::SetSequence(string sequence)  // 사용자 DNA서열초기화
{
  user_dna_sequence = sequence;
}
void User::FrameSetting()  // reading frame을 바꿔서 서열을 저장하는 함수
{
  string temp_dna =
      user_dna_sequence;  // DNA서열의 frame 별 저장을 위한 임시변수;

  dna1 = user_dna_sequence;  // 유저의 정방향, reading frame DNA 저장

  // 맨 앞 하나를 지워서 frameshift
  dna2 = temp_dna.erase(0, 1);  // 유저의 정방향, +1 reding frame DNA 저장

  // 맨 앞 둘을 지워서 frameshift
  dna3 = temp_dna.erase(0, 1);  // 유저의 정방향, +2 reding frame DNA 저장

  for (int i = (int)user_dna_sequence.length() - 1; i >= 0; i--) {
    reverse_dna1 +=
        user_dna_sequence[i];  // 유저의 역방향, reading frame DNA 저장
  }

  temp_dna = reverse_dna1;  // 역방향을 위한 저장

  // 맨 앞 하나를 지워서 frameshift
  reverse_dna2 =
      temp_dna.erase(0, 1);  // 유저의 역방향, +1 reding frame DNA 저장

  // 맨 앞 둘을 지워서 frameshift
  reverse_dna3 =
      temp_dna.erase(0, 1);  // 유저의 역방향, +2 reding frame DNA 저장
}
string User::GetName() {  // 사용자 이름반환
  return user_name;
}
string User::GetSequence() {  // 사용자 DNA서열반환
  return user_dna_sequence;
}
// 저장된 reading frame 변환된 서열을 반환하는 함수들
string User::GetDna1() { return dna1; }
string User::GetDna2() { return dna2; }
string User::GetDna3() { return dna3; }
string User::GetRDna1() { return reverse_dna1; }
string User::GetRDna2() { return reverse_dna2; }
string User::GetRDna3() { return reverse_dna3; }

string Orf::GetOriSeq() { return original_seq; }

// 스트링형태의 서열을 세개씩 나눠서 스트링 벡터에 저장하여
// 추후에 가공을 쉽게한다.
void Orf::TransferSeq()  // 즉, 스트링 벡터화하는 함수
{
  for (int i = 0; i < original_seq.length(); i += 3) {
    string triplet = original_seq.substr(i, 3);  // 아래의 코드를 사용해도 된다.
    // to_string(original_seq[i])+to_string(original_seq[i+1])+to_string(original_seq[i+2])
    // 하지만 이렇게 사용하려면 string의 길이보다 큰 인덱스 값을 사용하게 된다.
    // 그 경우에는 3의 배수가 되도록 서열을 다듬어야 해서 편의상 string에
    // 원래있는 스트링 클래스 멤버함수를 사용하였다. substr은 범위를 벗어나면
    // 더이상 반환하지 않는다.
    orf1.push_back(triplet);
  }
}
// 가공한 서열의 start codon과 stop codon의 인덱스를 저장하는 함수
void Orf::IndexFinder() {  // ATG(start codon)인 스트링을 찾아서 그때의 인덱스
                           // 값을 벡터에 저장
  for (int i = 0; i < orf1.size(); i++) {
    if (orf1[i] == "ATG") atg_index.push_back(i);
  }
  // TGA(stop codon)인 스트링을 찾아서 그때의 인덱스 값을 벡터에 저장
  for (int i = 0; i < orf1.size(); i++) {
    if (orf1[i] == "TGA") tga_index.push_back(i);
  }
  // TAA(stop codon)인 스트링을 찾아서 그때의 인덱스 값을 벡터에 저장
  for (int i = 0; i < orf1.size(); i++) {
    if (orf1[i] == "TAA") taa_index.push_back(i);
  }
  // TAG(stop codon)인 스트링을 찾아서 그때의 인덱스 값을 벡터에 저장
  for (int i = 0; i < orf1.size(); i++) {
    if (orf1[i] == "TAG") tag_index.push_back(i);
  }
}

// 저장된 인덱스의 값에 따라 시작(ATG)부터 끝(TAA,TAG,TGA)까지가 사용자의 범위에
// 들어가면 그 서열을 스트링 이중벡터에 저장한다.
void Orf::OrfFinder(User user) {
  int range1 = user.range1;  // 사용자가 입력한 범위를 가져온다.
  int range2 = user.range2;
  for (int i = 0; i < atg_index.size();
       i++)  // 저장된 시작코돈이 있는 인덱스를 사용하기위함
  {  // start codon은 동일하기 때문에 세개의 stop codon이 공유. 하위의 stop
     // codon을 세개의 for문으로 탐색
    for (int j = 0; j < tga_index.size();
         j++)  // stop codon이 저장된 인덱스를 사용하기위함
    {          // stop codon이 TGA인 경우
      if (atg_index[i] <
          tga_index[j])  // 시작 코돈이 stop codon보다 앞에 있을 때
      {  // 시작부터 끝까지의 서열의 크기가 사용자가 입력한 범위와 일치하는 경우
        if (tga_index[j] - atg_index[i] + 1 >= range1 &&
            tga_index[j] - atg_index[i] + 1 <= range2
           ) {  // 그 시작부터 끝까지(ORF)를 스트링 이중 벡터에
                            // 저장(push_back())한다.
          vector<string> temp_com_orf(orf1.begin() + atg_index[i],
                                      orf1.begin() + tga_index[j] + 1);
          complete_orf.push_back(temp_com_orf);

          // ORF가 추출 될 때마다의 인덱스 값을 저장
          SavedIndex temp_saved(0, atg_index[i], tga_index[j]);
          complete_index.push_back(temp_saved);
        }
      }
    }
    for (int j = 0; j < taa_index.size();
         j++)  // stop codon의 인덱스를 사용하기 위함
    {  // 위의 for문과 동일한 원리, 하지만 stop codon이 TAA인 경우
      if (atg_index[i] < taa_index[j]) {
        if (taa_index[j] - atg_index[i] + 1 >= range1 &&
            taa_index[j] - atg_index[i] + 1 <= range2) {
          vector<string> temp_com_orf(orf1.begin() + atg_index[i],
                                      orf1.begin() + taa_index[j] + 1);
          complete_orf.push_back(temp_com_orf);

          SavedIndex temp_saved(1, atg_index[i], taa_index[j]);
          complete_index.push_back(temp_saved);
        }
      }
    }
    for (int j = 0; j < tag_index.size();
         j++)  // stop codon의 인덱스를 사용하기 위함
    {  // 위의 for문과 동일한 원리, 하지만 stop codon이 TAG인 경우
      if (atg_index[i] < tag_index[j]) {
        if (tag_index[j] - atg_index[i] + 1 >= range1 &&
            tag_index[j] - atg_index[i] + 1 <= range2) {
          vector<string> temp_com_orf(orf1.begin() + atg_index[i],
                                      orf1.begin() + tag_index[j] + 1);
          complete_orf.push_back(temp_com_orf);

          SavedIndex temp_saved(2, atg_index[i], tag_index[j]);
          complete_index.push_back(temp_saved);
        }
      }
    }
  }
}
// 인트론을 제거하는 함수
void Orf::IntronFinder() {
  intron_removed =
      complete_orf;  // 찾은 orf를 인트론이 잘린 변수가 담길 곳에 복사

  for (int a = 0; a < intron_removed.size(); a++)  // 찾은 ORF의 수만큼 반복
  {
    for (int b = 0; b < intron_removed[a].size() - 5;
         b++)  // 끝까지 못가게 하여 없는 인덱스를 넣지 않도록
    {          // 전체 서열을 탐색
      // 탐색하다가 인트론의 시작 부분을 만나면
      if (intron_removed[a][b] == "GTA" && intron_removed[a][b + 1] == "AGT") {
        for (int c = b + 1; c < intron_removed[a].size();
             c++)  // 그 이후의 서열을 탐색
        {          // 인트론의 끝부분에 해당하는 부분을 찾으면
          if (intron_removed[a][c - 3] == "TTT" &&
              intron_removed[a][c - 2] == "TTT") {
            if (intron_removed[a][c - 1][0] == 'T' &&
                intron_removed[a][c - 1][1] == 'T') {
              if (intron_removed[a][c] == "CAG") {
                for (int i = b; i <= c; i++)  // 찾아낸 인트론의 처음과 끝을,
                {
                    if (intron_removed[a][i] == "TGA" ||
                        intron_removed[a][i] == "TAA" ||
                        intron_removed[a][i] == "TAG")
                        {}//stop codon일 경우 그대로 남겨놓음.
                    else
                        {
                        intron_removed[a][i] =
                      "itr";  // 해당하는 범위를 모두 인트론으로 바꿈
                    }
                
                }
                break;  // 완료한 후 for문을 탈출 하여 가장 처음으로 만난
                        // 인트론의 끝부분까지만
                // 인트론으로 바꿈.
              }
            }
          }
        }
      }
    }
  }
}

// 인트론 제거된 ORF를 단백질화 시키는 함수
void Orf::CodonDecipher() {
  protein = intron_removed;
  for (int a = 0; a < protein.size(); a++)  // 찾은 ORF에 접근
  {
    for (int b = 0; b < protein[a].size();
         b++)  // 각각 저장된 orf 안의 string들을 탐색
    {  // 그에 맞는 코돈을 만나면 protein 변수에 단백질로 바꿔서 저장.
      string pro_cpy = protein[a][b];
      if (pro_cpy == "TTT" || pro_cpy == "TTC")
        protein[a][b] = " F ";
      else if (pro_cpy == "TTA" || pro_cpy == "TTG" || pro_cpy == "CTT" ||
               pro_cpy == "CTC" || pro_cpy == "CTA" || pro_cpy == "CTG")
        protein[a][b] = " L ";
      else if (pro_cpy == "ATT" || pro_cpy == "ATC" || pro_cpy == "ATA")
        protein[a][b] = " L ";
      else if (pro_cpy == "ATG")
        protein[a][b] = " M ";
      else if (pro_cpy == "GTT" || pro_cpy == "GTC" || pro_cpy == "GTA" ||
               pro_cpy == "GTG")
        protein[a][b] = " V ";
      else if (pro_cpy == "TCT" || pro_cpy == "TCC" || pro_cpy == "TCA" ||
               pro_cpy == "TCG" || pro_cpy == "AGT" || pro_cpy == "AGC")
        protein[a][b] = " S ";
      else if (pro_cpy == "CCT" || pro_cpy == "CCC" || pro_cpy == "CCA" ||
               pro_cpy == "CCG")
        protein[a][b] = " P ";
      else if (pro_cpy == "ACT" || pro_cpy == "ACC" || pro_cpy == "ACA" ||
               pro_cpy == "ACG")
        protein[a][b] = " T ";
      else if (pro_cpy == "GCT" || pro_cpy == "GCC" || pro_cpy == "GCA" ||
               pro_cpy == "GCG")
        protein[a][b] = " A ";
      else if (pro_cpy == "TAT" || pro_cpy == "TAC")
        protein[a][b] = " Y ";
      else if (pro_cpy == "CAT" || pro_cpy == "CAC")
        protein[a][b] = " H ";
      else if (pro_cpy == "CAA" || pro_cpy == "CAG")
        protein[a][b] = " Q ";
      else if (pro_cpy == "AAT" || pro_cpy == "AAC")
        protein[a][b] = " N ";
      else if (pro_cpy == "AAA" || pro_cpy == "AAG")
        protein[a][b] = " K ";
      else if (pro_cpy == "GAT" || pro_cpy == "GAC")
        protein[a][b] = " D ";
      else if (pro_cpy == "GAA" || pro_cpy == "GAG")
        protein[a][b] = " E ";
      else if (pro_cpy == "TGT" || pro_cpy == "TGC")
        protein[a][b] = " C ";
      else if (pro_cpy == "TGG")
        protein[a][b] = " W ";
      else if (pro_cpy == "CGT" || pro_cpy == "CGC" || pro_cpy == "CGA" ||
               pro_cpy == "CGG" || pro_cpy == "AGA" || pro_cpy == "AGG")
        protein[a][b] = " R ";
      else if (pro_cpy == "GGT" || pro_cpy == "GGC" || pro_cpy == "GGA" ||
               pro_cpy == "GGG")
        protein[a][b] = " G ";
      else if (pro_cpy == "TAA" || pro_cpy == "TAG" || pro_cpy == "TGA")
        protein[a][b] = "stp";
    }
  }
}

void Orf::KozakCalculator() {
  for (int i = 0; i < complete_index.size(); i++) {
    int num = complete_index[i]
                  .start_index;  // start인덱스는 저장된 스트링에서의 번호니까
    int sequence_num = num * 3 + 1;  // 실제 서열에서의 순번을 저장.
    int score = 0;                   // 점수를 저장하기 위한 변수
    // 만약 시작코돈 주변을 접근하려는데 그곳에 서열이 존재하지 않을 경우.
    // 오류상황.
    if (sequence_num - 3 < 0 || sequence_num + 3 > original_seq.length()) {
      score = 7777;
    }
    // 오류가 아닌경우
    else {
      // 염기서열의 +4위치에 C또는 G가 있으면 1점 가산.
      if (original_seq[sequence_num + 3] == 'C' ||
          original_seq[sequence_num + 3] == 'G')
        score += 1;
      // 염기서열의 -2위치에 A 또는 G가 있으면 1점 가산.
      if (original_seq[sequence_num - 2] == 'A' ||
          original_seq[sequence_num - 2] == 'G')
        score += 1;
      // 염기서열의 -3위치에 A 또는 G가 있으면 1점 가산.
      if (original_seq[sequence_num - 3] == 'A' ||
          original_seq[sequence_num - 3] == 'G')
        score += 3;
    }
    complete_score.push_back(score);  // 합산한 최종 스코어를 저장.
  }
}
