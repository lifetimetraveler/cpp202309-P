#include "user.h"
using namespace std;

void User::SetName(string name)  // ����� �̸��ʱ�ȭ
{
  user_name = name;
}
void User::SetSequence(string sequence)  // ����� DNA�����ʱ�ȭ
{
  user_dna_sequence = sequence;
}
void User::FrameSetting()  // reading frame�� �ٲ㼭 ������ �����ϴ� �Լ�
{
  string temp_dna =
      user_dna_sequence;  // DNA������ frame �� ������ ���� �ӽú���;

  dna1 = user_dna_sequence;  // ������ ������, reading frame DNA ����

  // �� �� �ϳ��� ������ frameshift
  dna2 = temp_dna.erase(0, 1);  // ������ ������, +1 reding frame DNA ����

  // �� �� ���� ������ frameshift
  dna3 = temp_dna.erase(0, 1);  // ������ ������, +2 reding frame DNA ����

  for (int i = (int)user_dna_sequence.length() - 1; i >= 0; i--) {
    reverse_dna1 +=
        user_dna_sequence[i];  // ������ ������, reading frame DNA ����
  }

  temp_dna = reverse_dna1;  // �������� ���� ����

  // �� �� �ϳ��� ������ frameshift
  reverse_dna2 =
      temp_dna.erase(0, 1);  // ������ ������, +1 reding frame DNA ����

  // �� �� ���� ������ frameshift
  reverse_dna3 =
      temp_dna.erase(0, 1);  // ������ ������, +2 reding frame DNA ����
}
string User::GetName() {  // ����� �̸���ȯ
  return user_name;
}
string User::GetSequence() {  // ����� DNA������ȯ
  return user_dna_sequence;
}
// ����� reading frame ��ȯ�� ������ ��ȯ�ϴ� �Լ���
string User::GetDna1() { return dna1; }
string User::GetDna2() { return dna2; }
string User::GetDna3() { return dna3; }
string User::GetRDna1() { return reverse_dna1; }
string User::GetRDna2() { return reverse_dna2; }
string User::GetRDna3() { return reverse_dna3; }

string Orf::GetOriSeq() { return original_seq; }

// ��Ʈ�������� ������ ������ ������ ��Ʈ�� ���Ϳ� �����Ͽ�
// ���Ŀ� ������ �����Ѵ�.
void Orf::TransferSeq()  // ��, ��Ʈ�� ����ȭ�ϴ� �Լ�
{
  for (int i = 0; i < original_seq.length(); i += 3) {
    string triplet = original_seq.substr(i, 3);  // �Ʒ��� �ڵ带 ����ص� �ȴ�.
    // to_string(original_seq[i])+to_string(original_seq[i+1])+to_string(original_seq[i+2])
    // ������ �̷��� ����Ϸ��� string�� ���̺��� ū �ε��� ���� ����ϰ� �ȴ�.
    // �� ��쿡�� 3�� ����� �ǵ��� ������ �ٵ��� �ؼ� ���ǻ� string��
    // �����ִ� ��Ʈ�� Ŭ���� ����Լ��� ����Ͽ���. substr�� ������ �����
    // ���̻� ��ȯ���� �ʴ´�.
    orf1.push_back(triplet);
  }
}
// ������ ������ start codon�� stop codon�� �ε����� �����ϴ� �Լ�
void Orf::IndexFinder() {  // ATG(start codon)�� ��Ʈ���� ã�Ƽ� �׶��� �ε���
                           // ���� ���Ϳ� ����
  for (int i = 0; i < orf1.size(); i++) {
    if (orf1[i] == "ATG") atg_index.push_back(i);
  }
  // TGA(stop codon)�� ��Ʈ���� ã�Ƽ� �׶��� �ε��� ���� ���Ϳ� ����
  for (int i = 0; i < orf1.size(); i++) {
    if (orf1[i] == "TGA") tga_index.push_back(i);
  }
  // TAA(stop codon)�� ��Ʈ���� ã�Ƽ� �׶��� �ε��� ���� ���Ϳ� ����
  for (int i = 0; i < orf1.size(); i++) {
    if (orf1[i] == "TAA") taa_index.push_back(i);
  }
  // TAG(stop codon)�� ��Ʈ���� ã�Ƽ� �׶��� �ε��� ���� ���Ϳ� ����
  for (int i = 0; i < orf1.size(); i++) {
    if (orf1[i] == "TAG") tag_index.push_back(i);
  }
}

// ����� �ε����� ���� ���� ����(ATG)���� ��(TAA,TAG,TGA)������ ������� ������
// ���� �� ������ ��Ʈ�� ���ߺ��Ϳ� �����Ѵ�.
void Orf::OrfFinder(User user) {
  int range1 = user.range1;  // ����ڰ� �Է��� ������ �����´�.
  int range2 = user.range2;
  for (int i = 0; i < atg_index.size();
       i++)  // ����� �����ڵ��� �ִ� �ε����� ����ϱ�����
  {  // start codon�� �����ϱ� ������ ������ stop codon�� ����. ������ stop
     // codon�� ������ for������ Ž��
    for (int j = 0; j < tga_index.size();
         j++)  // stop codon�� ����� �ε����� ����ϱ�����
    {          // stop codon�� TGA�� ���
      if (atg_index[i] <
          tga_index[j])  // ���� �ڵ��� stop codon���� �տ� ���� ��
      {  // ���ۺ��� �������� ������ ũ�Ⱑ ����ڰ� �Է��� ������ ��ġ�ϴ� ���
        if (tga_index[j] - atg_index[i] + 1 >= range1 &&
            tga_index[j] - atg_index[i] + 1 <= range2
           ) {  // �� ���ۺ��� ������(ORF)�� ��Ʈ�� ���� ���Ϳ�
                            // ����(push_back())�Ѵ�.
          vector<string> temp_com_orf(orf1.begin() + atg_index[i],
                                      orf1.begin() + tga_index[j] + 1);
          complete_orf.push_back(temp_com_orf);

          // ORF�� ���� �� �������� �ε��� ���� ����
          SavedIndex temp_saved(0, atg_index[i], tga_index[j]);
          complete_index.push_back(temp_saved);
        }
      }
    }
    for (int j = 0; j < taa_index.size();
         j++)  // stop codon�� �ε����� ����ϱ� ����
    {  // ���� for���� ������ ����, ������ stop codon�� TAA�� ���
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
         j++)  // stop codon�� �ε����� ����ϱ� ����
    {  // ���� for���� ������ ����, ������ stop codon�� TAG�� ���
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
// ��Ʈ���� �����ϴ� �Լ�
void Orf::IntronFinder() {
  intron_removed =
      complete_orf;  // ã�� orf�� ��Ʈ���� �߸� ������ ��� ���� ����

  for (int a = 0; a < intron_removed.size(); a++)  // ã�� ORF�� ����ŭ �ݺ�
  {
    for (int b = 0; b < intron_removed[a].size() - 5;
         b++)  // ������ ������ �Ͽ� ���� �ε����� ���� �ʵ���
    {          // ��ü ������ Ž��
      // Ž���ϴٰ� ��Ʈ���� ���� �κ��� ������
      if (intron_removed[a][b] == "GTA" && intron_removed[a][b + 1] == "AGT") {
        for (int c = b + 1; c < intron_removed[a].size();
             c++)  // �� ������ ������ Ž��
        {          // ��Ʈ���� ���κп� �ش��ϴ� �κ��� ã����
          if (intron_removed[a][c - 3] == "TTT" &&
              intron_removed[a][c - 2] == "TTT") {
            if (intron_removed[a][c - 1][0] == 'T' &&
                intron_removed[a][c - 1][1] == 'T') {
              if (intron_removed[a][c] == "CAG") {
                for (int i = b; i <= c; i++)  // ã�Ƴ� ��Ʈ���� ó���� ����,
                {
                    if (intron_removed[a][i] == "TGA" ||
                        intron_removed[a][i] == "TAA" ||
                        intron_removed[a][i] == "TAG")
                        {}//stop codon�� ��� �״�� ���ܳ���.
                    else
                        {
                        intron_removed[a][i] =
                      "itr";  // �ش��ϴ� ������ ��� ��Ʈ������ �ٲ�
                    }
                
                }
                break;  // �Ϸ��� �� for���� Ż�� �Ͽ� ���� ó������ ����
                        // ��Ʈ���� ���κб�����
                // ��Ʈ������ �ٲ�.
              }
            }
          }
        }
      }
    }
  }
}

// ��Ʈ�� ���ŵ� ORF�� �ܹ���ȭ ��Ű�� �Լ�
void Orf::CodonDecipher() {
  protein = intron_removed;
  for (int a = 0; a < protein.size(); a++)  // ã�� ORF�� ����
  {
    for (int b = 0; b < protein[a].size();
         b++)  // ���� ����� orf ���� string���� Ž��
    {  // �׿� �´� �ڵ��� ������ protein ������ �ܹ����� �ٲ㼭 ����.
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
                  .start_index;  // start�ε����� ����� ��Ʈ�������� ��ȣ�ϱ�
    int sequence_num = num * 3 + 1;  // ���� ���������� ������ ����.
    int score = 0;                   // ������ �����ϱ� ���� ����
    // ���� �����ڵ� �ֺ��� �����Ϸ��µ� �װ��� ������ �������� ���� ���.
    // ������Ȳ.
    if (sequence_num - 3 < 0 || sequence_num + 3 > original_seq.length()) {
      score = 7777;
    }
    // ������ �ƴѰ��
    else {
      // ���⼭���� +4��ġ�� C�Ǵ� G�� ������ 1�� ����.
      if (original_seq[sequence_num + 3] == 'C' ||
          original_seq[sequence_num + 3] == 'G')
        score += 1;
      // ���⼭���� -2��ġ�� A �Ǵ� G�� ������ 1�� ����.
      if (original_seq[sequence_num - 2] == 'A' ||
          original_seq[sequence_num - 2] == 'G')
        score += 1;
      // ���⼭���� -3��ġ�� A �Ǵ� G�� ������ 1�� ����.
      if (original_seq[sequence_num - 3] == 'A' ||
          original_seq[sequence_num - 3] == 'G')
        score += 3;
    }
    complete_score.push_back(score);  // �ջ��� ���� ���ھ ����.
  }
}
