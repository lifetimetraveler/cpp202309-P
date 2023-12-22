#include <fstream>

#include "user.h"
using namespace std;

int main() {
  string name;                             // ������ �̸�
  string sequence;                         // ������ ����
  ifstream sequence_file("sequence.txt");  // ������ ���� ����
  User user;  // User�� �̸�, ����, ORF����, �ٸ� reading frame������
  // �����ϴ� Ŭ����

  // �̸��Է�, �׽�Ʈ�� ���� ������ �̸��� ����, ������ ����߽��ϴ�.
  cout << "�̸�: " << endl;
  //cin >> name;
  name = "��ö��";//�׽�Ʈ��

  // ���� �Է�
  cout << "DNA����:" << endl;
  try {
    if (!sequence_file)  // ���� �߻� �� �޼���//�̻��� ������ ���� �Է�
    {
      throw invalid_argument("������ �ҷ����� ���߽��ϴ�.");
    } else {
      sequence_file >> sequence;
    }
    sequence_file.close();
  } catch (invalid_argument& e) {
    cout << "����: " << e.what() << endl;
  }

  // ORF ���� �Է�
  cout << "ã�����ϴ� ORF�� ũ�� : " << endl;
  //cin >> user.range1 >> user.range2;  // ORF�� ������ �Է�
   user.range1 = 1;
   user.range2 = 100;
  //  range1�� �� ū ���� �Է��� ���, ���� �ٲ��ش�.
  if (user.range2 < user.range1) {
    int temp_length = user.range1;
    user.range1 = user.range2;
    user.range2 = temp_length;
  }
  user.SetName(name);          // ��ü���� �̸��ʱ�ȭ
  user.SetSequence(sequence);  // ��ü���� ���� �ʱ�ȭ
  // ������ ������ ���
  cout << "�̸� :" << user.GetName() << endl;
  cout << "���� : " << user.GetSequence() << endl;
  cout << "���� : " << user.range1 << "," << user.range2 << endl;

  // 6������ ������ readingframe���� ������ �����ϴ� �Լ�
  user.FrameSetting();
  // ������ �Ǿ����� Ȯ���ϰ� ���� ����ϴ� �ڵ�
  cout << "dna1: " << user.GetDna1() << endl;
  cout << "dna2: " << user.GetDna2() << endl;
  cout << "dna3: " << user.GetDna3() << endl;
  cout << "dna4: " << user.GetRDna1() << endl;
  cout << "dna5: " << user.GetRDna2() << endl;
  cout << "dna6: " << user.GetRDna3() << endl << endl;
 
  //������� ���� ��� �������� ������ ������ 6�� �ݺ��ϵ� reading frame�� �ٸ��� ����� ��ü�� ����մϴ�.
  // �� 6���� ������ ���� ���� ������� �м�. 6�� �ݺ�
  Orf user_seq1(user.GetDna1());  // Orf��ü�� �м��� ���� ���޹� Orf ��ü ����
  user_seq1.TransferSeq();  // ������ ��Ʈ�� ����ȭ
  cout << "����ȭdna1: " << endl;
  for (int i = 0; i < user_seq1.orf1.size();
       i++)  // ��Ʈ�� ����ȭ�� ������ ���.
  {
    cout << user_seq1.orf1[i] << " ";
  }
  cout << endl;

  user_seq1.IndexFinder();  // start codon�� stop codon�� �ش��ϴ� �ε��� ã��
  user_seq1.OrfFinder(user);  // orf ã��

  // �Ϸ��� orf ���
  cout << "ã�� ORF: " << endl;
  if (user_seq1.complete_orf.empty())  // ���ͺ�������� ���� �ȳ� ���
  {
    cout << "ORF�� Ž������ �ʾҽ��ϴ�." << endl;
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  {       // ��Ʈ�� ���ߺ��Ϳ� �����ϴ� ������ ���.
      //ã�� ORF�� ��� ���
    for (int k = 0; k < user_seq1.complete_orf.size(); k++) {
      for (int i = 0; i < user_seq1.complete_orf[k].size(); i++) {
        cout << user_seq1.complete_orf[k][i] << " ";
      }
      cout << endl;
    }
  }
  cout << endl << endl;

  // ��Ʈ�� ����
  cout << "��Ʈ�� ���� ��: " << endl;
  user_seq1.IntronFinder();
  if (user_seq1.intron_removed.empty())  // ���ͺ�������� ���� �ȳ� ���
  {
    cout << "ORF�� �������� �ʽ��ϴ�." << endl;
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  // ��Ʈ�� ���ߺ��Ϳ� �����ϴ� ������ ���.
  {
    // ã�� ORF ������ ����ŭ �ݺ�
    for (int k = 0; k < user_seq1.intron_removed.size();
         k++) {  // �� ORF���� string ���͸� ���
      for (int i = 0; i < user_seq1.intron_removed[k].size(); i++)  //
      {
        cout << user_seq1.intron_removed[k][i] << " ";
      }
      cout << endl;
    }
  }

  // codon ��ȯ �� Ȯ��
  user_seq1.CodonDecipher();  // ��Ʈ�� ���ŵ� ORF�� �ܹ���ȭ ��Ű�� �Լ�
  cout << "codon �ص� ��: " << endl;
  if (user_seq1.protein.empty())  // ���ͺ�������� ���� �ȳ� ���
  {
    cout << "�ܹ��� ������ �������� �ʽ��ϴ�.";
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  // ��Ʈ�� ���ߺ��Ϳ� �����ϴ� ������ ���.
  {
    // ã�� ORF ������ ����ŭ �ݺ�
    for (int k = 0; k < user_seq1.protein.size();
         k++) {  // �� ORF���� string ���͸� ���
      for (int i = 0; i < user_seq1.protein[k].size(); i++)  //
      {
        cout << user_seq1.protein[k][i] << " ";
      }
      cout << endl;
    }
  }
  cout << endl;
  // kozak score ���
  user_seq1.KozakCalculator();  // �� ORF�� kozak score�� ����Ͽ� ���Ϳ� ����
  // ����� ���͸� ����Ͽ� Kozak score�� ���
  cout << "Kozak score : ";
  for (int i = 0; i < user_seq1.complete_index.size(); i++) {
    // ORF�� �ʹ� ª�ٴ� ���� ������ �ʹ� ��꿡 �ʿ��� �ε����� ���� �Ұ��� ��
    if (user_seq1.complete_score[i] == 7777) {//���� �޽��� ���
      cout << "���Ұ�"
           << " ";
    } else {// �װ��� �ƴϸ� ����(score)�� ���
      cout << user_seq1.complete_score[i] << " ";
    }
  }
  cout << endl << endl;
  //��������� �����ϰ� 5�� �ٸ� ��ü(�ٸ� reading frame)�� �ݺ�

  Orf user_seq2(user.GetDna2());  // Orf��ü�� �м��� ���� ����
  user_seq2.TransferSeq();        // ������ ��Ʈ�� ����ȭ
  cout << "����ȭdna2: " << endl;
  for (int i = 0; i < user_seq2.orf1.size();
       i++)  // ��Ʈ�� ����ȭ�� ������ ���.
  {
    cout << user_seq2.orf1[i] << " ";
  }
  cout << endl;

  user_seq2.IndexFinder();  // start codon�� stop codon�� �ش��ϴ� �ε��� ã��
  user_seq2.OrfFinder(user);  // orf ã��
  cout << "ã�� ORF: " << endl;
  if (user_seq2.complete_orf.empty()) {
    cout << "ORF�� Ž������ �ʾҽ��ϴ�." << endl;
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  {
    for (int k = 0; k < user_seq2.complete_orf.size(); k++) {
      for (int i = 0; i < user_seq2.complete_orf[k].size(); i++) {
        cout << user_seq2.complete_orf[k][i] << " ";
      }
      cout << endl;
    }
  }
  cout << endl << endl;
  // ��Ʈ�� ����
  cout << "��Ʈ�� ���� ��: " << endl;
  user_seq2.IntronFinder();
  if (user_seq2.intron_removed.empty())  // ���ͺ�������� ���� �ȳ� ���
  {
    cout << "ORF�� �������� �ʽ��ϴ�." << endl;
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  {       // ��Ʈ�� ���ߺ��Ϳ� �����ϴ� ������ ���.

    // ã�� ORF ������ ����ŭ �ݺ�
    for (int k = 0; k < user_seq2.intron_removed.size();
         k++) {  // �� ORF���� string ���͸� ���
      for (int i = 0; i < user_seq2.intron_removed[k].size(); i++)  //
      {
        cout << user_seq2.intron_removed[k][i] << " ";
      }
      cout << endl;
    }
  }
  // codon ��ȯ �� Ȯ��
  user_seq2.CodonDecipher();  // ��Ʈ�� ���ŵ� ORF�� �ܹ���ȭ ��Ű�� �Լ�
  cout << "codon �ص� ��: " << endl;
  if (user_seq2.protein.empty())  // ���ͺ�������� ���� �ȳ� ���
  {
    cout << "�ܹ��� ������ �������� �ʽ��ϴ�.";
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  // ��Ʈ�� ���ߺ��Ϳ� �����ϴ� ������ ���.
  {
    // ã�� ORF ������ ����ŭ �ݺ�
    for (int k = 0; k < user_seq2.protein.size();
         k++) {  // �� ORF���� string ���͸� ���
      for (int i = 0; i < user_seq2.protein[k].size(); i++)  //
      {
        cout << user_seq2.protein[k][i] << " ";
      }
      cout << endl;
    }
  }
  cout << endl;
  //Kozak score ���. ������ ���� user_seq1�� ����
  user_seq2.KozakCalculator();
  cout << "Kozak score : ";
  for (int i = 0; i < user_seq2.complete_index.size(); i++) {
    if (user_seq2.complete_score[i] == 7777) {
      cout << "���Ұ�"
           << " ";
    } else {
      cout << user_seq2.complete_score[i] << " ";
    }
  }
  cout << endl << endl;

  Orf user_seq3(user.GetDna3());  // Orf��ü�� �м��� ���� ����
  user_seq3.TransferSeq();        // ������ ��Ʈ�� ����ȭ
  cout << "����ȭdna3: " << endl;
  for (int i = 0; i < user_seq3.orf1.size();
       i++)  // ��Ʈ�� ����ȭ�� ������ ���.
  {
    cout << user_seq3.orf1[i] << " ";
  }
  cout << endl;

  user_seq3.IndexFinder();  // start codon�� stop codon�� �ش��ϴ� �ε��� ã��
  user_seq3.OrfFinder(user);  // orf ã��
  cout << "ã�� ORF: " << endl;
  if (user_seq3.complete_orf.empty()) {
    cout << "ORF�� Ž������ �ʾҽ��ϴ�." << endl;
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  {
    for (int k = 0; k < user_seq3.complete_orf.size(); k++) {
      for (int i = 0; i < user_seq3.complete_orf[k].size(); i++) {
        cout << user_seq3.complete_orf[k][i] << " ";
      }
      cout << endl;
    }
  }
  cout << endl << endl;
  // ��Ʈ�� ����
  cout << "��Ʈ�� ���� ��: " << endl;
  user_seq3.IntronFinder();
  if (user_seq3.intron_removed.empty())  // ���ͺ�������� ���� �ȳ� ���
  {
    cout << "ORF�� �������� �ʽ��ϴ�." << endl;
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  {       // ��Ʈ�� ���ߺ��Ϳ� �����ϴ� ������ ���.

    // ã�� ORF ������ ����ŭ �ݺ�
    for (int k = 0; k < user_seq3.intron_removed.size();
         k++) {  // �� ORF���� string ���͸� ���
      for (int i = 0; i < user_seq3.intron_removed[k].size(); i++) {
        cout << user_seq3.intron_removed[k][i] << " ";
      }
      cout << endl;
    }
  }
  // codon ��ȯ �� Ȯ��
  user_seq3.CodonDecipher();  // ��Ʈ�� ���ŵ� ORF�� �ܹ���ȭ ��Ű�� �Լ�
  cout << "codon �ص� ��: " << endl;
  if (user_seq3.protein.empty())  // ���ͺ�������� ���� �ȳ� ���
  {
    cout << "�ܹ��� ������ �������� �ʽ��ϴ�.";
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  // ��Ʈ�� ���ߺ��Ϳ� �����ϴ� ������ ���.
  {
    // ã�� ORF ������ ����ŭ �ݺ�
    for (int k = 0; k < user_seq3.protein.size();
         k++) {  // �� ORF���� string ���͸� ���
      for (int i = 0; i < user_seq3.protein[k].size(); i++)  //
      {
        cout << user_seq3.protein[k][i] << " ";
      }
      cout << endl;
    }
  }
  cout << endl;
  // Kozak score ���. ������ ���� user_seq1�� ����
  user_seq3.KozakCalculator();
  cout << "Kozak score : ";
  for (int i = 0; i < user_seq3.complete_index.size(); i++) {
    if (user_seq3.complete_score[i] == 7777) {
      cout << "���Ұ�"
           << " ";
    } else {
      cout << user_seq3.complete_score[i] << " ";
    }
  }
  cout << endl << endl;

  Orf user_seq4(user.GetRDna1());  // Orf��ü�� �м��� ���� ����
  user_seq4.TransferSeq();         // ������ ��Ʈ�� ����ȭ
  cout << "����ȭdna4: " << endl;
  for (int i = 0; i < user_seq4.orf1.size();
       i++)  // ��Ʈ�� ����ȭ�� ������ ���.
  {
    cout << user_seq4.orf1[i] << " ";
  }
  cout << endl;
  user_seq4.IndexFinder();  // start codon�� stop codon�� �ش��ϴ� �ε��� ã��
  user_seq4.OrfFinder(user);  // orf ã��
  cout << "ã�� ORF: " << endl;
  if (user_seq4.complete_orf.empty()) {
    cout << "ORF�� Ž������ �ʾҽ��ϴ�." << endl;
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  {
    for (int k = 0; k < user_seq4.complete_orf.size(); k++) {
      for (int i = 0; i < user_seq4.complete_orf[k].size(); i++) {
        cout << user_seq4.complete_orf[k][i] << " ";
      }
      cout << endl;
    }
  }
  cout << endl << endl;
  // ��Ʈ�� ����
  cout << "��Ʈ�� ���� ��: " << endl;
  user_seq4.IntronFinder();
  if (user_seq4.intron_removed.empty())  // ���ͺ�������� ���� �ȳ� ���
  {
    cout << "ORF�� �������� �ʽ��ϴ�." << endl;
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  {       // ��Ʈ�� ���ߺ��Ϳ� �����ϴ� ������ ���.

    // ã�� ORF ������ ����ŭ �ݺ�
    for (int k = 0; k < user_seq4.intron_removed.size();
         k++) {  // �� ORF���� string ���͸� ���
      for (int i = 0; i < user_seq4.intron_removed[k].size(); i++)  //
      {
        cout << user_seq4.intron_removed[k][i] << " ";
      }
      cout << endl;
    }
  }
  // codon ��ȯ �� Ȯ��
  user_seq4.CodonDecipher();  // ��Ʈ�� ���ŵ� ORF�� �ܹ���ȭ ��Ű�� �Լ�
  cout << "codon �ص� ��: " << endl;
  if (user_seq4.protein.empty())  // ���ͺ�������� ���� �ȳ� ���
  {
    cout << "�ܹ��� ������ �������� �ʽ��ϴ�.";
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  // ��Ʈ�� ���ߺ��Ϳ� �����ϴ� ������ ���.
  {
    // ã�� ORF ������ ����ŭ �ݺ�
    for (int k = 0; k < user_seq4.protein.size();
         k++) {  // �� ORF���� string ���͸� ���
      for (int i = 0; i < user_seq4.protein[k].size(); i++)  //
      {
        cout << user_seq4.protein[k][i] << " ";
      }
      cout << endl;
    }
  }
  cout << endl;
  // Kozak score ���. ������ ���� user_seq1�� ����
  user_seq4.KozakCalculator();
  cout << "Kozak score : ";
  for (int i = 0; i < user_seq4.complete_index.size(); i++) {
    if (user_seq4.complete_score[i] == 7777) {
      cout << "���Ұ�"
           << " ";
    } else {
      cout << user_seq4.complete_score[i] << " ";
    }
  }
  cout << endl << endl;

  Orf user_seq5(user.GetRDna2());  // Orf��ü�� �м��� ���� ����
  user_seq5.TransferSeq();         // ������ ��Ʈ�� ����ȭ
  cout << "����ȭdna5: " << endl;
  for (int i = 0; i < user_seq5.orf1.size();
       i++)  // ��Ʈ�� ����ȭ�� ������ ���.
  {
    cout << user_seq5.orf1[i] << " ";
  }
  cout << endl;
  user_seq5.IndexFinder();  // start codon�� stop codon�� �ش��ϴ� �ε��� ã��
  user_seq5.OrfFinder(user);  // orf ã��
  cout << "ã�� ORF: " << endl;
  if (user_seq5.complete_orf.empty()) {
    cout << "ORF�� Ž������ �ʾҽ��ϴ�." << endl;
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  {
    for (int k = 0; k < user_seq5.complete_orf.size(); k++) {
      for (int i = 0; i < user_seq5.complete_orf[k].size(); i++) {
        cout << user_seq5.complete_orf[k][i] << " ";
      }
      cout << endl;
    }
  }
  cout << endl << endl;
  // ��Ʈ�� ����
  cout << "��Ʈ�� ���� ��: " << endl;
  user_seq5.IntronFinder();
  if (user_seq5.intron_removed.empty())  // ���ͺ�������� ���� �ȳ� ���
  {
    cout << "ORF�� �������� �ʽ��ϴ�." << endl;
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  {       // ��Ʈ�� ���ߺ��Ϳ� �����ϴ� ������ ���.

    // ã�� ORF ������ ����ŭ �ݺ�
    for (int k = 0; k < user_seq5.intron_removed.size();
         k++) {  // �� ORF���� string ���͸� ���
      for (int i = 0; i < user_seq5.intron_removed[k].size(); i++)  //
      {
        cout << user_seq5.intron_removed[k][i] << " ";
      }
      cout << endl;
    }
  }
  // codon ��ȯ �� Ȯ��
  user_seq5.CodonDecipher();  // ��Ʈ�� ���ŵ� ORF�� �ܹ���ȭ ��Ű�� �Լ�
  cout << "codon �ص� ��: " << endl;
  if (user_seq5.protein.empty())  // ���ͺ�������� ���� �ȳ� ���
  {
    cout << "�ܹ��� ������ �������� �ʽ��ϴ�.";
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  // ��Ʈ�� ���ߺ��Ϳ� �����ϴ� ������ ���.
  {
    // ã�� ORF ������ ����ŭ �ݺ�
    for (int k = 0; k < user_seq5.protein.size();
         k++) {  // �� ORF���� string ���͸� ���
      for (int i = 0; i < user_seq5.protein[k].size(); i++)  //
      {
        cout << user_seq5.protein[k][i] << " ";
      }
      cout << endl;
    }
  }
  cout << endl;
  // Kozak score ���. ������ ���� user_seq1�� ����
  user_seq5.KozakCalculator();
  cout << "Kozak score : ";
  for (int i = 0; i < user_seq5.complete_index.size(); i++) {
    if (user_seq5.complete_score[i] == 7777) {
      cout << "���Ұ�"
           << " ";
    } else {
      cout << user_seq5.complete_score[i] << " ";
    }
  }
  cout << endl << endl;

  Orf user_seq6(user.GetRDna3());  // Orf��ü�� �м��� ���� ����
  user_seq6.TransferSeq();         // ������ ��Ʈ�� ����ȭ
  cout << "����ȭdna6: " << endl;
  for (int i = 0; i < user_seq6.orf1.size();
       i++)  // ��Ʈ�� ����ȭ�� ������ ���.
  {
    cout << user_seq6.orf1[i] << " ";
  }
  cout << endl;
  user_seq6.IndexFinder();  // start codon�� stop codon�� �ش��ϴ� �ε��� ã��
  user_seq6.OrfFinder(user);  // orf ã��
  cout << "ã�� ORF: " << endl;
  if (user_seq6.complete_orf.empty()) {
    cout << "ORF�� Ž������ �ʾҽ��ϴ�." << endl;
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  {
    for (int k = 0; k < user_seq6.complete_orf.size(); k++) {
      for (int i = 0; i < user_seq6.complete_orf[k].size(); i++) {
        cout << user_seq6.complete_orf[k][i] << " ";
      }
      cout << endl;
    }
  }
  cout << endl << endl;

  // ��Ʈ�� ����
  cout << "��Ʈ�� ���� ��: " << endl;
  user_seq6.IntronFinder();
  if (user_seq6.intron_removed.empty())  // ���ͺ�������� ���� �ȳ� ���
  {
    cout << "ORF�� �������� �ʽ��ϴ�." << endl;
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  {       // ��Ʈ�� ���ߺ��Ϳ� �����ϴ� ������ ���.

    // ã�� ORF ������ ����ŭ �ݺ�
    for (int k = 0; k < user_seq6.intron_removed.size();
         k++) {  // �� ORF���� string ���͸� ���
      for (int i = 0; i < user_seq6.intron_removed[k].size(); i++)  //
      {
        cout << user_seq6.intron_removed[k][i] << " ";
      }
      cout << endl;
    }
  }
  // codon ��ȯ �� Ȯ��
  user_seq6.CodonDecipher();  // ��Ʈ�� ���ŵ� ORF�� �ܹ���ȭ ��Ű�� �Լ�
  cout << "codon �ص� ��: " << endl;
  if (user_seq6.protein.empty())  // ���ͺ�������� ���� �ȳ� ���
  {
    cout << "�ܹ��� ������ �������� �ʽ��ϴ�.";
  } else  // ���Ͱ� ���ִٸ� ORF�� ���
  // ��Ʈ�� ���ߺ��Ϳ� �����ϴ� ������ ���.
  {
    // ã�� ORF ������ ����ŭ �ݺ�
    for (int k = 0; k < user_seq6.protein.size();
         k++) {  // �� ORF���� string ���͸� ���
      for (int i = 0; i < user_seq6.protein[k].size(); i++)  //
      {
        cout << user_seq6.protein[k][i] << " ";
      }
      cout << endl;
    }
  }
  cout << endl;
  // Kozak score ���. ������ ���� user_seq1�� ����
  user_seq6.KozakCalculator();
  cout << "Kozak score : ";
  for (int i = 0; i < user_seq6.complete_index.size(); i++) {
    if (user_seq6.complete_score[i] == 7777) {
      cout << "���Ұ�"
           << " ";
    } else {
      cout << user_seq6.complete_score[i] << " ";
    }
  }
  cout << endl << endl;

  // ���˿��� �ε��� �׽�Ʈ
 // for (int i = 0; i < user_seq1.complete_index.size(); i++) {
   // cout << user_seq1.complete_index[i].case_num << " ";
    //cout << user_seq1.complete_index[i].start_index << " ";
    //cout << user_seq1.complete_index[i].stop_index << " ";
  //}

  // ���Ϸ� ��� ���
  string file_name =
      user.GetName() +
      ".txt";  // user�� �̸����� �ؽ�Ʈ���� �����ϱ� ���� �̸� ���� ����
  ofstream result(file_name,
                  ios::out | ios::trunc);  // �м������ ����ϱ� ����.
  string out_sequence = user.GetSequence();
  for (char ori : out_sequence)  // �������� ���� ���� ���
  {
    result << ori;
  }
  result << endl;
  // ��ȯ�� �ܹ��� ������ ���Ϸ� ���, ����ڰ� ������ ������ ��ġ�� ��ġ�ϵ���
  // ������ �߰�.
  // orf�� �����ϴ� ���� ������������ ���踦 ��Ÿ���� ����, �ո�ŭ ������
  // �ֱ����� �������� ó��
  // user_seq1�� ����� ���Ͽ� ���
  for (int a = 0; a < user_seq1.protein.size();
       a++)  // ã�� ORF�� ����ŭ �ݺ�
  {
    for (int i = 0; i < user_seq1.complete_index[a].start_index; i++) {
      result << "   ";//�����߰�
    }
    for (int b = 0; b < user_seq1.protein[a].size(); b++) {
      result << user_seq1.protein[a][b];//�ܹ��� �������
    }
    result << endl;
  }
  // user_seq2�� ����� ���Ͽ� ���
  for (int a = 0; a < user_seq2.protein.size(); a++) {
    result << " ";//�����߰�
    for (int i = 0; i < user_seq2.complete_index[a].start_index; i++) {
      result << "   ";//�����߰�
    }
    for (int b = 0; b < user_seq2.protein[a].size(); b++) {
      result << user_seq2.protein[a][b];//�ܹ��� �������
    }
    result << endl;
  }
  // user_seq3�� ����� ���Ͽ� ���
  for (int a = 0; a < user_seq3.protein.size(); a++) {
    result << "  ";//�����߰�
    for (int i = 0; i < user_seq3.complete_index[a].start_index; i++) {
      result << "   ";//�����߰�
    }
    for (int b = 0; b < user_seq3.protein[a].size(); b++) {
      result << user_seq3.protein[a][b];//�ܹ��� �������
    }
    result << endl;
  }
  // user_seq4�� ����� ���Ͽ� ���
  for (int a = 0; a < user_seq4.protein.size(); a++) {
    for (int i = 0; i < out_sequence.length() -
                            3 * (user_seq4.complete_index[a].stop_index + 1);
         i++) {
      result << " ";//�����߰�
    }
    for (int b = user_seq4.protein[a].size() - 1; b >= 0; b--) {
      result << user_seq4.protein[a][b];//�ܹ��� �������
    }
    result << endl;
  }
  // user_seq5�� ����� ���Ͽ� ���
  for (int a = 0; a < user_seq5.protein.size(); a++) {
    for (int i = 0;
         i < out_sequence.length() -
                 3 * (user_seq5.complete_index[a].stop_index + 1) - 1;
         i++) {
      result << " ";//�����߰�
    }
    for (int b = user_seq5.protein[a].size() - 1; b >= 0; b--) {
      result << user_seq5.protein[a][b];//�ܹ��� �������
    }
    result << endl;
  }
  // user_seq6�� ����� ���Ͽ� ���
  for (int a = 0; a < user_seq6.protein.size(); a++) {
    for (int i = 0;
         i < out_sequence.length() -
                 3 * (user_seq6.complete_index[a].stop_index + 1) - 2;
         i++) {
      result << " ";//�����߰�
    }
    for (int b = user_seq6.protein[a].size() - 1; b >= 0; b--) {
      result << user_seq6.protein[a][b];//�ܹ��� �������
    }
    result << endl;
  }
  result.close();

  return 0;
}