#include"user.h"
using namespace std;



void User::SetName(string name) {
	user_name = name;
	}
void User::SetSequence(string sequence) {
	user_dna_sequence = sequence;
	}
string User::GetName() {
	return user_name;
}
string User::GetSequence() {
	return user_dna_sequence;
}

