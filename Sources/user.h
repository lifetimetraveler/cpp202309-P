#include<iostream>
#include<string>
using namespace std;

class User {
	string user_name;
	string user_dna_sequence;
public:
	void SetName(string name);
	void SetSequence(string sequence);
	string GetName();
	string GetSequence();
};