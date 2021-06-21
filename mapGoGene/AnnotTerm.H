#ifndef _TISSUE
#define _TISSUE

#include <fstream>
#include <string>
#include <map>

using namespace std;

class AnnotTerm
{
	public:
		AnnotTerm();
		~AnnotTerm();
		//Although we are saying tissue, really we mean tissue and stage development
		int setAnnotTermName(const char*);
		const char* getAnnotTermName();
		
		int addMemberGene(const char*);
		map<string,int>& getMemberGenes();
		bool isMember(const string&);
		
		
	private:
		string termName;
		map<string,int> memberGenes;
};
#endif
