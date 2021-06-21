#ifndef _TISSUE_MGR
#define _TISSUE_MGR

#include <map>
#include <string>

using namespace std;

class AnnotTerm;
class AnnotMgr
{
	public:
		AnnotMgr();
		~AnnotMgr();
		int setTermType(const char*);
		int readAnnotTerms(const char*);
		map<string,AnnotTerm*>& getAllAnnotations();
		AnnotTerm* getAnnotTerm(const char*);
		int dumpAnnotTermGeneMatrix(const char*);
	private:
  int populateAnnotTermData(char*,int);
		map<string,string> annotNametoID;
		map<string,string> annotIDtoName;
		map<string,AnnotTerm*> annotSet;

		map<string,int> geneSet;
		int annotTermCnt;
		string termtype;
};
#endif
