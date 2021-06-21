#ifndef _IMAGO_TERMMANAGER
#define _IMAGO_TERMMANAGER

#include <map>
#include <string>

//#define OBO_FILE "/Users/sroy/mlresearch/sgddata/go/gene_ontology_edit_102607.obo"
//#define OBO_FILE "/Users/sroy/mlresearch/sgddata/go/gene_ontology.1_2.obo"
//#define OBO_FILE "inputdata/gene_ontology.1_2.obo_Jul2011"
using namespace std;

class GoTermManager
{
	public:
		GoTermManager();
		~GoTermManager();
		int setTermType(const char*);
		int readGoTerms(const char*);
		int readGoTermMap(const char*);
		int showHierarchy(const char*, int level);
        int saveNameIDMap();
		
		GoTerm* getTermForName(const char*);
		GoTerm* getTermForID(const char*);
		map<string,GoTerm*>& getAllGoTerms();
		int assignLevelFromRoot();
	private:
		map<string,GoTerm*> idToTermMap;
		map<string,GoTerm*> nameToTermMap;
		map<int,string> relationMap;
		string reqTermType;
};

#endif
