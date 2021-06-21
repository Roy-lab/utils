#ifndef _IMAGO_TERM
#define _IMAGO_TERM

#include <map>
#include <string>

using namespace std;

class AnnotTerm;

class GoTerm
{
	public:
		GoTerm();
		~GoTerm();
		int setID(int);
		int setAnnotID(const char*);
		int setName(const char*);
		int setChild(GoTerm*,bool partof=false);
		int setParent(GoTerm*);
		
		int getID();
		string& getAnnotID();
		string& getName();
		int getParentCnt();
		int getChildCnt();
		
		bool isChild(const char*);
		bool isDescendent(const char*);
		map<string,GoTerm*>& getParents();
		map<string,GoTerm*>& getChildren();

		map<string,int>& getWordSet();

		int setAnnotTerm(AnnotTerm*);
		AnnotTerm* getAnnotTerm();
		
		int createGeneList();
		map<string,int>& getGeneList();
		double getGeneOverlap(map<string,int>&);
		int getGeneCnt();

		int showAncestory(int level);
		int showChildTree(int currLevel, int finalLevel);
		int getChildren(map<string,GoTerm*>&,int,int);
		int getAncestors(map<string,GoTerm*>&);
		int assignLevelFromRoot(int);
		int getLevelFromRoot();
	private:
		int id;
		string annotId;
		//value of the name field
		string goName;
		map<string,int> relatedSynonyms;
		map<string,GoTerm*> children;
		map<string,GoTerm*> parents;
		map<string,int> partOfChild;

		map<string,int> wordSet;

		AnnotTerm* annotterm;
		map<string,int> geneList;
		int levelFromRoot;
};

#endif
