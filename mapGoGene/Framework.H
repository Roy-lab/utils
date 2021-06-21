#ifndef _FRAMEWORK
#define _FRAMEWORK

class GoTermManager;
class AnnotMgr;

class Framework
{
	public:
		Framework();
		~Framework();
  int init(const char* tissueFile,const char* ontFile,const char* stagename);
		int assignGoToTissue();
		int generateGeneTissueMap(const char*);
		int getCurrentNodes_Manual(const char*);
	private:
		double findOverlap(map<string,GoTerm*>&);
		GoTermManager imMgr;
		AnnotMgr tMgr;
		map<string,GoTerm*> termGoSet;
		map<string,GoTerm*> currNodes;
		map<string,double> wordDist; 
};
#endif
