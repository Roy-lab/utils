#include <iostream>
using namespace std;

#include "GoTerm.H"
#include "GoTermManager.H"
#include "AnnotTerm.H"
#include "AnnotMgr.H"
#include "Framework.H"



Framework::Framework()
{
}

Framework::~Framework()
{
}

int 
Framework::init(const char* genetermfile,const char* geneontfile,const char* gotype)
{
	imMgr.setTermType(gotype);
	imMgr.readGoTerms(geneontfile);
	imMgr.readGoTermMap(geneontfile);
	imMgr.assignLevelFromRoot();
	tMgr.setTermType(gotype);
	tMgr.readAnnotTerms(genetermfile);
	return 0;
}

int 
Framework::assignGoToTissue()
{
	map<string,AnnotTerm*>& allTerms=tMgr.getAllAnnotations();
	int miss=0;
	int hit=0;
	for(map<string,AnnotTerm*>::iterator tIter=allTerms.begin();tIter!=allTerms.end();tIter++)
	{
		AnnotTerm* aterm=tIter->second;
		GoTerm* imTerm=imMgr.getTermForID(aterm->getAnnotTermName());
		if(imTerm==NULL)
		{
			cout << "No term id " << aterm->getAnnotTermName()<< " in GO OBO file " << endl;
			miss++;
		}
		else
		{
			hit++;
		}
		if(imTerm!=NULL)
		{
			imTerm->setAnnotTerm(aterm);
			termGoSet[imTerm->getName()]=imTerm;
		}
	}
	cout <<"Found " << hit<< " with go info " << endl;
	cout <<"Missed " << miss<< " with no go info  " << endl;
	return 0;
}

//So here, we are going to get all nodes 
int
Framework::getCurrentNodes_Manual(const char* suff)
{
	map<string,GoTerm*>& imTerms=imMgr.getAllGoTerms();
	for(map<string,GoTerm*>::iterator aIter=imTerms.begin();aIter!=imTerms.end();aIter++)
	{
		GoTerm* img=aIter->second;
		img->createGeneList();
	}
	int otherRootNodeCnt=0;
	for(map<string,GoTerm*>::iterator aIter=imTerms.begin();aIter!=imTerms.end();aIter++)
	{
		GoTerm* img=aIter->second;
		currNodes[img->getName()]=img;
	}
	char aFName[256];
	sprintf(aFName,"%sgenecnt.txt",suff);
	ofstream oFile(aFName);
	for(map<string,GoTerm*>::iterator aIter=currNodes.begin();aIter!=currNodes.end();aIter++)
	{
		oFile << aIter->second->getName() << "\t" << aIter->second->getGeneCnt();
		if(aIter->second->getParentCnt()==0)
		{
			oFile <<"\tRoot\t"<<aIter->second->getLevelFromRoot() << endl;
		}
		else
		{
			oFile <<"\tChild\t"<< aIter->second->getLevelFromRoot() << endl;
		}
	}
	oFile.close();
	char gtMapFName[256];
	sprintf(gtMapFName,"%sgotermap.txt",suff);
	generateGeneTissueMap(gtMapFName);
	
	return 0;
}

int
Framework::generateGeneTissueMap(const char* aFName)
{
	map<string,map<string,int>*> geneTissueMap;
	//Make the gene set, where each gene is associated with a set of motifs
	for(map<string,GoTerm*>::iterator aIter=currNodes.begin();aIter!=currNodes.end();aIter++)
	{
		map<string,int>& genelist=aIter->second->getGeneList();
		for(map<string,int>::iterator gIter=genelist.begin();gIter!=genelist.end();gIter++)
		{
			map<string,int>* geneTissues=NULL;
			if(geneTissueMap.find(gIter->first)==geneTissueMap.end())
			{
				geneTissues=new map<string,int>;
				geneTissueMap[gIter->first]=geneTissues;
			}
			else
			{
				geneTissues=geneTissueMap[gIter->first];
			}
			(*geneTissues)[aIter->second->getName()]=gIter->second;
		}
	}

	cout<<"Total number of genes with go information "<< geneTissueMap.size() << endl;
	
	ofstream oFile(aFName);
	oFile << "GeneName\tGOTerm\tTermLevel" <<endl;
	for(map<string,map<string,int>*>::iterator gIter=geneTissueMap.begin();gIter!=geneTissueMap.end();gIter++)
	{
		map<string,int>* geneTissues=gIter->second;
		for(map<string,int>::iterator tIter=geneTissues->begin();tIter!=geneTissues->end();tIter++)
		{
			oFile << gIter->first.c_str()<<"\t"<<tIter->first.c_str() << "\t" << tIter->second << endl;
		}
	}

	oFile.close();
	return 0;
}


int
main(int argc, const char** argv)
{
	if(argc!=5)
	{
		cout <<"Usage: mapGoGene geneassociationfile geneontologyfile gotermtype outputsuff" << endl;
		return 0;
	}
	Framework fw;
	fw.init(argv[1],argv[2],argv[3]);
	fw.assignGoToTissue();
	fw.getCurrentNodes_Manual(argv[4]);
	return 0;
}
