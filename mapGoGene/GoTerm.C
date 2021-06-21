#include <iostream>
#include <string.h>
#include "AnnotTerm.H"
#include "GoTerm.H"

GoTerm::GoTerm()
{
	levelFromRoot=-1;
	annotterm=NULL;
}

GoTerm::~GoTerm()
{
}


int
GoTerm::setID(int aID)
{
	id=aID;
	return 0;
}

int 
GoTerm::setAnnotID(const char* anID)
{
	annotId.append(anID);
	return 0;
}

int 
GoTerm::setName(const char* aName)
{
	goName.append(aName);
	char* tempBuff=new char[strlen(aName)+1];
	strcpy(tempBuff,aName);
	char* tok=strtok(tempBuff," ");
	while(tok!=NULL)
	{
		string aKey(tok);
		wordSet[aKey]=0;
		tok=strtok(NULL," ");
	}
	delete[] tempBuff;
	return 0;
}

int 
GoTerm::setChild(GoTerm* childTerm,bool partof)
{
	if(children.find(childTerm->getAnnotID())==children.end())
	{
		children[childTerm->getAnnotID()]=childTerm;
	}
	else
	{
		cout << childTerm->getAnnotID().c_str() << " is already a child of " << goName.c_str() << endl;
	}
	if(partof)
	{
		partOfChild[childTerm->getAnnotID()]=0;
	}
	return 0;
}

int 
GoTerm::setParent(GoTerm* parentTerm)
{
	if(parents.find(parentTerm->getAnnotID())==parents.end())
	{
		parents[parentTerm->getAnnotID()]=parentTerm;
	}
	else
	{
		cout << parentTerm->getAnnotID() << " is already a parent of " << goName.c_str() << endl;
	}
	return 0;
}

int
GoTerm::getID()
{
	return id;
}

string&
GoTerm::getAnnotID()
{
	return annotId;
}

string& 
GoTerm::getName()
{
	return goName;
}


int 
GoTerm::getParentCnt()
{
	return parents.size();
}

int 
GoTerm::getChildCnt()
{
	return children.size();
}

bool 
GoTerm::isChild(const char* aName)
{
	string key(aName);
	if(children.find(key)!=children.end())
	{
		return true;
	}
	return false;
}

bool 
GoTerm::isDescendent(const char* aName)
{
	if(children.size()==0)
	{
		return false;
	}
	string key(aName);
	if(children.find(key)!=children.end())
	{
		return true;
	}
	bool found=false;
	map<string,GoTerm*>::iterator aIter=children.begin();
	while((!found) &&(aIter!=children.end()))
	{
		GoTerm* iTerm=aIter->second;
		found=iTerm->isDescendent(aName);
		aIter++;
	}
	return found;
}

map<string,GoTerm*>&
GoTerm::getParents()
{
	return parents;
}

map<string,GoTerm*>&
GoTerm::getChildren()
{
	return children;
}

map<string,int>&
GoTerm::getWordSet()
{
	return wordSet;
}


int 
GoTerm::setAnnotTerm(AnnotTerm* tPtr)
{
	if(annotterm==NULL)
	{
		annotterm=tPtr;
	}
	else
	{
		cout <<"Term " << annotterm->getAnnotTermName() << " already assigned to term " << goName.c_str() << endl;
	}
	return 0;
}

AnnotTerm* 
GoTerm::getAnnotTerm()
{
	return annotterm;
}

int
GoTerm::createGeneList()
{
	//Go down the hierachy and create the gene lists
	if(annotterm!=NULL)
	{
		map<string,int>& alist= annotterm->getMemberGenes();
		for(map<string,int>::iterator aIter=alist.begin();aIter!=alist.end();aIter++)
		{
			geneList[aIter->first]=0;
		}
	}
	//For each child that has tissue infor, get their genes
	for(map<string,GoTerm*>::iterator cIter=children.begin();cIter!=children.end();cIter++)
	{
		GoTerm* achild=cIter->second;
		if(achild->getAnnotTerm()==NULL)
		{
		//	continue;
		}
		achild->createGeneList();
		map<string,int>& childGeneList=achild->getGeneList();
		int newChildGenes=0;
		for(map<string,int>::iterator aIter=childGeneList.begin();aIter!=childGeneList.end();aIter++)
		{
			if(geneList.find(aIter->first)==geneList.end())
			{
				newChildGenes++;
			}
			int level=aIter->second;
			if(partOfChild.find(aIter->first)==partOfChild.end())
			{
				level=aIter->second+1;
			}
			if(geneList.find(aIter->first)!=geneList.end())
			{
				if(level>geneList[aIter->first])
				{
					geneList[aIter->first]=level;
				}
			}
			else
			{
				geneList[aIter->first]=level;
			}
		}
		if(newChildGenes>0)
		{
			//cout <<"Found " << newChildGenes <<" new genes in child " << achild->getName().c_str() 
			 //<< " absent in " << goName.c_str() << endl;
		}
	}
	return 0;
}

map<string,int>&
GoTerm::getGeneList()
{
	return geneList;
}

double
GoTerm::getGeneOverlap(map<string,int>& aList)
{
	double overlapGene=0;
	for(map<string,int>::iterator aIter=aList.begin();aIter!=aList.end();aIter++)
	{
		if(geneList.find(aIter->first)!=geneList.end())
		{
			overlapGene++;
		}
	}
	if(geneList.size()<aList.size())
	{
		overlapGene=overlapGene/(double)geneList.size();
	}
	else
	{
		overlapGene=overlapGene/(double)aList.size();
	}
	return overlapGene;
}

int
GoTerm::getGeneCnt()
{
	return geneList.size();
}

int
GoTerm::showAncestory(int level)
{
	for(map<string,GoTerm*>::iterator aIter=parents.begin();aIter!=parents.end();aIter++)
	{
		cout << goName.c_str() << " 's Parent at level " << level << " :" << aIter->second->getName().c_str() << endl;
		aIter->second->showAncestory(level+1);
	}
	return 0;
}

int
GoTerm::showChildTree(int currentLevel,int finalLevel)
{
	if(currentLevel> finalLevel)
	{
		return 0;
	}
	cout <<goName.c_str() << "'s children at level " << currentLevel << endl;
	for(map<string,GoTerm*>::iterator aIter=children.begin();aIter!=children.end();aIter++)
	{
		cout <<" : " << aIter->second->getName().c_str();
	}
	cout << endl;
	for(map<string,GoTerm*>::iterator aIter=children.begin();aIter!=children.end();aIter++)
	{
		aIter->second->showChildTree(currentLevel+1,finalLevel);
	}
	return 0;
}

int
GoTerm::getChildren(map<string,GoTerm*>& currNodeSet, int currLevel, int requiredLevel)
{
	//Only children at the requiredLevel are required.
	for(map<string,GoTerm*>::iterator aIter=children.begin();aIter!=children.end();aIter++)
	{
		if(currLevel<requiredLevel)
		{
			aIter->second->getChildren(currNodeSet,currLevel+1,requiredLevel);
		}
		else
		{
			if(aIter->second->getGeneCnt()>0)
			{
				currNodeSet[aIter->first]=aIter->second;
			}
		}
	}
	return 0;
}

int
GoTerm::getAncestors(map<string,GoTerm*>& ancestors)
{
	for(map<string,GoTerm*>::iterator aIter=parents.begin();aIter!=parents.end();aIter++)
	{
		ancestors[aIter->first]=aIter->second;
	}
	for(map<string,GoTerm*>::iterator aIter=parents.begin();aIter!=parents.end();aIter++)
	{
		aIter->second->getAncestors(ancestors);
	}
	return 0;
}

int
GoTerm::assignLevelFromRoot(int alevel)
{
	if((levelFromRoot==-1)||(levelFromRoot>alevel))
	{
		levelFromRoot=alevel;
	}
	for(map<string,GoTerm*>::iterator cIter=children.begin();cIter!=children.end();cIter++)
	{
		GoTerm* child=cIter->second;
		child->assignLevelFromRoot(alevel+1);
	}
	return 0;
}

int
GoTerm::getLevelFromRoot()
{
	return levelFromRoot;
}
