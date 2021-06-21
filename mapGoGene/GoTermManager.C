#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>

#include "AnnotTerm.H"
#include "GoTerm.H"
#include "GoTermManager.H"


GoTermManager::GoTermManager()
{
}


GoTermManager::~GoTermManager()
{
}


int 
GoTermManager::setTermType(const char* gotype)
{
	reqTermType.append(gotype);
	return 0;
}

int 
GoTermManager::readGoTerms(const char* geneontfile)
{
        //ifstream inFile(OBO_FILE);
        ifstream inFile(geneontfile);
	char buffer[4096];
	int lineCnt=0;
	int obsoleteTerm=0;
	string goid;
	string goname;
	string gotermtype;
	bool isObsolete=false;
	while(inFile.good())
	{
		inFile.getline(buffer,4095);
		if(strlen(buffer)==0)
		{
			continue;
		}
		if(lineCnt==0)
		{
			lineCnt++;
			continue;
		}
		if( (strstr(buffer,"[Term]")!=NULL) || (strstr(buffer,"[Typedef]")!=NULL))
		{
			//This marks the end of the existing term (if any) and the start of a new one
			if((strcmp(gotermtype.c_str(),reqTermType.c_str())==0) &&(!isObsolete))
			{
				if(strcmp(goid.c_str(),"GO:0065010")==0)
				{
					cout <<"Found term " << goid.c_str() << " " << goname.c_str() << endl;
				}
				GoTerm* goTerm=new GoTerm;
				goTerm->setAnnotID(goid.c_str());
				goTerm->setName(goname.c_str());
				idToTermMap[goid]=goTerm;
				nameToTermMap[goname]=goTerm;
			}
			if(isObsolete)
			{
				obsoleteTerm++;
			}
			goid.clear();
			goname.clear();
			gotermtype.clear();
			isObsolete=false;
		}
		else if(strstr(buffer,"id: GO:")!=NULL)
		{
			if(strstr(buffer,"alt_id:")==NULL)
			{
				char* pos=strstr(buffer,"GO:");
				goid.append(pos);
			}
		}
		else if(strstr(buffer,"name:")!=NULL)
		{
			char* pos=strchr(buffer,' ');
			pos++;
			while((*pos)==' ')
			{
				pos++;
			}
			goname.append(pos);
		}
		else if(strstr(buffer,"namespace:")!=NULL)
		{
			char* pos=strchr(buffer,' ');
			pos++;
			while((*pos)==' ')
			{
				pos++;
			}
			gotermtype.append(pos);
		}
		else if(strstr(buffer,"is_obsolete:")!=NULL)
		{
			if(strstr(buffer,"true")!=NULL)
			{
				isObsolete=true;
			}
		}
		
	}
	if(isObsolete)
	{
		obsoleteTerm++;
	}
	
	cout <<"Found " << obsoleteTerm <<  " obsolete terms " << endl;
	inFile.close();
    saveNameIDMap();
	return 0;
}


int 
GoTermManager::readGoTermMap(const char* geneontfile)
{
        ifstream inFile(geneontfile);
        //ifstream inFile(OBO_FILE);
	char buffer[4096];
	int lineCnt=0;
	int obsoleteTerm=0;
	string goid;
	string gotermtype;
	map<string,int> isaparents;
	bool isObsolete=false;
	while(inFile.good())
	{
		inFile.getline(buffer,4095);
		if(strlen(buffer)==0)
		{
			continue;
		}
		if(lineCnt==0)
		{
			lineCnt++;
			continue;
		}
		if((strstr(buffer,"[Term]")!=NULL) || (strstr(buffer,"[Typedef]")!=NULL))
		{
			//This marks the end of the existing term (if any) and the start of a new one
			if((strcmp(gotermtype.c_str(),reqTermType.c_str())==0) &&(!isObsolete))
			{
				GoTerm* childTerm=idToTermMap[goid];
				for(map<string,int>::iterator aIter=isaparents.begin();aIter!=isaparents.end();aIter++)
				{
					if(idToTermMap.find(aIter->first)==idToTermMap.end())
					{
						cout <<"No isa parent " << aIter->first.c_str() << " for " << goid.c_str() << endl;
						return -1;
					}
					GoTerm* parentTerm=idToTermMap[aIter->first.c_str()];
					if(aIter->second==0)
					{
						parentTerm->setChild(childTerm);
					}
					else
					{
						parentTerm->setChild(childTerm,true);
					}
					childTerm->setParent(parentTerm);
				}
			}
			if(isObsolete)
			{
				obsoleteTerm++;
			}
			goid.clear();
			gotermtype.clear();
			isObsolete=false;
			isaparents.clear();
		}
		else if(strstr(buffer,"id: GO:")!=NULL)
		{
			if(strstr(buffer,"alt_id:")==NULL)
			{
				char* pos=strstr(buffer,"GO:");
				goid.append(pos);
			}
		}
		else if(strstr(buffer,"is_a:")!=NULL)
		{
			char* pos=strchr(buffer,' ');
			pos++;
			while((*pos)==' ')
			{
				pos++;
			}
			char* endpos=strchr(pos,'!');
			if(endpos==NULL)
			{
				cout <<"No ! found with is_a string" << endl;
				return -1;
			}
			endpos--;
			while((*endpos)==' ')
			{
				endpos--;
			}
			*(endpos+1)='\0';
			string aparent(pos);
			isaparents[aparent]=0;
		}
		else if(strstr(buffer,"relationship: part_of")!=NULL)
		{
			char* ppos=strstr(buffer,"part_of");
			char* pos=strchr(ppos,' ');
			pos++;
			while((*pos)==' ')
			{
				pos++;
			}
			char* endpos=strchr(pos,'!');
			if(endpos==NULL)
			{
				cout <<"No ! found with is_a string" << endl;
				return -1;
			}
			endpos--;
			while((*endpos)==' ')
			{
				endpos--;
			}
			*(endpos+1)='\0';
			string aparent(pos);
			isaparents[aparent]=1;
		}
		else if(strstr(buffer,"namespace:")!=NULL)
		{
			char* pos=strchr(buffer,' ');
			pos++;
			while((*pos)==' ')
			{
				pos++;
			}
			gotermtype.append(pos);
		}
		else if(strstr(buffer,"is_obsolete:")!=NULL)
		{
			if(strstr(buffer,"true")!=NULL)
			{
				isObsolete=true;
			}
		}
	}
	return 0;
}

int GoTermManager::saveNameIDMap(){
        ofstream outmap("nameTerms.txt");
        for (map<string, GoTerm*>::iterator aIter=nameToTermMap.begin();aIter!=nameToTermMap.end();aIter++){
            GoTerm* aterm=aIter->second;
            outmap << aIter->first << "\t" << aterm->getAnnotID()<<endl;
        }
        outmap.close();    
        return 0;
}

int 
GoTermManager::showHierarchy(const char* aFName, int level)
{
	ofstream oFile(aFName);
	//Get all parents: This is the first level
	map<string,GoTerm*> currNodes;
	int currLevel=0;
	for(map<string,GoTerm*>::iterator aIter=nameToTermMap.begin();aIter!=nameToTermMap.end();aIter++)
	{
		GoTerm* aterm=aIter->second;
		if(aterm->getParentCnt()==0)
		{
			currNodes[aterm->getName()]=aterm;
		}
	}
	cout <<"Found " << currNodes.size() << " at level " << currLevel << endl;
	while(currLevel<=level)
	{
		map<string,GoTerm*> tempNodes;
		oFile <<"Level " << currLevel << " : ";
		for(map<string,GoTerm*>::iterator aIter=currNodes.begin();aIter!=currNodes.end();aIter++)
		{
			if(aIter!=currNodes.begin())
			{
				oFile <<",";
			}
			oFile <<" " << aIter->second->getName().c_str();
			map<string,GoTerm*>& children=aIter->second->getChildren();
			if(children.size()==0)
			{
				tempNodes[aIter->first]=aIter->second;
			}
			for(map<string,GoTerm*>::iterator bIter=children.begin();bIter!=children.end();bIter++)
			{
				tempNodes[bIter->first]=bIter->second;
			}
		}
		oFile << endl;

		//Clear out the currNodes
		currNodes.clear();
		for(map<string,GoTerm*>::iterator aIter=tempNodes.begin();aIter!=tempNodes.end();aIter++)
		{
			currNodes[aIter->first]=aIter->second;
		}
		currLevel++;
		cout <<"Found " << tempNodes.size() << " at level " << currLevel << endl;
		tempNodes.clear();
	}
	oFile.close();
	return 0;
}

GoTerm*
GoTermManager::getTermForName(const char* aName)
{
	string key(aName);
	if(nameToTermMap.find(key)==nameToTermMap.end())
	{
		return NULL;
	}
	return nameToTermMap[key];
}


GoTerm*
GoTermManager::getTermForID(const char* aid)
{
	string key(aid);
	if(idToTermMap.find(key)==idToTermMap.end())
	{
		return NULL;
	}
	return idToTermMap[key];
}




map<string,GoTerm*>& 
GoTermManager::getAllGoTerms()
{
	return nameToTermMap;
}

int
GoTermManager::assignLevelFromRoot()
{
	//First identify root node
	GoTerm* rootNode=NULL;
	map<string,GoTerm*>::iterator gIter=nameToTermMap.begin();
	while((rootNode==NULL) && (gIter!=nameToTermMap.end()))
	{
		GoTerm* aNode=gIter->second;
		if(aNode->getParents().size()==0)
		{
			rootNode=aNode;
		}
		gIter++;
	}
	
	//Do a traversal from the rootnode of the levels
	rootNode->assignLevelFromRoot(0);
	return 0;
}
