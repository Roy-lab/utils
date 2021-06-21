#include "AnnotTerm.H"

AnnotTerm::AnnotTerm()
{
}

AnnotTerm::~AnnotTerm()
{
}

//Although we are saying tissue, really we mean tissue and stage development
int 
AnnotTerm::setAnnotTermName(const char* aName)
{
	termName.append(aName);
	return 0;
}

const char*
AnnotTerm::getAnnotTermName()
{
	return termName.c_str();
}

int 
AnnotTerm::addMemberGene(const char* aName)
{
	string geneName(aName);
	memberGenes[geneName]=0;
	return 0;
}

map<string,int>& 
AnnotTerm::getMemberGenes()
{
	return memberGenes;
}

bool
AnnotTerm::isMember(const string& geneName)
{
	if(memberGenes.find(geneName)==memberGenes.end())
	{
		return false;
	}
	return true;
}
