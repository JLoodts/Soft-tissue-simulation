//////////////////////////////////////////////////////////////////
#ifndef INPUTFILE_H												//
#define INPUTFILE_H												//
// begin of InputFile.h											//
//////////////////////////////////////////////////////////////////

// The debugger can't handle symbols more than 255 characters long.
// STL often creates symbols longer than that.
// When symbols are longer than 255 characters, the warning is disabled.
#pragma warning(disable:4786)

#include <algorithm>// for Find

#include <string>	// for string
#include <vector>	// for vector
#include <fstream>	// for inlezen file

struct Scope
{
	std::string *begin;
	std::string *end;
};

typedef std::vector<std::string> STRINGVECTOR;

class InputFile
{
public:  
	InputFile(const char* fileName);
	void setScopeBegin(std::string *location){scope.begin = location;}
	void setScopeBegin(){scope.begin = location;}
	void setScope(std::string str);
//	void get(std::string str);
	// getInt returns -1 if the variable with name str could not be found
	int		getInt(std::string str);
	bool	getBool(std::string str);
	double	getDouble(std::string str);
	long	getLong(std::string str);
	int		getInt();
	bool	getBool();
	double	getDouble();
	long	getLong();
	std::string				*location;		
private:
	Scope					scope;
	STRINGVECTOR			text;
	STRINGVECTOR::iterator	textIt;
};


// exit(-200) : InputFile::setScope Error in determining the scope in inputfile
// exit(-201) : InputFile::setScope Error in searching for setting the scope
// exit(-202) : InputFile::getInt 
// exit(-203) : InputFile::getDouble
// exit(-204) : InputFile::getLong


InputFile::InputFile(const char* fileName)
{
	std::ifstream infile(fileName);
	std::istream_iterator<std::string> ifile(infile);
	std::istream_iterator<std::string> eos;

	std::copy(ifile, eos, std::inserter(text, text.begin()));
	scope.begin = text.begin();
	scope.end	= text.end();  
	location	= text.begin();
}

void InputFile::setScope(std::string str1)

	/* checks wether string str1 is part of text
	if so, then the scope is put on that section
	if not, the scope is kept on the entire text
	*/
{
	location = std::find(text.begin(), text.end(), str1);
	if (*location == str1) 
	{
		scope.begin = location;
		std::string str2("-end");
		str2.append(str1); 
		std::string *endLocation;
		endLocation = std::find(scope.begin, text.end(), str2);
		if (*endLocation == str2) scope.end = endLocation;
		else exit(-200);
	}
	else
	{
		exit(-201);
	}
}

int InputFile::getInt(std::string str)
{
	location = std::find(scope.begin, scope.end, str);
	if (*location == str) 
	{
		location++; // to move the location from the string to its numeric value in the file
		int value = atoi((*location).c_str());
		location++;
		return value;
	}
	else 
	{
		exit(-202);
		return 0;
	}
}

bool InputFile::getBool(std::string str)
{ // read an integer and return false if it equals 0
	location = std::find(scope.begin, scope.end, str);
	if (*location == str) 
	{
		location++; // to move the location from the string to its numeric value in the file
		int value = atoi((*location).c_str());
		location++;
		if(value==0) {
			return false;
		} else {
			return true;
		}
	}
	else 
	{
		exit(-202);
		return false;
	}
}

double InputFile::getDouble(std::string str)
{
	location = std::find(scope.begin, scope.end, str);
	if (*location == str) 
	{
		location++;
		double value = atof((*location).c_str());
		location++;
		return value;
	}
	else 
	{
		exit(-203);
		return 0;
	}
}

long InputFile::getLong(std::string str)
{
	location = std::find(scope.begin, scope.end, str);
	if (*location == str) 
	{
		location++;
		long value = atol((*location).c_str());
		location++;
		return value;
	}
	else 
	{
		exit(-204);
		return 0;
	}
}

int InputFile::getInt()
{
	int value = atoi((*location).c_str());
	location++;
	return value;
}

bool InputFile::getBool()
{// read an integer and return false if it equals 0
	int value = atoi((*location).c_str());
	location++;
	if(value==0) {
		return false;
	} else {
		return true;
	}
}

double InputFile::getDouble()
{
	double value = atof((*location).c_str());
	location++;
	return value;
}

long InputFile::getLong()
{
	long value = atol((*location).c_str());
	location++;
	return value;
}


//////////////////////////////////////////////////////////////////
#endif															//
// end of InputFile.h											//
//////////////////////////////////////////////////////////////////
