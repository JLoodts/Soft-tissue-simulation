#include<iostream>

int main()
{
	try
	{
		throw 5;
	}
	catch( int i)
	{
		std::cout << i << std::endl;
	}
	return 0;
}