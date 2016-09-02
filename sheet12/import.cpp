#include <fstream>
#include <string>
#include <iostream>
#include "helfer.h"

std::ifstream::pos_type filesize(const char* filename)
{
//	std::ifstream in(filename, std::ifstream::in | std::ifstream::binary);
	std::ifstream in(filename, std::ifstream::binary | std::ifstream::ate);
	//in.seekg(0,std::ifstream::end);
	return in.tellg();
}

int main ()
{
	std::ifstream import("msq.txt", std::ios::in);

// 	std::ifstream ifs("msq.txt");
// 	std::string content;
// 	content.assign( (std::istreambuf_iterator<char>(ifs) ), (std::istreambuf_iterator<char>() ) );
	
	int N = 25000;
	
	VecDoub a(N);
	
	int i = 0;
	
	while ( import.good() && !import.eof() && i < N)
		import >> a[i++];
	
	std::cout.precision(16);
	//std::cout << a[a.size()-1] << std::endl;
	

}
