#pragma once
#include<exception>
#include<iostream>
using namespace std;

class mException : public exception
{
public:
	mException(string);
	~mException();

	const char* what() const throw();

private:
	string msg;
};

