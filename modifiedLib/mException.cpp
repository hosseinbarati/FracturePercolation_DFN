#include "mException.h"



mException::mException(string errorMessage) : msg(errorMessage){}

mException::~mException(){

}


const char* mException::what() const throw() {
	return msg.c_str();
}