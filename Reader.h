#pragma once
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>

#include "Model.h"

class Reader
{
public:
	Reader() {}
	void read(string fileName, Model& model);
};


