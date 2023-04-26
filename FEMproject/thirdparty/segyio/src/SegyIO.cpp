/***********************************************************************************
Author  : yang.deng
Email   : dyc5810911@126.com
Version : v1.0
Date    : 2016.09.09
************************************************************************************/

#include "SegyIO.h"

CSegyIO::CSegyIO(std::string strSEGYFile, bool _bigEdit) : segyFile(strSEGYFile),
	bigEdit(_bigEdit),isOpen(false),id(0),sampleInterval(0),traceNumber(0),
	sampleNumber(0), format(0), trace(nullptr), tHeaders(nullptr){
	memset(fHeader, 0, sizeof(fHeader));
	memset(tHeader, 0, sizeof(tHeader));
}

void CSegyIO::switchByte(char* const c,const int &len)
{
	int index;char tmp;
	for (index = 0; index < len / 2; index++)
	{
		tmp = c[index];
		c[index] = c[len - index - 1];
		c[len - index - 1] = tmp;
	}
}

bool CSegyIO::open(){ return true; }

void CSegyIO::close(){}

CSegyIO::~CSegyIO()
{
	if (trace){ delete[]trace; trace = nullptr; }
	if (tHeaders){ delete[]tHeaders; tHeaders = nullptr; }
}
