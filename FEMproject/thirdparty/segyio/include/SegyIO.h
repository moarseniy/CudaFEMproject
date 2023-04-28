/***********************************************************************************
Author  : yang.deng
Email   :dyc5810911@126.com
Version : v1.0
Date    : 2016.09.09
************************************************************************************/
#ifndef SEGYIO_H
#define SEGYIO_H
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
class CSegyIO
{
	/* private copy-constructor function */
	CSegyIO(const CSegyIO &segyio);
protected:
	
	enum       { _FILE_HEAD_LEN = 27 };
	enum       { _TRACE_HEAD_LEN = 71 };
	enum       { _FILE_HEAD_SKIP_LEN = 3200 };
	enum       { _FILE_HEAD_BYTE_LEN = 400 };
	enum       { _TRACE_HEAD_BYTE_LEN = 240 };
	enum       { _LONGSIZE = 4 };
	enum       { _SHORTSIZE = 2 };
	enum       { _MAX_SAMPLE_NUMBER = 65535 };

	int        id;
	int        sampleNumber;                     /* segy trace sample number   */
	int        traceNumber;                      /* segy trace number          */
	float      sampleInterval;                   /* segy trace sample interval */
	int        format;                           /* segy file format           */
	int        tHeader[_TRACE_HEAD_LEN];         /* trace header     */
	int        fHeader[_FILE_HEAD_LEN];          /* file  header     */
	float      *trace;                           /* trace data       */
	int        *tHeaders;                        /* trace header all */
	
	long        size;
	std::ifstream fin;
	std::ofstream fout;
	std::string errorMessage;
	std::string segyFile;
	bool bigEdit;
	bool isOpen;

	inline void setErrorMessage(const std::string &strErrorMessage){
		errorMessage = strErrorMessage; 
	};
	void switchByte(char* const c,const int &len);
	typedef unsigned short _short;

	template <typename Type>
	bool isEffectiveValue(Type _value,const std::string &type)
	{
		if (_value <= 0)
		{
			setErrorMessage("Error : the value of " + type + " should be more than zero");
			return false;
		}
		else
			return true;
	}

public:
	CSegyIO(std::string strSEGYfile,bool bigEdit = false);
	inline std::string getErrorMessage() const { return errorMessage; }
	virtual bool open();
	virtual void close();
	virtual ~CSegyIO() = 0;
};
#endif
