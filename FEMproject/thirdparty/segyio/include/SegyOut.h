/***********************************************************************************
Author  : yang.deng
Email   : dyc5810911@126.com
Version : v1.0
Date    : 2016.09.09
************************************************************************************/
#ifndef SEGYOUT_H
#define SEGYOUT_H
#include "SegyIO.h"
class CSegyOut :
	public CSegyIO
{
	bool isWrite;
	bool writeSEGYTrace();
	void writeSEGYFileHeader();
	void writeSEGYTraceHeader();
	void setBuffer(char* const &src,
		char* const &buf,
		const int &offset,
		const int &len);
	void setSEGYFileHeader();
	bool setSEGYTraceHeaders();
	CSegyOut(const CSegyOut &segyout);
public:
	CSegyOut(const std::string &strSegyFile,
		const int &traceNumber = 0,
		const bool &bigEdit = false);
	~CSegyOut();
	bool open();
	void close();

	bool write();
	bool write(const int* const &ptHeaders,
		const int* const &pfHeader,
		const float* const &ptrace,
		const int &mtraceNumber);

	inline void setSegyID(const int &mID) { id = mID; }
	inline void setSegyFormat(const int &mFormat){ format = mFormat; }
	inline void setSampleInterval(const float &mSampleInterval){ sampleInterval = mSampleInterval; }
	inline void setSampleNumber(const int &mSampleNumber){ sampleNumber = mSampleNumber; }
	inline void setTraceNumber(const int &mTraceNumber){ traceNumber = mTraceNumber; }

	bool setSegyTrace(const float* const &trace,const int &length);

};
#endif
