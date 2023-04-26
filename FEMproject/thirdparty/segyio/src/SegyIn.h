/**********************************************************
Author : yang.deng
Email  : dyc5810911@126.com
***********************************************************/

#ifndef SEGYIN_H
#define SEGYIN_H
#include "SegyIO.h"
class CSegyIn :
	public CSegyIO
{
	CSegyIn(const CSegyIn & segyin);
	void readSEGYFileHeader();
	bool readSEGYTraceHeader();
	bool readSEGYTrace();
	void getBuf(const char* const src,
		char* const &buf,
		const int &offset,
		const int &len);
	long getFileSize();

public:
	~CSegyIn();
	CSegyIn(const std::string &strSegyFile,bool bigEdit = false);

	bool open();
	bool read();
	void close();

	inline int getSegyID() const{ return id; }
	inline int getSegyFormat() const { return format; }
	inline int getSampleNumber() const { return sampleNumber; }
	inline int getTraceNumber() const { return traceNumber; }
	inline float getSampleInterval() const { return sampleInterval; }

	void getSegyTrace(float* const &p_trace) const;
	void getSegyTraceHeaders(int* const &p_tHeaders) const;
	void getSegyFileHeader(int* const &p_fHeader) const;
};
#endif
