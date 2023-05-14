/***********************************************************************
Author : yang.deng
Email  : dyc5810911@126.com
***********************************************************************/

#include "SegyOut.h"

CSegyOut::CSegyOut(const std::string &strSegyFile,
	const int &mTraceNumber,const bool &bigEdian):
	CSegyIO(strSegyFile,bigEdian),isWrite(false){
	traceNumber = mTraceNumber;
}

CSegyOut::~CSegyOut(){}

bool CSegyOut::open()
{
	fout.open(segyFile.c_str(), std::ios::binary);
	if (!fout.is_open())
	{
		setErrorMessage("failed to open segy file ");
		return isOpen = false;
	}
	else
		return isOpen = true;
}

void CSegyOut::close()
{
	fout.close();
	isOpen = false;
	isWrite = false;
}

bool CSegyOut::write()
{
	if (isWrite)
		return true;

	int ir(0);
	if (!isEffectiveValue<int>(traceNumber, "trace number"))
		return false;
	if (!isEffectiveValue<int>(sampleNumber, "samples"))
		return false;
//	if (!isEffectiveValue<int>(id, "segy id"))
//		return false;
	if (!isEffectiveValue<float>(sampleInterval, "sample rate"))
		return false;
//	if (format != 5 && format != 2)
//	{
//		setErrorMessage("do not surpport this format of segy file");
//		return false;
//	}
	setSEGYFileHeader();

	setSEGYTraceHeaders();

	writeSEGYFileHeader();
	for (;ir < traceNumber; ir++)
	{
		writeSEGYTraceHeader();
		writeSEGYTrace();
	}
	isWrite = true;
	return true;
}

void CSegyOut::writeSEGYFileHeader()
{
	int index(0);
	long _offset(0);
	_short tmp;

	char _buf1[_FILE_HEAD_SKIP_LEN] = { 0 };
	fout.write(_buf1, 3200);

	char _buf2[_FILE_HEAD_BYTE_LEN] = { 0 };
	char _char2int[_LONGSIZE];
	char _char2short[_SHORTSIZE];

	for (; index < 3; index++)
	{
		memcpy(_char2int, &fHeader[index], _LONGSIZE);
		setBuffer(_char2int, _buf2, _offset, _LONGSIZE);
		_offset += _LONGSIZE;
	}
	for (; index < 27; index++)
	{
		tmp = static_cast<_short>(fHeader[index]);
		memcpy(_char2short, &tmp, _SHORTSIZE);
		setBuffer(_char2short, _buf2, _offset, _SHORTSIZE);
		_offset += _SHORTSIZE;
	}

	fout.write(_buf2, _FILE_HEAD_BYTE_LEN);
	return;
}

void CSegyOut::writeSEGYTraceHeader()
{
	int index(0);
	_short tmp;
	static int count = 0;
	long _offset = 0;
	long indx = count * _TRACE_HEAD_LEN;
	char _char2int[_LONGSIZE];
	char _char2short[_SHORTSIZE];
	char _strtHeader[_TRACE_HEAD_BYTE_LEN] = { 0 };
	for (; index < 7; index++)
	{
		memcpy(_char2int, &tHeaders[index+indx], _LONGSIZE);
		setBuffer(_char2int, _strtHeader, _offset, _LONGSIZE);
		_offset += _LONGSIZE;
	}
	for (; index < 11; index++)
	{
		tmp = static_cast<_short>(tHeaders[index+indx]);
		memcpy(_char2short, &tmp,_SHORTSIZE);
		setBuffer(_char2short, _strtHeader, _offset, _SHORTSIZE);
		_offset += _SHORTSIZE;
	}
	for (; index < 19;index++)
	{
		memcpy(_char2int, &tHeaders[index + indx], _LONGSIZE);
		setBuffer(_char2int, _strtHeader, _offset, _LONGSIZE);
		_offset += _LONGSIZE;
	}
	for (; index < 21; index++)
	{
		tmp = static_cast<_short>(tHeaders[index+indx]);
		memcpy(_char2short, &tmp, _SHORTSIZE);
		setBuffer(_char2short, _strtHeader, _offset, _SHORTSIZE);
		_offset += _SHORTSIZE;
	}
	for (; index < 25; index++)
	{
		memcpy(_char2int, &tHeaders[index + indx], _LONGSIZE);
		setBuffer(_char2int, _strtHeader, _offset, _LONGSIZE);
		_offset += _LONGSIZE;
	}
	for (; index < 67; index++)
	{
		tmp = static_cast<_short>(tHeaders[index + indx]);
		memcpy(_char2short, &tmp, _SHORTSIZE);
		setBuffer(_char2short, _strtHeader, _offset, _SHORTSIZE);
		_offset += _SHORTSIZE;
	}
	for (; index < 71; index++)
	{
		memcpy(_char2int, &tHeaders[index + indx], _LONGSIZE);
		setBuffer(_char2int, _strtHeader, _offset, _LONGSIZE);
		_offset += _LONGSIZE;
	}
	fout.write(_strtHeader, _TRACE_HEAD_BYTE_LEN);
	count++;
	return;
}

bool CSegyOut::writeSEGYTrace()
{
	static int count = 0;
	char _char2float[_LONGSIZE];
	char *strTrace = nullptr;
	long _offset = 0;
	try{
		strTrace = new char[sampleNumber*_LONGSIZE];
	}catch (std::bad_alloc)
	{
		setErrorMessage("failed to allocate memory for strTrace");
		return false;
	}
	int i(0);
	int index = count * sampleNumber;
	int tmp;
	for (; i < sampleNumber; i++)
	{
		if (format == 5)
			memcpy(_char2float, &trace[index + i], _LONGSIZE);
		else{
			tmp = static_cast<int>(trace[index + i]);
			memcpy(_char2float, &tmp, _LONGSIZE);
		}
		setBuffer(_char2float, strTrace, _offset, _LONGSIZE);
		_offset += _LONGSIZE;
	}
	count++;
	fout.write(strTrace, sampleNumber*_LONGSIZE);
	delete []strTrace;
	strTrace = nullptr;
	return true;
}

bool CSegyOut::write(const int* const &ptHeaders,
	const int* const &pfHeader,
	const float* const &ptrace,
	const int &mtraceNumber)
{
	int ir(0);
	if (isWrite)
		return true;
	traceNumber = mtraceNumber;
	sampleNumber = ptHeaders[38];
	format = pfHeader[9];
	if (!isEffectiveValue<int>(traceNumber, "trace number"))
		return false;
	if (!isEffectiveValue<int>(sampleNumber, "samples"))
		return false;
	memcpy(fHeader, pfHeader, sizeof(float)*_FILE_HEAD_LEN);
	if (tHeaders){
		delete[]tHeaders;
		tHeaders = nullptr;
	}try{
		tHeaders = new int[traceNumber*_TRACE_HEAD_LEN];
	}catch (std::bad_alloc)
	{
		setErrorMessage("failed to allocate memory for tHeaders");
		return false;
	}
	memcpy(tHeaders, ptHeaders, traceNumber*_TRACE_HEAD_LEN*sizeof(int));
	if (trace){
		delete[]trace;
		trace = nullptr;
	}try{
		trace = new float[traceNumber*sampleNumber];
	}catch (std::bad_alloc)
	{
		setErrorMessage("failed to allocate memory for trace");
		return false;
	}
	memcpy(trace, ptrace, traceNumber*sampleNumber*sizeof(float));
	
	writeSEGYFileHeader();
	
	for (; ir < traceNumber; ir++)
	{
		writeSEGYTraceHeader();
		writeSEGYTrace();
	}
	isWrite = true;
	return true;
}

void CSegyOut::setSEGYFileHeader()
{
	fHeader[3] = traceNumber;
	fHeader[5] = static_cast<int>(sampleInterval * 1000000L);
	fHeader[7] = sampleNumber;
	fHeader[9] = format;
}

bool CSegyOut::setSEGYTraceHeaders()
{
	int ir,index;
	try{
		tHeaders = new int[traceNumber*_TRACE_HEAD_LEN];
	}catch (std::bad_alloc)
	{
		setErrorMessage("failed to allocate memory for tHeaders");
		return false;
	}
	for (ir = 0; ir < traceNumber; ir++)
	{
		index = ir * _TRACE_HEAD_LEN;
		memset(tHeader, 0, sizeof(tHeader));

		tHeader[2] = id;
		tHeader[38] = sampleNumber;
		/* convert second to micro-second */
		tHeader[39] = static_cast<int>(sampleInterval*1000000L);

		memcpy(tHeaders + index, tHeader, sizeof(int)*_TRACE_HEAD_LEN);
	}
	return true;
}

bool CSegyOut::setSegyTrace(const float* const &p_trace,const int &length)
{
	/* avoid repeat-setting to cause memory leak */
	if (trace)
	{
		delete[]trace;
		trace = nullptr;
	}
	try{
		trace = new float[length];
	}catch (std::bad_alloc)
	{
		setErrorMessage("failed to allocate memory for trace");
		return false;
	}
  memcpy(trace,p_trace,length*sizeof(float));
	return true;
}

void CSegyOut::setBuffer(char* const &src,
	char* const &buf,
	const int &offset,
	const int &len)
{
	int index(0);
	if (!bigEdit)
		switchByte(src, len);
	for (; index < len; index++)
		buf[offset + index] = src[index];
	return;
}
