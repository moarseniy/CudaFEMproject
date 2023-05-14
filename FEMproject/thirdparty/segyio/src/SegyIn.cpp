/***********************************************************************************
Author  : yang.deng
Email   : dyc5810911@126.com
Version : v1.0
Date    : 2016.09.09
************************************************************************************/
#include "SegyIn.h"

CSegyIn::CSegyIn(const std::string &strSegyFile,bool bigEdit) :
	CSegyIO(strSegyFile,bigEdit){}

CSegyIn::~CSegyIn(){}

bool CSegyIn::open()
{
	fin.open(segyFile.c_str(), std::ios::binary);
	if (!fin.is_open())
	{
		setErrorMessage("failed to open segy file");
		return isOpen = false;
	}
	else
		return isOpen = true;
}

bool CSegyIn::read()
{
	if (!isOpen)
	{
		setErrorMessage("failed to open segy file");
		return false;
	}
	/* only read one time */
  static bool isRead = false;
  if (isRead)
    return true;
	isRead = true;
	int index;
	size = getFileSize();   /* file size : total bytes */
	readSEGYFileHeader();   /* read segy file header */
	readSEGYTraceHeader();  /* read fisrt trace header to get basic information */
	for (index = 0; index < traceNumber; index++)
	{
		if (index != 0)
			readSEGYTraceHeader();
		readSEGYTrace();
  }
	return true;
}

void CSegyIn::close()
{
	fin.close();
}

long CSegyIn::getFileSize()
{
	long size;
	fin.seekg(0L, std::ios::end);
	size = static_cast<long>(fin.tellg());
	fin.seekg(0L, std::ios::beg);
	return size;
}

bool CSegyIn::readSEGYTrace()
{
	static int count = 0;
	if (count == 0)
	{   /* only alloacate one time */
		try{/* allocate trace memory */
			trace = new float[traceNumber*sampleNumber];
		}catch (std::bad_alloc)
		{
			setErrorMessage("failed to allocate memory for segy trace");
			return false;
		}
	}
	/* allocate tmp trace memory */
	float *_trace(nullptr); // single trace
	try{
		_trace = new float[sampleNumber];
	}catch (std::bad_alloc)
	{
		setErrorMessage("failed to allocate memory for segy trace(single)");
		return false;
	}
	/* allocate char buffer memeory */
	char *buf(nullptr);
	try{ 
		buf = new char[sampleNumber * _LONGSIZE];
	}catch (std::bad_alloc)
	{
		setErrorMessage("failed to allocate memory for buffer(readSEGYTrace)");
		return false;
	}
	/* read buffer data */
	fin.read(buf, sampleNumber*_LONGSIZE);
	/* only surport format = 5 or 2 two types */
	int index; int tmp;
	char _char2int[_LONGSIZE];
	for (index = 0; index < sampleNumber; index++)
	{
		getBuf(buf, _char2int, index*_LONGSIZE, _LONGSIZE);
		if (format == 2)
		{
			memcpy(&tmp, _char2int, _LONGSIZE);
			_trace[index] = static_cast<float>(tmp);
		}
		else
			memcpy(&_trace[index], _char2int, _LONGSIZE);
	}
	/* copy data to trace pointer */
	memcpy(trace + count*sampleNumber, _trace, sizeof(float)*sampleNumber);
	count++;
	delete[]_trace;
	delete[]buf;
	return true;
}

void CSegyIn::readSEGYFileHeader()
{
	char buffer[_FILE_HEAD_BYTE_LEN];
	fin.seekg(_FILE_HEAD_SKIP_LEN, std::ios::beg);
	fin.read(buffer, _FILE_HEAD_BYTE_LEN);

	int index;
	_short tmp;
	char _char2int[_LONGSIZE];
	char _char2short[_SHORTSIZE];
	
	long _offset = 0; // stream point offset
	for (index = 0; index < 3; index++)
	{
		getBuf(buffer, _char2int, _offset, _LONGSIZE);
		memcpy(&fHeader[index], _char2int, _LONGSIZE);
		_offset += _LONGSIZE;
	}
	for (; index < 27; index++)
	{
		getBuf(buffer, _char2short, _offset, _SHORTSIZE);
		memcpy(&tmp, _char2short, _SHORTSIZE);
		fHeader[index] = static_cast<int>(tmp);
		_offset += _SHORTSIZE;
	}
	format = fHeader[9];
	return;
}

bool CSegyIn::readSEGYTraceHeader()
{
	int index;
	long _offset = 0;
	char buffer[_TRACE_HEAD_BYTE_LEN] = {0};
	fin.read(buffer, _TRACE_HEAD_BYTE_LEN);
	char _char2short[_SHORTSIZE];
	char _char2int[_LONGSIZE];
	_short tmp;
/************************************************************
	[0] trace seq #
	[1] trace seq #
	[2] Original FID
	[3] Trace # within FID
	[4] Source ID
	[5] CDP #
	[6] Trace # within CDP 
**************************************************************/

	for(index = 0; index < 7; index++)
	{
		getBuf(buffer,_char2int,_offset,_LONGSIZE);
		memcpy(&tHeader[index], _char2int, _LONGSIZE);
		_offset += _LONGSIZE;
	}

	id = tHeader[2];

/**************************************************************
	[7] Trace ID code
	[8]
	[9]
	[10] Data Use 
***************************************************************/

	for (; index < 11; index++)
	{
		getBuf(buffer, _char2short,_offset, _SHORTSIZE);
		memcpy(&tmp, _char2short, _SHORTSIZE);
		tHeader[index] = static_cast<int>(tmp);
		_offset += _SHORTSIZE;
	}

/***************************************************************
	[11] Distance
	[12] Receiver Elevation
	[13] Source elevation
	[14] Shot Depth below surface
	[15] Datum elevation at receivers
	[16] Datum elevation at source
	[17] water depth at source
	[18] water depth at receiver  
**************************************************************/

	for (; index < 19; index++)
	{
		getBuf(buffer, _char2int, _offset, _LONGSIZE);
		memcpy(&tHeader[index], _char2int, _LONGSIZE);
		_offset += _LONGSIZE;
	}
		
/*************************************************************
	[19] elevation scaler
	[20] coordinate scaler
*************************************************************/

	for (; index < 21; index++)
	{
		getBuf(buffer, _char2short, _offset, _SHORTSIZE);
		memcpy(&tmp, _char2short, _SHORTSIZE);
		tHeader[index] = static_cast<int>(tmp);
		_offset += _SHORTSIZE;
	}
		
/****************************************************************
	[21] source x
	[22] source y
	[23] receiver x
	[24] receiver y
*****************************************************************/

	for (; index < 25; index++)
	{
		getBuf(buffer, _char2int, _offset, _LONGSIZE);
		memcpy(&tHeader[index], _char2int, _LONGSIZE);
		_offset += _LONGSIZE;
	}
		
/****************************************************************
	[25] Coordinate unit
	[38] num of samples
	[39] sample interval
*****************************************************************/

	for (; index < 67; index++)
	{
		getBuf(buffer, _char2short, _offset, _SHORTSIZE);
		memcpy(&tmp, _char2short, _SHORTSIZE);
		tHeader[index] = static_cast<int>(tmp);
		_offset += _SHORTSIZE;
	}

	for (; index < 71; index++)
	{
		getBuf(buffer, _char2int, _offset, _LONGSIZE);
		memcpy(&tHeader[index], _char2int, _LONGSIZE);
		_offset += _LONGSIZE;
	}

	if (tHeader[38] > _MAX_SAMPLE_NUMBER)
	{
		setErrorMessage("error : the number of samples exceed the maxinum limit 65535");
		return false;
	}
	if (!isEffectiveValue<int>(tHeader[38], "sample number"))
		return false;

	sampleNumber = tHeader[38];
	sampleInterval = static_cast<float>(tHeader[39]);
	sampleInterval /= 1000000.0;        /* convert micro-second to second */

	/* calculate traceNumber   */
	traceNumber = (size - 3600) / (240 + sampleNumber * _LONGSIZE);

	static int count = 0;
	if (count == 0)
	{
		try{/* allocate all trace header */
			tHeaders = new int[traceNumber*_TRACE_HEAD_LEN];
		}
		catch (std::bad_alloc)
		{
			setErrorMessage("failed to allocate memory for tHeader(all)");
			return false;
		}
	}
	
	memcpy(tHeaders + count*_TRACE_HEAD_LEN, tHeader, sizeof(int)*_TRACE_HEAD_LEN);
	count++;

	return true;
}

void CSegyIn::getBuf(const char* const src,
	char* const &buf,
	const int &offset,
	const int &len)
{
	int index;
	for (index = 0; index < len; index++)
		buf[index] = src[offset + index];
	if (!bigEdit)
		switchByte(buf, len);
	return;
}

void CSegyIn::getSegyTrace(float* const &p_trace) const
{
	if (!trace) return;
	memcpy(p_trace, trace, sizeof(float)*traceNumber*sampleNumber);
}

void CSegyIn::getSegyFileHeader(int* const &p_fHeader) const
{
	memcpy(p_fHeader, fHeader, sizeof(int)*_FILE_HEAD_LEN);
}

void CSegyIn::getSegyTraceHeaders(int* const &p_tHeaders) const
{
	memcpy(p_tHeaders, tHeaders, sizeof(int)*traceNumber*_TRACE_HEAD_LEN);
}
