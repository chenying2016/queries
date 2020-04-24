#ifndef MACROS_H
#define MACROS_H

#ifndef BEGIN_NCBI_SCOPE
#define BEGIN_NCBI_SCOPE
#endif

#ifndef END_NCBI_SCOPE
#define END_NCBI_SCOPE
#endif

#ifndef NCBI_XALGOWINMASK_EXPORT
#define NCBI_XALGOWINMASK_EXPORT
#endif

#ifndef NCBI_XALGODUSTMASK_EXPORT
#define NCBI_XALGODUSTMASK_EXPORT
#endif

#ifndef NCBI_XNCBI_EXPORT
#define NCBI_XNCBI_EXPORT
#endif

#ifndef BEGIN_NCBI_NAMESPACE
#define BEGIN_NCBI_NAMESPACE
#define END_NCBI_NAMESPACE
#define BEGIN_STD_NAMESPACE
#define NCBI_NS_NCBI
#define END_STD_NAMESPACE
#endif

#define USING_SCOPE(sc) 

#undef  __DEPRECATED

#include "def.h"

#include <vector>
using std::vector;
//typedef vector<char> CSeqVector;

#include<iostream>
#define Warning std::cerr 
#ifndef ERR_POST
#define ERR_POST
#define LOG_POST
#endif

#ifndef TSeqPos
typedef Uint4 TSeqPos;
#endif

using std::pair;
using std::make_pair;

#include <string>
using std::string;
extern std::string kEmptyStr;

#include <list>
using std::list;

#include <memory>
using std::auto_ptr;

#  define _TRACE(message)                               ((void)0)

#include <fstream>
#include <iostream>
#include <strstream>

#define CNcbiIfstream std::ifstream 
#define CNcbiOstrstream std::ostrstream
#define CNcbiIstream std::istream
#define CNcbiOstream std::ostream
#define NcbiCout std::cout
#define CNcbiOfstream std::ofstream
#define Error std::cerr

using std::istringstream;
using std::ostringstream;
using std::endl;
using std::hex;
using std::dec;
using std::flush;

#include <algorithm>
using std::fill;
using std::max_element;
using std::min_element;
using std::less;

#define kMax_UI8 0xFFFFFFFFFFFFFFFFULL
#define kMax_UI4 UINT4_MAX
#define kMax_UI1 255

typedef std::size_t SIZE_TYPE;
const SIZE_TYPE NPOS = static_cast<SIZE_TYPE>(-1);

class NCBI_XNCBI_EXPORT CNcbiOstrstreamToString
{
    CNcbiOstrstreamToString(const CNcbiOstrstreamToString&);
    CNcbiOstrstreamToString& operator= (const CNcbiOstrstreamToString&);
public:
    CNcbiOstrstreamToString(CNcbiOstrstream& out)
        : m_Out(out)
	{
	}
    operator string(void) const
	{
		SIZE_TYPE length = (size_t)m_Out.pcount();
		if ( length == 0 )
			return string();
		const char* str = m_Out.str();
		m_Out.freeze(false);
		return string(str, length);
	}
private:
    CNcbiOstrstream& m_Out;
};

unsigned int
StringToUInt(const string& str, int flags = 0, int base = 10);

using std::streamsize;

#include <cstdlib>

//----------------------------------------------------------
class CSeqVector
{
public:
    typedef char value_type;
    typedef TSeqPos size_type;
    typedef Int4 difference_type;
    typedef std::random_access_iterator_tag iterator_category;
    typedef const value_type* pointer;
    typedef const value_type& reference;
    typedef const value_type* const_iterator;
    
private:
    pointer sequence;
    size_type sequence_length;
    
public:
    CSeqVector(pointer arg_sequence, size_type arg_sequence_length)
    {
        sequence = arg_sequence;
        sequence_length = arg_sequence_length;
    }
    
    ~CSeqVector() {}
    
    value_type operator[](TSeqPos pos) const
    {
        return sequence[pos];
    }
    
    bool empty() const
    {
        return sequence_length == 0;
    }
    
    TSeqPos size() const
    {
        return sequence_length;
    }
    
    const_iterator begin() const
    {
        return sequence;
    }
    
    const_iterator end() const
    {
        return sequence + sequence_length;
    }    
};

class CObject
{	
};

#endif
