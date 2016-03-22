/*  $Id: version.cpp 343034 2011-11-02 18:37:06Z ivanov $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Denis Vakatov, Vladimir Ivanov
 *
 * File Description:
 *   CVersionInfo -- a version info storage class
 *
 */

#include "ncbi_pch.hpp"
#include "version.hpp"
#include "ncbimisc.hpp"

#include <cstdlib>
#include <cstring>

BEGIN_NCBI_SCOPE

#define NCBI_PACKAGE                       1
#define NCBI_PACKAGE_NAME                  "blast"
#define NCBI_PACKAGE_VERSION_MAJOR         2
#define NCBI_PACKAGE_VERSION_MINOR         3
#define NCBI_PACKAGE_VERSION_PATCH         0
#define NCBI_PACKAGE_CONFIG                ""


CVersionInfo::CVersionInfo(void) 
    : m_Major(-1),
      m_Minor(-1),
      m_PatchLevel(-1),
      m_Name(kEmptyStr)

{
}

CVersionInfo::CVersionInfo(int ver_major,
                           int  ver_minor,
                           int  patch_level, 
                           const string& name) 
    : m_Major(ver_major),
      m_Minor(ver_minor),
      m_PatchLevel(patch_level),
      m_Name(name)

{
}


CVersionInfo::CVersionInfo(const string& version,
                           const string& name)
{
    FromStr(version);
    if (!name.empty()) {
        m_Name = name;
    }
}

CVersionInfo::CVersionInfo(EVersionFlags flags)
{
    ASSERT( flags == kAny || flags == kLatest);
    m_Major = m_Minor = m_PatchLevel = (flags == kAny) ? 0 : -1;
}

static
void s_ConvertVersionInfo(CVersionInfo* vi, const char* str)
{
    int major, minor, patch = 0;
	const char* err = "Invalid version format";
    if (!isdigit((unsigned char)(*str))) {
        //NCBI_THROW2(CStringException, eFormat, "Invalid version format", 0);
		error_and_exit(err);
    }
    major = atoi(str);
    if (major < 0) {
        //NCBI_THROW2(CStringException, eFormat, "Invalid version format", 0);
		error_and_exit(err);
    }
    for (; *str && isdigit((unsigned char)(*str)); ++str) {}
    if (*str != '.') {
        //NCBI_THROW2(CStringException, eFormat, "Invalid version format", 0);
		error_and_exit(err);
    }
    ++str;
    if (!isdigit((unsigned char)(*str))) {
        //NCBI_THROW2(CStringException, eFormat, "Invalid version format", 0);
		error_and_exit(err);
    }

    minor = atoi(str);
    if (minor < 0) {
        //NCBI_THROW2(CStringException, eFormat, "Invalid version format", 0);
		error_and_exit(err);
    }
    for (; *str && isdigit((unsigned char)(*str)); ++str) {}

    if (*str != 0) {
        if (*str != '.') {
            //NCBI_THROW2(CStringException, eFormat, "Invalid version format", 0);
			error_and_exit(err);
        }
        ++str;
        patch = atoi(str);
        if (patch < 0) {
            //NCBI_THROW2(CStringException, eFormat, "Invalid version format", 0);
			error_and_exit(err);
        }
    }

    vi->SetVersion(major, minor, patch);

}


void CVersionInfo::FromStr(const string& version)
{
    s_ConvertVersionInfo(this, version.c_str());
/*
    vector<string> lst;
    NStr::Tokenize(version, ".", lst, NStr::eNoMergeDelims);

    if (lst.size() == 0) {
        NCBI_THROW2(CStringException, eFormat,
                            "Invalid version format", 0);
    }

    for (unsigned i = 0; i < 3; ++i) {
        string tmp;
        if (i < lst.size()) {
            tmp = lst[i];
        }
        int value = tmp.empty() ? 0 : NStr::StringToInt(tmp);
        switch (i) {
        case 0: 
            if (value == 0) {
                NCBI_THROW2(CStringException, eFormat,
                            "Invalid version format (major is 0)", 0);
            }
            m_Major = value;
            break;
        case 1:
            m_Minor = value;
            break;
        case 2:
            m_PatchLevel = value;
            break;
        } 
    } // for
*/
}


CVersionInfo::CVersionInfo(const CVersionInfo& version)
    : m_Major(version.m_Major),
      m_Minor(version.m_Minor),
      m_PatchLevel(version.m_PatchLevel),
      m_Name(version.m_Name)
{
}

CVersionInfo& CVersionInfo::operator=(const CVersionInfo& version)
{
    m_Major = version.m_Major;
    m_Minor = version.m_Minor;
    m_PatchLevel = version.m_PatchLevel;
    return *this;
}


string CVersionInfo::Print(void) const
{
    if (m_Major < 0) {
        return kEmptyStr;
    }
    CNcbiOstrstream os;
    os << m_Major << "." << (m_Minor >= 0 ? m_Minor : 0);
    if (m_PatchLevel >= 0) {
        os << "." << m_PatchLevel;
    }
    if ( !m_Name.empty() ) {
        os << " (" << m_Name << ")";
    }
    return CNcbiOstrstreamToString(os);
}


CVersionInfo::EMatch 
CVersionInfo::Match(const CVersionInfo& version_info) const
{
    if (GetMajor() != version_info.GetMajor())
        return eNonCompatible;

    if (GetMinor() < version_info.GetMinor())
        return eNonCompatible;

    if (GetMinor() > version_info.GetMinor())
        return eBackwardCompatible;

    // Minor versions are equal.
    
    if (GetPatchLevel() == version_info.GetPatchLevel()) {
        return eFullyCompatible;
    }

    if (GetPatchLevel() > version_info.GetPatchLevel()) {
        return eBackwardCompatible;
    }

    return eConditionallyCompatible;

}

void CVersionInfo::SetVersion(int  ver_major,
                              int  ver_minor,
                              int  patch_level)
{
    m_Major      = ver_major;
    m_Minor      = ver_minor;
    m_PatchLevel = patch_level;
}



bool IsBetterVersion(const CVersionInfo& info, 
                     const CVersionInfo& cinfo,
                     int&  best_major, 
                     int&  best_minor,
                     int&  best_patch_level)
{
    int major = cinfo.GetMajor();
    int minor = cinfo.GetMinor();
    int patch_level = cinfo.GetPatchLevel();

    if (info.GetMajor() == -1) {  // best major search
        if (major > best_major) { 
            best_major = major;
            best_minor = minor;
            best_patch_level = patch_level;
            return true;
        }
    } else { // searching for the specific major version
        // Do not chose between major versions.
        // If they are not equal then they are not compatible.
        if (info.GetMajor() != major) {
            return false;
        }
    }

    if (info.GetMinor() == -1) {  // best minor search
        if (minor > best_minor) {
            best_major = major;
            best_minor = minor;
            best_patch_level = patch_level;
            return true;
        }
    } else { 
        if (info.GetMinor() > minor) {
            return false;
        }
        if (info.GetMinor() < minor) {
            best_major = major;
            best_minor = minor;
            best_patch_level = patch_level;
            return true;
        }
    }

    // Major and minor versions are equal.
    // always looking for the best patch
    if (patch_level > best_patch_level) {
            best_major = major;
            best_minor = minor;
            best_patch_level = patch_level;
            return true;
    }
    return false;    
}

/// Which end to truncate a string.
enum ETrunc {
	eTrunc_Begin,  ///< Truncate leading spaces only
	eTrunc_End,    ///< Truncate trailing spaces only
	eTrunc_Both    ///< Truncate spaces at both begin and end of string
};

void TruncateSpacesInPlace(string& str, ETrunc where = eTrunc_Both)
{
    SIZE_TYPE length = str.length();
    if (length == 0) {
        return;
    }
    SIZE_TYPE beg = 0;
    if ( where == eTrunc_Begin  ||  where == eTrunc_Both ) {
        // It's better to use str.data()[] to check string characters
        // to avoid implicit modification of the string by non-const operator[]
        ASSERT(beg < length);
        while ( isspace((unsigned char) str.data()[beg]) ) {
            if (++beg == length) {
                str.erase();
                return;
            }
        }
    }

    SIZE_TYPE end = length;
    if ( where == eTrunc_End  ||  where == eTrunc_Both ) {
        // It's better to use str.data()[] to check string characters
        // to avoid implicit modification of the string by non-const operator[]
        ASSERT(beg < end);
        while (isspace((unsigned char) str.data()[--end])) {
            if (beg == end) {
                str.erase();
                return;
            }
        }
        ASSERT(beg <= end  &&  !isspace((unsigned char) str.data()[end]));
        ++end;
    }
    ASSERT(beg < end  &&  end <= length);

#if defined(NCBI_COMPILER_GCC)  &&  (NCBI_COMPILER_VERSION == 304)
    // work around a library bug
    str.replace(end, length, kEmptyStr);
    str.replace(0, beg, kEmptyStr);
#else
    if ( beg | (end - length) ) { // if either beg != 0 or end != length
        str.replace(0, length, str, beg, end - beg);
    }
#endif
}

void ParseVersionString(const string&  vstr, 
                        string*        program_name, 
                        CVersionInfo*  ver)
{
    ASSERT(program_name);
    ASSERT(ver);

    if (vstr.empty()) {
        //NCBI_THROW2(CStringException, eFormat, "Version string is empty", 0);
		const char* err = "Version string is empty";
		error_and_exit(err);
    }

    program_name->erase();


    string lo_vstr(vstr); // NStr::ToLower(lo_vstr);
	for (string::iterator sit = lo_vstr.begin(); sit != lo_vstr.end(); ++sit)
	{
		int c = *sit;
		if (c >= 'A' && c <= 'Z') c = c - 'A' + 'a';
		*sit = c;
	}
    string::size_type pos;

    const char* vstr_str = vstr.c_str();

    // 2.3.4 (program)

    pos = lo_vstr.find("(");
    if (pos != string::npos) {
        string::size_type pos2 = lo_vstr.find(")", pos);
        if (pos2 == string::npos) { // not found
            //NCBI_THROW2(CStringException, 
            //            eFormat, "Version string format error", 0);
			const char* err = "Version string format error";
			error_and_exit(err);
        }
        for (++pos; pos < pos2; ++pos) {
            program_name->push_back(vstr.at(pos));
        }
        TruncateSpacesInPlace(*program_name);

        s_ConvertVersionInfo(ver, vstr.c_str());
        return;

    }


    // all other normal formats

    const char* version_pattern = "version";

    pos = lo_vstr.find(version_pattern);
    if (pos == string::npos) {
        version_pattern = "v.";
        pos = lo_vstr.find(version_pattern);

        if (pos == string::npos) {
            version_pattern = "ver";
            pos = lo_vstr.find(version_pattern);

            if (pos == string::npos) {
                version_pattern = "";
                // find the first space-digit and assume it's version
                const char* ch = vstr_str;
                for (; *ch; ++ch) {
                    if (isdigit((unsigned char)(*ch))) {
                        if (ch == vstr_str) {
                            // check if it's version
                            const char* ch2 = ch + 1;
                            for (;*ch2; ++ch2) {
                                if (!isdigit((unsigned char)(*ch2))) {
                                    break;
                                }
                            } // for
                            if (*ch2 == '.') {
                                pos = ch - vstr_str;
                                break;
                            } else {
                                continue;
                            }
                        } else {
                            if (isspace((unsigned char) ch[-1])) {
                                pos = ch - vstr_str;
                                break;
                            }
                        }
                    } // if digit
                    
                } // for
            }
        }
    }


    if (pos != string::npos) {
        int pname_end = (int)(pos - 1);
        for (; pname_end >= 0; --pname_end) {
            char ch = vstr[pname_end];
            if (!isspace((unsigned char) ch)) 
                break;
        } // for
        if (pname_end <= 0) {
        } else {
            program_name->append(vstr.c_str(), pname_end + 1);
        }

        pos += strlen(version_pattern);
        for(; pos < vstr.length(); ++pos) {
            char ch = vstr[pos];
            if (ch == '.') 
                continue;
            if (!isspace((unsigned char) ch)) 
                break;            
        } // for

        const char* ver_str = vstr_str + pos;
        s_ConvertVersionInfo(ver, ver_str);
        return;
    } else {
        *ver = CVersionInfo::kAny;
        *program_name = vstr;
        TruncateSpacesInPlace(*program_name);
        if (program_name->empty()) {
            //NCBI_THROW2(CStringException, eFormat, "Version string is empty", 0);
			const char* err = "Version string is empty";
			error_and_exit(err);
        }
    }


}

/////////////////////////////////////////////////////////////////////////////
//  CComponentVersionInfo

CComponentVersionInfo::CComponentVersionInfo( const string& component_name,
                                              int  ver_major,
                                              int  ver_minor,
                                              int  patch_level,
                                              const string& name)
    : CVersionInfo(ver_major, ver_minor, patch_level, name),
      m_ComponentName( component_name )
{
}

CComponentVersionInfo::CComponentVersionInfo( const string& component_name,
                                              const string& version,
                                              const string& name)
    : CVersionInfo( version, name),
      m_ComponentName( component_name )
      
{
}

CComponentVersionInfo::CComponentVersionInfo(
    const CComponentVersionInfo& version)
    : CVersionInfo( version),
      m_ComponentName( version.m_ComponentName )
      
{
}

CComponentVersionInfo& CComponentVersionInfo::operator=(
    const CComponentVersionInfo& version)
{
    m_ComponentName = version.m_ComponentName;
    CVersionInfo::operator=( version );
    return *this;
}

string CComponentVersionInfo::Print(void) const
{
    CNcbiOstrstream os;
    os << GetComponentName() << ": " << CVersionInfo::Print();
    return CNcbiOstrstreamToString(os);
}


END_NCBI_SCOPE
