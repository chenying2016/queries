/*  $Id: seq_masker_istat_ascii.cpp 122478 2008-03-19 19:14:23Z morgulis $
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
 * Author:  Aleksandr Morgulis
 *
 * File Description:
 *   Implementation for CSeqMaskerIstatAscii class.
 *
 */

#include <ncbi_pch.hpp>

#include "seq_masker_istat_ascii.hpp"

#include <cstring>

BEGIN_NCBI_SCOPE

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//CSeqMaskerIstatAscii::CSeqMaskerIstatAscii( const string & name,
//                                            Uint4 arg_threshold,
//                                            Uint4 arg_textend,
//                                            Uint4 arg_max_count,
//                                            Uint4 arg_use_max_count,
//                                            Uint4 arg_min_count,
//                                            Uint4 arg_use_min_count )
//    :   CSeqMaskerIstat(    arg_threshold, arg_textend, 
//                            arg_max_count, arg_use_max_count,
//                            arg_min_count, arg_use_min_count )
//{
//    CNcbiIfstream input_stream( name.c_str() );
//
//    //if( !input_stream )
//    //    NCBI_THROW( Exception, eStreamOpenFail,
//    //                string( "could not open " ) + name );
//	ASSERT(input_stream);
//
//    bool start = true;
//    Uint4 linenum = 0UL;
//    Uint4 ambig_len = kMax_UI4;
//    string line;
//
//    while( input_stream )
//    {
//        line.erase();
//        getline( input_stream, line );
//        ++linenum;
//
//        if( !line.length() || line[0] == '#' ) continue;
//
//        // Check if we have a precomputed parameter.
//        if( line[0] == '>' )
//        {
//            SIZE_TYPE name_end = line.find_first_of( " \t", 0 );
//            SIZE_TYPE val_start = line.find_first_not_of( " \t", name_end );
//
//            if( name_end == NPOS || val_start == NPOS )
//            {
//                CNcbiOstrstream str;
//                str << "at line " << linenum;
//                string msg = CNcbiOstrstreamToString(str);
//                //NCBI_THROW( Exception, eSyntax, msg);
//            }
//
//            string name = line.substr( 1, name_end - 1 );
//
//            if( name == "t_threshold" && get_threshold() == 0 )
//                set_threshold( 
//                    StringToUInt(line.substr(val_start, NPOS), 0, 0));
//
//            if( name == "t_extend" && get_textend() == 0 )
//                set_textend(
//                    StringToUInt(line.substr(val_start, NPOS), 0, 0));
//
//            if( name == "t_low" )
//                set_min_count(
//                    StringToUInt(line.substr(val_start, NPOS), 0, 0));
//
//            if( name == "t_high" && get_max_count() == 0 )
//                set_max_count(
//                    StringToUInt(line.substr(val_start, NPOS), 0, 0));
//
//            continue;
//        }
//
//        if( start )
//        {
//            start = false;
//            uset.set_unit_size( 
//                static_cast< Uint1 >( StringToUInt( line ) ) );
//            continue;
//        }
//
//        SIZE_TYPE unit_start = line.find_first_not_of( " \t", 0 );
//        SIZE_TYPE unit_end   = line.find_first_of( " \t", unit_start );
//        SIZE_TYPE cnt_start  = line.find_first_not_of( " \t", unit_end );
//
//        if( unit_start == NPOS || unit_end == NPOS || cnt_start == NPOS )
//        {
//            CNcbiOstrstream str;
//            str << "at line " << linenum;
//            string msg = CNcbiOstrstreamToString( str );
//            //NCBI_THROW( Exception, eSyntax, msg );
//        }
//
//        Uint4 unit = StringToUInt(line.substr(unit_start, 
//                                                    unit_end - unit_start),
//                                        0, 16);
//        Uint4 cnt = StringToUInt(line.substr(cnt_start));
//
//        if( cnt < ambig_len ) {
//            ambig_len = cnt;
//            set_ambig_unit( unit );
//        }
//
//        if( cnt >= get_min_count() ) 
//            uset.add_info( unit, cnt );
//    }
//
//    string bad_param;
//
//    if( get_threshold() == 0 )
//        bad_param += "t_threhold ";
//
//    if( get_textend() == 0 )
//        bad_param += "t_extend ";
//
//    if( get_max_count() == 0 )
//        bad_param += "t_high ";
//
//    if( get_min_count() == 0 )
//        bad_param += "t_low ";
//
//    //if( !bad_param.empty() )
//    //    NCBI_THROW( Exception, eParam, bad_param );
//	ASSERT(bad_param.empty());
//
//    if( get_use_min_count() == 0 )
//      set_use_min_count( (get_min_count() + 1)/2 );
//
//    if( get_use_max_count() == 0 )
//      set_use_max_count( get_max_count() );
//}

struct FindFirstSpaceOrTab
{
    std::size_t operator()(const char* str, std::size_t start_pos)
    {
        register const char* start = str + start_pos;
        while ( (*start) != '\0')
        {
            if ((*start) == ' ' || (*start) == '\t')
                return start - str;
            ++start;
        }
        return std::string::npos;
    }
};

struct FindFirstNotSpaceOrTab
{
    std::size_t operator()(const char* str, std::size_t start_pos)
    {
        register const char* start = str + start_pos;
        while ((*start) != '\0')
        {
            if((*start) != ' ' && (*start) != '\t')
                return start - str;
            ++start;
        }
        return std::string::npos;
    }
};

struct Substring
{
    void operator()(char* dst, const char* src, std::size_t start_pos, std::size_t num)
    {
        memcpy(dst, src + start_pos, num);
        dst[num] = '\0';
    }
};

struct StringEqual
{
    bool operator()(const char* a, const char* b, std::size_t n)
    {
        const char* aend = a + n;
        while ((a != aend) && (*a++ == *b++));
        if (a == aend) return true;
        return false;
    }
};

struct stringtouint
{
    Uint8 operator()(const char* str, char** endptr = NULL, int base = 10)
    {
        Uint8 ret = strtoull(str, endptr, base);
        return ret;
    }
};

struct MKFileBuffer
{
    char* _buf;
    std::size_t _buf_sz;
    std::size_t _cur;

    operator bool()
    {
        return _cur < _buf_sz;
    }

    const char* getline()
    {
        if (_cur >= _buf_sz) return NULL;

        const char* ret = _buf + _cur;

        while (_buf[_cur] != '\n') ++_cur;
        while (isspace(_buf[_cur]))++_cur;

        return ret;
    }

    MKFileBuffer(const char* name)
    {
        FILE* file = fopen(name, "r");
        if (file == NULL)
        {
            fprintf(stderr, "cannot open file %s for reading.\n", name);
            abort();
        }
        std::size_t flen;
        fseek(file, 0ul, SEEK_END);
        flen = ftell(file);
        fseek(file, 0ul, 0ul);
        std::size_t r;
        _buf = new char[flen + 1];
        r = fread(_buf, 1, flen, file);
        assert(r == flen);

        _buf[flen] = '\n';

        fclose(file);

        _buf_sz = flen;
        _cur = 0;
    }

    ~MKFileBuffer()
    {
        if (_buf)
        {
            delete[] _buf;
        }
    }
};

CSeqMaskerIstatAscii::CSeqMaskerIstatAscii(const string& name,
                                           Uint4 arg_threshold,
                                           Uint4 arg_textend,
                                           Uint4 arg_max_count,
                                           Uint4 arg_use_max_count,
                                           Uint4 arg_min_count,
                                           Uint4 arg_use_min_count)
    : CSeqMaskerIstat(arg_threshold, arg_textend,
                      arg_max_count, arg_use_max_count,
                      arg_min_count, arg_use_min_count)
{
   MKFileBuffer fb(name.c_str());
   const char* line;

   bool start = true;
   Uint4 linenum = 0UL;
   Uint4 ambig_len = UINT4_MAX;
   
   FindFirstNotSpaceOrTab ffnst;
   FindFirstSpaceOrTab ffst;
   stringtouint stui;
   StringEqual se;

   while (fb)
   {
       line = fb.getline();
      ++linenum;

      if (line[0] == '#') continue;

      /* Check if we have a precomputed parameter */
      if (line[0] == '>')
      {
          SIZE_TYPE name_end = ffst(line, 0);
          SIZE_TYPE val_start = ffnst(line, name_end);

          if (name_end == NPOS || val_start == NPOS)
          {
              fprintf(stderr, "[%s, %ul] syntax error.\n", __func__, linenum);
              abort();
          }

          if ( se(line, ">t_threshold", strlen(">t_threshold")) && get_threshold() == 0)
              set_threshold(stui(line + val_start, 0, 0));

          if ( se(line, ">t_extend", strlen(">t_extend")) && get_textend() == 0 )
              set_textend(stui(line + val_start, 0, 0));

          if ( se(line, ">t_low", strlen(">t_low")) )
              set_min_count(stui(line + val_start, 0, 0));

          if ( se(line, ">t_high", strlen(">t_high")) && get_max_count() == 0 )
              set_max_count(stui(line + val_start, 0, 0));

          continue;
      }

      if (start)
      {
          start = false;
          uset.set_unit_size(static_cast<Uint1>(stui(line, 0, 0)));
          continue;
      }

      char* pEnd;
      Uint4 unit = strtoull(line, &pEnd, 16);
      Uint4 cnt = strtoull(pEnd, 0, 10);

      if (cnt < ambig_len)
      {
          ambig_len = cnt;
          set_ambig_unit(unit);
      }

      if (cnt >= get_min_count())
      {
          uset.add_info(unit, cnt);
      }
   }

   std::size_t line_len = 0;
   char param[1024];
   if (get_threshold() == 0)
   {
      strcpy(param, "t_threshold") ;
      line_len += strlen("t_threshold");
   }

   if (get_textend() == 0)
   {
       strcpy(param + line_len, "t_extend");
       line_len += strlen("t_extend");
   }

   if (get_max_count() == 0)
   {
       strcpy(param + line_len, "t_high");
       line_len += strlen("t_high");
   }

   if (get_min_count() == 0)
   {
       strcpy(param + line_len, "t_low");
       line_len += strlen("t_low");
   }

   if (line_len > 0)
   {
       fprintf(stderr, "[%s] bad parameter: %s\n", __func__, param);
       abort();
   }

   if (get_use_min_count() == 0)
       set_use_min_count((get_min_count() + 1) / 2);

   if (get_use_max_count() == 0)
       set_use_max_count(get_max_count());
}

//------------------------------------------------------------------------------
Uint4 CSeqMaskerIstatAscii::trueat( Uint4 unit ) const
{ return uset.get_info( unit ); }

//------------------------------------------------------------------------------
Uint4 CSeqMaskerIstatAscii::at( Uint4 unit ) const
{
  Uint4 res = uset.get_info( unit );

  if( res == 0 || res < get_min_count() )
    return get_use_min_count();

  return (res > get_max_count()) ? get_use_max_count() : res;
}

END_NCBI_SCOPE
