/*  $Id: seq_masker_uset_simple.hpp 103491 2007-05-04 17:18:18Z kazimird $
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
 *   Definition for CSeqMaskerUsetSimple class.
 *
 */

#ifndef C_SEQ_MASKER_USET_SIMPLE_H
#define C_SEQ_MASKER_USET_SIMPLE_H

#include <vector>

#include <win_masker_macros.h>

//#include <ncbitype.h>
//#include <ncbistr.h>
//#include <ncbiobj.hpp>

BEGIN_NCBI_SCOPE

/**
 **\brief This class implements simple, unoptimized unit counts container.
 **/
class NCBI_XALGOWINMASK_EXPORT CSeqMaskerUsetSimple
{
    public:

        /**
         **\brief Object constructor.
         **\param arg_unit_size the unit size
         **/
        CSeqMaskerUsetSimple( Uint1 arg_unit_size = 15 ) 
            : unit_size( arg_unit_size )
        {}

        /**
         **\brief Get the unit size.
         **\return the unit size
         **/
        Uint1 get_unit_size() const { return unit_size; }

        /**
         **\brief Set the unit size.
         **\param arg_unit_size the new value of unit size
         **/
        void set_unit_size( Uint1 arg_unit_size )
        { unit_size = arg_unit_size; }

        /**
         **\brief Add information about the given unit.
         **
         ** When this method is called multiple times, units should be
         ** in ascending order.
         **
         **\param unit the target unit
         **\param count number of times unit and its reverse complement
         **             appear in the genome
         **/
        void add_info( Uint4 unit, Uint4 count );

        /**
         **\brief Lookup the count value for a given unit.
         **\param unit the target unit
         **\return the count of the given unit, or 0 if unit was not
         **        found
         **/
        Uint4 get_info( Uint4 unit ) const;	

    public:
        void DumpCSeqMaskerUsetSimple(FILE* file)
        {
            vector<Uint4>::size_type n;
            
            Uint4 t = unit_size;
            n = fwrite(&t, sizeof(Uint4), 1, file);
            ASSERT(n == 1);
            
            Uint4 size = units.size();
            
            n = fwrite(&size, sizeof(Uint4), 1, file);
            ASSERT(n == 1);
            
            Uint4* arr = &units[0];
            n = fwrite(arr, sizeof(Uint4), size, file);
            ASSERT(n == size);
            
            arr = &counts[0];
            n = fwrite(arr, sizeof(Uint4), size, file);
            ASSERT(n == size);
        }
        
        void LoadCSeqMaskerUsetSimple(FILE* file)
        {
            vector<Uint4>::size_type n;
            Uint4 size;
            
            using std::cout;
            using std::endl;
            
            Uint4 t;
            n = fread(&t, sizeof(Uint4), 1, file);
            ASSERT(n == 1);
            unit_size = t;
            
            n = fread(&size, sizeof(Uint4), 1, file);
            ASSERT(n == 1);
            
            units.resize(size);
            
            Uint4* arr = units.data();
            n = fread(arr, sizeof(Uint4), size, file);
            ASSERT(n == size);
            
            counts.resize(size);
            arr = counts.data();
            n = fread(arr, sizeof(Uint4), size, file);
            ASSERT(n == size);
        }
        
        void Print()
        {
            using std::endl;
            using std::cout;
            cout << endl;
            Uint4 size = units.size();
            Uint4* arr;
            
            arr = units.data();
            cout << "unit_size = " << (Uint4)unit_size << endl;
            cout << "units size = " << units.size() << endl;
            cout << "units[0] = " << arr[0] << endl;
            cout << "units last = " << arr[size-1] << endl;
            cout << "counts size = " << counts.size() << endl;
            
            arr = counts.data();
            cout << "counts[0] = " << arr[0] << endl;
            cout << "counts last = " << arr[size-1] << endl;
            cout << endl;
        }

    private:

        Uint1 unit_size;        /**<\internal The unit size. */

        vector< Uint4 > units;  /**<\internal The vector of units. */
        vector< Uint4 > counts; /**<\internal The vector of counts corresponding to units. */
};

END_NCBI_SCOPE

#endif
