/*  
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
 *   Header file for CSymDustMasker class.
 *
 */

#ifndef SYMDUST_H
#define	SYMDUST_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <memory>
#include <deque>
#include <list>

#include "def.h"

class CSymDustMasker
{
private:
    struct CIupac2Ncbi2na_converter
    {
        Uint1 operator()(Uint1 r) const
        {
            switch (r)
            {
                case 67: return 1;
                case 71: return 2;
                case 84: return 3;
                default: return 0;
            }
        }
    };

    //typedef std::string seq_t;
    typedef char* seq_t;
    typedef CIupac2Ncbi2na_converter convert_t;

public:
    typedef seq_t sequence_type;
    //typedef sequence_type::size_type size_type;
    typedef Uint4 size_type;
    typedef std::pair<size_type, size_type> TMaskedInterval;
    typedef std::vector<TMaskedInterval> TMaskList;

    static const Uint4 DEFAULT_LEVEL = 20;
    static const Uint4 DEFAULT_WINDOW = 64;
    static const Uint4 DEFAULT_LINKER = 1;

    struct perfect
    {
        TMaskedInterval bounds_;
        Uint4 score_;
        size_type len_;

        perfect(size_type start, size_type stop, Uint4 score, size_type len)
            : bounds_(start, stop), score_(score), len_(len) {}
    };

    typedef std::list<perfect> perfect_list_type;
    typedef std::vector<Uint4> thres_table_type;
    typedef Uint1 triplet_type;

    static const triplet_type TRIPLET_MASK = 0x3f;

    CSymDustMasker(Uint4 level = DEFAULT_LEVEL, size_type window = DEFAULT_WINDOW,
                   size_type linker = DEFAULT_LINKER);

    std::auto_ptr<TMaskList> operator()(const sequence_type& seq,
                                        size_type seq_len);

    std::auto_ptr<TMaskList> operator()(const sequence_type& seq, 
                                        size_type seq_len,
                                        size_type start, size_type stop);

    void GetMaskedLocs(const sequence_type&seq, size_type seq_len,
                       std::vector<std::pair<size_type, size_type> >& locs);
    void BLAST_ComplementMaskLocations(Uint4 len, 
                                       std::vector<std::pair<size_type, size_type> >& masked_locs,
                                       std::vector<std::pair<size_type, size_type> >& unmasked_locs);

private:
    //typedef sequence_type::const_iterator seq_citer_type;
    typedef const char* seq_citer_type;

    class triplets
    {
    public:
    	triplets(size_type window, Uint1 low_k,
    			 perfect_list_type& perfect_list,
    			 thres_table_type& thresholds);

    	size_type start() const { return start_; }
    	size_type stop() const { return stop_; }
    	size_type size() const { return triplet_list_.size(); }

    	void find_perfect();
    	bool shift_window(triplet_type t);
    	bool shift_high(triplet_type t);
    	bool needs_processing() const
    	{
    		Uint4 count = stop_ - L;
    		return count < triplet_list_.size() &&
    				10 * r_w > thresholds_[count];
    	}
    private:
    	typedef std::deque< triplet_type > impl_type;
    	typedef impl_type::const_iterator impl_citer_type;
    	typedef std::vector<Uint1> counts_type;

    	void add_triplet_info(Uint4& r, counts_type& c, triplet_type t)
    	{ r += c[t]; ++c[t]; }

    	void rem_triplet_info(Uint4& r, counts_type& c, triplet_type t)
    	{ --c[t]; r -= c[t]; }

    	impl_type triplet_list_;
    	size_type start_;
    	size_type stop_;
    	size_type max_size_;

    	Uint1 low_k_;
    	Uint4 L;

    	perfect_list_type& P;
    	thres_table_type& thresholds_;

    	counts_type c_w;
    	counts_type c_v;
    	Uint4 r_w;
    	Uint4 r_v;
    	Uint4 num_diff;
    };

    void save_masked_regions(TMaskList& res, size_type w, size_type start);

    Uint4 level_;
    size_type window_;
    size_type linker_;
    Uint1 low_k_;
    perfect_list_type P;
    thres_table_type thresholds_;
    convert_t converter_;
};

#endif	/* SYMDUST_H */

