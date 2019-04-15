/*  $Id$
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
 * Authors:  Aleksandr Morgulis
 *
 */

/** \file libmatchhits/fast_seeds.hpp
    \brief Internal seeder.
*/

#ifndef LIBREADFINDER_FAST_SEEDS_HPP
#define LIBREADFINDER_FAST_SEEDS_HPP

#include <boost/dynamic_bitset.hpp>

#include <config.h>

#include <libreadfinder/batch.hpp>
#include <libreadfinder/fsdefs.hpp>
#include <libreadfinder/fsidx.hpp>
// #include <libreadfinder/seeddata.hpp>
#include <libreadfinder/defs.hpp>

READFINDER_NS_BEGIN

//==============================================================================
class CFastSeeds : public CFastSeedsDefs
{
private:

    static size_t const WMAP_SIZE = 1ULL + (1ULL<<(ANCHOR_BITS + LB));
    static size_t const ANCHOR_TBL_SIZE = (1ULL<<ANCHOR_BITS);

    //--------------------------------------------------------------------------
    struct WordCountingJobData;
    struct WordCountingJob;

    //--------------------------------------------------------------------------
    typedef boost::dynamic_bitset< uint64_t > BitSet;

    static TSeqLen const WORD_STRIDE = 3;

    //--------------------------------------------------------------------------
    struct WordTableGeneratorJob;
    struct WordTableSortJob;

    struct HashWord
    {
        WordData wd;
        OrdId readid;
        TReadOff hashoff;
        EStrand strand;
        EMate mate;
    };

    friend std::ostream & operator<<( std::ostream & os, HashWord const & hw );

    static_assert( sizeof( HashWord ) == 12, "" );

    struct ExtHashWord
    {
        HashWord hw;
        bool exact;
    };

    // typedef std::vector< HashWord, TAlloc > WordTable;
    typedef std::vector< HashWord > WordTable;
    typedef std::vector< uint32_t > WordMap;

    class HashWordSource;
    class HashMaskSource;

    //--------------------------------------------------------------------------
    struct AnchorListEntry
    {
        union
        {
            struct
            {
                uint32_t wanchor : ANCHOR_BITS + LB;
                uint32_t ins : 1;
                uint32_t del : 1;
                bool exact : 1;
            } f;

            uint32_t w;
        };
    };

    static_assert( sizeof( AnchorListEntry ) == 4, "" );

    static size_t const ANCHOR_LIST_MAX_LEN = 292;

    typedef AnchorListEntry AnchorList[ANCHOR_LIST_MAX_LEN];

    struct AnchorTableEntry
    {
        size_t len = 0;
        uint8_t n_subst = 0,
                n_del   = 0,
                n_ins   = 0;
        AnchorList data;
    };

    typedef std::vector< AnchorTableEntry > AnchorTable;

    struct SeedSearchJobData;
    struct SeedSearchJob;

    //--------------------------------------------------------------------------
    struct Hit
    {
        uint32_t refpos;
        OrdId read;
        TReadOff readpos;
        uint8_t strand : 1;
        uint8_t mate : 1;
        bool exact : 1;
        bool mid : 1;

        TSeqOff GetDiag() const { return refpos - readpos; }

        friend std::ostream & operator<<( std::ostream & os, Hit const & h )
        {
            os << h.refpos << ':' << h.read << ':' << h.readpos << ':'
               << (int)h.strand << ':' << (int)h.mate;
            return os;
        }

        friend bool operator<( Hit const & x, Hit const & y )
        {
            // ASSUMPTION: x and y are from the same reference
            //
            auto tx( (x.strand + x.mate)%2 ),
                 ty( (y.strand + y.mate)%2 );
            int64_t dx( x.refpos - x.readpos ),
                    dy( y.refpos - y.readpos );
            return x.read == y.read ?
                   tx == ty ?
                   x.strand == y.strand ?
                   dx < dy :
                   x.strand < y.strand :
                   tx < ty :
                   x.read < y.read;
        }
    };

    static_assert( sizeof( Hit ) == 12, "" );

    typedef std::vector< Hit > Hits;
    typedef HashWord TaskEntry;
    // typedef std::vector< CSeedData > SeedData;

    struct HitFilteringJob;
    struct ReadMarkingJob;

public:

    CFastSeeds( CBatch & bctx, bool seeder_mode = false );
    void Run();
    // SeedData * GetSeeds() { return seeds_.release(); }

private:

    void CreateAnchorTable();
    void UpdateAnchorLists( uint32_t wanchor );
    void UpdateAnchorUseMap( uint32_t wa );
    void CreateWordTable();
    void ComputeSeeds();

    CRefData const & GetRefs() const { return *bctx_.GetSearchCtx().refs; }
    CReadData const & GetReads() const { return bctx_.GetReads(); }

    CBatch & bctx_;
    CFastSeedsIndex fsidx_;
    // std::unique_ptr< SeedData > seeds_;
    WordTable wt_;
    WordMap wmap_ = WordMap( WMAP_SIZE, 0 );
    AnchorTable atbl_ = AnchorTable( ANCHOR_TBL_SIZE );
    BitSet anchor_use_map_;
    bool seeder_mode_ = false;
};

READFINDER_NS_END

#endif

