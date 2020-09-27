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

#ifndef LIBREADFINDER_FSDEFS_HPP
#define LIBREADFINDER_FSDEFS_HPP

#include <libreadfinder/refdata.hpp>
#include <libreadfinder/defs.hpp>

READFINDER_NS_BEGIN

//==============================================================================
struct CFastSeedsDefs
{
    typedef CRefData::TRefOId TRefOId;
    typedef std::vector< uint32_t > WordMap;

    static ssize_t const ANCHOR_BASES = 9;
    static ssize_t const HALF_WORD_BASES = 6;
    static ssize_t const WORD_BASES = HALF_WORD_BASES*2;
    static ssize_t const NMER_BASES = ANCHOR_BASES + WORD_BASES;
    static ssize_t const MIN_INEXACT_BASES = 10;

protected:

    static ssize_t const ANCHOR_BITS  = LB*ANCHOR_BASES;
    static ssize_t const WORD_BITS = LB*WORD_BASES;
    static ssize_t const NMER_BITS = LB*NMER_BASES;

    union WordData
    {
        struct
        {
            uint32_t word     : WORD_BITS;
            bool erepeat      : 1;
            bool repeat       : 1;
        } w;

        uint32_t data;

        friend std::ostream & operator<<(
                std::ostream & os, WordData const & wd )
        {
            return os << std::hex
                      << "{ w: " << wd.w.word
                      << "; erep: " << wd.w.erepeat
                      << "; rep: " << wd.w.repeat
                      << " }" << std::dec;
        }
    };

public:

    struct FreqTableEntry
    {
        union
        {
            struct
            {
                /*
                uint64_t word   : WORD_BITS;
                uint64_t anchor : ANCHOR_BITS;
                */
                uint64_t nmer   : NMER_BITS;
                uint64_t freq   : 8;
            } f;

            uint64_t d;
        } data;

        FreqTableEntry() {}

        FreqTableEntry( uint64_t anchor, uint64_t word, uint64_t f = 0 )
        {
            data.f.nmer = (anchor<<WORD_BITS) + word;
            data.f.freq = f;
        }

        friend bool operator<(
            FreqTableEntry const & x, FreqTableEntry const & y )
        {
            return x.data.f.nmer < y.data.f.nmer;
        }
    };
};

READFINDER_NS_END

#endif

