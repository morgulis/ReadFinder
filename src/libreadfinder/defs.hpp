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

#ifndef LIBREADFINDER_DEFS_HPP
#define LIBREADFINDER_DEFS_HPP

#include <cmath>
#include <vector>

#include <libseq/coding.hpp>
#include <libseq/scores.hpp>
#include <libseq/defs.hpp>

#include <config.h>

#define READFINDER_NS_LCL readfinder

#ifdef OUTER_NS
#   define READFINDER_NS OUTER_NS::ATB_NS::READFINDER_NS_LCL
#   define READFINDER_NS_BEGIN namespace OUTER_NS { \
                              namespace ATB_NS { \
                              namespace READFINDER_NS_LCL {
#   define READFINDER_NS_END }}}
#else
#   define READFINDER_NS ATB_NS::READFINDER_NS_LCL
#   define READFINDER_NS_BEGIN namespace ATB_NS { namespace READFINDER_NS_LCL {
#   define READFINDER_NS_END }}
#endif

READFINDER_NS_BEGIN

using namespace TOOLS_NS;
using namespace SEQ_NS;

//------------------------------------------------------------------------------
typedef uint64_t TWord;     ///< Common word size used to store sequence data.

/** Number of bits in a NCBI2NA letter. */
static size_t const LB = CCode< eNCBI2NA >::LBITS;

/** NCBI2NA letter mask */
static TWord const LMASK = (1<<LB) - 1;

/** Shift to apply when converting NCBI2NA letters to bits. */
static size_t const LSHIFT = 1;

/** Number of NCBI2NA letters in a word. */
static TSeqLen const WL = sizeof( TWord )*CBL< LB >::V;

typedef int16_t TReadLen;   ///< Type to represent short read lengths.
typedef int16_t TReadOff;   ///< Type to represent offsets into short reads.

typedef uint32_t OrdId; ///< Read ordinal id.

/// Storage for subject letters in spliced alignments.
typedef std::vector< uint8_t > TLetters;

/// Splice site representation.
typedef std::pair< TSeqOff, TSeqOff > SpliceSite;

/// No splice indicator.
static SpliceSite const NO_SPLICE = SpliceSite { -1, -1 };

/// Max supported mate length.
static TReadLen const MAX_READ_LEN = std::numeric_limits< TReadLen >::max();

/// Largest (worst) score to start with.
static TScore const BIG_SCORE = std::numeric_limits< TScore >::max()/2.0;

/// Compute intron penalty.
inline TScore GetIntronPenalty(
        double scale, double exp_scale, TSeqLen ilen, TReadLen rlen )
{
    return rlen*scale*( exp( exp_scale*ilen ) - 1.0 );
}

inline uint8_t GetMateIdx( uint8_t sidx, uint8_t csidx )
{
    return (sidx + csidx)%2;
}

READFINDER_NS_END

#endif

