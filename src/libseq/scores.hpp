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
/*! \file libseq/scores.hpp
    \brief Scoring for non-spliced alignments.
*/

#ifndef LIBSEQ_SCORES_HPP
#define LIBSEQ_SCORES_HPP

#include <libseq/defs.hpp>

SEQ_NS_BEGIN

//------------------------------------------------------------------------------
typedef float TScore;   ///< Type to represent alignment scores/penalties.

//------------------------------------------------------------------------------
/** Penalties for simple (non-spliced) alignments. */
struct CPenalties
{
    TScore match         = -1.0,
           mismatch      = 1.0,
           gap_open      = 0.0,
           gap_extend    = 1.0,
           clip_len_base = 0.1; /**< Linear coefficient for computation of
                                     penalties of soft-clipped regions of a
                                     read. The penalty is proportional to
                                     the length of the region with this
                                     coefficient of proportionality. */

    friend std::ostream & operator<<( std::ostream & os, CPenalties const & x )
    {
        return os <<    "   match penalty: " << x.match << '\n'
                  <<    "   mismatch penalty: " << x.mismatch << '\n'
                  <<    "   gap open penalty: " << x.gap_open << '\n'
                  <<    "   gap extension penalty: " << x.gap_extend << '\n'
                  <<    "   softmask base penalty: "
                        << x.clip_len_base << '\n';
    }
};

SEQ_NS_END

#endif

