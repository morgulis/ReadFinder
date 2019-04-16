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
    static ssize_t const SFX_BASES  = 3;
    static ssize_t const NMER_BASES = ANCHOR_BASES + WORD_BASES;
    static ssize_t const EXT_NMER_BASES = NMER_BASES + SFX_BASES;
    static ssize_t const MIN_INEXACT_BASES = 10;

protected:

    static ssize_t const ANCHOR_BITS  = LB*ANCHOR_BASES;
    static ssize_t const WORD_BITS = LB*WORD_BASES;
    static ssize_t const SFX_BITS  = LB*SFX_BASES;
    static ssize_t const NMER_BITS = LB*NMER_BASES;
    static ssize_t const EXT_NMER_BITS = LB*EXT_NMER_BASES;

    union WordData
    {
        struct
        {
            uint32_t word     : WORD_BITS;
            uint32_t sfx      : SFX_BITS;
            bool erepeat      : 1;
            bool repeat       : 1;
        } w;

        uint32_t data;
    };
};

READFINDER_NS_END

#endif

