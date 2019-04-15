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
//------------------------------------------------------------------------------
/*! \file libseq/defs.hpp 
    \brief High level definitions for sequence tools library. 
*/

#ifndef LIBSEQ_DEFS_HPP
#define LIBSEQ_DEFS_HPP

#include <cstdint>
#include <string>

#include <libtools/defs.hpp>

#include <config.h>

//------------------------------------------------------------------------------
/*! \brief Internal namespace name for libseq entities.

    If macro \c OUTER_NS is defined, then its value becomes the enclosing
    namespace for \c ATB_NS::SEQ_NS_LCL.
*/
#define SEQ_NS_LCL seq

//------------------------------------------------------------------------------
/*! \def SEQ_NS
    \brief Fully qualified namespace name for libseq entities.

    The value should be used when a qualifier for a libseq entity 
    is needed.
*/
/*! \def SEQ_NS_BEGIN
    \brief Open libseq namespace.
*/
/*! \def SEQ_NS_END
    \brief Close libseq namespace.
*/
#ifdef OUTER_NS
#   define SEQ_NS OUTER_NS::ATB_NS::SEQ_NS_LCL
#   define SEQ_NS_BEGIN namespace OUTER_NS { \
                             namespace ATB_NS { \
                             namespace SEQ_NS_LCL {
#   define SEQ_NS_END }}}
#else
#   define SEQ_NS ATB_NS::SEQ_NS_LCL
#   define SEQ_NS_BEGIN namespace ATB_NS { namespace SEQ_NS_LCL {
#   define SEQ_NS_END }}
#endif

SEQ_NS_BEGIN

//------------------------------------------------------------------------------
typedef int32_t TSeqLen;    ///< Type to represent sequence lengths.
typedef int32_t TSeqOff;    ///< Type to represent offsets into sequences.
typedef TSeqOff TDiag;      ///< Type to represent hit diagonals.

/// Type to use for textual sequence ids/accessions.
typedef std::string TSeqId; 

/// Constants to represent sequence strands.
typedef uint8_t EStrand;
static EStrand const eNONE = 0;  ///< No strand.
static EStrand const eFWD  = 1;  ///< Forward strand.
static EStrand const eREV  = 2;  ///< Reverse strand.
static EStrand const eANY  = 3;  ///< Any (undetermined but present).
static EStrand const MAX_STRANDS = 3;   ///< Number of strand values.

/// Constants to represent read mates.
typedef uint8_t EMate;
static EMate const eBOTH = 0;   ///< Both mates.
static EMate const eFIRST = 1;  ///< First mate.
static EMate const eSECOND = 2; ///< Second mate.
static EMate const MAX_MATES = 3;   ///< Number of mate type constants.

SEQ_NS_END

#endif

