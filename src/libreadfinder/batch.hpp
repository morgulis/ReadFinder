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

/** \file libreadfinder/batch.hpp
    \brief Single batch of reads.
*/

#ifndef LIBREADFINDER_BATCH_HPP
#define LIBREADFINDER_BATCH_HPP

#include <condition_variable>
#include <memory>
#include <mutex>

#include <libreadfinder/rf_ctx.hpp>
#include <libreadfinder/readdata.hpp>
#include <libreadfinder/defs.hpp>

READFINDER_NS_BEGIN

//==============================================================================
class CBatch
{
public:

    CBatch( CSearchContext & ctx, size_t batch_num = 0 );
    bool Run();
    bool RunSeeder();

    CSearchContext const & GetSearchCtx() const { return ctx_; }
    CSearchContext & GetSearchCtx() { return ctx_; }
    CReadData const & GetReads() const { return *reads_; }
    size_t GetJobSize() const { return reads_per_job_; }

private:

    typedef std::vector< std::string * > OutStr;

    static std::vector< std::string > const BatchStatDescriptions;

    void OutputThread();

    CSearchContext & ctx_;
    std::unique_ptr< CReadData > reads_;
    OutStr out_str_;
    std::mutex out_mtx_;
    std::condition_variable out_cvar_;
    size_t reads_per_job_;
    size_t batch_num_ = 0;
};

//------------------------------------------------------------------------------
inline std::unique_ptr< CBatch > MakeBatch(
        CSearchContext & ctx, size_t batch_num )
{
    return std::unique_ptr< CBatch >( new CBatch( ctx, batch_num ) );
}

READFINDER_NS_END

#endif

