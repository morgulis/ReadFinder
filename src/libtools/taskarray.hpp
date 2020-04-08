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
/*! \file libtools/taskarray.hpp
    \brief Support for a array of parallel tasks.
*/

#ifndef LIBTOOLS_TASKARRAY_HPP
#define LIBTOOLS_TASKARRAY_HPP

#include <memory>
#include <vector>
#include <thread>

#include <libtools/exception.hpp>
#include <libtools/defs.hpp>

TOOLS_NS_BEGIN

//------------------------------------------------------------------------------
/** Parallel execution of several instances of the same task.

    T_Worker must implement the following interface.

        void operator()();

    \tparam T_Worker Worker task.
*/
template< typename T_Worker >
class CTaskArray
{
public:

    /** Instance constructor.

        Creates the necessary number of T_Worker instances.

        \tparam T_Args Type vector for T_Worker constructor parameters.

        \param [in] n_threads   Task pool to use.
        \param [in] args        Arguments forwarded to worker constructors.
    */
    template< typename ... T_Args > 
    CTaskArray( size_t n_threads, T_Args && ... args );

    /** Starts the task array.
        
        Starts all worker tasks. Waits for all tasks to complete.
    */
    void Start();

private:

    void WorkerProc( size_t idx ) throw();

    /** Type by which to refer to individual worker task instances. */
    typedef std::unique_ptr< T_Worker > Worker;

    std::vector< Worker > workers_;         ///< Worker task instances.
    std::vector< std::string > err_msg_;    ///< Error messages from workers.
};

//------------------------------------------------------------------------------
template< typename T_Worker >
template< typename ... T_Args >
inline CTaskArray< T_Worker >::CTaskArray(
        size_t n_threads, T_Args && ... args )
    : workers_( n_threads ), err_msg_( n_threads )
{
    for( auto & worker : workers_ )
    {
        worker.reset( new T_Worker( args ... ) );
    }
}

//------------------------------------------------------------------------------
template< typename T_Worker >
void CTaskArray< T_Worker >::WorkerProc( size_t idx ) throw()
{
    try
    {
        (*workers_[idx])();
    }
    catch( std::exception const & e )
    {
        err_msg_[idx] = e.what();
    }
    catch( ... )
    {
        err_msg_[idx] = "unknown exception";
    }
}

//------------------------------------------------------------------------------
template< typename T_Worker >
inline void CTaskArray< T_Worker >::Start()
{
    std::vector< std::thread * > worker_threads( workers_.size(), nullptr );

    for( size_t i( 0 ); i < workers_.size(); ++i )
    {
        worker_threads[i] = new std::thread( 
                [this]( size_t i ){ WorkerProc( i ); }, i );
    }

    std::string em;

    for( size_t i( 0 ); i < workers_.size(); ++i )
    {
        worker_threads[i]->join();

        if( !err_msg_[i].empty() )
        {
            em += std::to_string( i ) + ":" + err_msg_[i] + ";";
        }

        delete worker_threads[i];
        worker_threads[i] = nullptr;
    }

    if( !em.empty() )
    {
        M_THROW( em );
    }
}

TOOLS_NS_END

#endif

