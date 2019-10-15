cmake_minimum_required( VERSION 2.8.10 )
project( "readfinder" )

option( NOASSERT "disable asserts in Debug build" )

#set( NGS_BASE_LCL $ENV{NGS_BASE} )

#if( "${NGS_BASE_LCL}" STREQUAL "" )
#    set( NGS_BASE_LCL "/panfs/traces01.be-md.ncbi.nlm.nih.gov/trace_software/ngs" )
#endif()


if( "${ENABLE_LCL_BOOST}" STREQUAL "Y" )
    if( "${BOOST_TAG}" STREQUAL "" )
        set( BOOST_TAG "-gcc73-mt-1_62-64" )
    endif()
    set( BOOST_HOME "/netopt/ncbi_tools64/boost-1.62.0-ncbi1" )
    set( BOOST_INC "${BOOST_HOME}/include" )
    set( BOOST_LIB "${BOOST_HOME}/lib" )
endif()

if( "${ENABLE_NGS}" STREQUAL "Y" )
set( USE_NGS "1" )
include_directories( ${CMAKE_CURRENT_BINARY_DIR}
                     ${CMAKE_CURRENT_SOURCE_DIR} 
                     ${BOOST_INC}
                     ${NGS_HOME}/include
                     ${VDB_HOME}/include )
link_directories( ${BOOST_LIB} ${NGS_HOME}/lib64 ${VDB_HOME}/lib64 )
else()
include_directories( ${CMAKE_CURRENT_BINARY_DIR}
                     ${CMAKE_CURRENT_SOURCE_DIR} 
                     ${BOOST_INC} )
link_directories( ${BOOST_LIB} )
endif()
#link_directories( ${NGS_BASE_LCL}/linux/release/x86_64/lib )

if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
    set( ATB_GCC "1" )
    set( CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -std=gnu++11 -Wall -pthread " )
    set( CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS_DEBUG} -O0 " )

    if( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
        if( ${NOASSERT} )
            set( CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS_DEBUG} -DNDEBUG " )
            message( "assertions are disabled" )
        else()
            message( "assertions are enabled" )
        endif()
    endif()
endif()

configure_file( config.h.cmk config.h @ONLY )

add_subdirectory( libtools )
add_subdirectory( libseq )
add_subdirectory( libreadfinder )
add_subdirectory( readfinder )

