cmake_minimum_required( VERSION 2.8.10 )
project( "readfinder" )

option( NOASSERT "disable asserts in Debug build" )

set( NGS_BASE_LCL $ENV{NGS_BASE} )

if( "${NGS_BASE_LCL}" STREQUAL "" )
    set( NGS_BASE_LCL "/panfs/traces01.be-md.ncbi.nlm.nih.gov/trace_software/ngs" )
endif()

include_directories( ${CMAKE_CURRENT_BINARY_DIR}
                     ${CMAKE_CURRENT_SOURCE_DIR} 
                     ${NGS_BASE_LCL}/include )

link_directories( ${NGS_BASE_LCL}/linux/release/x86_64/lib )

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

