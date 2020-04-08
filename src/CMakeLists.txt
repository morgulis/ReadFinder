cmake_minimum_required( VERSION 2.8.10 )
project( "readfinder" )

option( NOASSERT "disable asserts in Debug build" )

set( rf_includes ${CMAKE_CURRENT_BINARY_DIR}
                 ${CMAKE_CURRENT_SOURCE_DIR} )
set( rf_libdir "" )
set( rf_libs "" )

if( "$ENV{NCBI_BOOST}" STREQUAL "Y" )
    set( rf_includes ${rf_includes}
                     $ENV{BOOST_INCLUDE} )
    set( rf_libdir ${rf_libdir}
                   $ENV{BOOST_LIBDIR} )
    set( rf_libs ${rf_libs}
                 boost_program_options$ENV{BOOST_TAG}
                 boost_iostreams$ENV{BOOST_TAG} )
else()
    find_package( Boost 1.62 REQUIRED
                  COMPONENTS program_options
                             iostreams )
    set( rf_includes ${rf_includes}
                     ${Boost_INCLUDE_DIRS} )
    set( rf_libs ${rf_libs}
                 ${Boost_LIBRARIES} )
endif()

if( "${ENABLE_NGS}" STREQUAL "Y" )
    set( USE_NGS "1" )
    set( rf_includes ${rf_includes}
                     $ENV{NGS_INCLUDE}
                     $ENV{NGS_ROOT}/include
                     $ENV{VDB_ROOT}/include )
    set( rf_libdir ${rf_libdir}
                   $ENV{NGS_LIBDIR}
                   $ENV{NGS_ROOT}/lib64
                   $ENV{VDB_ROOT}/lib64 )
    set( rf_libs ${rf_libs}
                 ncbi-ngs-c++ ngs-c++ ncbi-vdb )
endif()

include_directories( ${rf_includes} )
link_directories( ${rf_libdir} )

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

