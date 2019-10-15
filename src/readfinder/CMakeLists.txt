add_executable( rf readfinder.cpp )

set_target_properties( rf PROPERTIES OUTPUT_NAME readfinder )

if( "${ENABLE_LCL_BOOST}" STREQUAL "Y" )
    set( BOOST_LIBS boost_program_options${BOOST_TAG}-static boost_iostreams${BOOST_TAG}-static )
else()
    set( BOOST_LIBS boost_program_options boost_iostreams )
endif()

if( "${ENABLE_NGS}" STREQUAL "Y" )
target_link_libraries( rf readfinder seq tools
                          ncbi-ngs-c++ ngs-c++ ncbi-vdb-static
                          ${BOOST_LIBS}
                          dl z bz2 )
else()
target_link_libraries( rf readfinder seq tools
                          ${BOOST_LIBS}
                          dl z bz2 )
endif()

