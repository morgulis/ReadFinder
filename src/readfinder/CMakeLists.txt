add_executable( rf readfinder.cpp )

set_target_properties( rf PROPERTIES OUTPUT_NAME readfinder )
target_link_libraries( rf readfinder seq tools
                          ncbi-ngs-c++ ngs-c++ ncbi-vdb-static
                          boost_program_options
                          boost_iostreams
                          dl z bz2 )

