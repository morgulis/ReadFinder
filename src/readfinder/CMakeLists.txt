add_executable( rf readfinder.cpp )
set_target_properties( rf PROPERTIES OUTPUT_NAME readfinder )
target_link_libraries( rf readfinder seq tools ${rf_libs} dl z )

