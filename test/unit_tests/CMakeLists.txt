add_executable(geodetic geodetic.cpp)
add_executable(spherical spherical.cpp)
add_executable(typeWrappers type_wrappers.cpp)
add_executable(typeWrappersCartesian type_wrappers_cartesian.cpp)

target_link_libraries(geodetic PRIVATE geodesy)
target_link_libraries(spherical PRIVATE geodesy)
target_link_libraries(typeWrappers PRIVATE geodesy)
target_link_libraries(typeWrappersCartesian PRIVATE geodesy)

add_test(NAME geodetic COMMAND geodetic)
add_test(NAME spherical COMMAND spherical)
add_test(NAME typeWrappers COMMAND typeWrappers)
add_test(NAME typeWrappersCartesian COMMAND typeWrappers)
