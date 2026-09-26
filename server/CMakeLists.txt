# Pure classification helpers and report wrapper, shared with the unit tests
add_library(kraken2_server_core
    classify_core.cc
    report_server.cc)

target_include_directories(kraken2_server_core PUBLIC .)

target_link_libraries(kraken2_server_core
    classify)

# Create executable for the server
add_executable(kraken2_server
    kraken2_server.cc
    classify_server.cc)

target_include_directories(kraken2_server PUBLIC .)

# Add dependencies / links for this executable
target_link_libraries(kraken2_server
    kraken2_server_core
    kraken2_proto # Proto files lib
    server_client_utils # Additional utils lib
    classify
    ${_PROTOBUF_LIBPROTOBUF}
    ${_REFLECTION}
    ${_GRPC_GRPCPP})
