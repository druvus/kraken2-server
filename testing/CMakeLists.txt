# Test tools, built only with -DBUILD_TEST_TOOLS=ON
add_executable(raw_client raw_client.cc)

target_link_libraries(raw_client
    kraken2_proto
    ${_PROTOBUF_LIBPROTOBUF}
    ${_REFLECTION}
    ${_GRPC_GRPCPP})
