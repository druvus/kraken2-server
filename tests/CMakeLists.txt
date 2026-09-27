# Unit tests with doctest (conda-forge package `doctest`). Run with ctest or
# ./tests/unit_tests after configuring with -DBUILD_TESTS=ON.
find_package(doctest CONFIG REQUIRED)

add_executable(unit_tests
    main.cc
    test_classify_core.cc
    test_report.cc
    test_thread_safe_queue.cc
    test_messages.cc
    test_cli.cc
    test_worker_pool.cc
    test_translate.cc
    test_merged_taxonomy.cc)

target_link_libraries(unit_tests
    doctest::doctest
    kraken2_server_core
    server_client_utils
    kraken2_proto
    classify
    ${_PROTOBUF_LIBPROTOBUF}
    ${_GRPC_GRPCPP}
    ${_ABSL_SYNC})

add_test(NAME unit_tests COMMAND unit_tests)

project_warnings(unit_tests)
