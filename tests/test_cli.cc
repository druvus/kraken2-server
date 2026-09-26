#include <doctest/doctest.h>

#include <string>
#include <vector>

#include "cli.h"

namespace
{
    struct Args
    {
        std::vector<std::string> store;
        std::vector<char *> ptrs;
        Args(std::initializer_list<const char *> a) : store(a.begin(), a.end())
        {
            for (auto &s : store) ptrs.push_back(&s[0]);
            ptrs.push_back(nullptr);
        }
        int argc() { return (int) store.size(); }
        char **argv() { return ptrs.data(); }
    };

    struct Values
    {
        std::string db;
        int port = 8080;
        double confidence = 0.0;
        bool flag = false;
    };

    cli::Parser MakeParser(Values &v)
    {
        cli::Parser p("prog");
        p.add({"db", 'd', true, "[path]", "Database", [&](const std::string &s) { v.db = s; }, true});
        p.add({"port", 'p', true, "[int]", "Port",
               [&](const std::string &s) { v.port = cli::Parser::ParseInt(s, "--port", 0, 65535); }});
        p.add({"confidence", 'c', true, "[double]", "Confidence",
               [&](const std::string &s) { v.confidence = cli::Parser::ParseDouble(s, "--confidence", 0, 1); }});
        p.add({"flag", 'f', false, "", "A flag", [&](const std::string &) { v.flag = true; }});
        return p;
    }
}

TEST_CASE("cli parses long options, short options and upper case aliases")
{
    Values v;
    auto p = MakeParser(v);
    Args a{"prog", "--db", "/x", "-P", "9000", "-c", "0.5", "-F"};
    auto r = p.parse(a.argc(), a.argv());
    CHECK(r.status == cli::Status::Ok);
    CHECK(v.db == "/x");
    CHECK(v.port == 9000);
    CHECK(v.confidence == doctest::Approx(0.5));
    CHECK(v.flag);
}

TEST_CASE("cli accepts --opt=value")
{
    Values v;
    auto p = MakeParser(v);
    Args a{"prog", "--db=/y", "--port=1"};
    CHECK(p.parse(a.argc(), a.argv()).status == cli::Status::Ok);
    CHECK(v.db == "/y");
    CHECK(v.port == 1);
}

TEST_CASE("cli reports help for -h, -H, -? and --help")
{
    Values v;
    auto p = MakeParser(v);
    for (const char *h : {"-h", "-H", "-?", "--help"})
    {
        Args a{"prog", h};
        CHECK(p.parse(a.argc(), a.argv()).status == cli::Status::Help);
    }
}

TEST_CASE("cli errors")
{
    Values v;
    auto p = MakeParser(v);

    SUBCASE("missing required option")
    {
        Args a{"prog", "-p", "80"};
        auto r = p.parse(a.argc(), a.argv());
        CHECK(r.status == cli::Status::Error);
        CHECK(r.message == "--db is required.");
    }
    SUBCASE("missing value")
    {
        Args a{"prog", "--db"};
        auto r = p.parse(a.argc(), a.argv());
        CHECK(r.status == cli::Status::Error);
        CHECK(r.message.find("requires a value") != std::string::npos);
    }
    SUBCASE("unknown option")
    {
        Args a{"prog", "--db", "x", "--bogus"};
        auto r = p.parse(a.argc(), a.argv());
        CHECK(r.status == cli::Status::Error);
        CHECK(r.message.find("Unknown option") != std::string::npos);
    }
    SUBCASE("value out of range")
    {
        Args a{"prog", "--db", "x", "--port", "70000"};
        auto r = p.parse(a.argc(), a.argv());
        CHECK(r.status == cli::Status::Error);
        CHECK(r.message.find("--port") != std::string::npos);
    }
    SUBCASE("value not a number")
    {
        Args a{"prog", "--db", "x", "--confidence", "high"};
        auto r = p.parse(a.argc(), a.argv());
        CHECK(r.status == cli::Status::Error);
        CHECK(r.message.find("--confidence") != std::string::npos);
    }
}

TEST_CASE("cli usage lists every option once with both short forms")
{
    Values v;
    auto p = MakeParser(v);
    std::string u = p.usage();
    CHECK(u.find("Usage: prog [options]") == 0);
    CHECK(u.find("-d, -D, --db [path]") != std::string::npos);
    CHECK(u.find("-f, -F, --flag") != std::string::npos);
    CHECK(u.find("*  -d") != std::string::npos);  // required marker
}

TEST_CASE("cli value helpers")
{
    CHECK(cli::Parser::ParseInt("42", "x", 0, 100) == 42);
    CHECK_THROWS_AS(cli::Parser::ParseInt("-1", "x", 0, 100), std::invalid_argument);
    CHECK_THROWS_AS(cli::Parser::ParseInt("4x", "x", 0, 100), std::invalid_argument);
    CHECK_THROWS_AS(cli::Parser::ParseInt("", "x", 0, 100), std::invalid_argument);
    CHECK(cli::Parser::ParseDouble("0.25", "x", 0, 1) == doctest::Approx(0.25));
    CHECK_THROWS_AS(cli::Parser::ParseDouble("1.5", "x", 0, 1), std::invalid_argument);
}
