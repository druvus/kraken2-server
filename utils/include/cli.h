#pragma once

// A small table-driven command line parser over getopt_long, shared by the
// server and the client. Each option is declared once; both the lower and
// upper case forms of the short flag are accepted, as the programs have
// always done. Parsing does not exit, so it can be unit tested; callers
// act on the returned status.

#include <functional>
#include <string>
#include <vector>

namespace cli
{
    struct Option
    {
        std::string long_name;   // e.g. "port" for --port
        char short_name;         // e.g. 'p'; 0 for none
        bool has_arg;            // takes a value
        std::string arg_name;    // shown in help, e.g. "[int]"
        std::string help;
        // Called with the value, or an empty string for flags. Throw
        // std::invalid_argument to report a bad value.
        std::function<void(const std::string &)> apply;
        bool required = false;
    };

    enum class Status { Ok, Help, Error };

    struct Result
    {
        Status status;
        std::string message;  // set for Error
    };

    class Parser
    {
    public:
        explicit Parser(std::string program);

        void add(Option option);

        // Parse argv. --help / -h / -H / -? give Status::Help. Unknown
        // options, missing values, values rejected by apply(), and missing
        // required options give Status::Error with a message.
        Result parse(int argc, char **argv) const;

        std::string usage() const;

        // Value helpers for apply() callbacks; throw std::invalid_argument
        // with a message naming the option when out of range.
        static int ParseInt(const std::string &value, const std::string &name, long min, long max);
        static double ParseDouble(const std::string &value, const std::string &name, double min, double max);

    private:
        std::string program_;
        std::vector<Option> options_;
    };
}
