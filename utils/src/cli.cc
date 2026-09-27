#include "cli.h"

#include <cctype>
#include <cstdlib>
#include <getopt.h>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace cli
{

Parser::Parser(std::string program) : program_(std::move(program)) {}

void Parser::add(Option option)
{
    options_.push_back(std::move(option));
}

int Parser::ParseInt(const std::string &value, const std::string &name, long min, long max)
{
    char *end = nullptr;
    long v = strtol(value.c_str(), &end, 10);
    if (value.empty() || *end != '\0' || v < min || v > max)
    {
        throw std::invalid_argument(
            name + " must be an integer in [" + std::to_string(min) + ", " + std::to_string(max) + "], got '" + value + "'");
    }
    return (int) v;
}

double Parser::ParseDouble(const std::string &value, const std::string &name, double min, double max)
{
    char *end = nullptr;
    double v = strtod(value.c_str(), &end);
    if (value.empty() || *end != '\0' || v < min || v > max)
    {
        std::ostringstream ss;
        ss << name << " must be a number in [" << min << ", " << max << "], got '" << value << "'";
        throw std::invalid_argument(ss.str());
    }
    return v;
}

Result Parser::parse(int argc, char **argv) const
{
    // Long options get values above the character range so they map back
    // to an option index without clashing with short flags.
    const int kLongBase = 1000;
    std::vector<struct option> longopts;
    std::string shortopts = "hH?";
    for (size_t i = 0; i < options_.size(); i++)
    {
        const Option &o = options_[i];
        longopts.push_back({o.long_name.c_str(), o.has_arg ? required_argument : no_argument,
                            nullptr, kLongBase + (int) i});
        if (o.short_name != 0)
        {
            shortopts += o.short_name;
            if (o.has_arg) shortopts += ':';
            if (isalpha(o.short_name))
            {
                shortopts += (char) toupper(o.short_name);
                if (o.has_arg) shortopts += ':';
            }
        }
    }
    longopts.push_back({"help", no_argument, nullptr, 'h'});
    longopts.push_back({nullptr, 0, nullptr, 0});
    // Leading ':' makes getopt report a missing value as ':' rather than '?'.
    shortopts = ":" + shortopts;

    std::vector<bool> seen(options_.size(), false);
    // Reset getopt state so parse() can be called more than once per process
    // (unit tests). BSD libc needs optreset as well as optind.
    optind = 1;
    opterr = 0;
    optopt = 0;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset = 1;
#endif
    int opt;
    while ((opt = getopt_long(argc, argv, shortopts.c_str(), longopts.data(), nullptr)) != -1)
    {
        if (opt == 'h' || opt == 'H')
        {
            return {Status::Help, ""};
        }
        if (opt == ':')
        {
            return {Status::Error, std::string("Option ") + argv[optind - 1] + " requires a value."};
        }
        if (opt == '?')
        {
            // '?' is both getopt's unknown-option code and a declared help
            // flag, and optopt is not set for a declared flag. Look at the
            // token itself: "-?" is help, anything else is an error.
            std::string token = argv[optind - 1];
            if (token == "-?") return {Status::Help, ""};
            std::string bad = optopt ? std::string("-") + (char) optopt : token;
            return {Status::Error, "Unknown option " + bad + "."};
        }
        size_t idx = options_.size();
        if (opt >= kLongBase)
        {
            idx = opt - kLongBase;
        }
        else
        {
            for (size_t i = 0; i < options_.size(); i++)
            {
                if (options_[i].short_name != 0 && tolower(options_[i].short_name) == tolower(opt))
                {
                    idx = i;
                    break;
                }
            }
        }
        if (idx >= options_.size())
        {
            return {Status::Error, "Unknown option."};
        }
        try
        {
            options_[idx].apply(options_[idx].has_arg ? std::string(optarg) : std::string());
            seen[idx] = true;
        }
        catch (const std::invalid_argument &ex)
        {
            return {Status::Error, ex.what()};
        }
    }
    for (size_t i = 0; i < options_.size(); i++)
    {
        if (options_[i].required && !seen[i])
        {
            return {Status::Error, "--" + options_[i].long_name + " is required."};
        }
    }
    return {Status::Ok, ""};
}

std::string Parser::usage() const
{
    std::ostringstream ss;
    ss << "Usage: " << program_ << " [options]\n\n"
       << "Options: (* mandatory)\n";
    auto line = [&](bool required, const std::string &flags, const std::string &help) {
        ss << (required ? "*" : " ") << "  " << std::left << std::setw(36) << flags << help << "\n";
    };
    line(false, "-h, -H, -?, --help", "Show this help");
    for (const Option &o : options_)
    {
        std::string flags;
        if (o.short_name != 0)
        {
            flags += std::string("-") + o.short_name;
            if (isalpha(o.short_name)) flags += std::string(", -") + (char) toupper(o.short_name);
            flags += ", ";
        }
        flags += "--" + o.long_name;
        if (o.has_arg) flags += " " + o.arg_name;
        line(o.required, flags, o.help);
    }
    return ss.str();
}

} // namespace cli
