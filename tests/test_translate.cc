#include <doctest/doctest.h>

#include <set>
#include <string>
#include <vector>

#include "aa_translate.h"

// Guards against the translation tables not being initialised: upstream
// kraken2 requires an explicit initLookUpTables() call, and without it every
// codon translates to 'K'. kraken2 indexes frames by the position of each
// codon's last base, so frame +0 is aa_seqs[2]; the checks below compare
// frame sets rather than positions.
TEST_CASE("TranslateToAllFrames gives the standard genetic code after initLookUpTables")
{
    kraken2::initLookUpTables();
    // ATG GCT AAA TGA : M A K *
    std::string dna = "ATGGCTAAATGA";
    std::vector<std::string> frames(6);
    kraken2::TranslateToAllFrames(dna, frames);

    std::set<std::string> forward(frames.begin(), frames.begin() + 3);
    std::set<std::string> reverse(frames.begin() + 3, frames.end());
    // +0: MAK*  +1: TGG CTA AAT = WLN  +2: GGC TAA ATG = G*M
    CHECK(forward == std::set<std::string>{"MAK*", "WLN", "G*M"});
    // reverse complement TCATTTAGCCAT: SFSH, HLA (CAT TTA GCC), I*P (ATT TAG CCA)
    CHECK(reverse == std::set<std::string>{"SFSH", "HLA", "I*P"});
}

TEST_CASE("TranslateToAllFrames marks ambiguous codons with X")
{
    kraken2::initLookUpTables();
    std::string dna = "ATGNCTAAA";
    std::vector<std::string> frames(6);
    kraken2::TranslateToAllFrames(dna, frames);
    std::set<std::string> forward(frames.begin(), frames.begin() + 3);
    CHECK(forward.count("MXK") == 1);
}
