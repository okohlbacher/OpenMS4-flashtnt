// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong, Ayesha Feroz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/CHEMISTRY/ProForma.h> // added for ProForma
#include <OpenMS/FORMAT/FLASHTnTFile.h>

namespace OpenMS
{
/**
  @brief FLASHTnT output *.tsv file format
   @ingroup FileIO
**/
void FLASHTnTFile::writeTagHeader(std::fstream& fs)
{
  fs << "TagIndex\tScan\tRetentionTime\tProteoformIndex\tProteinAccession\tProteinDescription\tTagSequence\tNmass\tCmass\tStartPosition\tDeltaMass\tL"
        "ength\tDeNovoScore\tMasses\tMassScores\n";
}

/// write header line for PrSM file
void FLASHTnTFile::writePrSMHeader(std::fstream& fs)
{
  fs << "PrSMIndex\tScan\tRetentionTime\tNumMass\tProteinAccession\tProteinDescription\tPrecursorMass\tProteoformMass\tProteoformMassFromComplementaryFragmentIonPairs\tDatabaseSequence\tProteinSequence\tProf"
        "orma\tMatchingFragments\tCoverage(%)\tStartPosition\tEndPosition"
        "\tTagCount\tTagIndices\tModCount\tModMass\tModID\tModAccession\tModStart\tModEnd\tPrecursorQscore\tPrecursorSNR\tScore\tPrSMLevelQvalue\tProteoformLevelQvalue\n";
}

/// write header line for Proteoform file
void FLASHTnTFile::writeProHeader(std::fstream& fs)
{
  fs << "ProteoformIndex\tPrSMIndices\tScan\tRetentionTime\tNumMass\tProteinAccession\tProteinDescription\tPrecursorMass\tProteoformMass\tProteoformMassFromComplementaryFragmentIonPairs\tDatabaseSequence\tProteinSequence\tProf"
        "orma\tMatchingFragments\tCoverage(%)\tStartPosition\tEndPosition"
        "\tTagCount\tTagIndices\tModCount\tModMass\tModID\tModAccession\tModStart\tModEnd\tPrecursorQscore\tPrecursorSNR\tScore\tProteoformLevelQvalue\n";
}

/// write the features in regular file output
void FLASHTnTFile::writeTags(const FLASHTnTAlgorithm& tnt, double flanking_mass_tol, std::fstream& fs)
{
  std::stringstream ss;
  auto tags = std::vector<FLASHTnTHelpers::Tag>();
  tnt.getTags(tags);

  for (int c = 0; c < 2; c++)
  {
    for (const auto& tag : tags)
    {
      auto hits = std::vector<ProteinHit>();
      tnt.getProteoformHitsMatchedBy(tag, hits);
      if (c == 0 && hits.empty()) continue;
      if (c == 1 && ! hits.empty()) continue;

      std::string acc = "";
      std::string description = "";
      std::string hitindices = "";
      std::string positions = "";
      std::string delta_masses = "";
      for (const auto& hit : hits)
      {
        if (! hitindices.empty()) {
          acc += ";";
          description += ";";
          hitindices += ";";
          positions += ";";
          delta_masses += ";";
        }

        acc += hit.getAccession();

        std::string proteindescription = hit.getDescription();
        if (proteindescription.empty()) { proteindescription = " "; }
        description += proteindescription;
        hitindices += (std::string)hit.getMetaValue("Index");

        auto pos = std::vector<int>();
        auto masses = std::vector<double>();
        auto pos_in_truncated = std::vector<int>();
        auto masses_in_truncated = std::vector<double>();

        int protein_start_position = hit.getMetaValue("StartPosition");
        int protein_end_position = hit.getMetaValue("EndPosition");
        protein_start_position--;
        std::string seq = hit.getSequence();
        FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs(pos, masses, -1, seq, tag);

        if (protein_end_position >= 0) seq = seq.substr(0, protein_end_position);
        if (protein_start_position >= 0) seq = seq.substr(protein_start_position);
        FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs(pos_in_truncated, masses_in_truncated, flanking_mass_tol, seq, tag);
        if (!pos_in_truncated.empty())
        {
          for (int i = 0; i < pos.size(); i++)
          {
            if (pos[i] - (protein_start_position < 0 ? 0 : protein_start_position) != pos_in_truncated[0]) continue;
            positions += std::to_string(pos[i] + 1);
            delta_masses += std::to_string(masses[i]);
            break;
          }
        }
      }

      ss << tag.getIndex() << "\t" << tag.getScan() << "\t" << tag.getRetentionTime() << "\t" << hitindices << "\t" << acc << "\t" << description
         << "\t" << tag.getSequence() << "\t" << std::to_string(tag.getNtermMass()) << "\t" << std::to_string(tag.getCtermMass()) << "\t" << positions
         << "\t" << delta_masses << "\t" << tag.getLength() << "\t" << tag.getScore() << "\t";

      for (const auto& mz : tag.getMzs())
      {
        ss << std::to_string(mz) << ",";
      }
      ss << "\t";
      for (size_t i = 0; i <= tag.getLength(); i++) // Fixed signed/unsigned comparison issue
      {
        ss << std::to_string(tag.getScore(i)) << ",";
      }
      ss << "\n";
    }
  }
  fs << ss.str();
}

std::string FLASHTnTFile::generateProFormaString_(const std::string& sequence,
                                             int seq_start,
                                             int seq_end,
                                             const std::vector<double>& mod_masses,
                                             const std::vector<int>& mod_starts,
                                             const std::vector<int>& mod_ends,
                                             const std::vector<std::string>& mod_ids)
{
  if (seq_start < 0) seq_start = 0;
  if (seq_end < 0) seq_end = sequence.length();
  const auto truncated_seq = sequence.substr(seq_start, seq_end - seq_start);
  auto proforma = ProForma::fromAASequence(AASequence::fromString(truncated_seq));
  if (mod_starts.size() != mod_masses.size() || mod_ends.size() != mod_masses.size() || mod_ids.size() != mod_masses.size())
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "Inconsistent FLASHTnT modification records", sequence);
  }
  std::map<int, std::pair<int, ProForma::Modification>> modifications;
  for (size_t i = 0; i < mod_masses.size(); ++i)
  {
    ProForma::Modification modification;
    if (!mod_ids[i].empty())
    {
      modification.alternatives.emplace_back(ProForma::NamedMod{std::nullopt, mod_ids[i]}, std::nullopt);
    }
    else
    {
      modification.alternatives.emplace_back(ProForma::MassDelta{ProForma::MassDelta::Source::NONE, mod_masses[i], {}}, std::nullopt);
    }
    const int start = mod_starts[i] - seq_start;
    const int end = mod_ends[i] - seq_start;
    if (start == -1 && end == -1)
    {
      proforma.n_term_mods.push_back(modification);
    }
    else if (start == static_cast<int>(truncated_seq.size()) && end == start)
    {
      proforma.c_term_mods.push_back(modification);
    }
    else if (start < 0 || end < start || end >= static_cast<int>(truncated_seq.size()) || !modifications.emplace(start, std::make_pair(end, modification)).second)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Invalid FLASHTnT modification position", std::to_string(start));
    }
  }
  std::vector<ProForma::SequenceSection> sections;
  for (int pos = 0; pos < static_cast<int>(proforma.sequence.size()); ++pos)
  {
    auto element = std::get<ProForma::SequenceElement>(proforma.sequence[pos]);
    auto found = modifications.find(pos);
    if (found == modifications.end())
    {
      sections.emplace_back(element);
    }
    else if (found->second.first == pos)
    {
      element.modifications.push_back(found->second.second);
      sections.emplace_back(element);
    }
    else
    {
      ProForma::ModifiedRange range;
      const int end = found->second.first;
      for (int index = pos; index <= end; ++index)
      {
        if (index > pos && modifications.contains(index))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Overlapping FLASHTnT modification ranges", std::to_string(index));
        }
        range.elements.push_back(std::get<ProForma::SequenceElement>(proforma.sequence[index]));
      }
      range.modifications.push_back(found->second.second);
      sections.emplace_back(range);
      pos = end;
    }
  }
  proforma.sequence = std::move(sections);
  return ProForma::toString(proforma);
}


void OpenMS::FLASHTnTFile::writePrSMs(const std::vector<ProteinHit>& hits, std::fstream& fs)
{
  std::stringstream ss;
  for (const auto& hit : hits)
  {
    if (! hit.metaValueExists("Index")) continue;
    std::string tagindices = "", modmasses = "", modstarts = "", modends = "", modids = "", modaccs = "";

    int cntr = 0;
    std::vector<FLASHTnTHelpers::Tag> tags;
    std::vector<int> indices;
    if (hit.metaValueExists("TagIndices")) indices = (std::vector<int>)hit.getMetaValue("TagIndices").toIntList();
    for (const int& index : indices)
    {
      if (! tagindices.empty()) tagindices += ";";
      tagindices += std::to_string(index);
      cntr++;
    }

    std::vector<double> mod_masses = hit.getMetaValue("Modifications");
    std::vector<int> mod_starts = hit.getMetaValue("ModificationStarts");
    std::vector<int> mod_ends = hit.getMetaValue("ModificationEnds");
    std::vector<std::string> mod_ids = hit.getMetaValue("ModificationIDs");
    std::vector<std::string> mod_accs = hit.getMetaValue("ModificationACCs");

    for (size_t i = 0; i < mod_masses.size(); i++) // Fixed signed/unsigned comparison issue
    {
      if (i > 0)
      {
        modmasses += ";";
        modstarts += ";";
        modends += ";";
        modids += ";";
        modaccs += ";";
      }

      modmasses += std::to_string(mod_masses[i]);
      modstarts += std::to_string(mod_starts[i] + 1);
      modends += std::to_string(mod_ends[i] + 1);
      modids += mod_ids[i];
      modaccs += mod_accs[i];
    }

    int start = hit.getMetaValue("StartPosition");
    int end = hit.getMetaValue("EndPosition");

    int start_in_seq = start < 0 ? 0 : (start - 1);
    int end_in_seq = end < 0 ? hit.getSequence().size() : end;

    // Use ProForma for sequence generation
    std::string proformaStr = generateProFormaString_(hit.getSequence(), start_in_seq, end_in_seq, mod_masses, mod_starts, mod_ends, mod_ids);
    ss << hit.getMetaValue("Index") << "\t" << hit.getMetaValue("Scan") << "\t" << hit.getMetaValue("RT") << "\t" << hit.getMetaValue("NumMass")
       << "\t" << hit.getAccession() << "\t" << hit.getDescription() << "\t" << hit.getMetaValue("GivenMass")  << "\t" << hit.getMetaValue("Mass") << "\t" << ((int)hit.getMetaValue("ProteoformMassByFragmentMass") > 0) << "\t" <<  hit.getSequence() << "\t"
       << hit.getSequence().substr(start_in_seq, end_in_seq - start_in_seq) << "\t" << proformaStr << "\t" << hit.getMetaValue("MatchedAA") << "\t"
       << 100.0 * hit.getCoverage() << "\t" << start << "\t" << end << "\t" << cntr << "\t" << tagindices << "\t" << mod_masses.size() << "\t"
       << modmasses << "\t" << modids << "\t" << modaccs << "\t" << modstarts << "\t" << modends << "\t" << hit.getMetaValue("PrecursorScore")  << "\t" << hit.getMetaValue("PrecursorSNR")  << "\t" << hit.getScore() << "\t"
       << std::to_string((hit.metaValueExists("qvalue") ? (double)hit.getMetaValue("qvalue") : -1)) << "\t"
       << std::to_string((hit.metaValueExists("proqvalue") ? (double)hit.getMetaValue("proqvalue") : -1)) << "\n";
  }
  fs << ss.str();
}

void OpenMS::FLASHTnTFile::writeProteoforms(const std::vector<ProteinHit>& hits, std::fstream& fs, double pro_fdr)
{
  std::stringstream ss;
  int pro_index = 0;
  for (const auto& hit : hits)
  {
    if (! hit.metaValueExists("Index")) continue;
    if ((int)hit.getMetaValue("Representative") == 0) continue;
    if (hit.metaValueExists("proqvalue") && (double)hit.getMetaValue("proqvalue") > pro_fdr) continue;
    std::string prsmindices = "", tagindices = "", modmasses = "", modstarts = "", modends = "", modids = "", modaccs = "";

    int cntr = 0;
    //std::vector<FLASHTnTHelpers::Tag> tags;
    std::vector<int> tag_indices;
    if (hit.metaValueExists("TagIndices")) tag_indices = (std::vector<int>)hit.getMetaValue("TagIndices").toIntList();
    for (const int& index : tag_indices)
    {
      if (! tagindices.empty()) tagindices += ";";
      tagindices += std::to_string(index);
      cntr++;
    }

    std::vector<double> mod_masses = hit.getMetaValue("Modifications");
    std::vector<int> mod_starts = hit.getMetaValue("ModificationStarts");
    std::vector<int> mod_ends = hit.getMetaValue("ModificationEnds");
    std::vector<std::string> mod_ids = hit.getMetaValue("ModificationIDs");
    std::vector<std::string> mod_accs = hit.getMetaValue("ModificationACCs");

    for (size_t i = 0; i < mod_masses.size(); i++) // Fixed signed/unsigned comparison issue
    {
      if (i > 0)
      {
        modmasses += ";";
        modstarts += ";";
        modends += ";";
        modids += ";";
        modaccs += ";";
      }

      modmasses += std::to_string(mod_masses[i]);
      modstarts += std::to_string(mod_starts[i] + 1);
      modends += std::to_string(mod_ends[i] + 1);
      modids += mod_ids[i];
      modaccs += mod_accs[i];
    }
    if (hit.metaValueExists("PrSMIndices"))
    {
      auto prsm_indices = hit.getMetaValue("PrSMIndices").toIntList();
      for (const int& index : prsm_indices)
      {
        if (! prsmindices.empty()) prsmindices += ";";
        prsmindices += std::to_string(index);
      }
    }
    int start = hit.getMetaValue("StartPosition");
    int end = hit.getMetaValue("EndPosition");

    int start_in_seq = start < 0 ? 0 : (start - 1);
    int end_in_seq = end < 0 ? hit.getSequence().size() : end;

    // Use ProForma to generate ProForma string
    std::string proformaStr = generateProFormaString_(hit.getSequence(), start_in_seq, end_in_seq, mod_masses, mod_starts, mod_ends, mod_ids);

    ss << pro_index++ << "\t" << prsmindices << "\t" << hit.getMetaValue("Scan") << "\t" << hit.getMetaValue("RT") << "\t" << hit.getMetaValue("NumMass")
       << "\t" << hit.getAccession() << "\t" << hit.getDescription() << "\t" << hit.getMetaValue("GivenMass") << "\t" << hit.getMetaValue("Mass") << "\t" << ((int)hit.getMetaValue("ProteoformMassByFragmentMass") > 0) << "\t" << hit.getSequence() << "\t"
       << hit.getSequence().substr(start_in_seq, end_in_seq - start_in_seq) << "\t" << proformaStr << "\t" << hit.getMetaValue("MatchedAA") << "\t"
       << 100.0 * hit.getCoverage() << "\t" << start << "\t" << end << "\t" << cntr << "\t" << tagindices << "\t" << mod_masses.size() << "\t"
       << modmasses << "\t" << modids << "\t" << modaccs << "\t" << modstarts << "\t" << modends << "\t" << hit.getMetaValue("PrecursorScore") << "\t" << hit.getMetaValue("PrecursorSNR") << "\t" << hit.getScore() << "\t"
       << std::to_string((hit.metaValueExists("proqvalue") ? (double)hit.getMetaValue("proqvalue") : -1)) << "\n";
  }
  fs << ss.str();
}

} // namespace OpenMS
// namespace OpenMS
