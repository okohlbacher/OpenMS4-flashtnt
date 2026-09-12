// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Oliver Kohlbacher $
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTHelpers.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/FORMAT/FLASHTnTFile.h>
#include <stdexcept>

int main()
{
  using OpenMS::FLASHTnTFile;
  auto require = [](bool condition) { if (!condition) throw std::runtime_error("FLASHTnT port contract failed"); };
  require(FLASHTnTFile::generateProFormaString_("PEPTIDE", 0, 7, {}, {}, {}, {}) == "PEPTIDE");
  require(FLASHTnTFile::generateProFormaString_("MPEPTIDE", 1, 8, {15.9949}, {2}, {2}, {"Oxidation"}) == "PE[Oxidation]PTIDE");
  require(FLASHTnTFile::generateProFormaString_("MPEPTIDE", 1, 8, {15.9949}, {2}, {2}, {"Oxidation,"}) == "PE[Oxidation]PTIDE");
  const auto ambiguous = FLASHTnTFile::generateProFormaString_("MPEPTIDE", 1, 8, {15.9949}, {2}, {2}, {"Oxidation,Hydroxylation,"});
  const auto parsed_ambiguous = OpenMS::ProForma::parse(ambiguous);
  const auto& ambiguous_residue = std::get<OpenMS::ProForma::SequenceElement>(parsed_ambiguous.sequence[1]);
  const auto& delta = std::get<OpenMS::ProForma::MassDelta>(ambiguous_residue.modifications[0].alternatives[0].first);
  require(std::abs(delta.mass - 15.9949) < 0.00001);
  const auto range = FLASHTnTFile::generateProFormaString_("PEPTIDE", 0, 7, {79.9663}, {1}, {3}, {"Phospho"});
  require(range == "P(EPT)[Phospho]IDE");
  require(OpenMS::ProForma::toString(OpenMS::ProForma::parse(range)) == range);
  require(FLASHTnTFile::generateProFormaString_("PEPTIDE", 0, 7, {42.0106}, {-1}, {-1}, {"Acetyl"}) == "[Acetyl]-PEPTIDE");
  bool rejected = false;
  try { FLASHTnTFile::generateProFormaString_("PEPTIDE", 0, 7, {15.0}, {}, {}, {}); }
  catch (const OpenMS::Exception::InvalidValue&) { rejected = true; }
  require(rejected);

  std::vector<double> masses{1.0, 2.0, 3.0};
  std::vector<int> scores{1, 2, 3};
  OpenMS::FLASHTnTHelpers::Tag tag("Pe", 0.0, -1.0, masses, scores, 10);
  require(tag.getSequence() == "PE" && tag.getUppercaseSequence() == "PE");
  require(tag.getScore() == 6 && tag.getIndex() == -1);
  OpenMS::FLASHTnTHelpers::DAG graph(3);
  boost::dynamic_bitset<> visited(3);
  visited[2] = true;
  require(graph.addEdge(1, 2, visited));
  require(graph.addEdge(0, 1, visited));
  require(!graph.addEdge(0, 3, visited));
  std::vector<std::vector<OpenMS::Size>> paths;
  graph.findAllPaths(0, 2, paths, 10);
  require(paths == std::vector<std::vector<OpenMS::Size>>{{0, 1, 2}});

  // Two PET ladders 0.5 Da apart exceed 5 ppm and must remain distinct tags.
  // An integer abs overload incorrectly truncates their difference to zero.
  OpenMS::DeconvolvedSpectrum ladders(1);
  double mass = 100.0;
  for (const auto& residue : std::vector<std::string>{"", "P", "E", "T"})
  {
    if (!residue.empty()) mass += OpenMS::AASequence::fromString(residue).getMonoWeight(OpenMS::Residue::Internal);
    for (double shift : {0.0, 0.5})
    {
      OpenMS::PeakGroup peak;
      peak.setMonoisotopicMass(mass + shift);
      peak.setQscore(1.0);
      ladders.push_back(peak);
    }
  }
  OpenMS::PeakGroup sink;
  sink.setMonoisotopicMass(600.0);
  sink.setQscore(1.0);
  ladders.push_back(sink);
  OpenMS::FLASHTaggerAlgorithm tagger;
  auto parameters = tagger.getParameters();
  parameters.setValue("min_length", 3);
  parameters.setValue("max_length", 3);
  tagger.setParameters(parameters);
  tagger.run(ladders, 5.0);
  std::vector<OpenMS::FLASHTnTHelpers::Tag> ladder_tags;
  tagger.fillTags(ladder_tags);
  std::set<double> nterm_masses;
  for (const auto& candidate : ladder_tags)
  {
    if (candidate.getSequence() == "PET" && candidate.getNtermMass() >= 0)
    {
      nterm_masses.insert(candidate.getNtermMass());
    }
  }
  require(nterm_masses == std::set<double>{100.0, 100.5});

  OpenMS::MSExperiment experiment;
  OpenMS::MSSpectrum spectrum;
  spectrum.setMSLevel(2);
  spectrum.resize(5);
  experiment.addSpectrum(spectrum);
  const std::vector<OpenMS::FASTAFile::FASTAEntry> proteins{{"target", "fixture", "PEPTIDE"}};
  rejected = false;
  try { OpenMS::FLASHTnTAlgorithm().run(experiment, proteins); }
  catch (const OpenMS::Exception::MissingInformation&) { rejected = true; }
  require(rejected);
  experiment[0].setMetaValue("DeconvMassInfo", "tol=10;precursorscan=-1;qscore=1,;snr=1,2,3,4,5,");
  rejected = false;
  try { OpenMS::FLASHTnTAlgorithm().run(experiment, proteins); }
  catch (const OpenMS::Exception::InvalidValue&) { rejected = true; }
  require(rejected);
}
