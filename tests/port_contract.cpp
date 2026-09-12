// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Oliver Kohlbacher $
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTHelpers.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
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
  require(tag.getSequence() == "Pe" && tag.getUppercaseSequence() == "PE");
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
