// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Oliver Kohlbacher $
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTHelpers.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/FORMAT/FLASHTnTFile.h>
#include <stdexcept>

int main()
{
  using OpenMS::FLASHTnTFile;
  auto require = [](bool condition) { if (!condition) throw std::runtime_error("FLASHTnT port contract failed"); };
  require(FLASHTnTFile::generateProFormaString_("PEPTIDE", 0, 7, {}, {}, {}, {}) == "PEPTIDE");
  require(FLASHTnTFile::generateProFormaString_("MPEPTIDE", 1, 8, {15.9949}, {2}, {2}, {"Oxidation"}) == "PE[Oxidation]PTIDE");
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
}
