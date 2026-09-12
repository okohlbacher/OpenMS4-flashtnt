// Copyright (c) 2002-2025, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTHelpers.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <boost/dynamic_bitset.hpp>
#include <iomanip>
#include <iostream>

namespace OpenMS
{
/**
@brief Identify top-down proteoforms by matching sequence tags and extending candidate protein hits.
@ingroup Topdown
*/

class FLASHTnTAlgorithm : public DefaultParamHandler, public ProgressLogger
{
public:
  /// constructor
  FLASHTnTAlgorithm();

  /// destructor
  ~FLASHTnTAlgorithm() override = default;

  /// copy constructor
  FLASHTnTAlgorithm(const FLASHTnTAlgorithm&);

  /// move constructor
  FLASHTnTAlgorithm(FLASHTnTAlgorithm&& other) = default;

  /// assignment operator
  FLASHTnTAlgorithm& operator=(const FLASHTnTAlgorithm& other);

  /**
    @brief Generate sequence tags, match the protein database, and score extended proteoform candidates.

    @param[in] map Spectra deconvolved by FLASHDeconv, including DeconvMassInfo metadata.
    @param[in] fasta_entry Target and optional decoy protein sequences to search.

    @throws Exception::MissingInformation if FLASHDeconv metadata is missing.
    @throws Exception::InvalidValue if metadata dimensions are invalid or no target proteins are supplied.

  */
  void run(const MSExperiment& map, const std::vector<FASTAFile::FASTAEntry>& fasta_entry);
  void getProteoformHitsMatchedBy(const FLASHTnTHelpers::Tag& tag, std::vector<ProteinHit>& hits) const;
  void getTags(std::vector<FLASHTnTHelpers::Tag>& tags) const;
  void getProteoforms(std::vector<ProteinHit>& hits) const
  {
    for (const auto& hit : proteoform_hits_)
    {
      hits.push_back(hit);
    }
  }

protected:
  void updateMembers_() override;
  /// implemented for DefaultParamHandler
  void setDefaultParams_();

private:

  Param tagger_param_, extender_param_;
  double decoy_factor_ = 0, prsm_fdr_ = 1;
  /// how many nodes (massses) will be used
  Size max_node_cntr_ = 300;

  bool keep_decoy_ = false;
  bool keep_underdetermined_ = true;
  bool multiple_hits_per_spec_ = false;
  std::vector<ProteinHit> proteoform_hits_;
  std::vector<FLASHTnTHelpers::Tag> tags_;
  std::map<int, std::vector<int>> matching_hits_indices_;
  static bool areConsistent_(const ProteinHit& a, const ProteinHit& b, double tol) ;
  void markRepresentativeProteoformHits_(double tol);
  static void vectorizeProteinSequence_(const std::vector<std::string>& cleaned_protein_seqs,
                            std::vector<std::unordered_set<int>>& vec_pro,
                            std::vector<std::unordered_set<int>>& rev_vec_pro);
};
} // namespace OpenMS
