// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTHelpers.h>
namespace OpenMS
{
  FLASHTnTHelpers::Tag::Tag(std::string seq, double n_mass, double c_mass, std::vector<double>& mzs, std::vector<int>& scores,  int scan) :
      seq_(std::move(seq)), n_mass_(n_mass), c_mass_(c_mass), mzs_(mzs), scores_(scores), scan_(scan), length_(mzs.size() - 1)
  {
    // The upstream String::toUpper() mutated seq_; keep that normalization
    // because tag deduplication and extension also inspect getSequence().
    upper_seq_ = StringUtils::toUpper(seq_);
  }

  const std::string& FLASHTnTHelpers::Tag::getSequence() const
  {
    return seq_;
  }

  const std::string& FLASHTnTHelpers::Tag::getUppercaseSequence() const
  {
    return upper_seq_;
  }

  Size FLASHTnTHelpers::Tag::getLength() const
  {
    return length_;
  }


  const std::vector<double>& FLASHTnTHelpers::Tag::getMzs() const
  {
    return mzs_;
  }

  double FLASHTnTHelpers::Tag::getNtermMass() const
  {
    return n_mass_;
  }

  double FLASHTnTHelpers::Tag::getCtermMass() const
  {
    return c_mass_;
  }

  int FLASHTnTHelpers::Tag::getScore() const
  {
    return std::accumulate(scores_.begin(), scores_.end(), 0);
  }

  int FLASHTnTHelpers::Tag::getScore(int pos) const
  {
    if (pos < 0 || pos >= scores_.size()) return 0;
    return scores_[pos];
  }

  int FLASHTnTHelpers::Tag::getScan() const
  {
    return scan_;
  }

  bool FLASHTnTHelpers::Tag::operator<(const Tag& a) const
  {
    if (this->length_ == a.length_)
    {
      if (this->seq_ == a.seq_)
      {
        if (this->c_mass_ >= 0 && a.c_mass_ >= 0) // c term tag
        {
          return this->c_mass_ < a.c_mass_;
        }
        else if (this->n_mass_ >= 0 && a.n_mass_ >= 0)
        {
          return this->n_mass_ < a.n_mass_;
        }
        else
        {
          return this->n_mass_ < a.n_mass_;
        }
      }
      return this->seq_ < a.seq_;
    }
    return this->length_ < a.length_;
  }

  bool FLASHTnTHelpers::Tag::operator>(const Tag& a) const
  {
    if (this->length_ == a.length_)
    {
      if (this->seq_ == a.seq_)
      {
        if (this->c_mass_ >= 0 && a.c_mass_ >= 0) // c term tag
        {
          return this->c_mass_ > a.c_mass_;
        }
        else if (this->n_mass_ >= 0 && a.n_mass_ >= 0)
        {
          return this->n_mass_ > a.n_mass_;
        }
        else
        {
          return this->n_mass_ > a.n_mass_;
        }
      }
      return this->seq_ > a.seq_;
    }
    return this->length_ > a.length_;
  }

  bool FLASHTnTHelpers::Tag::operator==(const Tag& a) const
  {
    return this->seq_ == a.seq_ && this->n_mass_ == a.n_mass_ && this->c_mass_ == a.c_mass_;
  }


  std::string FLASHTnTHelpers::Tag::toString() const
  {
    std::string ret;
    if (n_mass_ >= 0)
      ret = '[' + std::to_string(n_mass_) + "]\t";
    ret += seq_;
    if (c_mass_ >= 0)
      ret += "\t[" + std::to_string(c_mass_) + "]";
    ret += "\tscore : " + std::to_string(getScore()) + "\tmzs : ";

    for (auto mz : mzs_)
    {
      ret += std::to_string(mz) + " ";
    }

    return ret;
  }
} // namespace OpenMS
