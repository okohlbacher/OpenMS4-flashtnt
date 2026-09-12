// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>
#include <map>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace OpenMS
{
  struct FLASHTnTHelpers
  {
    /// Sequence tag used for FLASHTaggerAlgorithm. No mass gap is allowed in the seq.
    /// The mass gap containing tag should be enumerated into multiple Tag instances from outside.
    class Tag
    {
    public:
      /// constructor
      explicit Tag(std::string seq, double n_mass, double c_mass, std::vector<double>& mzs, std::vector<int>& scores, int scan);

      /// copy constructor
      Tag(const Tag&) = default;

      /// destructor
      ~Tag() = default;

      bool operator<(const Tag& a) const;
      bool operator>(const Tag& a) const;
      bool operator==(const Tag& other) const;

      [[nodiscard]] const std::string& getSequence() const;
      [[nodiscard]] const std::string& getUppercaseSequence() const;
      [[nodiscard]] double getNtermMass() const;
      [[nodiscard]] double getCtermMass() const;
      Size getLength() const;
      int getScore() const;
      int getScore(int pos) const;
      int getScan() const;
      int getIndex() const
      {
        return index_;
      }

      void setIndex(int i)
      {
        index_ = i;
      }

      float getRetentionTime() const
      {
        return rt_;
      }

      void setRetentionTime(float i)
      {
        rt_ = i;
      }

      std::string toString() const;
      const std::vector<double>& getMzs() const;

    private:
      std::string seq_, upper_seq_;
      double n_mass_ = -1, c_mass_ = -1;
      float rt_ = -1;
      std::vector<double> mzs_;
      std::vector<int> scores_;
      int scan_;
      int index_ = -1;
      Size length_;
    };

    /// A generic Directed Acyclic Graph structure for FLASHTaggerAlgorithm and FLASHExtensionAlgorithm.
    /// It is a 1-D graph so multiple dimension graph should be changed into 1-D to use this structure.
    class DAG
    {
    public:
      DAG()
      {
      }
      explicit DAG(Size vertice_count): vertex_count_(vertice_count)
      {
      }

      Size size() const
      {
        return vertex_count_;
      }

      bool addEdge(Size vertex1, Size vertex2, boost::dynamic_bitset<>& visited)
      {
        if (vertex1 >= vertex_count_ || vertex2 >= vertex_count_) return false;
        if (! visited[vertex2]) return false;
        visited[vertex1] = true;
        adj_list_for_speed_up_[vertex2].insert(vertex1); //
        return true;
      }

      /**
       * Simple path finding algorithm on DAG
       * @param source source
       * @param sink sink
       * @param all_paths the resulting paths
       * @param max_count maximum path count
       */
      void findAllPaths(Size source, Size sink, std::vector<std::vector<Size>>& all_paths, Size max_count)
      {
        std::unordered_set<Size> visited;
        std::vector<Size> path;

        if (adj_list_.empty())
        {
          for (const auto& [v2, vertices] : adj_list_for_speed_up_)
          {
            for (const auto& v1 : vertices)
            {
              adj_list_[v1].push_back(v2);
            }
          }
          for (auto& v : adj_list_)
          {
            std::sort(v.second.begin(), v.second.end());
          }
        }
        findAllPaths_(source, sink, visited, path, all_paths, max_count); // reverse traveling
      }

    private:
      Size vertex_count_ = 0;
      // 0, 1, 2, ... ,vertex_count - 1
      std::map<Size, std::vector<Size>> adj_list_; //
      std::unordered_map<Size, std::unordered_set<Size>> adj_list_for_speed_up_; // vertex swapped for fast storing of edges
      /// recursive part for the findAllPaths(..)
      void findAllPaths_(Size current,
                         Size destination,
                         std::unordered_set<Size>& visited,
                         std::vector<Size>& path,
                         std::vector<std::vector<Size>>& all_paths,
                         Size max_count)
      {
        if (max_count > 0 && all_paths.size() >= max_count) return;

        visited.insert(current);
        path.push_back(current);

        if (current == destination) { all_paths.push_back(path);}
        else
        {
          for (Size i : adj_list_[current])
          {
            if (visited.find(i) == visited.end()) { findAllPaths_(i, destination, visited, path, all_paths, max_count); }
          }
        }

        // Backtrack
        path.pop_back();
        visited.erase(current);
      }
    };

  };
}
