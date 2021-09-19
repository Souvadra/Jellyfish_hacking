/*  This file is part of Jellyfish.

    This work is dual-licensed under 3-Clause BSD License or GPL 3.0.
    You can choose between one of them if you use this work.

`SPDX-License-Identifier: BSD-3-Clause OR  GPL-3.0`
*/


#ifndef __MER_ITERATOR_HPP__
#define __MER_ITERATOR_HPP__

#include <iterator>
#include <jellyfish/mer_dna.hpp>

namespace jellyfish {
template<typename SequencePool, typename MerType>
class mer_iterator : public std::iterator<std::input_iterator_tag,MerType> {
  typename SequencePool::job* job_;
  const char*                 cseq_;
  MerType                     m_; // mer
  MerType                     rcm_; // reverse complement mer
  unsigned int                filled_;
  const bool                  canonical_;
  uint32_t                    read_number = 0; //Souvadra's addition
  uint64_t                    kmer_int[2] = {0,0}; // Souvadra's addition
  uint64_t                    mask1 = (1ULL<<2 * m_.k()) - 1; // Souvadra's addition
  uint64_t                    shift1 = 2 * (m_.k()-1); // Souvadra's addition
  bool                        should_skip = false; // Souvadra's addition 
public:
  typedef MerType      mer_type;
  typedef SequencePool sequence_parser_type;

  mer_iterator(SequencePool& seq, bool canonical = false) :
    job_(new typename SequencePool::job(seq)), cseq_(0), filled_(0), canonical_(canonical)
  {
    if(job_->is_empty()) {
      delete job_;
      job_ = 0;
    } else {
      cseq_ = (*job_)->start;
      this->operator++();
    }
  }
  mer_iterator() : job_(0), cseq_(0), filled_(0), canonical_(false) { }
  //  mer_iterator(const mer_iterator& rhs) : job_(rhs.job_), cseq_(rhs.cseq_), m_(rhs.m_), filled_(rhs.filled_) { }
  ~mer_iterator() {
    delete job_;
  }

  bool operator==(const mer_iterator& rhs) const { return job_ == rhs.job_; }
  bool operator!=(const mer_iterator& rhs) const { return job_ != rhs.job_; }

  operator void*() const { return (void*)job_; }
  const mer_type& operator*() const { return !canonical_ || m_.get_kmer_int() < rcm_.get_kmer_int() ? m_ : rcm_; } // Souvadra's modification
  const mer_type* operator->() const { return &this->operator*(); }
  mer_iterator& operator++() {
    while(true) {
      while(cseq_ == (*job_)->end) {
        //std::cout << "line 59 @mer_iterator.hpp" <<std::endl; // Souvadra's addition
        job_->next();
        if(job_->is_empty()) {
          delete job_;
          job_  = 0;
          cseq_ = 0;
          return *this;
        }
        cseq_   = (*job_)->start;
        //std::cout << cseq_ << std::endl; // Souvadra's addition
        //std::cout << "line 67 @mer_iterator.hpp " << std::endl;
        should_skip = true; // Souvadra's addition
        filled_ = 0;
      }

      do {
        int code = m_.code(*cseq_++);
        //std::cout << code << std::endl; // Souvadra's addition
        if(code >= 0) {
          m_.shift_left(code);
          //std::cout << "filled = " << filled_ << ";  ";
          //if (filled_ >= m_.k()-1) std::cout << m_ << std::endl; // Souvadra's addition
          if(canonical_){
            rcm_.shift_right(rcm_.complement(code));
          }
          filled_ = std::min(filled_ + 1, mer_dna::k());

          kmer_int[0] = (kmer_int[0] << 2 | code) & mask1; // forward k-mer // Souvadra's addition
          if (canonical_) kmer_int[1] = (kmer_int[1] >> 2) | (3ULL^code) << shift1; // reverse k-mer // Souvadra's addition
        } else {
            //std::cout << "line 85 @mer_iterator.hpp" << std::endl;
            should_skip = false;
            read_number += 1; // // need to send this signal to count_main.cc somehow to update "rid" variable
            filled_ = 0;
        }
      } while(filled_ < m_.k() && cseq_ < (*job_)->end);
      
      if(filled_ >= m_.k())
      {
        if (should_skip == true) std::cout << "should skip" << "  --> "; // Just to check <-- Souvadra
        m_.set_skip(should_skip);
        should_skip = false;
        std::cout << m_ << std::endl;
        m_.set_rid(read_number); // Souvadra's addition
        m_.set_kmer_int(kmer_int[0]); // Souvadra's addition
        m_.set_strand(0); // Souvadra's addition
        if (canonical_) { // Souvadra's addition
          rcm_.set_rid(read_number); 
          rcm_.set_kmer_int(kmer_int[1]);
          rcm_.set_strand(1); 
        }
        break;
      }
      
    }

    #if 0
    m_.alter_it();
    #endif 
    return *this;

  }

  mer_iterator operator++(int) {
    mer_iterator res(*this);
    ++*this;
    return res;
  }
};

} // namespace jellyfish {

#endif /* __MER_ITERATOR_HPP__ */
