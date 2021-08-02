/*  This file is part of Jellyfish.

    This work is dual-licensed under 3-Clause BSD License or GPL 3.0.
    You can choose between one of them if you use this work.

`SPDX-License-Identifier: BSD-3-Clause OR  GPL-3.0`
*/


#ifndef __MER_ITERATOR_HPP__
#define __MER_ITERATOR_HPP__

#include <iterator>
#include <jellyfish/mer_dna.hpp>
#include <iostream>         // I added this

namespace jellyfish {
template<typename SequencePool, typename MerType>
class mer_iterator : public std::iterator<std::input_iterator_tag,MerType> {
  typename SequencePool::job* job_;
  const char*                 cseq_;
  MerType                     m_; // mer
  MerType                     rcm_; // reverse complement mer
  unsigned int                filled_;
  const bool                  canonical_;

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
  const mer_type& operator*() const { return !canonical_ || m_ < rcm_ ? m_ : rcm_; }
  const mer_type* operator->() const { return &this->operator*(); }
  mer_iterator& operator++() {
    while(true) {
      while(cseq_ == (*job_)->end) {
        job_->next();
        if(job_->is_empty()) {
          delete job_;
          job_  = 0;
          cseq_ = 0;
          return *this;
        }
        cseq_   = (*job_)->start;
        filled_ = 0;
      }

      do {
        int code = m_.code(*cseq_++);
          // std::cout << code << "   code on line 70 in mer_iterator.hpp\n";
          //std::cout << "I'm inside the do loop on line 71 in mer_iterator.hpp\n";

        if(code >= 0) {
          // std::cout << m_ << "   m_ before shift left on line 74 in mer_iterator.hpp\n";
          m_.shift_left(code);
          // std::cout << m_ << "   m_ after shift left on line 76 in mer_iterator.hpp\n";
          if(canonical_){
            // std::cout << rcm_ << "   rcm_ before shift right on line 78 in mer_iterator.hpp\n";
            rcm_.shift_right(rcm_.complement(code));
            // std::cout << rcm_ << "   rcm_ after shift right on line 80 in mer_iterator.hpp\n";
          }
          filled_ = std::min(filled_ + 1, mer_dna::k());
          // std::cout << filled_ << "     filled_ on line 82 in mer_iterator.hpp\n";
        } else
          filled_ = 0;
      } while(filled_ < m_.k() && cseq_ < (*job_)->end);
      // std::cout << filled_ << "     filled_ after exiting the loop line 87 in mer_iterator.hpp\n";
      if(filled_ >= m_.k())
      {
        // std::cout << filled_ << "     filled_ in the if condition on line 90 in mer_iterator.hpp\n";
        // std::cout << m_.k() << "     m_.k() in the if condition on line 91 in mer_iterator.hpp\n";
        break;
      }
    }
    // -------  std::cout << m_ << "     final m_ being returned on line 95 in mer_iterator.hpp\n";
    // -------  std::cout << rcm_ << "     final rcm_ being returned on line 96 in mer_iterator.hpp\n";
    // mer_type temporary = m_ ;
    // std::cout << temporary << "     temporary before modifying m on line 98 in mer_iterator.hpp\n";
    #if 0
    m_.alter_it();
    #endif 
    // std::cout << temporary << "     temporary after modifying m on line 98 in mer_iterator.hpp\n";
    // std::cout << m_.alter_it() << "     m_.alter_it() on line 97 in mer_iterator.hpp\n";
    // --------  std::cout << m_ << "     final m_ being returned on line 98 in mer_iterator.hpp\n\n";
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
