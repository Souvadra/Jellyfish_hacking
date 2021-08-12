/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <unistd.h>
#include <assert.h>
#include <signal.h>
#include <stdlib.h> // Souvadra's addition
#include <stdint.h> // Souvadra's addition
// #include "kvec.h"   // Souvadra's addition // don't need it actually, remove kvec.h later

#include <cctype>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <memory>
#include <chrono>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <jellyfish/err.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/locks_pthread.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>
#include <jellyfish/mer_qual_iterator.hpp>
#include <jellyfish/jellyfish.hpp>
#include <jellyfish/merge_files.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <jellyfish/generator_manager.hpp>
#include <sub_commands/count_main_cmdline.hpp>

static count_main_cmdline args; // Command line switches and arguments

namespace err = jellyfish::err;

using std::chrono::system_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

template<typename DtnType>
inline double as_seconds(DtnType dtn) { return duration_cast<duration<double>>(dtn).count(); }

using jellyfish::mer_dna;
using jellyfish::mer_dna_bloom_counter;
using jellyfish::mer_dna_bloom_filter;
typedef std::vector<const char*> file_vector;
typedef jellyfish::stream_manager<file_vector::const_iterator> stream_manager_type;

// Types for parsing arbitrary sequence ignoring quality scores
typedef jellyfish::mer_overlap_sequence_parser<stream_manager_type> sequence_parser;
typedef jellyfish::mer_iterator<sequence_parser, mer_dna> mer_iterator;

// Types for parsing reads with quality score. Interface match type
// above.
class sequence_qual_parser :
  public jellyfish::whole_sequence_parser<jellyfish::stream_manager<file_vector::const_iterator> >
{
  typedef jellyfish::stream_manager<file_vector::const_iterator> StreamIterator;
  typedef jellyfish::whole_sequence_parser<StreamIterator> super;
public:
  sequence_qual_parser(uint16_t mer_len, uint32_t max_producers, uint32_t size, size_t buf_size,
                       StreamIterator& streams) :
    super(size, 100, max_producers, streams)
  { }
};

class mer_qual_iterator : public jellyfish::mer_qual_iterator<sequence_qual_parser, mer_dna> {
  typedef jellyfish::mer_qual_iterator<sequence_qual_parser, mer_dna> super;
public:
  static char min_qual;
  mer_qual_iterator(sequence_qual_parser& parser, bool canonical = false) :
    super(parser, min_qual, canonical)
  { }
};
char mer_qual_iterator::min_qual = '!'; // Phred+33

// k-mer filters. Organized in a linked list, interpreted as a &&
// (logical and). I.e. all filter must return true for the result to
// be true. By default, filter returns true.
struct filter {
  filter* prev_;
  filter(filter* prev = 0) : prev_(prev) { }
  virtual ~filter() { }
  virtual bool operator()(const mer_dna& x) { return and_res(true, x); }
  bool and_res(bool r, const mer_dna& x) const {
    return r ? (prev_ ? (*prev_)(x) : true) : false;
  }
};

struct filter_bc : public filter {
  const mer_dna_bloom_counter& counter_;
  filter_bc(const mer_dna_bloom_counter& counter, filter* prev = 0) :
    filter(prev),
    counter_(counter)
  { }
  bool operator()(const mer_dna& m) {
    unsigned int c = counter_.check(m);
    return and_res(c > 1, m);
  }
};

struct filter_bf : public filter {
  mer_dna_bloom_filter& bf_;
  filter_bf(mer_dna_bloom_filter& bf, filter* prev = 0) :
    filter(prev),
    bf_(bf)
  { }
  bool operator()(const mer_dna& m) {
    unsigned int c = bf_.insert(m);
    return and_res(c > 0, m);
  }
};

// ****************************** Souvadra's addition starts ************************* //
unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

typedef struct { uint64_t x, y; } mm128_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;

static inline uint64_t hash64(uint64_t key, uint64_t mask) {
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

typedef struct { // a simplified version of kdq
	int front, count;
	int a[32];
} tiny_queue_t;

#define star_mers_type jellyfish::mer_dna_ns::mer_base_static<long unsigned int, 0>
class minimizer_factory {
private:
  int k, w;
  uint64_t shift1, mask;
  uint64_t kmer[2] = {0,0};
  int i, j, l, buf_pos, min_pos = 0;
  uint32_t rid = 0; //hardcoded right now, will need to change later
  int kmer_span;
  tiny_queue_t tq;
  mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
  std::vector<star_mers_type> buf_mer; //Souvadra's invention
  star_mers_type min_mer; // Souvadra's invention
public:
  bool signal; // tells Jellyfish if it should put the returned mer into an object or not
  std::vector<star_mers_type> return_mer; // Souvadra's  addition

  minimizer_factory(int k, int w) {
    assert((w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
    this->k = k;
    this->w = w;
    shift1 = 2 * (k - 1);
    mask = (1ULL<<2*k) - 1;
    kmer_span = k; // Do I even need this variable anymore ???
    memset(buf, 0xff, w * 16);
    memset(&tq, 0, sizeof(tiny_queue_t));
  }

  void buf_mer_add(star_mers_type mer, int buf_pos) {
    // a custom function that makes sure the size of the buf_mer_vector always remains
    // equal to this->w.
    if ((int)buf_mer.size() < this->w) { // Need to check if this typecast is working properly // Souvadra
      buf_mer.push_back(mer);
    } else {
      buf_mer[buf_pos] = mer;
    }
  }

  #if 1
  // new version of this code 
  void minimizer_helper(star_mers_type mer, int c) {
    mm128_t info = { UINT64_MAX, UINT64_MAX };
    star_mers_type info_mer; // Souvadra's addition
    if (c < 4) { // not an ambiguous base
      int z;
      kmer_span = l + 1 < k? l + 1 : k;
      kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
      kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
      
      z = kmer[0] <= kmer[1]? 0 : 1; // strand // Souvadra: convert the < to <= to skip dealing 
      //          with the situation where both the forward and the reverse k-mers are the same
      ++l;
      if (l >= k && kmer_span < 256) {
        info.x = hash64(kmer[z], mask) << 8 | kmer_span;
			  info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
        info_mer = mer; // Souvadra's addition
      }
    } else {
      std::cout << "HEEYYYYYY!!!!        AMBIGUOUS BASE FOUND          DO SOMETHING \n" << std::endl;
      l = 0, tq.count = tq.front = 0; kmer_span = 0; // THE CODE SHOULD NEVER COME HERE, NEVER !!
    }
    buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
    buf_mer_add(info_mer, buf_pos); // Souvadra's addition
    if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window -because identical k-mers are not stored yet
      for (j = buf_pos + 1; j < w; ++j)
        if (min.x == buf[j].x && buf[j].y != min.y) return_mer.push_back(buf_mer[j]);
      for (j = 0; j < buf_pos; ++j)
        if (min.x == buf[j].x && buf[j].y != min.y) return_mer.push_back(buf_mer[j]);
		}
    if (info.x <= min.x) { // a new minimum; then write the old min
        if (l >= w + k && min.x != UINT64_MAX) return_mer.push_back(min_mer);
        min = info, min_pos = buf_pos, min_mer = info_mer;
    } else if (buf_pos == min_pos) { // old min has moved outside the window
      if (l >= w + k - 1 && min.x != UINT64_MAX) return_mer.push_back(min_mer);
      for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
        if (min.x >= buf[j].x) min = buf[j], min_pos = j, min_mer = buf_mer[j]; //  >= is important s.t. min is always the closest k-mer
      for (j = 0; j <= buf_pos; ++j)
        if (min.x >= buf[j].x) min = buf[j], min_pos = j, min_mer = buf_mer[j];
      if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
        for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
          if (min.x == buf[j].x && min.y != buf[j].y) return_mer.push_back(buf_mer[j]);
        for (j = 0; j <= buf_pos; ++j)
          if (min.x == buf[j].x && min.y != buf[j].y) return_mer.push_back(buf_mer[j]);
			}
    }
    if (++buf_pos == w) buf_pos = 0;
  }
  #endif

  void select_minimizer(star_mers_type mer) {
    if (buf_mer.empty()) {
      std::string str = mer.to_str();
      for (i = l = 0; i < (int)str.length(); ++i) {
        int c = seq_nt4_table[(uint8_t)str[i]];
        minimizer_helper(mer, c);
      }
    } else {
      char str = mer.to_str().back();
      int c = seq_nt4_table[(uint8_t)str];
      minimizer_helper(mer,c);
    }
    signal = true;
  }

  star_mers_type last_minimizer() {
    return this->min_mer;
  }

};

// ****************************** Souvadra's addition ends ************************* //
enum OPERATION { COUNT, PRIME, UPDATE };
template<typename MerIteratorType, typename ParserType>
class mer_counter_base : public jellyfish::thread_exec {
  mer_hash&            ary_;
  ParserType           parser_;
  filter*              filter_;
  OPERATION            op_;

public:
  mer_counter_base(int nb_threads, mer_hash& ary, stream_manager_type& streams,
                   OPERATION op, filter* filter = new struct filter)
    : ary_(ary)
    , parser_(mer_dna::k(), streams.nb_streams(), 3 * nb_threads, 4096, streams)
    , filter_(filter)
    , op_(op)
  {
    ary_.reset_done();
  }
  
  virtual void start(int thid) {
    size_t count = 0;
    MerIteratorType mers(parser_, args.canonical_flag);
    minimizer_factory mmf(6,1); // k and w value hardcoded, NEET TO CHANGE
    switch(op_) {
     case COUNT:
      std::cout << "Counting Happening" << std::endl; // Souvadra's addition
      #if 1
      for (; mers; ++mers) {
        if((*filter_)(*mers)) {
          mmf.select_minimizer(*mers);
          //if((rand() % 100) / 100.0 <= (2.0 / (mers->k() + 1.0))) {                         
          if (!mmf.return_mer.empty()) {
            std::cout << "yes ";
            //std::cout << "count = " << count << ", " <<  *mers << " is being added to hash" << std::endl;          
            auto selected = mmf.return_mer.back();
            std::string ANS = selected.to_str();
            std::cout << "Received: " << ANS << std::endl; // souvadra's addition
            mmf.return_mer.pop_back();
            ary_.add(selected, 1);
          } else std::cout << "no" << std::endl;
        }
        ++count;
      }
      if (true) { // souvadra's addition
        std::cout << "yes ";
        auto last_mer = mmf.last_minimizer();
        std::string ANSWER = last_mer.to_str();
        std::cout << "Received: " << last_mer << std::endl;
        ary_.add(last_mer, 1);
      }
      #endif 
      break;

    case PRIME:
      std::cout << "Priming Happening" << std::endl; // Souvadra's addition
      for( ; mers; ++mers) {
        if((*filter_)(*mers))
          ary_.set(*mers);
        ++count;
      }
      break;

    case UPDATE:
      std::cout << "Update Happning" << std::endl; // Souvadra's addition
      mer_dna tmp;
      for( ; mers; ++mers) {
        if((*filter_)(*mers))
          ary_.update_add(*mers, 1, tmp);
        ++count;
      }
      break;
    }

    ary_.done();
  }
};

// Counter with and without quality value
typedef mer_counter_base<mer_iterator, sequence_parser> mer_counter;
typedef mer_counter_base<mer_qual_iterator, sequence_qual_parser> mer_qual_counter;

mer_dna_bloom_counter* load_bloom_filter(const char* path) {
  std::ifstream in(path, std::ios::in|std::ios::binary);
  jellyfish::file_header header(in);
  if(!in.good())
    err::die(err::msg() << "Failed to parse bloom filter file '" << path << "'");
  if(header.format() != "bloomcounter")
    err::die(err::msg() << "Invalid format '" << header.format() << "'. Expected 'bloomcounter'");
  if(header.key_len() != mer_dna::k() * 2)
    err::die("Invalid mer length in bloom filter");
  jellyfish::hash_pair<mer_dna> fns(header.matrix(1), header.matrix(2));
  auto res = new mer_dna_bloom_counter(header.size(), header.nb_hashes(), in, fns);
  if(!in.good())
    err::die("Bloom filter file is truncated");
  in.close();
  return res;
}

// If get a termination signal, kill the manager and then kill myself.
static pid_t manager_pid = 0;
static void signal_handler(int sig) {
  if(manager_pid)
    kill(manager_pid, SIGTERM);
  signal(sig, SIG_DFL);
  kill(getpid(), sig);
  _exit(EXIT_FAILURE); // Should not be reached
}

int count_main(int argc, char *argv[])
{
  auto start_time = system_clock::now();

  jellyfish::file_header header;
  header.fill_standard();
  header.set_cmdline(argc, argv);

  args.parse(argc, argv);

#ifndef HAVE_HTSLIB
  if(args.sam_given)
    count_main_cmdline::error() << "SAM/BAM/CRAM not supported (missing htslib).";
#endif


  if(args.min_qual_char_given) {
    if(args.min_qual_char_arg.size() != 1)
      count_main_cmdline::error("[-Q, --min-qual-char] must be one character.");
    const char min_qual = args.min_qual_char_arg[0];
    if(!isprint(min_qual))
      count_main_cmdline::error () << "Invalid non-printable quality character";
    if(min_qual < '!' || min_qual > '~')
      count_main_cmdline::error() << "Quality character '" << min_qual
                                  << "' is outside of the range [!, ~]";
    mer_qual_iterator::min_qual = min_qual;
  }
  if(args.min_quality_given) {
    if(args.quality_start_arg < '!' || args.quality_start_arg > '~')
      count_main_cmdline::error() << "Quality start " << args.quality_start_arg
                                  << " is outside the range [" << (int)'!' << ", "
                                  << (int)'~' << ']';
    int min_qual = args.quality_start_arg + args.min_quality_arg;
    if(min_qual < '!' || min_qual > '~')
      count_main_cmdline::error() << "Min quality " << args.min_quality_arg
                                  << " is outside the range [0, "
                                  << ((int)'~' - args.quality_start_arg) << ']';
    mer_qual_iterator::min_qual = min_qual;
  }

  mer_dna::k(args.mer_len_arg);

  std::unique_ptr<jellyfish::generator_manager> generator_manager;
  if(args.generator_given) {
    auto gm =
      new jellyfish::generator_manager(args.generator_arg, args.Generators_arg,
                                       args.shell_given ? args.shell_arg : (const char*)0);
    generator_manager.reset(gm);
    generator_manager->start();
    manager_pid = generator_manager->pid();
    struct sigaction act;
    memset(&act, '\0', sizeof(act));
    act.sa_handler = signal_handler;
    assert(sigaction(SIGTERM, &act, 0) == 0);
  }

  header.canonical(args.canonical_flag);
  mer_hash ary(args.size_arg, args.mer_len_arg * 2, args.counter_len_arg, args.threads_arg, args.reprobes_arg);
  if(args.disk_flag)
    ary.do_size_doubling(false);

  std::unique_ptr<jellyfish::dumper_t<mer_array> > dumper;
  if(args.text_flag)
    dumper.reset(new text_dumper(args.threads_arg, args.output_arg, &header));
  else
    dumper.reset(new binary_dumper(args.out_counter_len_arg, ary.key_len(), args.threads_arg, args.output_arg, &header));
  ary.dumper(dumper.get());
  
  // Counting is starting here ************************************************************************
  auto after_init_time = system_clock::now();
  
  OPERATION do_op = COUNT;
  if(args.if_given) {
    stream_manager_type streams(args.Files_arg);
    streams.paths(args.if_arg.begin(), args.if_arg.end());
    mer_counter counter(args.threads_arg, ary, streams, PRIME);
    counter.exec_join(args.threads_arg);
    do_op = UPDATE;
  }
  // Iterators to the multi pipe paths. If no generator manager,
  // generate an empty range.
  auto pipes_begin = generator_manager.get() ? generator_manager->pipes().begin() : args.file_arg.end();
  auto pipes_end = (bool)generator_manager ? generator_manager->pipes().end() : args.file_arg.end();

  stream_manager_type streams(args.Files_arg);
  streams.paths(args.file_arg.begin(), args.file_arg.end());
  streams.pipes(pipes_begin, pipes_end);
#ifdef HAVE_HTSLIB
  streams.sams(args.sam_arg.begin(), args.sam_arg.end());
#endif

  // Bloom counter read from file to filter out low frequency
  // k-mers. Two pass algorithm.
  std::unique_ptr<filter> mer_filter(new filter);
  std::unique_ptr<mer_dna_bloom_counter> bc;
  if(args.bc_given) {
    bc.reset(load_bloom_filter(args.bc_arg));
    mer_filter.reset(new filter_bc(*bc));
  }

  // Bloom filter to filter out low frequency k-mers. One pass
  // algorithm.
  std::unique_ptr<mer_dna_bloom_filter> bf;
  if(args.bf_size_given) {
    bf.reset(new mer_dna_bloom_filter(args.bf_fp_arg, args.bf_size_arg));
    mer_filter.reset(new filter_bf(*bf));
  }

  if(args.min_qual_char_given || args.min_quality_given) {
    mer_qual_counter counter(args.threads_arg, ary, streams,
                             do_op, mer_filter.get());
    counter.exec_join(args.threads_arg);
  } else {
    // std::cout << do_op << std::endl; // Souvadra's addition
    mer_counter counter(args.threads_arg, ary, streams,
                        do_op, mer_filter.get());
    std::cout << "count_main.cc @line 335" << std::endl; // Souvadra's addition
    // THIS IS WHERE THE mer_iterator IS BEING CALLED 
    counter.exec_join(args.threads_arg);
    // THE ABOVE LINE IS DOING ALL THE mer_iteration WORK
    std::cout << "count_main.cc @line 339" << std::endl; // Souvadra's addition
  }

  // If we have a manager, wait for it
  if(generator_manager) {
    signal(SIGTERM, SIG_DFL);
    manager_pid = 0;
    if(!generator_manager->wait())
      err::die("Some generator commands failed");
    generator_manager.reset();
  }

  auto after_count_time = system_clock::now();
  // Counting is stopping here ******************************************************************

  // If no intermediate files, dump directly into output file. If not, will do a round of merging
  if(!args.no_write_flag) {
    if(dumper->nb_files() == 0) {
      dumper->one_file(true);
      if(args.lower_count_given)
        dumper->min(args.lower_count_arg);
      if(args.upper_count_given)
        dumper->max(args.upper_count_arg);
      dumper->dump(ary.ary());
    } else {
      dumper->dump(ary.ary());
      if(!args.no_merge_flag) {
        std::vector<const char*> files = dumper->file_names_cstr();
        uint64_t min = args.lower_count_given ? args.lower_count_arg : 0;
        uint64_t max = args.upper_count_given ? args.upper_count_arg : std::numeric_limits<uint64_t>::max();
        try {
          merge_files(files, args.output_arg, header, min, max);
        } catch(MergeError &e) {
          err::die(err::msg() << e.what());
        }
        if(!args.no_unlink_flag) {
          for(int i =0; i < dumper->nb_files(); ++i)
            unlink(files[i]);
        }
      } // if(!args.no_merge_flag
    } // if(!args.no_merge_flag
  }

  auto after_dump_time = system_clock::now();

  if(args.timing_given) {
    std::ofstream timing_file(args.timing_arg);
    timing_file << "Init     " << as_seconds(after_init_time - start_time) << "\n"
                << "Counting " << as_seconds(after_count_time - after_init_time) << "\n"
                << "Writing  " << as_seconds(after_dump_time - after_count_time) << "\n";
  }

  return 0;
}
