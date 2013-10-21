#ifndef SEQPARSER_H
#define SEQPARSER_H

#include <string>
#include <list>
#include <fstream>
#include "genomeseq.h"

using std::string;
using std::list;
using std::ofstream;

class SeqParser
{
public:
    SeqParser();

    static list <GenomeSeq> * parse_fasta(const char * filename);
    static list <GenomeSeq> * parse_fastq(const char * filename);

    static void filter_fastq_by_gc(const char * in_file, ofstream & fout,
                                   double min_qc, double max_qc,
                                   const list <GenomeSeq> * seq_list);
};

#endif // SEQPARSER_H
