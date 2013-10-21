#include "seqparser.h"
#include <fstream>
#include <string>
#include <list>

#define FASTA_STR_LEN 70

using std::ifstream;
using std::ofstream;
using std::endl;
using std::string;
using std::list;

SeqParser::SeqParser()
{
}

list <GenomeSeq> * SeqParser::parse_fasta(const char * filename)
{
    return NULL;
}

list <GenomeSeq> * SeqParser::parse_fastq(const char * filename)
{
    list <GenomeSeq> * seq_list = new list <GenomeSeq>;
    ifstream fin(filename, ifstream::in);

    string seq_id;
    string seq;

    while (!fin.eof()) {
        string tmp;
        getline(fin, seq_id);
        getline(fin, seq);
        getline(fin, tmp);
        getline(fin, tmp);
        if (!seq_id.empty()) {
            seq_id.erase(seq_id.begin());
            seq_list->push_back(GenomeSeq(seq_id, seq));
        }
    }

    fin.close();

    return seq_list;
}

void SeqParser::filter_fastq_by_gc(const char * in_file, ofstream & fout,
                                   double min_gc, double max_gc,
                                   const list <GenomeSeq> * seq_list)
{
    ifstream fin(in_file, ifstream::in);

    string seq, seq_id;
    list <GenomeSeq>::const_iterator curr_seq = seq_list->begin();

    while (!fin.eof()) {
        string tmp;
        getline(fin, seq_id);
        getline(fin, seq);
        getline(fin, tmp);
        getline(fin, tmp);
        if (!seq_id.empty()) {
            seq_id.erase(seq_id.begin());
            if (curr_seq->get_seq_id().compare(seq_id) == 0 &&
                    curr_seq->get_gc_content() >= min_gc &&
                    curr_seq->get_gc_content() <= max_gc) {
                fout << '>' << seq_id << endl;
                size_t seq_len = seq.length();
                for (size_t i = 0; i < seq_len; i += FASTA_STR_LEN) {
                    tmp = seq.substr(i, FASTA_STR_LEN);
                    fout << tmp << endl;
                }
            }
            ++curr_seq;
        }
    }

    fin.close();
    fout.close();
}
