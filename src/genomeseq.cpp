#include "genomeseq.h"

GenomeSeq::GenomeSeq(const string & seq_id, const string & seq)
    :seq_id(seq_id)
{
    count_gc(seq);
}

GenomeSeq::~GenomeSeq() {}

void GenomeSeq::count_gc(const string & seq)
{
    int gc = 0;
    int len = 0;
    string::const_iterator c = seq.begin();
    for (; c != seq.end(); ++c) {
        if (*c != 'N') {
            ++len;
            if (*c == 'G' || *c == 'C') {
                ++gc;
            }
        }
    }
    if (len) {
        gc_content = (double)gc / len;
    } else {
        gc_content = 0;
    }
}

const string & GenomeSeq::get_seq_id() const
{
    return seq_id;
}

double GenomeSeq::get_gc_content() const
{
    return gc_content;
}
