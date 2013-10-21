#ifndef GENOMESEQ_H
#define GENOMESEQ_H

#include <string>
using std::string;

class GenomeSeq {
public:
    GenomeSeq(const string & seq_id, const string & seq);
    ~GenomeSeq();

    void count_gc(const string & seq);

    const string & get_seq_id() const;
    double get_gc_content() const;

private:
    string seq_id;
    double gc_content;
};

#endif // GENOMESEQ_H
