#include <iostream>
#include <vector>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <opencv2/ml/ml.hpp>

#include "genomeseq.h"
#include "seqparser.h"
#include "model.h"

using std::cout;
using std::endl;
using std::list;
using std::pair;

#define NUM 10000

int main(int argc, char ** argv)
{
    if (argc < 4) {
        cout << "USAGE: ./filter_qc reads.fastq gc_deviation out.fa" << std::endl;
        return -1;
    }

    list <GenomeSeq> * reads = SeqParser::parse_fastq(argv[1]);
    list <pair<double, double> > means = Model::find_means(reads, 10, 1000);
    auto max_mean = means.begin();
    for (auto it = means.begin(); it != means.end(); ++it) {
        if (it->second > max_mean->second) {
            max_mean = it;
        }
        std::cout << it->first << std::endl;
    }

    means.erase(max_mean);

    ofstream fout(argv[3]);
    for (auto it = means.begin(); it != means.end(); ++it) {
        double min_gc = it->first -
                boost::lexical_cast <double> (argv[2]) / 100;
        double max_gc = it->first +
                boost::lexical_cast <double> (argv[2]) / 100;

        SeqParser::filter_fastq_by_gc(argv[1], fout, min_gc, max_gc, reads);
    }

    delete reads;

    return 0;
}

