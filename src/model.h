#ifndef MODEL_H
#define MODEL_H
#include <vector>
#include <list>
#include "genomeseq.h"
using std::vector;
using std::list;
using std::pair;

class Model
{
public:
    Model();
    static list<pair <double, double> > find_means(
            const list <GenomeSeq> * reads,
            int clusters_sz = 10,
            int num = 1000);
    static void build_distribution(
            int num,
            const list <GenomeSeq> * reads,
            vector <double> & distribution);
    static void filter_residuals(
            vector <double> & distribution,
            int rad = 10);
    static void log_distribution(
            vector <double> & distribution);
};

#endif // MODEL_H
