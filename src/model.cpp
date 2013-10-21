#include <cmath>
#include "model.h"
#include <iostream>

Model::Model()
{
}
void Model::build_distribution(
        int num,
        const list <GenomeSeq> * reads,
        vector <double> & distribution)
{
    distribution = vector<double>(num, 0);
    for (auto it = reads->begin(); it != reads->end(); ++it) {
        if (it->get_gc_content() > 0.01) {
            distribution[int(it->get_gc_content()*num)] += 1;
        }
    }
}

void Model::filter_residuals(
        vector <double> & distribution, int rad)
{
    int num = distribution.size();
//    std::cout << num << std::endl;
    list <pair<double, double> > minmaxes;
    size_t sum = 0;
    for (size_t i = 0; i < num - 1; ++i) {
//        std::cout << i << std::endl;
        sum += distribution[i];
        for (auto it = minmaxes.rbegin();
             it != minmaxes.rend() &&
             (it->first < distribution[i] || it->second > distribution[i]);
             ++it) {
            if (it->first < distribution[i]) {
                it->first = distribution[i];
            }
            if (it->second > distribution[i]) {
                it->second = distribution[i];
            }
        }
        minmaxes.push_back(std::make_pair(distribution[i], distribution[i]));
        if (minmaxes.size() == rad + 1) {
            sum -= distribution[i - rad];
            minmaxes.pop_front();
            if (!(minmaxes.front().first == 0 && minmaxes.front().second == 0)) {
                if (distribution[i-rad] - minmaxes.front().first <= 0.01 ||
                        minmaxes.front().second - distribution[i-rad] <= 0.01) {
                    distribution[i-rad] = sum / rad;
                }
            }
        }
    }
}

void Model::log_distribution(
        vector <double> & distribution)
{
    size_t num = distribution.size();
    for (size_t i = 0; i < num; ++i) {
        if (distribution[i] >= 0.01) {
            distribution[i] = log(distribution[i]);
        }
    }
}

list <pair <double, double> > Model::find_means(
        const list <GenomeSeq> * reads,
        int clusters_sz, int num)
{
    vector <double> gcs;
    build_distribution(num, reads, gcs);
    filter_residuals(gcs);
    log_distribution(gcs);

    /*
    for (size_t i = 0; i < num; ++i) {
        std::cout << i << ", " << gcs[i] << std::endl;
    }
    */

    vector <double> partial_means(num/clusters_sz, 0);
    list <pair <double, double> > means;
    bool increase = true;
    for (size_t i = 0; i < num/clusters_sz; ++i) {
        for (int j = 0; j < clusters_sz; ++j) {
            partial_means[i] += gcs[i*clusters_sz + j];
        }
        partial_means[i] /= clusters_sz;
        if (i != 0) {
            if (increase && partial_means[i-1] > partial_means[i]) {
                increase = false;
                means.push_back(std::make_pair((double)(i-1)*clusters_sz/num, partial_means[i-1]));
            }
            if (!increase && partial_means[i-1] < partial_means[i]) {
                increase = true;
            }
        }
    }

    return means;
}
