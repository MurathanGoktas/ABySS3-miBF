#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdint.h>

#include "btl_bloomfilter/MIBloomFilter.hpp"
#include "btl_bloomfilter/MIBFConstructSupport.hpp"
//#include "btl_bloomfilter/vendor/stHashIterator.hpp"
//#include "Common/sntHashIterator.hpp"

#include "btl_bloomfilter/BloomFilter.hpp"

#include "Common/Options.h"

#include <tuple>
#include <google/dense_hash_map>
#include <google/sparse_hash_map>
#include <google/dense_hash_set>
#include <sdsl/int_vector.hpp>

#include "btl_bloomfilter/vendor/stHashIterator.hpp"
#include "Common/sntHashIterator.hpp"
#include "btl_bloomfilter/MIBFConstructSupport.hpp"



#include <zlib.h>
#include <stdio.h>

template<typename H>
H hashIterator(const string &seq,
        const vector<vector<unsigned> > &seedVal,
        unsigned hashNum,
        unsigned kmerSize) {
    return H(seq, seedVal, hashNum, 1, kmerSize);
}

template<   typename ID=uint32_t, 
            typename H=sntHashIterator>
void myFunc(){
    //ID = uint32_t; 
    //H = sntHashIterator;
    // -- if using spaced seeds
    // H = stHashIterator

    std::string filePrefix = "pref";
    size_t m_expectedEntries = 1000; 
    size_t m_kmerSize = 15;
    size_t hashNum = 5; 
    double occ = 0.8;
    vector<string> spaced_seeds;
    // -- if using spaced seeds
    // spaced_seeds = {"111100001","111100001"} // -- your spaced seed design

    std::string sequence = "";
    MIBFConstructSupport<ID, H> miBFCS(m_expectedEntries, m_kmerSize,
            hashNum, occ, spaced_seeds);
    vector<vector<unsigned> > ssVal;
    if (!spaced_seeds.empty()) {
        ssVal =	MIBloomFilter<ID>::parseSeedString(spaced_seeds);
    }
    // calculate and allocate the memory for bit vector memory
    
    // Stage 1 for bit vector insertion
    H itr = hashIterator<H>(sequence, ssVal, hashNum, m_kmerSize);
    miBFCS.insertBVColli(itr);

    //generateBV(miBFCS, ssVal);

    // initiate the BF with the created bit vector
    MIBloomFilter<ID> *miBF = miBFCS.getEmptyMIBF();

    // Stage 2 for insertion ------------

    // loop for seqeuences
    sequence = "";
    itr = hashIterator<H>(sequence, ssVal, hashNum, m_kmerSize);
    // define the id of insertion here
    uint id = 0;
    miBFCS.insertMIBF(*miBF, itr, id);

    // Stage 3 for saturation -------------

    // loop for every sequence in previous loop
    sequence = "";
    itr = hashIterator<H>(sequence, ssVal, hashNum, m_kmerSize);
    id = 0;
    miBFCS.insertSaturation(*miBF, itr, id);    

    // store Bloom filter
    miBF->store(filePrefix + ".bf");

    delete(miBF);
}

int main(){

    myFunc();
}