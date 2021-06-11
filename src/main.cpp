#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdint.h>

#include "btl_bloomfilter/MIBloomFilter.hpp"
#include "btl_bloomfilter/MIBFConstructSupport.hpp"
//#include "btl_bloomfilter/vendor/stHashIterator.hpp"
//#include "Common/sntHashIterator.hpp"

//#include "btl_bloomfilter/BloomFilter.hpp"

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

template<   typename ID = uint32_t, 
            typename H = stHashIterator> // H = sntHashIterator -- if multi hash kmer used
void myFunc(){
    //ID = uint32_t; 
    //H = sntHashIterator;
    // -- if using spaced seeds
    // H = stHashIterator

    std::string filePrefix = "pref";
    size_t m_expectedEntries = 2; 
    size_t m_kmerSize = 15;
    size_t hashNum = 2; 
    double occ = 0.5;
    vector<string> spaced_seeds;
    // -- if using spaced seeds
    spaced_seeds = {"111111111100001","100001111111111"}; // -- your spaced seed design

    // sequence for test purposes
    std::string sequence = "AAAAAAAAAAAAAAA";
    // create MIBFCS -- this class is essential and needed for most operations as responsible for
    // operations such as random sampling
    MIBFConstructSupport<ID, H> miBFCS(m_expectedEntries, m_kmerSize,
            hashNum, occ, spaced_seeds);
    // Parses spaced seeds
    vector<vector<unsigned> > ssVal;
    if (!spaced_seeds.empty()) {
        ssVal =	MIBloomFilter<ID>::parseSeedString(spaced_seeds);
    }
    
    // Stage 1 for bit vector insertion -------
    // calculate and allocate the memory for bit vector memory
    // Bit vector is first of 3 arrays of miBF data structure and acts as a bit bloom filter
    // miBFCS is responsible for creating and populating this vector 
    // later, this bv will be passed to miBF constructor implicitly

    // loop for sequences
    H itr = hashIterator<H>(sequence, ssVal, hashNum, m_kmerSize);
    miBFCS.insertBVColli(itr);
    // Stage 1 ends -------------



    // miBF must be created AFTER stage 1 BEFORE stage 2
    MIBloomFilter<ID> *miBF = miBFCS.getEmptyMIBF();



    // Stage 2 for ID insertion ------------
    
    // loop for seqeuences
    sequence = "AAAAAAAAAAAAAAA";
    itr = hashIterator<H>(sequence, ssVal, hashNum, m_kmerSize);
    // define the id of insertion here
    uint id = 1;
    miBFCS.insertMIBF(*miBF, itr, id);
    // Stage 2 ends -----------



    // Stage 3 for saturation -------------

    // loop for every sequence in previous loop
    sequence = "AAAAAAAAAAAAAAA";
    itr = hashIterator<H>(sequence, ssVal, hashNum, m_kmerSize);
    id = 1;
    miBFCS.insertSaturation(*miBF, itr, id);    
    // Stage 3 ends -------------




    // QUERYING STAGE --------------

    // Reusable vector for indexes
    vector<uint64_t> m_rank_pos(miBF->getHashNum());
    // Reusable vector for IDs in the indexes
	vector<ID> m_data(miBF->getHashNum());

    itr = hashIterator<H>(sequence, ssVal, hashNum, m_kmerSize);

    if(miBF->atRank(*itr,m_rank_pos)){                      // if its a hit
        m_data = miBF->getData(m_rank_pos);                 // m_data has ID's
        for(unsigned m = 0; m < miBF->getHashNum(); m++){   // iterate over ID's
            if(m_data[m] > miBF->s_mask){                     // if ID is saturated
			    // whatever the logic is for saturation
                continue;
			}
            else{
                // code here
                std::cout << "index: " << m_data[m] << std::endl;   // for test purposes
            }
        }
    }

    
    delete(miBF);
}

int main(){
    myFunc();
}