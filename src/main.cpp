#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdint.h>

#include "../btl_bloomfilter/MIBloomFilter.hpp"
#include "../btl_bloomfilter/MIBFConstructSupport.hpp"
#include "../btl_bloomfilter/vendor/stHashIterator.hpp"
#include "../Common/sntHashIterator.hpp"

#include "../btl_bloomfilter/BloomFilter.hpp"

#include "../Common/Options.h"

#include <tuple>
#include <google/dense_hash_map>
#include <google/sparse_hash_map>
#include <google/dense_hash_set>
#include <sdsl/int_vector.hpp>

#include <zlib.h>
#include <stdio.h>
#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
#include "../Common/kseq.h"


int main(){

    MIBFConstructSupport<ID, H> miBFCS(m_expectedEntries, m_kmerSize,
				opt::hashNum, occ, opt::sseeds);
		vector<vector<unsigned> > ssVal;
		if (!opt::sseeds.empty()) {
			ssVal =	MIBloomFilter<ID>::parseSeedString(opt::sseeds);
		}
        // calculate and allocate the memory for bit vector memory
		generateBV(miBFCS, ssVal);

        // initiate the BF with the created bit vector
        MIBloomFilter<ID> *miBF = miBFCS.getEmptyMIBF();

        // Stage 1 for insertion ------------

        // loop for seqeuences
        H itr = hashIterator<H>(sequence, ssVal);
        // define the id of insertion here
        uint id = 0;
        miBFCS.insertMIBF(*miBF, itr, id);

        // Stage 2 for saturation -------------

        H itr = hashIterator<H>(sequence, ssVal);
        id = 0;
        miBFCS.insertSaturation(*miBF, itr, id);    


        // store Bloom filter
        miBF->store(filePrefix + ".bf");

        delete(miBF);
    }