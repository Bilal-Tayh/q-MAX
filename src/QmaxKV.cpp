#include "QmaxKV.hpp"
#include <iostream>
#include <random>
#include <math.h> 

QMaxKV::QMaxKV(int q, float gamma) {
	_actualsize = q * (1+gamma);
	_actualsizeMinusOne = _actualsize - 1;
	_curIdx = _actualsize;
	_K = (key*) malloc(sizeof(key) * _actualsize);
	_V = (val*) malloc(sizeof(val) * _actualsize);
	if (!_V || !_K) {
		exit(1);
	}
	_gamma = gamma;
	_q = q;
	_qMinusOne = q - 1;
	_nminusq = _actualsize - q;
	_phi = 0;
	_k_phi = 0;
}

QMaxKV::~QMaxKV() {
	free(_K);
	free(_V);
}

val QMaxKV::getMinimalVal() {
	return _phi;
}

key QMaxKV::getMinimalKey() {
	return _k_phi;
}

void QMaxKV::update(key k, val v) {
	for(int i=0; i < _actualsize; ++i) {
		if(_K[i] == k) {
			_V[i] = v;
			return;
		}
	}
}

void QMaxKV::insert(key k, val v) {
	if (v < _phi) { // DIRECT COMPARSION: compare if new val is smaller than smallest val
		return;
	} else {
		_K[--_curIdx] = k;
		_V[_curIdx] = v;
		if (_curIdx){
			return;
		} else {
			maintenance();
		}
	}
}

void QMaxKV::maintenance() {
	_phi = findKthLargestAndPivot();
	_curIdx = _nminusq;
}

outputkv QMaxKV::largestQ() {
	outputkv out;
	maintenance();
	out.keyArr = _K + _nminusq;
	out.valArr = _V + _nminusq;
	return out;
}

inline void QMaxKV::swap(int a, int b) {
	if (a==b) return;
	key k = _K[a];
	_K[a] = _K[b];
	_K[b] = k;
	val v = _V[a];
	_V[a] = _V[b];
	_V[b] = v;
}



inline void QMaxKV::swap1(int a, int b, val* num) {
	if (a==b) return;
	val v = num[a];
	num[a] = num[b];
	num[b] = v;
}



int QMaxKV::PartitionAroundPivot(int left, int right, int pivot_idx, val* nums) {
	val pivot_value = nums[pivot_idx];
	int new_pivot_idx = right;
	swap(pivot_idx, right);
	for (int i = right-1; i >= left; --i) {
		if (nums[i] > pivot_value) { // DIRECT COMPARSION: compare if nums[i] is bigger than pivot_value
			swap(i, --new_pivot_idx);
		}
	}
	swap(right, new_pivot_idx);
//     std::cout<<_V[new_pivot_idx] <<" " << _V[right]<<" " << _V[left]<<std::endl;
	return new_pivot_idx;
}
    
    
    
    
    
    
    
    
    
    
int QMaxKV::PartitionAroundPivot1(int left, int right, int pivot_idx, val* nums) {
	val pivot_value = nums[pivot_idx];
	int new_pivot_idx = right;
	swap1(pivot_idx, right,nums);
	for (int i = right-1; i >= left; --i) {
		if (nums[i] > pivot_value) { // DIRECT COMPARSION: compare if nums[i] is bigger than pivot_value
			swap1(i, --new_pivot_idx,nums);
		}
	}
	swap1(right, new_pivot_idx,nums);
	return new_pivot_idx;
}




int QMaxKV::GenerateRandom(int min,int max){
  // uniform_real_distribution documentation
  // http://www.cplusplus.com/reference/random/uniform_real_distribution/
  std::uniform_real_distribution<double> distribution(min,max+1);
  double u = distribution(_generator);
  return int(u);
}





int QMaxKV::findValueIndex(val value){
    for(int i=0;i<_actualsize;i++){
        if(_V[i]==value){
            return i;
        }
    }
    return -1;
}





// check if the conditions holds for the possible pivot "value"
// return the pivot index in _V if it hold otherwise return -1
int QMaxKV::checkPivot(val value, double psi){
//     int index=0;
//     int bigger=0;
//     int smaller=0;
//     for(int i=0;i<_actualsize;i++){
//         if(_V[i]==value){
//             index=i;
//         }
//         else if(_V[i]>value){
//             bigger++;
//         }
//         else smaller++;
//     } 
//     
    
    
    
    
    int left = 0, right = _actualsizeMinusOne;
    int pivot_idx = findValueIndex(value);
    
     if(pivot_idx==-1){
        return -1;
    }
    
    int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _V);
    
    
//         std::cout<<"*****" <<new_pivot_idx <<"*******"<<std::endl;
//         std::cout<<"bigger= "<< bigger << ",   n-q= "<< _nminusq<<std::endl;
//         std::cout<<"should= "<< _q << ",   should = "<< _gamma*_q*psi<<std::endl;
    
    
    if (new_pivot_idx <= _nminusq) {
        if(new_pivot_idx >= _gamma*_q*psi) {  // new_pivot_idx < _q - 1.
            return _V[new_pivot_idx];
        }
    }
    return -1;
}








val QMaxKV::findKthLargestAndPivot() {
	int left = 0, right = _actualsizeMinusOne;
	while (left <= right) {
		int pivot_idx = _q;
		int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _V);
		if (new_pivot_idx == _nminusq) {
			_k_phi = _K[new_pivot_idx];
			return _V[new_pivot_idx];
		} else if (new_pivot_idx > _nminusq) {
			right = new_pivot_idx - 1;
		} else {  // new_pivot_idx < _q - 1.
			left = new_pivot_idx + 1;
		}
	}
}



/*

val QMaxKV::findKthLargestAndPivot() {
    double delta = 1.0-0.999;
    double alpha=0.83;
    double psi = 2.0/3.0;
    int k=ceil( ((alpha*_gamma*(2+_gamma - alpha*_gamma)) / (pow(_gamma -alpha*_gamma,2))) * log(1/delta));
    int Z = (int)  ( (k*(1+_gamma)) / (alpha*_gamma) );
    
    int tries=2;
    while(tries!=0){
        
        // B should contain Z random values from _V
        val *B = (val*) malloc(sizeof(val) * Z);
        for(int i=0;i<Z;i++){
            int j=GenerateRandom(0,_actualsize);
            B[i]=_V[j];
        }
        
        
        int left1 = 0, right1 = Z-1;
        int Kth_minimal_idx=0;
        //find kth minimal value in B
        while (left1 <= right1) {
            int pivot_idx = left1;
            int new_pivot_idx = PartitionAroundPivot1(left1, right1, pivot_idx, B);
            if (new_pivot_idx == k) {
                Kth_minimal_idx = new_pivot_idx;
                break;
            } else if (new_pivot_idx > k) {
                right1 = new_pivot_idx - 1;
            } else {  // new_pivot_idx < k - 1.
                left1 = new_pivot_idx + 1;
            }
        }


        
        // check if the conditions holds for the possible pivot B[Kth_minimal_idx]
        int idx = checkPivot(B[Kth_minimal_idx],psi);
        
        free(B);
        if(idx!=-1){
            return _V[idx];   
        }
        // if the conditions dont hold try sample Z elemnts from _V again...
        tries--;
    }
 int left = 0, right = _actualsizeMinusOne;
	while (left <= right) {
		int pivot_idx = left;
		int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _V);
		if (new_pivot_idx == _nminusq) {
			return _V[new_pivot_idx];
		} else if (new_pivot_idx > _nminusq) {
			right = new_pivot_idx - 1;
		} else {  // new_pivot_idx < _q - 1.
			left = new_pivot_idx + 1;
		}
	}
	
}*/
