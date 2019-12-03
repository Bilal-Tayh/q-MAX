#include "Qmax.hpp"
#include <iostream>
#include <random>
#include <math.h> 

void QMax::print(){
	for (int i = 0; i < _actualsize; ++i)
		std::cout << _A[i] << " ";
	std::cout << std::endl;
	std::cout << "_phi  = " << _phi << std::endl;
}

QMax::QMax(int q, float gamma){
	_actualsize = q * (1+gamma);
	_actualsizeMinusOne = _actualsize - 1;
	_curIdx = _actualsize;
	_A = (int*) malloc(sizeof(int) * _actualsize);
	if (!_A)
		exit(1);
	_gamma = gamma;
	_q = q;
	_qMinusOne = q - 1;
	_nminusq = _actualsize - q;
	_phi = -1;
}
void QMax::insert(int v){
	if (v < _phi){
		return;
	}
	else{
		_A[--_curIdx] = v;
		if (_curIdx){
			return;
		}
		else {
			maintenance();
		}
	}
}

void QMax::maintenance(){
	_phi = findKthLargestAndPivot();
	_curIdx = _nminusq;
}

int* QMax::largestQ(){
	maintenance();
	return _A + _nminusq;
}

inline void swap(int &x, int &y){
	int z = x;
	x = y;
	y = z;
}
int QMax::PartitionAroundPivot(int left, int right, int pivot_idx, int* nums) {
	int pivot_value = nums[pivot_idx];
	int new_pivot_idx = right;
	swap(nums[pivot_idx], nums[right]);
	for (int i = right-1; i >= left; --i) {
		if (nums[i] > pivot_value) {
			swap(nums[i], nums[--new_pivot_idx]);
		}
	}
	swap(nums[right], nums[new_pivot_idx]);
	return new_pivot_idx;
}

/*
int QMax::GenerateRandom(int min,int max){
  // uniform_real_distribution documentation
  // http://www.cplusplus.com/reference/random/uniform_real_distribution/
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  double u = distribution(_generator);
  int A = max-min+1;
  for(int i=0;i<A;i++){
      if(u<(double)(i+1)/(double)A){
         return min+i;
      }
  }
  
  return 0;
}*/



int QMax::GenerateRandom(int min,int max){
  // uniform_real_distribution documentation
  // http://www.cplusplus.com/reference/random/uniform_real_distribution/
  std::uniform_real_distribution<double> distribution(min,max+1);
  double u = distribution(_generator);
  return int(u);
}







// if value dont exist return -1
int QMax::findValueIndex(int value){
    for(int i=0;i<_actualsize;i++){
        if(_A[i]==value){
            return i;
        }
    }
    return -1;
}









int QMax::checkPivot(int value, double psi){
//     int index=0;
//     int bigger=0;
//     int smaller=0;
//     for(int i=0;i<_actualsize;i++){
//         if(_A[i]==value){
//             index=i;
//         }
//         else if(_A[i]>value){
//             bigger++;
//         }
//         else smaller++;
//     } 
    
    
    
    
    
    int left = 0, right = _actualsizeMinusOne;
    int pivot_idx = findValueIndex(value);
    
    if(pivot_idx==-1){
        return -1;
    }
    
    int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _A);
    
    
//         std::cout<<"*****" <<new_pivot_idx <<"*******"<<std::endl;
//         std::cout<<"bigger= "<< bigger << ",   n-q= "<< _nminusq<<std::endl;
//         std::cout<<"should= "<< _q << ",   should = "<< _gamma*_q*psi<<std::endl;
    
    
    if (new_pivot_idx <= _nminusq) {
        if(new_pivot_idx >= _gamma*_q*psi) {  // new_pivot_idx < _q - 1.
            return _A[new_pivot_idx];
        }
    }
    return -1;
}











int QMax::findKthLargestAndPivot() {
    double delta = 1.0-0.999;
    double alpha=0.83;
    double psi = 2.0/3.0;
    int k=ceil( ((alpha*_gamma*(2+_gamma - alpha*_gamma)) / (pow(_gamma -alpha*_gamma,2))) * log(1/delta));
    int Z = (int)  ( (k*(1+_gamma)) / (alpha*_gamma) );
    
    int tries=2;
    while(tries!=0){
        
        // B should contain Z random values from _A
        int *B = (int*) malloc(sizeof(int) * Z);
        for(int i=0;i<Z;i++){
            int j=GenerateRandom(0,_actualsize);
            B[i]=_A[j];
        }
        
        
        int left1 = 0, right1 = Z-1;
        int Kth_minimal_idx=0;
        //find kth minimal value in B
        while (left1 <= right1) {
            int pivot_idx = left1;
            int new_pivot_idx = PartitionAroundPivot(left1, right1, pivot_idx, B);
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
            return _A[idx];   
        }
        // if the conditions dont hold try sample Z elemnts from _A again...
        tries--;
    }
 int left = 0, right = _actualsizeMinusOne;
	while (left <= right) {
		int pivot_idx = left;
		int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _A);
		if (new_pivot_idx == _nminusq) {
			return _A[new_pivot_idx];
		} else if (new_pivot_idx > _nminusq) {
			right = new_pivot_idx - 1;
		} else {  // new_pivot_idx < _q - 1.
			left = new_pivot_idx + 1;
		}
	}
	
}

void QMax::reset(){
	_phi = -1;
	_curIdx = _actualsize;
}

