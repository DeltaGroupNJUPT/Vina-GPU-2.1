/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#ifndef VINA_CURL_H
#define VINA_CURL_H

#include "common.h"

#if 1 // prefer this to "hard curl"? 一定会执行

/*
	input——e、v：double；deriv：任意类型
	先判断e > 0并且v < 0.1*max_fl(当前处理器下double型变量最大值)
	如果v<2.22045e-16 ,e=0，deriv=0否则
	e=e*v / (v + e)，deriv=deriv*v / (v + e)
*/
template<typename T> // T = fl or vec
void curl(fl& e, T& deriv, fl v) {
	if(e > 0 && not_max(v)) { // FIXME authentic_v can be gotten rid of everywhere now
		fl tmp = (v < epsilon_fl) ? 0 : (v / (v + e));
		e *= tmp;
		deriv *= sqr(tmp);
	}   
}
/*
	input——e、v：double
	先判断e > 0并且v < 0.1*max_fl(当前处理器下double型变量最大值)
		如果v<2.22045e-16 ,e=0，否则
		e=e*v / (v + e)
*/
inline void curl(fl& e, fl v) {
	if(e > 0 && not_max(v)) {           //v<0.1*max_fl(当前处理器下double型变量最大值)
		fl tmp = (v < epsilon_fl) ? 0 : (v / (v + e));    //epsilon_fl=2.22045e-16
		e *= tmp;
	}
}

#else      //被注释掉
/*
	input——e、v：double；deriv：任意类型
	如果e<v ,e=v，deriv=0
*/
template<typename T> // T = fl or vec
void curl(fl& e, T& deriv, fl v) {
	if(e > v) {
		e = v;
		deriv = 0;
	}
}

/*
	input——e、v：double
	如果e<v ,e=v
*/
inline void curl(fl& e, fl v) {
	if(e > v) {
		e = v;
	}
}
#endif

#endif
