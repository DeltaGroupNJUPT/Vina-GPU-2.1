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


/*
  变量声明以及函数原型；
  evals：unsigned；initial_factor 、min_factor、up 、down：double；
  print()函数为了实现对  evals、initial_factor 、min_factor、up 、down值的输出；
  其中对这四个变量赋了初值， evals=300、initial_factor=1e-4 、min_factor=1e-6、up=1.6 、down=0.5；
  operator（）函数原型，在ssh中会使用到；
*/

#ifndef VINA_SSD_H
#define VINA_SSD_H

#include "model.h"


/*


evals：unsigned； initial_factor、min_factor、up、down：double；
初始化：evals=300，initial_factor=1e-4、min_factor=1e-6、up=1.6、down=0.5；
声明 operator()函数原型；
print函数用于输出evals、initial_factor、min_factor、up、down的值


*/
struct ssd {
	unsigned evals;
	fl initial_factor;
	fl min_factor;
	fl up;
	fl down;
	std::vector<int> bfgs_steps;
	void print() const { std::cout << "evals=" << evals << ", initial_factor=" << initial_factor << ", min_factor=" << min_factor << ", up=" << up << ", down=" << down; }
	ssd() : evals(300), initial_factor(1e-4), min_factor(1e-6), up(1.6), down(0.5) {}
	// clean up
	void operator()(model& m, const precalculate& p, const igrid& ig, output_type& out, change& g, const vec& v) const; // g must have correct size
};

#endif
