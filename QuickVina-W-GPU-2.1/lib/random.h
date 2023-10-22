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

#ifndef VINA_RANDOM_H
#define VINA_RANDOM_H

#include <boost/random.hpp>
#include "common.h"

typedef boost::mt19937 rng;

fl random_fl(fl a, fl b, rng& generator); // 在double类型的输入a、b之间生成均匀分布的随机数，并返回函数  double类型
fl random_normal(fl mean, fl sigma, rng& generator); // 在double类型的输入mean、sigma之间生成正态分布的随机数，并返回函数  double类型
int random_int(int a, int b, rng& generator); // 在int类型的输入a、b之间生成均匀分布的随机整数，并返回函数  int类型
sz random_sz(sz a, sz b, rng& generator); // 在unsigned int类型的输入a、b之间生成均匀分布的随机整数，并返回函数 unsigned int类型
vec random_inside_sphere(rng& generator); // 当信号ture到来的时候，返回一个数组，并且数组的三个数在以0为圆心，半径为1的球体内  double
vec random_in_box(const vec& corner1, const vec& corner2, rng& generator); // 输入两个数组，产生一个随机数组，要求生成在corner1[i]-corner2[i]之间   double
int auto_seed(); // 返回时间和当前进程的标示符(PID)的乘积   int类型

#endif
