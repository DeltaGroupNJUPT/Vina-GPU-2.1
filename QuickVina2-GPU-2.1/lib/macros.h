/*

   Copyright (c) 2004-2010, Dr. Oleg Trott <ot14@columbia.edu>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

*/

#ifndef VINA_MACROS_H
#define VINA_MACROS_H
/*   i从0循环到V的长度次  
input  V ：容器  unit类型   i：容器的下标  UINT类型       后接语句无输出   
size_t类型表示C中任何对象所能达到的最大长度，它是无符号整数。
*/
#define VINA_FOR_IN(i, v) for(std::size_t VINA_MACROS_TMP = (v).size(), (i) = 0; (i) < VINA_MACROS_TMP; ++(i))


//  i从0循环到n次        input  i 、nUINT类型        后接语句无输出    size_t类型表示C中任何对象所能达到的最大长度，它是无符号整数。
#define VINA_FOR(i, n)    for(std::size_t VINA_MACROS_TMP = (n),        (i) = 0; (i) < VINA_MACROS_TMP; ++(i))



//  i从0循环到n次            input  i、n UINT类型         后接语句无输出
#define VINA_U_FOR(i, n)  for(unsigned    VINA_MACROS_TMP = (n),        (i) = 0; (i) < VINA_MACROS_TMP; ++(i))



//   i从a循环到b次      input  i 、a 、b  UINT类型       后接语句无输出   size_t类型表示C中任何对象所能达到的最大长度，它是无符号整数。
#define VINA_RANGE(i, a, b)   for(std::size_t VINA_MACROS_TMP = (b), (i) = (a); (i) < VINA_MACROS_TMP; ++(i))


//   i从a循环到b次       input  i 、a 、b  UINT类型       后接语句无输出
#define VINA_U_RANGE(i, a, b) for(unsigned    VINA_MACROS_TMP = (b), (i) = (a); (i) < VINA_MACROS_TMP; ++(i))



//   i从a循环到b次       input  i 、a  、  b  int类型    后接语句无输出
#define VINA_I_RANGE(i, a, b) for(int         VINA_MACROS_TMP = (b), (i) = (a); (i) < VINA_MACROS_TMP; ++(i))



//  遍历容器的所有元素值   input  i 、v UINT类型  t:  vector<int>  迭代器不能改变    后接语句无输出
#define VINA_LOOP_CONST(t, i, v) for(t::const_iterator (i) = (v).begin(); (i) != (v).end(); ++(i))
//  遍历容器的所有元素值   input  i 、v UINT类型  t:  vector<int>  迭代器可以改变 
#define       VINA_LOOP(t, i, v) for(t::iterator       (i) = (v).begin(); (i) != (v).end(); ++(i))




//执行打印一次输出 x=        
#define VINA_SHOW(x)       do { std::cout << #x << " = " << (x) << std::endl; } while(false)
//执行打印一次输出 x=
#define VINA_SHOW_FAST(x)  do { std::cout << #x << " = " << (x) <<      '\n'; } while(false)
//执行打印一次输出 x=        用于迅速输出出错信息，是标准错误
#define VINA_ESHOW(x)      do { std::cerr << #x << " = " << (x) <<      '\n'; } while(false)

#endif 
