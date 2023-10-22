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

#include "mutate.h"

/*
* 函数：count_mutable_entities
* input:conf类对象
* output：unsigned int类型counter
* counter大小为：2倍conf类对象的成员ligands的大小+ligands容器内每个元素（ligand_conf类）的成员torsions容器的大小之和+flex容器内每个元素（residue_conf类）的成员torsions容器的大小之和
* ligands容器和flex容器都是conf类的数据成员
*/
sz count_mutable_entities(const conf& c) {
	//unsigned int类型计数值
	sz counter = 0;

	//遍历容器ligands，该容器内元素为ligand_conf结构体
	//循环体内依次执行如下操作：1.查找容器ligands的第i个元素，其数据成员torsions（含有double元素的容器）的大小
	//2.counter=counter+第1步得到的大小+2

	VINA_FOR_IN(i, c.ligands)
		counter += 2 + c.ligands[i].torsions.size();

	//遍历容器flex，该容器内元素为residue_conf结构体
    //循环体内依次执行如下操作：1.查找容器flex的第i个元素，其数据成员torsions（含有double元素的容器）的大小
    //2.counter=counter+第1步得到的大小

	VINA_FOR_IN(i, c.flex)
		counter += c.flex[i].torsions.size();
	//返回counter
	return counter;
}


/*
* 函数：mutate_conf
* input:conf类对象，model类对象，浮点数amplitude，rng类
* output:无
* 随机对输入的conf类对象的一个数据成员赋值，具体赋值见下面详细解释
* ligands容器和flex容器都是conf类的数据成员
*/
// does not set model
void mutate_conf(conf& c, const model& m, fl amplitude, rng& generator) { // ONE OF: 2A for position, similar amp for orientation, randomize torsion
	//计算conf类对象的可变实体数，返回给mutable_entities_num
	sz mutable_entities_num = count_mutable_entities(c);
	//可变实体数为0，则退出函数，否则继续执行余下
	if(mutable_entities_num == 0) return;
	//生成一个0到mutable_entities_num - 1之间的整型随机数which_int，该随机数符合均匀分布
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	//判断整形随机数which_int是否大于等于0，否则报错，终止程序执行
	VINA_CHECK(which_int >= 0);
	//将整型数which_int转换成unsigned int类型which
	sz which = sz(which_int);
	//判断which是否小于mutable_entities_num，否则报错，终止程序执行
	VINA_CHECK(which < mutable_entities_num);

	//遍历容器ligands，该容器内元素为ligand_conf结构体
	//循环体内依次执行如下操作：
	//1.如果which等于0，
	//令容器ligands的第i个元素，其数据成员rigid（rigid_conf类）的数据成员position（vec矢量，含三个double类型数据）执行公式1，并直接退出函数
	//否则进行第2步

	//2.如果which不等于0，which先进行自减1，再判断which是否等于0，
	//如果等于0，定义一个浮点数gr，根据循环次数i赋值（具体计算公式见model文件gyration_radius函数)
	//再比较gr是否大于运行编译程序的计算机所能识别的最小非零浮点数，如果是则
	//（1）创建一个含三个double类型元素的矢量rotation，初始化rotation（公式2）,再对容器ligands的第i个元素，其数据成员rigid（rigid_conf类）的数据成员orientation依据rotation矢量利用quaternion_increment重新赋新值（具体计算公式见quaternion.cpp文件）
	//最后直接退出函数，否则进行第3步（这个否则跟第2步开头的如果配套）




	//3.which先进行自减1，再判断which是否小于{容器ligands的第i个元素，其数据成员torsions（含有double元素的容器）的大小},如果小于，
	//令容器ligands的第i个元素，其数据成员torsions（含有double元素的容器）第which个元素等于
	//-3.1415926535897931到3.1415926535897931之间的double类型随机数，该随机数符合均匀分布，并直接退出函数，
	//否则令which等于which-{容器ligands的第i个元素，其数据成员torsions（含有double元素的容器）的大小}
	//ligands容器和flex容器都是conf类的数据成员


	VINA_FOR_IN(i, c.ligands) {
		if(which == 0) { c.ligands[i].rigid.position += amplitude * random_inside_sphere(generator); return; }
		--which;
		if(which == 0) { 
			fl gr = m.gyration_radius(i); 
			if(gr > epsilon_fl) { // FIXME? just doing nothing for 0-radius molecules. do some other mutation?
				vec rotation; 
				rotation = amplitude / gr * random_inside_sphere(generator); 
				quaternion_increment(c.ligands[i].rigid.orientation, rotation);
			}
			return; 
		}
		--which;
		if(which < c.ligands[i].torsions.size()) { c.ligands[i].torsions[which] = random_fl(-pi, pi, generator); return; }
		which -= c.ligands[i].torsions.size();
	}

	//遍历容器flex，该容器内元素为residue_conf结构体
	//循环体内依次执行如下操作：1.查找容器flex的第i个元素，其数据成员torsions（含有double元素的容器）的大小
	
	//2.如果which小于第1步找到的数，令容器flex的第i个元素，其数据成员torsions（含有double元素的容器）第which个元素等于
	//-3.1415926535897931到3.1415926535897931之间的double类型随机数，该随机数符合均匀分布，并直接进行退出函数
	
	//3.如果which小于第1步找到的数，令which等于which-第1步找到的数

	VINA_FOR_IN(i, c.flex) {
		if(which < c.flex[i].torsions.size()) { c.flex[i].torsions[which] = random_fl(-pi, pi, generator); return; }
		which -= c.flex[i].torsions.size();
	}
}
