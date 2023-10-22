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

#ifndef VINA_CONF_H
#define VINA_CONF_H

#include <boost/ptr_container/ptr_vector.hpp> // typedef output_container

#include "quaternion.h"
#include "random.h"

/* 定义配体变化中的位置、方向、扭矩变化
	position、orientation、torsion:double
	position=position_；orientation=orientation_；torsion=torsion_
*/
struct scale {
	fl position;
	fl orientation;
	fl torsion;
	scale(fl position_, fl orientation_, fl torsion_) : position(position_), orientation(orientation_), torsion(torsion_) {}
};



/*
	ligands、flex：unint向量；

*/
struct conf_size {
	szv ligands;
	szv flex;
	/*
		sum（）函数是将向量中的元素进行累加；
		demo：vertor<int> ligands=[1,1,1,1];
			  vertor<int> flex=[2,2,2,2];
			  返回4+8+4=16；
		
	*/
	sz num_degrees_of_freedom() const {
		return sum(ligands) + sum(flex) + 6 * ligands.size();
	}
};



/*	
	input——torsions：double向量；
	该函数是将torsions向量元素都赋值为0
*/
inline void torsions_set_to_null(flv& torsions) {
	VINA_FOR_IN(i, torsions)               //for(i=0;i<torsions.size;i++)
		torsions[i] = 0;
}

/*	对torsions向量进行标准化
	input——torsions：double向量；c：double向量；factor：double；

*/
inline void torsions_increment(flv& torsions, const flv& c, fl factor) { // new torsions are normalized
	VINA_FOR_IN(i, torsions) {                            //for(i=0;i<torsions.size;i++)
		torsions[i] += normalized_angle(factor * c[i]);   //normalized_angle（）对角度进行标准化，范围为[-pi,pi];
		normalize_angle(torsions[i]);                     //normalized_angle（）对角度进行标准化，范围为[-pi,pi];
	}
}

/*	
    该函数生成[-pi,pi]的随机向量torsions（扭矩）
	input——torsions：double向量，generator：随机数生成器
*/
inline void torsions_randomize(flv& torsions, rng& generator) {
	VINA_FOR_IN(i, torsions)                        //for(i=0;i<torsions.size;i++)
		torsions[i] = random_fl(-pi, pi, generator);//生成随机数
}

/*	判断torsions1和torsions2是否相似（扭矩十分相似）
	input——torsions1、torsions2：double向量，cutoff：double；
	output——bool类型
	在torsions1和torsions2向量元素个数相等的情况下，根据（torsions1[i] - torsions2[i]）的角度标准化
	的绝对值是否大于输入cutoff，是的话返回false，否则true
*/
inline bool torsions_too_close(const flv& torsions1, const flv& torsions2, fl cutoff) {
	assert(torsions1.size() == torsions2.size());
	VINA_FOR_IN(i, torsions1)                           //for(i=0;i<torsions1.size;i++)
		if(std::abs(normalized_angle(torsions1[i] - torsions2[i])) > cutoff) 
			return false;
	return true;
}

/*	生成配体扭矩
	input——torsions：double向量，spread、rp：double，rs：double向量指针，generator：随机数生成器
	首先要判断rs是否空指针或者rs指针向量是否与torsions向量元素个数相等，是的话继续否则终止；
	if（rs部为空指针并且生成[0,1]随机数小于rp）
		将rs向量赋值给torsions向量；
	else
		torsions[i]=torsions[i]+生成[-spread, spread]随机数	
*/
inline void torsions_generate(flv& torsions, fl spread, fl rp, const flv* rs, rng& generator) {
	assert(!rs || rs->size() == torsions.size());                  //if present, rs should be the same size as torsions
	VINA_FOR_IN(i, torsions)                                       //for(i=0;i<torsions1.size;i++)
		if(rs && random_fl(0, 1, generator) < rp)                  //random_fl生成[0,1]随机数
			torsions[i] = (*rs)[i];
		else
			torsions[i] += random_fl(-spread, spread, generator); //random_fl生成[-spread, spread]随机数
}


/*
	刚体变化
*/
struct rigid_change {
	vec position;                //vec结构体，用于存放四元数虚部
	vec orientation;		     //vec结构体，用于存放四元数虚部
	/*   参数初始化
	   position结构体中： data[0] = 0；data[1] = 0;data[2] = 0;
	   orientation结构体中： data[0] = 0；data[1] = 0;data[2] = 0;
	*/
	rigid_change() : position(0, 0, 0), orientation(0, 0, 0) {}
	/*
	   逐一打印position、orientation中内容,并以“,”分割
	*/
	void print() const {
		::print(position);   //逐一打印position中内容,并以“,”分割
		::print(orientation);
	}
};

/* 
	刚体配置
*/
struct rigid_conf {
	vec position;                    //vec结构体，用于存放四元数虚部
	qt orientation;					 //四元数
	/*   参数初始化
	   position结构体中： data[0] = 0；data[1] = 0;data[2] = 0;
	   orientation：一个实部为1，虚部都为0的四元数
	*/
	rigid_conf() : position(0, 0, 0), orientation(qt_identity) {}

	/*   参数设为null
		vec结构体中： data[0] = 0；data[1] = 0;data[2] = 0;
		orientation：一个实部为1，虚部都为0的四元数
	*/
	void set_to_null() {
		position = zero_vec;
		orientation = qt_identity;
	}

/*
	input——c：rigid_change结构体，factor：double；
	inter——rotation：vec 结构体，orientation：四元数；
	c.position、c.orientation:vec结构体；
	position=(position.data[0]+factor * c.position.data[0],position.data[1]+factor * c.position.data[1], position.data[2]+factor * c.position.data[2])
	rotation=（factor * c.orientation.data[0]，factor * c.orientation.data[1]，factor * c.orientation.data[2]）
	quaternion_increment(orientation, rotation)四元数乘法以及近似标准化：
	orientation（w1,x1,y1,z1）;rotation（w2,x2,y2,z2）
	orientation * rotation =
	(w1*w2 - x1*x2 - y1*y2 - z1*z2) +
	(w1*x2 + x1*w2 + y1*z2 - z1*y2) i +
	(w1*y2 - x1*z2 + y1*w2 + z1*x2) j +
	(w1*z2 + x1*y2 - y1*x2 + z1*w2) k
*/
	void increment(const rigid_change& c, fl factor) {
		position += factor * c.position;       //factor * data[0], factor * data[1], factor * data[2]
		                                       //值依次存入vec结构体的data[0]~data[2]
		vec rotation;                          //vec结构体，用于存放四元数虚部
		rotation = factor * c.orientation;
		quaternion_increment(orientation, rotation); // orientation does not get normalized; tests show rounding errors growing very slowly
	}

/*生成随机配体的变化位置、方向
	input——corner1、corner2：vec结构体，generator：随机数生成器
*/
	void randomize(const vec& corner1, const vec& corner2, rng& generator) {
		position = random_in_box(corner1, corner2, generator);//输入两个vec结构体，产生一个[corner1[i], corner2[i]]随机vec结构体
		orientation = random_orientation(generator);
	}

	/*
		input——c：rigid_conf结构体，position_cutoff、orientation_cutoff：double；
		output——bool；
		if  vec_distance_sqr(position, c.position)计算这两个vec中逐个data数组元素差的平方大于
		sqr(position_cutoff)=position_cutoff*position_cutoff，返回false；
		if  计算quaternion_difference(orientation, c.orientation)返回的vec中逐个data数组元素的平方和 大于
		orientation_cutoff*orientation_cutoff，返回false；
		else
		返回true；
	*/
	bool too_close(const rigid_conf& c, fl position_cutoff, fl orientation_cutoff) const {
		if(vec_distance_sqr(position, c.position) > sqr(position_cutoff)) return false;
		if(sqr(quaternion_difference(orientation, c.orientation)) > sqr(orientation_cutoff)) return false;
		return true;
	}

	/*	配体位置变化
		input——spread：double，generator：随机数生成器；
		position：vec结构体，用于存放四元数虚部
		random_inside_sphere(generator)：在单位球内生成随机vec结构体tmp
	position=(position.data[0]+spread * tmp.data[0],position.data[1]+spread * tmp.data[1], position.data[2]+spread * tmp.data[2])
	*/
	void mutate_position(fl spread, rng& generator) {
		position += spread * random_inside_sphere(generator);//在单位球内生成随机vec结构体
	}

/*    配体方向变化
	input——spread：double，generator：随机数生成器；
	inter——orientation：四元数
	tmp：vec结构体，用于存放四元数虚部
	random_inside_sphere(generator)：在单位球内生成随机vec结构体tmp，这个tmp结构体和这个函数中的tmp不是一个
	tmp=(spread * tmp.data[0],spread * tmp.data[1],spread * tmp.data[2])；
	最后进行两个四元数乘法以及近似标准化：quaternion_increment(orientation, tmp)
*/
	void mutate_orientation(fl spread, rng& generator) {
		vec tmp;
		tmp = spread * random_inside_sphere(generator);
		quaternion_increment(orientation, tmp);
	}


/*  生成配体变化：位置、方向、扭矩
	input——position_spread、orientation_spread、rp：double，rs：rigid_conf结构体指针，generator：随机数生成器；
	inter—— position：vec结构体,orientation:四元数  
	if rs不为空指针并且在[0,1]之间生成的double随机数小于rp
		position=rs结构体中的position；
		orientation=rs结构体中的orientation；
	else
		position=(position.data[0]+position_spread * tmp.data[0],position.data[1]+position_spread * tmp.data[1], position.data[2]+position_spread * tmp.data[2])；
		quaternion_increment(orientation, (orientation_spread * tmp.data[0],orientation_spread * tmp.data[1],orientation_spread * tmp.data[2]))；
		（对两个四元数进行乘法以及近似标准化）
		tmp：在单位球内生成随机vec结构体；
*/
	void generate(fl position_spread, fl orientation_spread, fl rp, const rigid_conf* rs, rng& generator) {
		if(rs && random_fl(0, 1, generator) < rp)
			position = rs->position;
		else
			mutate_position(position_spread, generator);
		if(rs && random_fl(0, 1, generator) < rp)
			orientation = rs->orientation;
		else
			mutate_orientation(orientation_spread, generator);
	}
	
/*
	input——in、out：vec结构体向量，begin、end：unint
	inter——m：mat结构体,position:vec结构体
	orientation（a,b,c,d）----------->
	aa = a*a;ab = a*b;ac = a*c;ad = a*d;bb = b*b
	bc = b*c;bd = b*d;cc = c*c;cd = c*d;dd = d*d 
	----------->m结构体中data
	data[0]=(aa + bb - cc - dd)；data[3]=2 * (-ad + bc)     ；data[6]=2 * (ac + bd)      ;
	data[1]=( 2 * (ad + bc)    ；data[4]=(aa - bb + cc - dd)；data[7]= 2 * (-ab + cd)    ;
	data[2]=2 * (-ac + bd)     ；data[5]=2 * (ab + cd)      ；data[8]=(aa - bb - cc + dd);
	----------->out[i]结构体中data
	(m.data[0]*in[i][0] + m.data[3]*in[i][1] + m.data[6]*in[i][2]+position.data[0],
	 m.data[1]*in[i][0] + m.data[4]*in[i][1] + m.data[7]*in[i][2]+position.data[1],
	 m.data[2]*in[i][0] + m.data[5]*in[i][1] + m.data[8]*in[i][2]+position.data[2])
*/	
	void apply(const vecv& in, vecv& out, sz begin, sz end) const {
		assert(in.size() == out.size());   
		const mat m = quaternion_to_r3(orientation);
		VINA_RANGE(i, begin, end)                     //for(i=begin;i<end;i++)
			out[i] = m * in[i] + position;
	}

/*  打印配体的位置、方向
	逐一打印vec position中内容,并以“,”分割
	打印出orientation四元数
*/
	void print() const {
		::print(position);
		::print(orientation);
	}

private:
	friend class boost::serialization::access; //友元函数
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & position;
		ar & orientation;
	}
};

/* 
	改变配体
*/
struct ligand_change {
	rigid_change rigid;      //rigid_change结构体
	flv torsions;            //double向量                
	/*
	逐一打印vec结构体position、orientation中内容,并以“,”分割
	以及打印torsions变量
	*/
	void print() const {
		rigid.print();
		printnl(torsions);
	}
};

/*
	配体配置
*/
struct ligand_conf {
	rigid_conf rigid;		//rigid_change结构体
	flv torsions;           //double向量  
	/*   参数设为null
		rigid结构体中：
		data[0] = 0；data[1] = 0;data[2] = 0;
		orientation：一个实部为1，虚部都为0的四元数
		torsions向量元素都赋值为0
	*/
	void set_to_null() {
		rigid.set_to_null();
		torsions_set_to_null(torsions);
	}

	/*
		input——c、rigid：ligand_change结构体，factor：double，c.rigid：rigid_change结构体；
		c.rigid.position、c.rigid.orientation:vec结构体；
	  （rigid_conf结构体中）position=(position.data[0]+factor *c.rigid.position.data[0],position.data[1]+factor * c.rigid.position.data[1], position.data[2]+factor * c.rigid.position.data[2])
		rotation=（factor * c.rigid.orientation.data[0]，factor * c.rigid.orientation.data[1]，factor * c.rigid.orientation.data[2]）
		quaternion_increment(orientation（rigid_conf结构体中）, rotation)四元数乘法以及近似标准化：
		orientation（w1,x1,y1,z1）;rotation（w2,x2,y2,z2）
		orientation * rotation =
		(w1*w2 - x1*x2 - y1*y2 - z1*z2) +
		(w1*x2 + x1*w2 + y1*z2 - z1*y2) i +
		(w1*y2 - x1*z2 + y1*w2 + z1*x2) j +
		(w1*z2 + x1*y2 - y1*x2 + z1*w2) k
		最后对torsions向量进行标准化
	*/
	void increment(const ligand_change& c, fl factor) {
		rigid.increment(c.rigid, factor);   //c.rigid：rigid_change结构体
		torsions_increment(torsions, c.torsions, factor);//对torsions向量进行标准化
	}

	/*	生成随机的position结构体和orientation四元数以及[-pi,pi]的随机向量torsions
		input——rigid：ligand_change结构体，corner1、corner2：vec结构体，generator：随机数生成器
		inter——torsions：double向量，rigid：rigid_change结构体；
		
	*/
	void randomize(const vec& corner1, const vec& corner2, rng& generator) {
		rigid.randomize(corner1, corner2, generator);
		torsions_randomize(torsions, generator);
	}
	/* 打印刚性变化中的位置、方向、扭矩信息
		inter——rigid：rigid_change结构体；
		逐一打印vec position中内容,并以“,”分割和打印出orientation四元数以及打印torsions变量
	*/
	void print() const {
		rigid.print();
		printnl(torsions);
	}
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & rigid;
		ar & torsions;
	}
};

/*   残基变化
	torsions：double 向量
*/
struct residue_change {
	flv torsions;
	/*
	打印torsions变量
	*/
	void print() const {
		printnl(torsions);
	}
};

/*   残基配置
	torsions：double 向量

*/
struct residue_conf {
	flv torsions;

	/*
	input——torsions：double向量；
	该函数是将torsions向量元素都赋值为0
	*/
	void set_to_null() {
		torsions_set_to_null(torsions);
	}


	/*	对torsions向量进行标准化
	input——c：residue_change结构体；factor：double；
	*/
	void increment(const residue_change& c, fl factor) {
		torsions_increment(torsions, c.torsions, factor);
	}

	/*该函数生成[-pi,pi]的随机向量torsions（扭矩）
	 input——generator：随机数生成器
	*/
	void randomize(rng& generator) {
		torsions_randomize(torsions, generator);
	}

	/*
		打印扭矩信息
	*/
	void print() const {
		printnl(torsions);
	}
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & torsions;
	}
};

/*   
	s:conf_size结构体
	ligands：ligand_change结构体向量；
	flex：residue_change结构体向量；
	ligands向量的元素个数为conf_size结构体中的ligands uint向量的元素个数;
	flex向量的元素个数为conf_size结构体中的flex uint向量的元素个数;

	ligands[i].torsions double向量中元素个数为s.ligands[ligands.size]个，其值都为0；
	flex[i].torsions double向量中元素个数为s.flex[flex.size]个，其值都为0；

*/
struct change {
	std::vector<ligand_change> ligands;    //ligand_change结构体向量
	std::vector<residue_change> flex;
	change(const conf_size& s) : ligands(s.ligands.size()), flex(s.flex.size()) {

		VINA_FOR_IN(i, ligands)     //for (i = 0; i<ligands.size; i++)
			ligands[i].torsions.resize(s.ligands[i], 0);
		VINA_FOR_IN(i, flex)	  //for (i = 0; i<flex.size; i++)
			flex[i].torsions.resize(s.flex[i], 0);
	}

	/*
		重载change（sz），返回对象
		input——index：unint ；

		ligands：ligand_change结构体向量；
		flex：residue_change结构体向量；
		lig为ligands[i]结构体的别名，改变lig就是改变ligands[i]；
		res为flex[i]结构体的别名，改变res就是改变flex[i]


		for (i = 0; i<ligands.size; i++)
		if（index<3）
			return      rigid_change结构体中position结构体的第index个元素
		index=index-3；（上面if没有执行就一定会执行该语句，注意这里的index为无符号整型）
		if（index<3）
		return      rigid_change结构体中orientation结构体的第index个元素
		index=index-3；（上面if没有执行就一定会执行该语句，注意这里的index为无符号整型）
		if（ligands[i]结构体中torsions向量的元素个数<3）
		return      ligand_change结构体中torsions向量的第index个元素
		index=ligands[i]结构体中torsions向量的元素个数-3；（上面if没有执行就一定会执行该语句 ，注意这里的index为无符号整型）
	
	    for (i = 0; i<flex.size; i++)
		if（flex[i]结构体中torsions向量的元素个数<3）
		return      residue_change结构体中torsions向量的第index个元素
		index=flex[i]结构体中torsions向量的元素个数-3；（上面if没有执行就一定会执行该语句 ，注意这里的index为无符号整型）
		上面的return有却只会执行一遍
	*/

	fl operator()(sz index) const { // returns by value
		VINA_FOR_IN(i, ligands) {  //for (i = 0; i<ligands.size; i++)
			const ligand_change& lig = ligands[i];
			if(index < 3) return lig.rigid.position[index];
			index -= 3;
			if(index < 3) return lig.rigid.orientation[index];
			index -= 3;
			if(index < lig.torsions.size()) return lig.torsions[index];
			index -= lig.torsions.size();
		}
		VINA_FOR_IN(i, flex) {    //for (i = 0; i<flex.size; i++)
			const residue_change& res = flex[i];
			if(index < res.torsions.size()) return res.torsions[index];
			index -= res.torsions.size();
		}
		VINA_CHECK(false); 
		return 0; // shouldn't happen, placating the compiler
	}

/* 重载change（sz），返回引用
	与上面重载函数实现的内容一致，只是返回不同；
*/
	fl& operator()(sz index) {
		VINA_FOR_IN(i, ligands) {
			ligand_change& lig = ligands[i];
			if(index < 3) return lig.rigid.position[index];
			index -= 3;
			if(index < 3) return lig.rigid.orientation[index];
			index -= 3;
			if(index < lig.torsions.size()) return lig.torsions[index];
			index -= lig.torsions.size();
		}
		VINA_FOR_IN(i, flex) {
			residue_change& res = flex[i];
			if(index < res.torsions.size()) return res.torsions[index];
			index -= res.torsions.size();
		}
		VINA_CHECK(false); 
		return ligands[0].rigid.position[0]; // shouldn't happen, placating the compiler
	}


/*	计算ligand_change、residue_change结构体向量中torsions元素个数总和加上6*ligands.size的总和
	ligands：ligand_change结构体向量；
	flex：residue_change结构体向量；
*/
	sz num_floats() const {
		sz tmp = 0;
		VINA_FOR_IN(i, ligands)                       //for (i = 0; i<ligands.size; i++)
			tmp += 6 + ligands[i].torsions.size();
		VINA_FOR_IN(i, flex)                          //for (i = 0; i<flex.size; i++)
			tmp += flex[i].torsions.size();
		return tmp;
	}

	/*
	逐一打印ligands结构体向量中所有vec结构体position、orientation中内容,并以“,”分割
	以及打印其中的torsions变量。
	打印flex结构体向量其中的torsions变量。
	ligands：ligand_change结构体向量；
	flex：residue_change结构体向量；
	
	*/

	void print() const {
		VINA_FOR_IN(i, ligands)
			ligands[i].print();
		VINA_FOR_IN(i, flex)
			flex[i].print();
	}
};

/*  配置
		ligands：ligand_conf结构体向量；刚性配置
	    flex：residue_conf结构体向量；  柔性配置
		ligands向量的元素个数为conf_size结构体中的ligands unint向量的元素个数;
		flex向量的元素个数为conf_size结构体中的flex uint向量的元素个数;

		ligands[i].torsions double向量中元素个数为s.ligands[ligands.size]个，其值都为0；
		flex[i].torsions double向量中元素个数为s.flex[flex.size]个，其值都为0；
*/
struct conf {
	std::vector<ligand_conf> ligands;
	std::vector<residue_conf> flex;
	conf() {}
	conf(const conf_size& s) : ligands(s.ligands.size()), flex(s.flex.size()) {
		VINA_FOR_IN(i, ligands)
			ligands[i].torsions.resize(s.ligands[i], 0); // FIXME?
		VINA_FOR_IN(i, flex)
			flex[i].torsions.resize(s.flex[i], 0); // FIXME?
	}

	/*  参数设为null
		ligands结构体向量中：
		data[0] = 0；data[1] = 0;data[2] = 0;
		orientation：一个实部为1，虚部都为0的四元数
		torsions向量元素都赋值为0；
		flex结构体向量中：
		torsions向量元素都赋值为0
	*/
	void set_to_null() {
		VINA_FOR_IN(i, ligands)
			ligands[i].set_to_null();
		VINA_FOR_IN(i, flex)
			flex[i].set_to_null();
	}

	/*	配体刚性变化中位置、方向、扭矩以及柔性中扭矩进行标准化
		生成c.ligands[i]结构体中的position、orientation 结构体以及对ligand_change结构体中torsions向量进行标准化
	以及对flex结构体向量中的torsions向量进行标准化
    */
	void increment(const change& c, fl factor) { // torsions get normalized, orientations do not
		VINA_FOR_IN(i, ligands)
			ligands[i].increment(c.ligands[i], factor);
		VINA_FOR_IN(i, flex)
			flex[i]   .increment(c.flex[i],    factor);
	}
  /* 内部配体变化是否相似
    判断ligands结构体向量中的torsions是否和conf.结构体向量中的torsions相似
	是的话返回true，否则返回false
	input——c：conf结构体，torsions_cutoff：double
	ligands：ligand_conf结构体向量；
	flex：residue_conf结构体向量；
   */
	bool internal_too_close(const conf& c, fl torsions_cutoff) const {
		assert(ligands.size() == c.ligands.size());               //判断ligands结构体向量的元素个数是否和conf结构体中的ligands结构体向量的元素个数相等
		VINA_FOR_IN(i, ligands)                                   //for (i = 0; i<ligands.size; i++)
			if(!torsions_too_close(ligands[i].torsions, c.ligands[i].torsions, torsions_cutoff))
				return false;   //判断torsions1和torsions2是否相似
		return true;
	}
/*	外部配体变化是否相似
	input——c：conf结构体；cutoff：scale结构体
	先判断ligands[i].rigid中的position、orientation是否和c.ligands[i].rigid中的position、orientation相似，不是的话返回false；
	否则判断flex[i]中的torsions向量和c.flex[i]中的torsions向量是否相似，不是的话返回false；
	要不最后返回true
*/
	bool external_too_close(const conf& c, const scale& cutoff) const {
		assert(ligands.size() == c.ligands.size());
		VINA_FOR_IN(i, ligands)
			if(!ligands[i].rigid.too_close(c.ligands[i].rigid, cutoff.position, cutoff.orientation))
				return false;
		assert(flex.size() == c.flex.size());
		VINA_FOR_IN(i, flex)
			if(!torsions_too_close(flex[i].torsions, c.flex[i].torsions, cutoff.torsion))
				return false;
		return true;
	}
/*	外部配体和内部配体变化是否都相似
	input——c：conf结构体；cutoff：scale结构体;
	

*/
	bool too_close(const conf& c, const scale& cutoff) const {
		return internal_too_close(c, cutoff.torsion) &&
			   external_too_close(c, cutoff); // a more efficient implementation is possible, probably
	}

/*  生成配体内部变化：位置、方向、扭矩。
	input——torsion_spread、rp：double，rs:conf结构体指针，generator：随机数生成器
    位置： position（0，0，0）；
	方向： orientation（1，0，0，0）；
	扭矩：ligands[i].torsions=ligands[i].torsions+生成[-torsion_spread, torsion_spread]随机数
*/
	void generate_internal(fl torsion_spread, fl rp, const conf* rs, rng& generator) { // torsions are not normalized after this
		VINA_FOR_IN(i, ligands) {
			ligands[i].rigid.position.assign(0);
			ligands[i].rigid.orientation = qt_identity;
			const flv* torsions_rs = rs ? (&rs->ligands[i].torsions) : NULL;//rs结构体指针部为空时赋值torsions向量
			torsions_generate(ligands[i].torsions, torsion_spread, rp, torsions_rs, generator);
		}
	}


   /*  生成外部变化：刚性（位置、方向、扭矩）、柔性（扭矩）。
   input——spread：scale结构体，rp：double，rs:conf结构体指针，generator：随机数生成器
	*/

	void generate_external(const scale& spread, fl rp, const conf* rs, rng& generator) { // torsions are not normalized after this
		VINA_FOR_IN(i, ligands) {
			const rigid_conf* rigid_conf_rs = rs ? (&rs->ligands[i].rigid) : NULL;
			ligands[i].rigid.generate(spread.position, spread.orientation, rp, rigid_conf_rs, generator);
		}
		VINA_FOR_IN(i, flex) {
			const flv* torsions_rs = rs ? (&rs->flex[i].torsions) : NULL;
			torsions_generate(flex[i].torsions, spread.torsion, rp, torsions_rs, generator);
		}
	}

	/*生成刚性的随机位置、方向、扭矩和柔性中的随机扭矩；
	ligands结构体向量中生成随机的position结构体和orientation四元数以及[-pi,pi]的随机向量torsions，
	flex结构体向量中生成[-pi,pi]的随机向量torsions
	*/
	void randomize(const vec& corner1, const vec& corner2, rng& generator) {
		VINA_FOR_IN(i, ligands)
			ligands[i].randomize(corner1, corner2, generator);
		VINA_FOR_IN(i, flex)
			flex[i].randomize(generator);
	}

	/*
	 打印刚性变化中的位置、方向、扭矩信息以及柔性中的扭矩信息
	*/
	void print() const {
		VINA_FOR_IN(i, ligands)
			ligands[i].print();
		VINA_FOR_IN(i, flex)
			flex[i].print();
	}
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & ligands;
		ar & flex;
	}
};


/* 输出类型
     c——conf 结构体，别名；      e——double；   coords——vec结构体向量   
	 初始化

 */
struct output_type {
	conf c;
	fl e;
	vecv coords;
	output_type(const conf& c_, fl e_) : c(c_), e(e_) {}
};

typedef boost::ptr_vector<output_type> output_container;   //output_type结构体向量


/*
	重载<
*/
inline bool operator<(const output_type& a, const output_type& b) { // for sorting output_container
	return a.e < b.e;
}

#endif
