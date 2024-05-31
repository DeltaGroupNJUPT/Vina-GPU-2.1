#include "kernel2.h"

//inline float dist2_cl(individual_container* list, output_type_cl* now_x, int neighbour) {
//	float out = 0;
//	for (int i = 0; i < 3; i++) {
//		float d = list->list_cl[neighbour].position[i] - now_x->position[i];
//		out += d * d;
//	}
//	for (int i = 0; i < 4; i++) {
//		float d = list->list_cl[neighbour].orientation[i] - now_x->orientation[i];
//		out += d * d;
//	}
//	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++) {
//		float d = list->list_cl[neighbour].lig_torsion[i] - now_x->lig_torsion[i];
//		out += d * d;
//	}
//	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++) {
//		float d = list->list_cl[neighbour].flex_torsion[i] - now_x->flex_torsion[i];
//		out += d * d;
//	}
//	return out;
//}

inline float dist2_g(ele_cl* list, output_type_cl* now_x, int torsion_size) {
	float out = 0;
	for (int i = 0; i < 3; i++) {
		float d = list->position[i] - now_x->position[i];
		out += d * d;
	}
	for (int i = 0; i < 4; i++) {
		float d = list->orientation[i] - now_x->orientation[i];
		out += d * d;
	}
	for (int i = 0; i < torsion_size; i++) {
		float d = list->lig_torsion[i] - now_x->lig_torsion[i];
		out += d * d;
	}
	/*for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++) {
		float d = list->flex_torsion[i] - now_x->flex_torsion[i];
		out += d * d;
	}*/
	return out;
}

inline float dist2_3D_cl(__global ele_cl* list_cl, output_type_cl* now_x) {
	float d, out = 0;
	for (int i = 0; i < 3; i++) {
		d = list_cl->position[i] - now_x->position[i];
		out += d * d;
	}
	return out;
}

int get_n(int thread, int search_depth) {
	for (int i = 0;; i += 3) {
		if (thread * search_depth < (2<<i))
			return i;
	}
}

inline void ele_init_cl(ele_cl* present_ele, output_type_cl* x_cl, float f_cl, change_cl* d_cl, int torsion_size) {
	
	for (int i = 0; i < 3; i++)present_ele->position[i] = x_cl->position[i];
	for (int i = 0; i < 4; i++)present_ele->orientation[i] = x_cl->orientation[i];
	for (int i = 0; i < torsion_size; i++)present_ele->lig_torsion[i] = x_cl->lig_torsion[i];
	/*for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++)present_ele->flex_torsion[i] = x_cl->flex_torsion[i];*/
	present_ele->energy = f_cl;
	present_ele->d_zero = 0;
	present_ele->d_positive = 0;

	const long ONE  = 1;
	long bitMask    = 0;

	for (int i = 0; i < 3; i++) {
		bitMask = ONE << i;
		if (d_cl->position[i] == 0)
			present_ele->d_zero |= bitMask;
		else if (d_cl->position[i] > 0)
			present_ele->d_positive |= bitMask;
	}
	for (int i = 0; i < 3; i++) {
		bitMask = ONE << i;
		if (d_cl->orientation[i] == 0)
			present_ele->d_zero |= bitMask;
		else if (d_cl->orientation[i] > 0)
			present_ele->d_positive |= bitMask;
	}
	for (int i = 0; i < torsion_size; i++) {
		bitMask = ONE << i;
		if (d_cl->lig_torsion[i] == 0)
			present_ele->d_zero |= bitMask;
		else if (d_cl->lig_torsion[i] > 0)
			present_ele->d_positive |= bitMask;
	}
	/*for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++) {
		bitMask = ONE << i;
		if (d_cl->flex_torsion[i] == 0)
			present_ele->d_zero |= bitMask;
		else if (d_cl->flex_torsion[i] > 0)
			present_ele->d_positive |= bitMask;
	}*/
}


inline vec3_cl origin_init(float center_x, float center_y, float center_z) {
	vec3_cl origin;
	origin.x = center_x;
	origin.y = center_y;
	origin.z = center_z;
	return origin;
}

inline vec3_cl boxsize_init(float size_x, float size_y, float size_z) {
	vec3_cl boxsize;
	boxsize.x = size_x;
	boxsize.y = size_y;
	boxsize.z = size_z;
	return boxsize;
}


inline int set_Binary_Id_cl(int thread, int search_depth, output_type_cl* x, vec3_cl* origin, vec3_cl* dimension) {
	int n = get_n(thread, search_depth) - 3;
	int binary_id = 0;
	int oct[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	// int oct[10];
	for (int i = 0; i < n / 3; i++) {
		if (x->position[0] >= origin->x) oct[i] |= 1;
		if (x->position[1] >= origin->y) oct[i] |= 2;
		if (x->position[2] >= origin->z) oct[i] |= 4;
		switch (oct[i])
		{
		case 0: origin->x = origin->x - dimension->x / (4 << (2 * i));
			origin->y = origin->y - dimension->y / (4 << (2 * i));
			origin->z = origin->z - dimension->z / (4 << (2 * i));
			break;
		case 1: origin->x = origin->x - dimension->x / (4 << (2 * i));
			origin->y = origin->y - dimension->y / (4 << (2 * i));
			origin->z = origin->z + dimension->z / (4 << (2 * i));
			break;
		case 2: origin->x = origin->x - dimension->x / (4 << (2 * i));
			origin->y = origin->y + dimension->y / (4 << (2 * i));
			origin->z = origin->z - dimension->z / (4 << (2 * i));
			break;
		case 3: origin->x = origin->x - dimension->x / (4 << (2 * i));
			origin->y = origin->y + dimension->y / (4 << (2 * i));
			origin->z = origin->z + dimension->z / (4 << (2 * i));
			break;
		case 4: origin->x = origin->x + dimension->x / (4 << (2 * i));
			origin->y = origin->y - dimension->y / (4 << (2 * i));
			origin->z = origin->z - dimension->z / (4 << (2 * i));
			break;
		case 5: origin->x = origin->x + dimension->x / (4 << (2 * i));
			origin->y = origin->y - dimension->y / (4 << (2 * i));
			origin->z = origin->z + dimension->z / (4 << (2 * i));
			break;
		case 6: origin->x = origin->x + dimension->x / (4 << (2 * i));
			origin->y = origin->y + dimension->y / (4 << (2 * i));
			origin->z = origin->z - dimension->z / (4 << (2 * i));
			break;
		case 7: origin->x = origin->x + dimension->x / (4 << (2 * i));
			origin->y = origin->y + dimension->y / (4 << (2 * i));
			origin->z = origin->z + dimension->z / (4 << (2 * i));
			break;
		default:
		    printf("Error at function set_Binary_Id_cl \n");
		}
	}
	for (int i = 0; i < n / 3; i++) {
		binary_id += (oct[i] << (n - 3 - 3 * i));
	}
	return binary_id;
}


inline void circularvisited_init_cl(individual_container* list) {
	list->n_variable = 0;
	list->p = 0;
	list->full = false;
	list->counter = 0;
}

inline void global_init_cl(global_container* g_container) {
	g_container->binary_ID = 0;
	g_container->n_variable = 0;
	g_container->counter = 0;
}

//inline void global_buffer_init_cl_with_ele(__global ele_cl* ptr, ele_cl* e) {
//	
//	/*for (int i = 0; i < 3; i++) AtomicAdd(&ptr->position[i], e->position[i]);
//	for (int i = 0; i < 4; i++) AtomicAdd(&ptr->orientation[i], e->orientation[i]);
//	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++) AtomicAdd(&ptr->lig_torsion[i], e->lig_torsion[i]);
//	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++) AtomicAdd(&ptr->flex_torsion[i], e->flex_torsion[i]);
//	AtomicAdd_long(&ptr->energy, e->energy);
//	AtomicAdd_long(&ptr->d_zero, e->d_zero);
//	AtomicAdd_long(&ptr->d_zero, e->d_zero);*/
//	for (int i = 0; i < 3; i++) ptr->position[i] = e->position[i];
//	for (int i = 0; i < 4; i++) ptr->orientation[i] = e->orientation[i];
//	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++) ptr->lig_torsion[i] = e->lig_torsion[i];
//	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++) ptr->flex_torsion[i] = e->flex_torsion[i];
//	ptr->energy = e->energy;
//	ptr->d_zero = e->d_zero;
//	ptr->d_positive = e->d_positive;
//}


inline void add_to_individual_buffer(	individual_container*		list, 
										output_type_cl*				conf_v, 
										float						f, 
										change_cl*					change_v,
										int							torsion_size
) 
{
	list->tempf = f;
	ele_cl e;
	ele_init_cl(&e, conf_v, f, change_v, torsion_size);
	if (list->counter == 0)
	{
		list->n_variable = torsion_size + 6;
	}
	if (!list->full) {
		list->list_cl[list->counter] = e;
		list->counter++;
		if (list->counter >= 5 * list->n_variable) {
			list->full = true;
			list->p = 0;
		}
	}
	else {
		list->list_cl[list->p] = e;
		list->p = (list->p + 1) % (5 * list->n_variable);
	}
}



inline void add_to_global_buffer(	global_container*		global_c, 
									int						thread, 
									int						search_depth, 
									output_type_cl*			conf_v, 
									float					f, 
									change_cl*				change_v, 
									vec3_cl*				origin, 
									vec3_cl*				dimension, 
						__global	ele_cl*					ptr,
						__global    int*					count_id,
									int						torsion_size
) 
{
	//global_c->tempf = f;
	ele_cl e;
	ele_init_cl(&e, conf_v, f, change_v, torsion_size);
	global_c->binary_ID = set_Binary_Id_cl(thread, search_depth, conf_v, origin, dimension);

	//add e to global buffer according to binary_ID
	for (int i = 0; i < 3; i++) (*(ptr + global_c->binary_ID * 8 + count_id[global_c->binary_ID])).position[i] = e.position[i];
	for (int i = 0; i < 4; i++) (*(ptr + global_c->binary_ID * 8 + count_id[global_c->binary_ID])).orientation[i] = e.orientation[i];
	for (int i = 0; i < torsion_size; i++) (*(ptr + global_c->binary_ID * 8 + count_id[global_c->binary_ID])).lig_torsion[i] = e.lig_torsion[i];
	//for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++) (*(ptr + global_c->binary_ID * 8 + count_id[global_c->binary_ID])).flex_torsion[i] = e.flex_torsion[i];
	(*(ptr + global_c->binary_ID * 8 + count_id[global_c->binary_ID])).energy = e.energy;
	(*(ptr + global_c->binary_ID * 8 + count_id[global_c->binary_ID])).d_zero = e.d_zero;
	(*(ptr + global_c->binary_ID * 8 + count_id[global_c->binary_ID])).d_positive = e.d_positive;
	count_id[global_c->binary_ID]++;
	if (count_id[global_c->binary_ID] >= 8) {
		count_id[global_c->binary_ID] %= 8;
	}
}

inline void getPointsWithinCutoff_cl(		global_container*	global_c,
											output_type_cl*		point,
											ele_cl*				results,
							__global		ele_cl*				ptr,
							__global		int*				count_id,
											int 				torsion_size
)
{
	global_c->counter = count_id[global_c->binary_ID];
	for (int m = 0; m < count_id[global_c->binary_ID]; m++) {
		for (int j = 0; j < 3; j++) results[m].position[j] = (*(ptr + global_c->binary_ID * 8 + m)).position[j];
		for (int j = 0; j < 4; j++) results[m].orientation[j] = (*(ptr + global_c->binary_ID * 8 + m)).orientation[j];
		for (int j = 0; j < torsion_size; j++) results[m].lig_torsion[j] = (*(ptr + global_c->binary_ID * 8 + m)).lig_torsion[j];
		/*for (int j = 0; j < MAX_NUM_OF_FLEX_TORSION; j++) results[m].flex_torsion[j] = (*(ptr + global_c->binary_ID * 8 + m)).flex_torsion[j];*/
		results[m].energy = (*(ptr + global_c->binary_ID * 8 + m)).energy;
		results[m].d_zero = (*(ptr + global_c->binary_ID * 8 + m)).d_zero;
		results[m].d_positive = (*(ptr + global_c->binary_ID * 8 + m)).d_positive;
	//	distances[m] = dist2_g(&results[m], point);
	}

}

inline bool check_cl(ele_cl* list_cl, output_type_cl* now_x, float now_f, change_cl* now_d,int neighbour, int torsion_size) {
	bool newXBigger, newYBigger;
	const long ONE = 1;
	long bitMask;
	newYBigger = (now_f - list_cl[neighbour].energy) > 0;

	for (int i = 0; i < 3; i++) {
		bitMask = ONE << i;
		if ((list_cl[neighbour].d_zero & bitMask) || !(now_d->position[i])) {
			continue;
		}
		else {
			const bool nowPositive = now_d->position[i] > 0;
			const bool dPositive = list_cl[neighbour].d_positive & bitMask;
			if (nowPositive ^ dPositive) {
				continue;
			}
			else {
				newXBigger = (now_x->position[i] - list_cl[neighbour].position[i]) > 0;
				if (!(nowPositive ^ newXBigger ^ newYBigger)) {
					return true;
				}
				else {
					return false;
				}
			}
		}
	}
	for (int i = 0; i < 3; i++) {
		bitMask = ONE << i;
		if ((list_cl[neighbour].d_zero & bitMask) || !(now_d->orientation[i])) {
			continue;
		}
		else {
			const bool nowPositive = now_d->orientation[i] > 0;
			const bool dPositive = list_cl[neighbour].d_positive & bitMask;
			if (nowPositive ^ dPositive) {
				continue;
			}
			else {
				newXBigger = (now_x->orientation[i] - list_cl[neighbour].orientation[i]) > 0;
				if (!(nowPositive ^ newXBigger ^ newYBigger)) {
					return true;
				}
				else {
					return false;
				}
			}
		}
	}
	for (int i = 0; i < torsion_size; i++) {
		bitMask = ONE << i;
		if ((list_cl[neighbour].d_zero & bitMask) || !(now_d->lig_torsion[i])) {
			continue;
		}
		else {
			const bool nowPositive = now_d->lig_torsion[i] > 0;
			const bool dPositive = list_cl[neighbour].d_positive & bitMask;
			if (nowPositive ^ dPositive) {
				continue;
			}
			else {
				newXBigger = (now_x->lig_torsion[i] - list_cl[neighbour].lig_torsion[i]) > 0;
				if (!(nowPositive ^ newXBigger ^ newYBigger)) {
					return true;
				}
				else {
					return false;
				}
			}
		}
	}
	/*for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++) {
		if ((list_cl[neighbour].d_zero & bitMask) || !(now_d->flex_torsion[i])) {
			continue;
		}
		else {
			const bool nowPositive = now_d->flex_torsion[i] > 0;
			const bool dPositive = list_cl[neighbour].d_positive & bitMask;
			if (nowPositive ^ dPositive) {
				continue;
			}
			else {
				newXBigger = (now_x->flex_torsion[i] - list_cl[neighbour].flex_torsion[i]) > 0;
				if (!(nowPositive ^ newXBigger ^ newYBigger)) {
					return true;
				}
				else {
					return false;
				}
			}
		}
	}*/
	return false;
}



inline int global_interesting_cl(	global_container*		global_mem,
									output_type_cl*			now_x, 
									float					now_f, 
									change_cl*				now_d, 
									int						excluded, 
									int						thread, 
									int						search_depth, 
									vec3_cl*				origin, 
									vec3_cl*				dimension, 
					__global		ele_cl*					ptr,
					__global		int*					count_id,
									int 					torsion_size

)
{
	global_mem->n_variable = torsion_size + 6;
	bool notYetChecked[MAX_NUM_OF_RESULTS];
	ele_cl nearbyPoints[MAX_NUM_OF_RESULTS];
	float distances[MAX_NUM_OF_RESULTS];
	for (int i = 0; i < MAX_NUM_OF_RESULTS; i++) {
		notYetChecked[i] = true;
	}
//	const float CUTOFF = 5.0;
	getPointsWithinCutoff_cl(global_mem,now_x, nearbyPoints, ptr, count_id, torsion_size);
//	getPointsWithinCutoff_cl(CUTOFF, global_mem, thread, search_depth, now_x, nearbyPoints, distances, origin, dimension, ptr, count_id);
	int nearbyPoints_size = global_mem->counter;
	for (int j = 0; j < nearbyPoints_size; j++) {
		distances[j] = dist2_g(&nearbyPoints[j], now_x, torsion_size);
	}
	const int grandMaxCheck = 1 * global_mem->n_variable; //1N in this case
	const int maxCheck = (nearbyPoints_size <= grandMaxCheck) ? nearbyPoints_size : grandMaxCheck;

	double min = 1e10;
	int i = 0; //counts checked done so far
	int p = 0; //pointer to current nearest point
	for (; i < maxCheck; i++) {
		min = 1e10;
		for (int j = 0; j < nearbyPoints_size; j++) {
			if (notYetChecked[j] && (distances[j] <= min)) {
				p = j;
				min = distances[p];
			}
		}
		notYetChecked[p] = false;

		if (check_cl(nearbyPoints, now_x, now_f, now_d, p, torsion_size)) {
			return -1; //i.e. return success
		}
	}
	return i;
}

inline int individual_interesting_cl(individual_container* list, output_type_cl* now_x, float now_f, change_cl* now_d, int excluded, int torsion_size) {
	int len = list->counter;
	float dist[MAX_SIZE_OF_LIST];
	bool notPicked[MAX_SIZE_OF_LIST];
	
	if (len == 0) {
		return -1; //i.e. interesting
	}
	else {
		if (!list->full) {
			return -1; //i.e. interesting
		}
		for (int i = 0; i < len; i++) {
			notPicked[i] = true;
		}
		for (int j = 0; j < len; j++) {
			dist[j] = dist2_g(&list->list_cl[j], now_x, torsion_size);
			//printf("distances:%f\n", dist[j]);
		}
		double min = 1e10;
		int p = 0;
		const int maxCheck = 4 * list->n_variable - excluded;
		int i = 0;
		for (; i < maxCheck; i++) {
			min = 1e10;
			for (int j = 0; j < len; j++) {
				if (notPicked[j] && (dist[j] < min)) {
					p = j;
					min = dist[j];
				}
			}
			notPicked[p] = false;
			if (check_cl(list->list_cl, now_x, now_f, now_d, p, torsion_size))
				return -1; //i.e. interesting
		}
		return i;
	}
}