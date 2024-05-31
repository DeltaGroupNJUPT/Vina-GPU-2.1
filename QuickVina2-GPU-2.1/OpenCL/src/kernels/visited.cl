//#include "kernel2.h"
inline void ele_init_cl(ele_cl* present, output_type_cl* x, float f_cl, change_cl* d, int lig_torsion_size) {
	for (int i = 0; i < 3; i++)present->x_cl.position[i] = x->position[i];
	for (int i = 0; i < 4; i++)present->x_cl.orientation[i] = x->orientation[i];
	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++)present->x_cl.lig_torsion[i] = x->lig_torsion[i];
	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++)present->x_cl.flex_torsion[i] = x->flex_torsion[i];
	present->energy = f_cl;
	for (int i = 0; i < 3; i++)present->d_cl.position[i] = d->position[i];
	for (int i = 0; i < 3; i++)present->d_cl.orientation[i] = d->orientation[i];
	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++)present->d_cl.lig_torsion[i] = d->lig_torsion[i];
	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++)present->d_cl.flex_torsion[i] = d->flex_torsion[i];
	present->d_zero = 0;
	present->d_positive = 0;
	const long ONE = 1;
	long bitMask = 0;
	for (int i = 0; i < 3; i++) {
		bitMask = ONE << i;
		if (d->position[i] == 0) present->d_zero |= bitMask;
		else if (d->position[i] > 0) present->d_positive |= bitMask;
	}
	for (int i = 0; i < 3; i++) {
		bitMask = ONE << i;
		if (d->orientation[i] == 0) present->d_zero |= bitMask;
		else if (d->orientation[i] > 0) present->d_positive |= bitMask;
	}
	for (int i = 0; i < lig_torsion_size; i++) {
			bitMask = ONE << i;
			if (d->lig_torsion[i] == 0) present->d_zero |= bitMask;
			else if (d->lig_torsion[i] > 0) present->d_positive |= bitMask;
	}
	//for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++) {
	//	bitMask = ONE << i;
	//	if (d->lig_torsion[i] == 0) present->d_zero |= bitMask;
	//	else if (d->orientation[i] > 0) present->d_positive |= bitMask;
	//}

	//for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++) {
	//	bitMask = ONE << i;
	//	if (d->flex_torsion[i] == 0) present->d_zero |= bitMask;
	//	else if (d->flex_torsion[i] > 0) present->d_positive |= bitMask;
	//}
}
void visited_init(visited_cl* visited)
{
	visited->n_variable = 0;
	visited->p = 0;
	visited->full = false;
	visited->index = 0;
}
inline void add(visited_cl* visited, output_type_cl* conf_v, float f, change_cl* change_v, int lig_torsion_size)
{
	//visited->tempf = f;
	//output_type_cl tempx;
	//for (int i = 0; i < 3; i++)tempx.position[i] = conf_v->position[i];
	//for (int i = 0; i < 4; i++)tempx.orientation[i] = conf_v->orientation[i];
	//for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++)tempx.lig_torsion[i] = conf_v->lig_torsion[i];
	//for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++)tempx.flex_torsion[i] = conf_v->flex_torsion[i];
	//output_type_cl_init_with_output(&tempx, conf_v);
	ele_cl e;
	ele_init_cl(&e, conf_v, f, change_v, lig_torsion_size);
	if (visited->index == 0)
	{
		visited->n_variable = lig_torsion_size+6;
	}
	if (!visited->full)
	{
		visited->list_cl[visited->index] = &e;
		visited->index++;
		if (visited->index >= 10 * visited->n_variable)
		{
			visited->full = true;
			visited->p = 0;
		}
	}
	else
	{
		visited->list_cl[visited->p] = &e;
		visited->p = (visited->p + 1) % (10 * visited->n_variable);
	}
	//return true;
}
bool check(visited_cl* visited, output_type_cl* now_x, float now_f, change_cl* now_d, int neighbor, int lig_torsion_size)
{
	bool newXBigger, newYBigger;
	const long ONE = 1;
	long bitMask = 0;
	newYBigger = (now_f - visited->list_cl[neighbor]->energy) > 0;
	//int nowDSize = now_d->lig_torsion_size;
	for (int i = 0; i < 3; i++)
	{
		bitMask = ONE << i;
		if ((visited->list_cl[neighbor]->d_zero & bitMask) || !(now_d->position[i])) { //�������жϵķ�ǰ�� 2022.02.24

		}
		else
		{
			const bool nowPositive = now_d->position[i] > 0;
			const bool dPositive = visited->list_cl[neighbor]->d_positive & bitMask;
			if (nowPositive ^ dPositive) {

			}
			else {
				newXBigger = (now_x->position[i] - visited->list_cl[neighbor]->x_cl.position[i]) > 0;
				if (nowPositive ? (newXBigger ^ newYBigger) : (!(newXBigger ^ newYBigger))) {

				}
				else {
					return false;
				}
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		bitMask = ONE << i;
		if ((visited->list_cl[neighbor]->d_zero & bitMask) || !(now_d->orientation[i])) {

		}
		else
		{
			const bool nowPositive = now_d->orientation[i] > 0;
			const bool dPositive = visited->list_cl[neighbor]->d_positive & bitMask;
			if (nowPositive ^ dPositive) {

			}
			else {
				newXBigger = (now_x->orientation[i] - visited->list_cl[neighbor]->x_cl.orientation[i]) > 0;
				if (nowPositive ? (newXBigger ^ newYBigger) : (!(newXBigger ^ newYBigger))) {

				}
				else {
					return false;
				}
			}
		}
	}
	for (int i = 0; i < lig_torsion_size; i++)
	{
		bitMask = ONE << i;
		if ((visited->list_cl[neighbor]->d_zero & bitMask) || !(now_d->lig_torsion[i])) {

		}
		else
		{
			const bool nowPositive = now_d->position[i] > 0;
			const bool dPositive = visited->list_cl[neighbor]->d_positive & bitMask;
			if (nowPositive ^ dPositive) {

			}
			else {
				newXBigger = (now_x->lig_torsion[i] - visited->list_cl[neighbor]->x_cl.lig_torsion[i]) > 0;
				if (nowPositive ? (newXBigger ^ newYBigger) : (!(newXBigger ^ newYBigger))) {

				}
				else {
					return false;
				}
			}
		}
	}
	return false;
}
inline float dist2_cl(output_type_cl* now, visited_cl* visited, int neighbor, int lig_torsion_size) {
	float out = 0;
	for (int i = 0; i < 3; i++) {
		float d = visited->list_cl[neighbor]->x_cl.position[i] - now->position[i];
		out += d * d;
	}
	for (int i = 0; i < 4; i++) {
		float d = visited->list_cl[neighbor]->x_cl.orientation[i] - now->orientation[i];
		out += d * d;
	}
	for (int i = 0; i < lig_torsion_size; i++) {
		float d = visited->list_cl[neighbor]->x_cl.lig_torsion[i] - now->lig_torsion[i];
		out += d * d;
	}
	return out;
}
bool interesting(output_type_cl* conf_v, float f, change_cl* g, visited_cl* visited, int lig_torsion_size)
{
	int len = visited->index;
	if (visited->index == 0) {
		return true;
	}
	else {
		if (!visited->full) {
			return true;
		}
		else {
			float dist[SIZE_OF_LIST];
			bool notPicked[SIZE_OF_LIST];


			for (int i = 0; i < len; i++)
			{
				dist[i] = dist2_cl(conf_v, visited, i, lig_torsion_size);
			}


			bool flag = false;
			float min = 1e10;
			int p = 0;
			const int maxCheck = visited->n_variable;//meiwenti
			for (int i = 0; i < maxCheck; i++)
			{
				min = 1e10;
				for (int j = 0; j < len; j++) {
					if (notPicked[j] && (dist[j] < min)) {
						p = j;
						min = dist[j];
					}
				}
				notPicked[p] = false;
				flag = check(visited, conf_v, f, g, p, lig_torsion_size);
				if (flag) break;
			}
			return flag;
		}
	}
	return true;
}