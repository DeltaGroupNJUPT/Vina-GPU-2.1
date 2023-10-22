#pragma once
void main_procedure_cl(cache& c, const std::vector<model>& ms, const precalculate& p, const parallel_mc par,
	const vec& corner1, const vec& corner2, const int seed, std::vector<output_container>& out, std::string opencl_binary_path,
	const std::vector<std::vector<std::string>> ligand_names, const int rilc_bfgs);
